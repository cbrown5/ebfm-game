#make a shiny app

#UP TO HERE
# make a performance indicators page
# allow to run multiple scenarios

library(shiny)
library(tidyverse)
source("functions.R")

theme_set(theme_classic())

#Parameters - fixed
fmort_init <- c(0, 0)
nloc <- 12
nfishery <- 3
N0 <- rep(200, nloc)
P0 <- rep(50, nloc)
R0 <- rep(1, nloc)
discount_rate <- 0.05
nonret <- 0.9 #non-retention
tmax <- 400
tmax2 <- 200
steps_per_year <- 10
irec <- paste0("V",1:3)
icomm <- 4:nloc
Dmat <- diag(1, nloc)
#Rows need to sum to 1
# Rows are from, columns are to, e.g. Dmat[1,2] is from loc 1 to loc 2
#predator has mostly retention
Dmatpred <- diag(nloc)
diag(Dmatpred) <- 1 - nonret
#set off diagonals to nonret/(nloc-1)
Dmatpred[!diag(nloc)] <- nonret/(nloc-1)

prey_params <- list(R0 = 0.8, #steepness of the stock-recruitment relationship
                    K = 200, #carrying capacity for recruitment
                    M = 0.1, #natural mortality
                    Fmort = fmort_init[1], #fishing mortality
                    Dmat = Dmat, #dispersal matrix
                    MPAs = NULL) #protected areas

pred_params <- list(R0 = 1, K = 10, M = 0.05,
 a = 0.005, Dmat = Dmatpred, 
  conversion = 0.05, 
  tmat = 3*steps_per_year, #delay between spawning and recruitment
  Fmort = fmort_init[2], MPAs = NULL)
sim <- simulate(N0, P0, R0, prey_params, pred_params, tmax, nloc)

#Create a 3x4 grid that has numbers 1-12 on each grid
loc <- expand.grid(x = 1:3, y = 1:4)
loc$iloc <- 1:nloc
loc$col <- c(rep("Recreational", 3), rep("Commercial", 9))
#Map the grid
loc_map <- ggplot(loc, aes(x = x, y = y, label = iloc)) +
    geom_tile(aes(fill = col), colour = "black", size = 3) +
    scale_fill_manual(values = c("blue2", "lightblue")) +
  geom_text(size = 9)
input <- NULL

#Create a UI that takes predator and prey fmort values and MPAvect

ui <- fluidPage(
  titlePanel("Marine ecosystem model"),
  sidebarLayout(
    sidebarPanel(
      sliderInput("prey_fmort", "Prey fishing mortality", min = 0, max = 1, value = 0.1),
      sliderInput("pred_fmort", "Predator fishing mortality", min = 0, max = 0.3, value = 0.1),
      #allow user to select MPA locations, including multiple locations
        selectInput("MPAs", "Marine protected areas", choices = 1:nloc, multiple = TRUE)
    ),
    mainPanel(
    tabsetPanel(
      tabPanel("Performance Indicators", plotOutput("perf_ind")),
      tabPanel("Abundance", plotOutput("abundance")),
      tabPanel("Catches", plotOutput("catches")),
      tabPanel("Location Map", plotOutput("loc_map"))
      
    )
    )
  )
)
server <- function(input, output) {
  
 xout <- reactive({
  pred_params2 <- pred_params
  prey_params2 <- prey_params

  prey_params2$MPAs <- as.numeric(input$MPAs)
  pred_params2$MPAs <- as.numeric(input$MPAs)
    prey_params2$Fmort <- input$prey_fmort
    pred_params2$Fmort <- input$pred_fmort
  
  sim2 <- simulate(sim$N[tmax,], sim$P[tmax,], sim$R[tmax,], prey_params2, 
                   pred_params2, tmax2, nloc)
  
  xout <- wrangle_output(sim, sim2)
  xout$catch_long$fishery <- ifelse(xout$catch_long$loc %in% irec, "Recreational", "Commercial")
  return(xout)
    })

  perf_ind <- reactive({
    xout <- xout()
    perf_ind <- xout$catch_long %>%
      filter(t > tmax) %>%
      mutate(t = t-tmax) %>%
      mutate(Sector = if_else(species == "prey" & fishery == "Recreational", "Commercial prey", paste(fishery, species))) %>%
      group_by(Sector, species) %>%
      summarize(Metric = discounted_catch(catch, discount_rate, t)) %>%
      select(-species)
  
    #function for gini coefficient
    gini <- function(x){
      n <- length(x)
      x <- sort(x)
      G <- sum((2*(1:n) - n - 1)*x)/(n*sum(x))
      return(1 - G)
    }

    #catch equity for commercial fishery
    catch_eq <- xout$catch_long %>%
      filter(t > tmax) %>%
      mutate(t = t-tmax) %>%
      mutate(Sector = if_else(species == "prey" & fishery == "Recreational", "Commercial prey", paste(fishery, species))) %>%
      filter(Sector == "Commercial predator") %>%
      group_by(loc) %>%
      summarize(NPV = discounted_catch(catch, discount_rate, t)) %>%
      summarize(Metric = gini(NPV)) %>%
      mutate(Sector = "Commercial equity")

    Perf_bio <- xout$sim_long %>%
      filter(t == max(t)) %>%
      mutate(Sector = paste("Biomass", species)) %>%
      group_by(Sector) %>%
      summarize(Metric = sum(abund)) %>%
      filter(!grepl("recruits", Sector))
      perfout <- rbind(perf_ind, catch_eq, Perf_bio)
      
      perfout$Metric <- perfout$Metric/c(500, 8000, 200, 1, 500, 8000)
    return(perfout)
  })


output$loc_map <- renderPlot({
    loc_map + 
    #change colour of MPA locations
    geom_tile(data = loc[loc$iloc %in% as.numeric(input$MPAs),], 
    aes(x = x, y = y), color = "green", size = 3, fill = NA)
  })



  output$abundance <- renderPlot({
    ggplot(xout()$sim_long, aes(t/steps_per_year, abund, colour = loc)) +
      geom_line() +
      geom_vline(xintercept = tmax/steps_per_year, linetype = "dashed") +
      facet_wrap(~species, scales = "free") +
      ylim(0, NA) +
      labs(x = "Time", y = "Abundance")
  })
  
  output$catches <- renderPlot({
    xtemp <- xout()$catch_long %>%
      group_by(t, fishery, species) %>% 
      summarize(catch = sum(catch))
    ggplot(xtemp, aes(t/steps_per_year, catch, colour = fishery)) +
      geom_line() +
            geom_vline(xintercept = tmax/steps_per_year, linetype = "dashed") +
      facet_wrap(~species, scales = "free") +
      labs(x = "Time", y = "Catches")
  })

  output$perf_ind <- renderPlot({
    ggplot(perf_ind(), aes(Sector, Metric, fill = Sector)) +
      geom_bar(stat = "identity", position = "dodge")
  })
  
}

shinyApp(ui, server)

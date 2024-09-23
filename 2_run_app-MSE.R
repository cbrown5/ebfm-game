#make a shiny app
# MSE application

library(shiny)
library(tidyverse)
source("functions.R")


#increase font size for all ggplots
theme_set(theme_classic(base_size = 25))

#Parameters - fixed
fmort_init <- c(0.25, 0.25)
nloc <- 1
nfishery <- 1
N0 <- rep(200, nloc)
P0 <- rep(50, nloc)
R0 <- rep(1, nloc)

B0 <- 1620#calculated given paramaters below. 
discount_rate <- 0.05
nonret <- 0.9 #non-retention
tmax <- 200
tmax2 <- 200
steps_per_year <- 10
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
                    recruit_sd = 0.2, #recruitment variability
                    MPAs = NULL) #protected areas

pred_params <- list(R0 = 1, K = 10, M = 0.05,
                    a = 0.001, Dmat = Dmatpred, 
                    conversion = 0.05, 
                    tmat = 3*steps_per_year, #delay between spawning and recruitment
                    Fmort = fmort_init[2], MPAs = NULL)

input <- NULL

#Create a UI that takes predator and prey fmort values and MPAvect

ui <- fluidPage(
  titlePanel("Marine ecosystem model"),
  sidebarLayout(
    sidebarPanel(
      sliderInput("prey_fmort", "Status quo fishing mortality", min = 0, max = 1, value = 0.25),
      sliderInput("prey_fmort_max", "Max fishing mortality", min = 0, max = 1, value = 0.25),
      numericInput("nsims", "Number of simulations", value = 100, min = 10, max = 1000, step = 50),
      sliderInput("TRP", "Target reference point", min = 0, max = 1, value = 0.4),
      sliderInput("LRP", "Limit reference point", min = 0, max = 1, value = 0.2),
      sliderInput("obs_err", "Observation error", min = 0, max = 1, value = 0.1)
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("Performance Indicators", tableOutput("perf_ind")),
        tabPanel("Broken stick plot", plotOutput("broken_stick_plot")),
        tabPanel("Depletion", plotOutput("Brel")),
        tabPanel("Catches", plotOutput("catches")),
        tabPanel("Biomass", plotOutput("abundance")),
        tabPanel("Biomass lines", plotOutput("abundance_lines"))
      )
    )
  )
)

server <- function(input, output) {
  
  xout <- reactive({
    pred_params2 <- pred_params
    prey_params2 <- prey_params
    prey_params$Fmort <- input$prey_fmort
    prey_params2$Fmort <- input$prey_fmort

    catch_all <- data.frame(t = rep(1:(tmax+tmax2), input$nsims),
                            isim = rep(1:input$nsims, each = tmax+tmax2),
                            catch = NA)
    abund_all <- data.frame(t = rep(1:(tmax+tmax2), input$nsims),
                             isim = rep(1:input$nsims, each = tmax+tmax2),
                             N = NA)
    for (i in 1:input$nsims){
      epsilon <- rep(0, tmax + tmax2)
      for (t in 2:(tmax + tmax2)){
        epsilon[t] <- epsilon[t-1]*0.8 + rnorm(1, 0, prey_params$recruit_sd)
      }
      sim <- simulate(N0, P0, R0, prey_params, pred_params, tmax, nloc, exp(epsilon[1:tmax]))
      sim2 <- simulate2(sim$N[tmax,], sim$P[tmax,], sim$R[tmax,], prey_params2, 
                        pred_params2, tmax2, nloc, input$obs_err, B0, input$prey_fmort_max, input$TRP, input$LRP,
                        exp(epsilon[(tmax+1):(tmax+tmax2)]))
      catch_all[catch_all$isim==i,"catch"] <- c(sim$Ncatch, sim2$Ncatch)
      abund_all[abund_all$isim==i,"N"] <- c(sim$N, sim2$N)
    }
    abund_all$Brel <- abund_all$N/B0
    xout <- list(catch_long = catch_all, sim_long = abund_all)
    return(xout)
  })
  
  perf_ind <- reactive({
    
    xout2 <- xout()
    catchquant <- xout2$catch_long %>%
      group_by(t) %>%
      summarize(med = median(catch),
                q25 = quantile(catch, 0.25),
                q75 = quantile(catch, 0.75))
    
    Nquant <- xout2$sim_long %>%
      group_by(t) %>%
      summarize(med = median(N),
                q25 = quantile(N, 0.25),
                q75 = quantile(N, 0.75))
    
    Brelquant <- xout2$sim_long %>%
      group_by(t) %>%
      summarize(med = median(Brel),
                q25 = quantile(Brel, 0.25),
                q75 = quantile(Brel, 0.75))
    final_Brel <- xout2$sim_long %>%
      filter(t == (tmax+tmax2)) %>%
      summarize(x = sum(Brel < 0.2)/input$nsims)
    
    #mean catch
    perfout <- data.frame(Metric = c("Median catch", "Median biomass", "Median B/B0", "Risk"), 
                          Value = c(pull(catchquant[nrow(catchquant), "med"]),
                                    pull(Nquant[nrow(Nquant), "med"]),
                                    pull(Brelquant[nrow(Brelquant), "med"]),
                                    pull(final_Brel, x)))
    
    return(list(perfout = perfout, catchquant = catchquant, Nquant= Nquant, Brelquant = Brelquant))
  })
  
  output$Brel <- renderPlot({
    ggplot(perf_ind()$Brelquant, aes(t/steps_per_year, med)) +
      geom_line() +
      geom_hline(yintercept = input$LRP, linetype = "dashed", color = "orange") +
      geom_hline(yintercept = input$TRP, linetype = "dashed", color = "blue") +
      geom_ribbon(aes(ymin = q25, ymax = q75), alpha = 0.3) +
      geom_vline(xintercept = tmax/steps_per_year, linetype = "dashed") +
      ylim(0, NA) +
      labs(x = "Time", y = "Depletion")
  })
  
  output$catches <- renderPlot({
    ggplot(perf_ind()$catchquant, aes(t/steps_per_year, med)) +
      geom_line() +
      geom_ribbon(aes(ymin = q25, ymax = q75), alpha = 0.3) +
      ylim(0, NA) +
      labs(x = "Time", y = "Catch")
  })
  
  output$abundance <- renderPlot({
    ggplot(perf_ind()$Nquant, aes(t/steps_per_year, med)) +
      geom_line() +
      geom_ribbon(aes(ymin = q25, ymax = q75), alpha = 0.3) +
      geom_vline(xintercept = tmax/steps_per_year, linetype = "dashed") +
      ylim(0, NA) +
      labs(x = "Time", y = "Biomass")
  })
  
  output$abundance_lines <- renderPlot({
    ggplot(filter(xout()$sim_long, isim <101), aes(t/steps_per_year, N, group = isim)) +
      geom_line(alpha = 0.3) +
      geom_vline(xintercept = tmax/steps_per_year, linetype = "dashed") +
      ylim(0, NA) +
      labs(x = "Time", y = "Biomass")
  })
  
  output$broken_stick_plot <- renderPlot({
   x <- seq(0, 1, length.out = 100)
   y <- rep(NA, 100)
   for (i in 1:100){
     y[i] <- broken_stick(x[i], input$prey_fmort_max, input$TRP, input$LRP)
   }
   ggplot(data.frame(x = x, y = y), aes(x, y)) +
     geom_line() +
     geom_hline(yintercept = input$prey_fmort_max, linetype = "dashed") +
     geom_vline(xintercept = input$LRP, linetype = "dashed", color = "orange") +
     geom_vline(xintercept = input$TRP, linetype = "dashed", color = "blue") +
     labs(x = "Depletion", y = "F mortality")
  })
  
  # 
  # output$catches <- renderPlot({
  #   xtemp <- xout()$catch_long %>%
  #     mutate(Sector = if_else(species == "prey" & fishery == "Recreational", "Commercial prey", paste(fishery, species))) %>%
  #     group_by(t, Sector) %>% 
  #     summarize(catch = sum(catch))
  #   ggplot(xtemp, aes(t/steps_per_year, catch, colour = Sector)) +
  #     geom_line() +
  #     geom_vline(xintercept = tmax/steps_per_year, linetype = "dashed") +
  #     facet_wrap(~Sector, scales = "free") +
  #     labs(x = "Time", y = "Catches") + 
  #     #change facet label size
  #     theme(strip.text = element_text(size = 15), legend.position = "none")
  # })
  
  # output$perf_ind <- renderPlot({
  #   ggplot(perf_ind()$perfout, aes(Metric, Value)) +
  #     geom_bar(stat = "identity", position = "dodge") +
  #     #rotate x-axis text
  #     facet_grid(Metric~., scales = "free") +
  #     theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  #     ylim(0, 1.1)
  # })
  
  #make a table of the performance indicators
  output$perf_ind <- renderTable({
    perf_ind()$perfout
  })
  
}

shinyApp(ui, server)

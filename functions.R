#Create a model of a population with two lifestages and beverton holt recrutiment function

# Function to calculate the population size at time t+1
vbh <- function(Nt, R0, K) {
  Nt1 <- (R0 * Nt) / (1 + (R0*Nt/K))
  return(Nt1)
}


#function for gini coefficient
gini <- function(x){
  n <- length(x)
  x <- sort(x)
  G <- sum((2*(1:n) - n - 1)*x)/(n*sum(x))
  return(1 - G)
}

#Function for discounted catch
# Function for discounted catch with time component
discounted_catch <- function(catch, discount_rate, t){
  # discounted catch with time component
  disc_catch <- sum(catch/((1+discount_rate)^t))
  return(disc_catch)
}

# population growth of the predator and prey
deltaN <- function(Nt, Pt, pred_recruits, prey_params, pred_params){
    Nt1 <- Nt
    Pt1 <- Pt
    # add recruits, includes larval dispersal
    Nt1 <- Nt + vbh(Nt%*%prey_params$Dmat, prey_params$R0, prey_params$K) 
    # Pt1 <- Pt + vbh(Pt%*%prey_params$Dmat, pred_params$R0, pred_params$K)

    # predation
    delta_pred <- apply(rbind(Pt*Nt1*pred_params$a, Nt), 2, min)
    Nt1 <- Nt1 - delta_pred
    #one option -predator pop grows with more food, 
    # then then has feedback to recruitment next step
    # Pt1 <- Pt1 + delta_pred*pred_params$conversion
    #Second option for predation, becomes food dependent
    Pspawn <-  delta_pred * pred_params$conversion
    recruits_plus <- vbh(Pspawn %*% pred_params$Dmat, pred_params$R0, pred_params$K)
    #catch of new recruits
    Rcatch <- pred_recruits - pred_params$F
    #add recruits to predator population
    Pt1 <- Pt + pred_recruits

    # natural mortality
    Nt1 <- Nt1 - Nt1*prey_params$M
    Pt1 <- Pt1 - Pt1*pred_params$M
   
    # fish catches
    Ncatch <- Nt1*prey_params$Fmort
    Pcatch <- Pt1*pred_params$Fmort

    # fishing mortality
    Nt1 <- Nt1 - Ncatch
    Pt1 <- Pt1 - Pcatch

    return(rbind(Nt1, Pt1, recruits_plus, Ncatch, Pcatch))
}

# Function to simulate population dynamics

simulate <- function(N0, P0, R0, prey_params, pred_params, tmax, nloc){
  N <- matrix(NA, nrow = tmax, ncol = nloc)
  P <- matrix(NA, nrow = tmax, ncol = nloc)
  R <- matrix(NA, nrow = tmax, ncol = nloc)
  Ncatch <- matrix(0, nrow = tmax, ncol = nloc)
  Pcatch <- matrix(0, nrow = tmax, ncol = nloc)
  N[1,] <- N0
  P[1,] <- P0
  R[1:pred_params$tmat,] <- R0

 # Modify Fmort to allow for protected areas
 prey_params$Fmort <- rep(prey_params$Fmort, nloc)
 if (!is.null(prey_params$MPAs)){
    prey_params$Fmort[prey_params$MPA] <- 0
    #redistribution fishing effort
    fact <- length(prey_params$MPAs)*max(prey_params$Fmort)/(nloc - length(prey_params$MPAs))
    #add displaced effort to fished areas
    prey_params$Fmort[-prey_params$MPA] <- prey_params$Fmort[-prey_params$MPA] + 
      fact - prey_params$Fmort[-prey_params$MPA] *fact
  }

  pred_params$Fmort <- rep(pred_params$Fmort, nloc)
  if (!is.null(pred_params$MPAs)){
    pred_params$Fmort[pred_params$MPA] <- 0
    #redistribution fishing effort
    fact <- length(pred_params$MPAs)*max(pred_params$Fmort)/(nloc - length(pred_params$MPAs))
    #add displaced effort to fished areas
    pred_params$Fmort[-pred_params$MPA] <- pred_params$Fmort[-pred_params$MPA] + 
      fact - pred_params$Fmort[-pred_params$MPA] *fact
  }

  for (t in 2:tmax){
    if (t <= pred_params$tmat){
      Rtemp <- R0
    } else {
      Rtemp <- R[t-pred_params$tmat,]
    }
    NP <- deltaN(N[t-1,], P[t-1,], Rtemp, prey_params, pred_params)
    N[t,] <- NP[1,]
    P[t,] <- NP[2,]
    R[t,] <- NP[3,]
    Ncatch[t,] <- (NP[4,])
    Pcatch[t,] <- (NP[5,])
  }
  return(list(N=N, P=P, R= R, Ncatch = Ncatch, Pcatch = Pcatch))
}

wrangle_output <- function(sim, sim2){
  tmax <- nrow(sim$N)
  tmax2 <- nrow(sim2$N)
  #wrangle data to long format
sim_long <- sim$N %>% 
  rbind(sim2$N) %>%
  as.data.frame() %>% 
  mutate(t = 1:(tmax+tmax2)) %>% 
  pivot_longer(-t, names_to = "loc", values_to = "N") %>% 
  mutate(species = "prey") %>% 
  rename(abund = N) %>%
  bind_rows(
    sim$P %>% 
    rbind(sim2$P) %>%
      as.data.frame() %>% 
        mutate(t = 1:(tmax+tmax2)) %>% 
      pivot_longer(-t, names_to = "loc", values_to = "P") %>% 
      mutate(species = "predator")%>% 
  rename(abund = P) 
  ) %>%
  bind_rows(
    sim$R %>% 
    rbind(sim2$R) %>%
      as.data.frame() %>% 
        mutate(t = 1:(tmax+tmax2)) %>% 
      pivot_longer(-t, names_to = "loc", values_to = "R") %>% 
      mutate(species = "recruits")%>%
  rename(abund = R)
  )

#wrangel catches
#reset cathes in first year of sim2
sim2$Ncatch[1,] <- sim$Ncatch[tmax,]
sim2$Pcatch[1,] <- sim$Pcatch[tmax,]
catch_long <- sim$Ncatch %>% 
  rbind(sim2$Ncatch) %>%
  as.data.frame() %>% 
  mutate(t = 1:(tmax+tmax2)) %>% 
  pivot_longer(-t, names_to = "loc", values_to = "N") %>% 
  mutate(species = "prey") %>% 
  rename(catch = N) %>%
  bind_rows(
    sim$Pcatch %>% 
    rbind(sim2$Pcatch) %>%
      as.data.frame() %>% 
        mutate(t = 1:(tmax+tmax2)) %>% 
      pivot_longer(-t, names_to = "loc", values_to = "P") %>% 
      mutate(species = "predator")%>% 
  rename(catch = P) 
  ) 
  return(list(sim_long = sim_long, catch_long = catch_long))
}
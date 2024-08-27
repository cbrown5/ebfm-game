
#test the functions 
library(tidyverse)

source("functions.R")


theme_set(theme_classic())
nloc <- 3
nfishery <- 3
Dmat <- diag(1, nloc)
#Rows need to sum to 1
# Rows are from, columns are to, e.g. Dmat[1,2] is from loc 1 to loc 2
#predator has mostly retention
Dmatpred <- diag(nloc)
nonret <- 0.9 #non-retention
diag(Dmatpred) <- 1 - nonret
#set off diagonals to nonret/(nloc-1)
Dmatpred[!diag(nloc)] <- nonret/(nloc-1)

MPAvect <- c(2)
prey_params <- list(R0 = 1, #steepness of the stock-recruitment relationship
                    K = 200, #carrying capacity for recruitment
                    M = 0.1, #natural mortality
                    Fmort = 0.1, #fishing mortality
                    Dmat = Dmat, #dispersal matrix
                    MPAs = NULL) #protected areas
prey_params2 <- prey_params
prey_params2$MPAs <- MPAvect

pred_params <- list(R0 = 0.5, K = 10, M = 0.1,
 a = 0.005, Dmat = Dmatpred, 
  conversion = 0.1, tmat = 3, 
  recruit_weight = 0.5, #relative weight of recruits
  Fmort = 0.1, MPAs = NULL)

pred_params2 <- pred_params
pred_params2$MPAs <- MPAvect



#plot vbh
# N <- 1:100
# plot(N, vbh(N, 500, 50), type = "l", 
# col = "blue", xlab = "N", ylab = "Nt+1", 
# ylim = c(0, 100))
# lines(N, vbh(N, 100, 50), col = "red")


N0 <- rep(10, nloc)
P0 <- rep(10, nloc)
R0 <- rep(1, nloc)

tmax <- 100
sim <- simulate(N0, P0, R0, prey_params, pred_params, tmax, nloc)
sim2 <- simulate(sim$N[tmax,], sim$P[tmax,], sim$R[tmax,], prey_params2, 
pred_params2, tmax, nloc)

xout <- wrangle_output(sim, sim2)
sim_long <- xout$sim_long
catch_long <- xout$catch_long

#plot abundance over time, colours by location, facet by prey/predator
ggplot(sim_long, aes(t, abund, colour = loc)) +
  geom_line() +
  geom_vline(xintercept = tmax, linetype = "dashed") +
  facet_wrap(~species, scales = "free") +
  ylim(0, NA) +
  labs(x = "Time", y = "Abundance")

ggplot(catch_long, aes(t, catch, colour = loc)) +
  geom_line() +
    facet_wrap(~species, scales = "free") +
  labs(x = "Time", y = "Catches")


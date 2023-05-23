#############################
#Author: Veronica Malizia
#R version: 4.2.2
#
#Initial conditions and initial population to start simulations
#Species: Schistosoma mansoni
#############################

## Load initial population 
load("Equilibrium_age_distribution.RData") 
prob_death <- read.csv("prob_death_Uganda_2019.csv")
#Initial human cohort
cohort  <-  cohort %>%
  mutate(jw1 = 1, 
         Ind_sus = rgamma(nrow(cohort), shape = parms$parasite$k_w, scale = 1/parms$parasite$k_w))

## Initial conditions
init <- list(# Human cohort
  humans = list(cohort = cohort,
                #Initial cumulative exposure
                cum_exp = sum(Age_profile_exp(cohort$age)$y * cohort$Ind_sus),
                #Index of SAC in the cohort
                SAC = which(cohort$age >= 5 & cohort$age <= 15)),
  #Environment
  environment = list(eggs0 = 0,
                     contributions0 = 0),
  #Snail population 
  snails = list(snail.pop=10000,
                E0=0,
                I0=0,
                C0=0))
init$snails$S0 = init$snails$snail.pop - sum(init$snails$E0 + init$snails$I0) 

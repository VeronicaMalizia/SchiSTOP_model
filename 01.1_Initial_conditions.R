#############################
#Author: Veronica Malizia
#Date: 20/06/2022
#R version: 4.2.2
#
# The script sets the initial conditions to start simulations with SchiSTOP, about:
# - human population
# - environment
# - snail population
#
# Input Data: 
# - "Equilibrium_age_distribution.RData" generated from script 00_Demography.R
# - "prob_death_Uganda_2019.csv" (age-specific death probabilities)
# Species of interest: Schistosoma mansoni
#############################

## Load initial population 
load("Equilibrium_age_distribution.RData") 
prob_death <- read.csv("prob_death_Uganda_2019.csv")
#Initial human cohort
# cohort  <-  cohort %>%
#   mutate(jw1 = 1, 
#          Ind_sus = rgamma(nrow(cohort), shape = 0.15, scale = 1/0.15))

## Initial conditions
init <- list(# Human cohort
              humans = list(pop = cbind(cohort, #%>%
                                        jw1 = 1,
                                        ID = c(1:nrow(cohort)),
                                        death_age = 100,
                                        age_group = as.numeric(cut(cohort$age, c(-1, prob_death$Age_hi))),
                                        rate = 0,
                                        wp1 = 0,
                                        wp2 = 0,
                                        wp3 = 0,
                                        jw2 = 0,
                                        jw3 = 0,
                                        cum_dwp = 0,
                                        mu = 0,
                                        ec = 0),
                            #Initial cumulative exposure
                            #cum_exp = sum(Age_profile_exp(cohort$age)$y * cohort$Ind_sus),
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

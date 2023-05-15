#############################
#Author: Veronica Malizia
#R version: 4.1.2

#Calibration code for the three transmission parameters:
## transmission on humans (zeta)
## aggregation of worms
## transmission on snails

#Figures to match: 
## Moderate endemicity: 30% prevalence in SAC, 6% HI prevalence in SAC, 6% snail prevalence
## Low endemicity: 10% prevalence in SAC, ~ 0% HI prevalence in SAC, 6% snail prevalence
## High endemicity: 60% prevalence in SAC, 20% HI prevalence in SAC, 6% snail prevalence
#############################
rm(list = ls())

#Loading packages
#.libPaths(c("C:/Program Files/R/R-4.1.2/library",.libPaths()))
library(tidyverse)
library(readxl)
library(ggridges)
library(foreach)
library(doParallel)
#library(profvis)
library(rstudioapi)
library(beepr)

###########
#Priors
###########

# zeta (> 0)
# worms_aggr (in (0,1))
# tr_snails (> 0, very small)



#######################
# OUTPUT
#######################
#Select the three prevalence indicators at the end of simulation, averaged by seed

## Load initial population 
source.dir <- dirname(getActiveDocumentContext()$path)
setwd(source.dir)
load("Equilibrium_age_distribution.RData") 
prob_death <- read.csv("prob_death_Uganda_2019.csv")
#Load functions
source("01_Handy_functions.R")
##Load parameters
source("02_Parameters_Smansoni.R")

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

# Input (user choices)
#In this file you can set the modelling scenario you want to simulate
source("Setting_simulation_scenario.R")

#OUTPUT FUNCTION
#cohort is the initial population (comes from line 54)
#prob_death contains the death probabilities for the demography (comes from line 47)
#init sets initial conditions (comes from line 59)
#parms includes the list of parameters (comes from line 51)
#T is the time horizon (comes from line 77)
#stoch_scenarios sets the modelling scenario (comes from line 77)

output <- function(cohort, prob_death, init, parms, T, endem, stoch_scenarios){
  #Run model
  # Will give results as object 'results'
  source("04_Model_specification.R")
  
  #Check and remove (if any) indexes with errors in the results 
  #Errors mean that the snail system goes to zero (not good combination of parameters)
  index <- which(sapply(results, length)>2)
  results <- results[index]
  
  #Collating, averaging (by seeds), and saving population-level output
  res <- bind_rows(results) %>%
    group_by(zeta, worms_aggr, tr_snails, Immunity, Snails, DDF) %>%
    summarise(prevSAC = mean(eggs_prev_SAC), 
              Hprev = mean(Heggs_prev_SAC),
              snailprev = mean(inf_snail/(susc_snail+exp_snail+inf_snail)))
  
  return(res)
  # Save if needed
  # saveRDS(res, file = file.path(pop.output.dir, 
  #                               paste(setting, ".RDS", sep = "")))
  #then read with readRDS
}

#Check if works
output(cohort, prob_death, init, parms, T, endem, stoch_scenarios)

###########################
#Desired output of vector 3
###########################
## Moderate endemicity: 30% prevalence in SAC, 6% HI prevalence in SAC, 6% snail prevalence
## Low endemicity: 10% prevalence in SAC, ~ 0% HI prevalence in SAC, 6% snail prevalence
## High endemicity: 60% prevalence in SAC, 20% HI prevalence in SAC, 6% snail prevalence

high <- c(0.6, 0.2, 0.06)
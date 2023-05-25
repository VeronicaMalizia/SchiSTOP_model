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
.libPaths(c("C:/Program Files/R/R-4.1.2/library",.libPaths()))
library(tidyverse)
library(readxl)
library(rstudioapi)
library(beepr)
library(EasyABC)
library(deSolve)
library(splitstackshape)

# Working directory
source.dir <- dirname(getActiveDocumentContext()$path)
setwd(source.dir)

# 'use_seed' must be turned to TRUE.
#n=10

###########
#Priors
###########
#theta = c(zeta, worms_aggr, tr_snails)

theta_prior = list(c("normal", -5, 0.5), # zeta (> 0) #-7.5
                   c("unif", 0, 0.5)) # worms_aggr (in (0,1)) or beta(che può essere skewed)
                   #c("lognormal", -23, 1)) # tr_snails (> 0, very small) 
theta_prior

#######################
# OUTPUT
#######################
#Select the three prevalence indicators at the end of simulation

#Load functions
source("01_Handy_functions.R")
##Load parameters
source("02_Parameters_Smansoni.R")
#Load initial conditions
source("01.1_Initial_conditions.R")

# Input (user choices)
#Through this file the user can set the desired modelling scenario
source("Setting_simulation_scenario.R")
#set scenario-specific parameters
scen <- stoch_scenarios
#parms$parasite$zeta = scen$zeta #overall exposure rate.  
parms$immunity$imm = case_when(scen$imm_strength== "Absent" ~ 0,
                               scen$imm_strength== "Mild" ~ 0.0005,
                               scen$imm_strength== "Strong" ~ 0.005) #immunity slope parameter
parms$snails$carrying.capacity = case_when(scen$snails == "Absent" ~ 1, #No module
                                           scen$snails == "Mild" ~ 20000,
                                           scen$snails == "Strong" ~ 10000)
parms$parasite$eggs$alpha = case_when(scen$DDF_strength== "Absent" ~ alpha_lin,
                                      scen$DDF_strength== "Mild" ~ alpha_mild,
                                      scen$DDF_strength== "Strong" ~ alpha_strong) #fecundity parameter
parms$parasite$eggs$z = case_when(scen$DDF_strength== "Absent" ~ 0,
                                  scen$DDF_strength== "Mild" ~ 0.0005,
                                  scen$DDF_strength== "Strong" ~ 0.0007) #severity of density dependency
# parms$parasite$k_w = scen$worms_aggr
# parms$snails$snail_transmission_rate = scen$tr_snails


#OUTPUT (model) FUNCTION
schisto_model <- function(theta = c(-5, 0.25, 10^(-10))){
  
  parms$parasite$zeta = exp(theta[1])
  parms$parasite$k_w = theta[2]
  parms$snails$snail_transmission_rate = theta[3]

  #Run model
  # Will give results as object 'results'
  source("04_Model_specification_simpleABC.R")
  
  return(res[1:2])
  # Save if needed
  # saveRDS(res, file = file.path(pop.output.dir, 
  #                               paste(setting, ".RDS", sep = "")))
  #then read with readRDS
}

#Check if it works
schisto_model()

###########################
#Desired output of vector 3
###########################
## Moderate endemicity: 30% prevalence in SAC, 6% HI prevalence in SAC, 6% snail prevalence
## Low endemicity: 10% prevalence in SAC, ~ 0% HI prevalence in SAC, 6% snail prevalence
## High endemicity: 60% prevalence in SAC, 20% HI prevalence in SAC, 6% snail prevalence

high_stat_obs = c(0.6, 0.2) #c(0.6, 0.2, 0.06)

###### Specific parameters #######
## to perform the Beaumont et al. (2009)'s method:
tolerance=3 #c(10, 0.5)


ABC_Delmoral<-ABC_sequential(method="Delmoral", model=schisto_model, prior=theta_prior, 
                             #n_cluster = parallel::detectCores(logical = FALSE),
                             #prior_test = "X3 < X1", #inside_prior = FALSE,
                             nb_simul=10, use_seed=TRUE, verbose = FALSE, progress_bar = TRUE,
                             summary_stat_target=high_stat_obs, tolerance_target=tolerance)
ABC_Beaumont



method="Beaumont"
model=schisto_model 
prior=theta_prior #inside_prior = FALSE,
nb_simul=10 
use_seed=TRUE 
verbose = TRUE 
progress_bar = TRUE
summary_stat_target=high_stat_obs 
tolerance_tab=tolerance

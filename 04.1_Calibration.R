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

theta_prior = list(c("lognormal", -2, 0.5), # zeta (> 0) #-7.5
                   c("unif", 0, 0.5), # worms_aggr (in (0,1)) or beta(che può essere skewed)
                   c("lognormal", -23, 1)) # tr_snails (> 0, very small) 
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

#OUTPUT (model) FUNCTION
schisto_model <- function(theta = c(0.15, 0.15, 10^(-10))){
  
  parms$parasite$zeta = theta[1]
  parms$parasite$k_w = theta[2]
  parms$snails$snail_transmission_rate = theta[3]

  #Run model
  # Will give results as object 'results'
  source("04_Model_specification_simpleABC.R")
  
  #Check and remove (if any) indexes with errors in the results 
  #Errors mean that the snail system goes to zero (not good combination of parameters)

  return(res)
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

high_stat_obs = c(0.6, 0.2, 0.06)

###### Specific parameters #######
## to perform the Beaumont et al. (2009)'s method:
tolerance=c(10, 1, 0.5)


ABC_Beaumont<-ABC_sequential(method="Beaumont", model=schisto_model, prior=theta_prior, #inside_prior = FALSE,
                             nb_simul=10, use_seed=TRUE, verbose = TRUE, progress_bar = TRUE,
                             summary_stat_target=high_stat_obs, tolerance_tab=tolerance)
ABC_Beaumont


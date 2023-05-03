#############################
#Author: Veronica Malizia
#R version: 4.1.2

#This script runs the model from R source, with human demography (aging/birth/deaths)
#and transmission dynamics including an infection reservoir where snail population dynamcis are explicitely modelled.

#Data required: 
# 1. Initial age distribution, after equilibrium is reached in 00_Demography_eq.R
#    The file is called "Equilibrium_age_distribution.RData"
# 2. Death_probabilities from WHO App - South Africa
#############################
rm(list = ls())

#Loading packages
.libPaths(c("C:/Program Files/R/R-4.1.2/library",.libPaths()))
library(tidyverse)
library(readxl)
library(ggridges)
library(foreach)
library(doParallel)
library(patchwork)
library(profvis)
library(rstudioapi)
library(beepr)

################
#Initial population
#Building up the initial cohort
################
source.dir <- dirname(getActiveDocumentContext()$path)
  #"C:/Users/Z541213/Documents/Project/Model/Schisto_model"
setwd(source.dir)
#Load age distribution at equilibrium and death rates
load("Equilibrium_age_distribution.RData") #the object is called "to.save"
#death_rates <- read.csv("death_rates_Uganda_2019.csv")
prob_death <- read.csv("prob_death_Uganda_2019.csv")

#0=male, 1=female
cohort <- c()
for(i in 1:nrow(to.save)){
  tmp <- tibble(age = rep(to.save$age[i], round(to.save$n[i])))
  cohort <- rbind(cohort, tmp) 
}
cohort <- cohort %>%
  mutate(sex = as.numeric(rbernoulli(nrow(cohort), 0.5)))
#Checks
hist(cohort$age, main = "Initial age distribution", xlab = "Age (ys)")
table(cohort$sex)

################
#Loadings
################

#Load functions
source("01_Handy_functions.R")
#Water contact data (duration)
age_groups <- c(0, 10, 20, 150)
exposure_rates <- c(0.62, 1, 0.51, 0.51) #Relative Age-specific exposure (minutes/person)
#(frequency)
exposure_rates <- c(0.75, 1, 0.50, 0.50) #Relative Age-specific exposure (activity/person)

#Checks age-exposure and contribution
plot(approxfun(x=age_groups, y=exposure_rates, method = "constant"), 
     xlim = c(0, 100), ylim = c(0, 1), 
     type = 's', xlab = "Age", ylab = "Relative exposure")
abline(v = c(5, 15), col = 'red')
#lines(approx(x=c(0, 5, 10, 16, 100), y=c(0.01, 0.61, 1, 0.12, 0), method = "linear"), type = 'l', col='red')
lines(approx(x=c(0, 10, 200), y=c(1, 1, 1), method = "linear"), col = "brown")

##Load parameters
source("02_Parameters_Smansoni.R")

################
#Initializing
################

#Worms are initialized in the cohort with an artificial FOI (1 worm pair per person per month)
cohort <- cohort %>%
  mutate(jw1 = 1, 
         Ind_sus = rgamma(nrow(cohort), shape = parms$parasite$k_w, scale = 1/parms$parasite$k_w))


#Environment is initialized as empty
eggs0 <- 0
contributions0 <- 0

#Initial cumulative exposure
cum_exp <- sum(Age_profile_exp(cohort$age)$y * cohort$Ind_sus) 
#Index of SAC in the cohort
SAC <- which(cohort$age >= 5 & cohort$age <= 15)

#Snail population is initialized
snail.pop=10000
E0=0
I0=0
S0=snail.pop - sum(E0, I0)
C0=0

################
#Simulation settings
################
source("Setting_simulation_scenario.R")

################
#Set output directory to save results
################
#This will be the directory where the individual output is automatically saved throughout the simulations
if(write.output == TRUE){
  ind.output.dir <- file.path(source.dir, paste("Output/Individual/", setting, sep = "")) 
  if(!file.exists(ind.output.dir)){
    dir.create(ind.output.dir) #Add check: this command to be run only if the directory is not existent
  }
  if(file.exists(ind.output.dir)){
  #Empty the Output folder (only if needed)
  unlink(file.path(ind.output.dir, "/*")) 
  }
}

#parms$snails$snail_transmission_rate = parms$snails$snail_transmission_rate/100

################
#Run the model
################
time.start <- Sys.time()
source("04_Model_specification.R")
time.end <- Sys.time()
time.end - time.start
beep()

################
#Collating and saving population-level results 
################
#Population-level results
#Set endemicity setting:
pop.output.dir <- file.path(source.dir, "Output/Population/")
if(!file.exists(pop.output.dir)){
  dir.create(pop.output.dir)
}
#Collating and saving population-level output
#Individual output is automatically saved through the simulations
res <- bind_rows(results)
saveRDS(res, file = file.path(pop.output.dir, 
                           paste(setting, ".RDS", sep = "")))

################
#Running plotting code
################
##Checks when running single scenario
# source("06_Plotting_single_scenario.R")
# 
# ################
# #Inspecting and saving plots (if desired)
# ################
# 
# #Prevalence timelines
# tiff("Plots/Modeling scenarios/IMM_mild/Prevalence_HIGHsetting.tif", width=7, height=6, units = "in", res = 300)
# Fig1
# dev.off()
# 
# #Prevalence of infected snails
# tiff("Plots/Prevalence_snail_low.tif", width=7, height=6, units = "in", res = 300)
# Fig2
# dev.off()
# 
# Fig3
# Fig4
# cerc
# mirac
# cercVSmirac

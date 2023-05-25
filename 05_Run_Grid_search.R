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
#.libPaths(c("C:/Program Files/R/R-4.1.2/library",.libPaths()))
library(tidyverse)
library(readxl)
library(ggridges)
library(foreach)
library(doParallel)
library(patchwork)
library(profvis)
library(rstudioapi)
library(deSolve)

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
#Checks
hist(cohort$age, main = "Initial age distribution", xlab = "Age (ys)")
table(cohort$sex)

################
#Loadings
################

#Load functions
source("01_Handy_functions.R")
age_groups <- c(0, 10, 20, 150)
exposure_rates <- c(0.62, 1, 0.51, 0.51) #Relative age-specific exposure (minutes/person)
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
#Preparing simulation setting
################
source("Setting_simulation_scenario.R")

################
#Run the model
################
setwd(source.dir)

time.start <- Sys.time()
source("04_Model_specification.R")
time.end <- Sys.time()
time.end - time.start

################
#Collating and saving population-level results 
################
#Population-level results
#Set endemicity setting:
pop.output.dir <- file.path(source.dir, "Grid search")
if(!file.exists(pop.output.dir)){
  dir.create(pop.output.dir)
}

#Check indexes without errors (for ex. when the snail system goes to zero)
index <- which(sapply(results, length)>2)
results <- results[index]

#Collating and saving population-level output
#Individual output is automatically saved through the simulations
res <- bind_rows(results)
saveRDS(res, file = file.path(pop.output.dir, 
                           paste(setting, ".RDS", sep = "")))
#then read with readRDS



# ################
# Inspection
# ################
sims <- res %>%
  group_by(zeta, worms_aggr, tr_snails, Immunity, Snails, DDF) %>%
  summarise(prevSAC = mean(eggs_prev_SAC), 
            Hprev = mean(Heggs_prev),
            snailprev = mean(inf_snail/(susc_snail+exp_snail+inf_snail))) %>%
  mutate(RMSE = sqrt((sum(0.6-prevSAC, 0.2-Hprev, 0.06-snailprev)^2)/3))

sims <- sims %>%
  group_by(Immunity, Snails, DDF) %>%
  slice_min(order_by = RMSE)

ggplot(sims, aes(x=zeta, y=worms_aggr))+
  geom_point(aes(color = RMSE)) + 
  facet_grid(~ DDF)

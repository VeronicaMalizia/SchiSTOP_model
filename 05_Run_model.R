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

################
#Initial population
#Building up the initial cohort
################
source.dir <- "C:/Users/Z541213/Documents/Project/Model/Schisto_model"
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
#Checks age-exposure and contribution
plot(approxfun(x=age_groups, y=exposure_rates, method = "constant"), xlim = c(0, 200), ylim = c(0, 1), 
     type = 'l', xlab = "Age", ylab = "Relative exposure rate")
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
T <- 300 #number of years simulated
seeds <- 1
fr <- 10 #frequency for printing to file the individual output [years]

################
#SETTING THE MODELLING SCENARIO: limiting mechanism
################
#For each scenario, choose level: 'No', 'Mild', 'Strong'
#Anti-reinfection immunity
imm_strength <- "Strong"
parms$immunity$imm = case_when(imm_strength== "Absent" ~ 0,
                               imm_strength== "Mild" ~ 0.0005,
                               imm_strength== "Strong" ~ 0.005) #immunity slope parameter

#Snails 
snails = "Strong" #Choose "yes" or "no"
#If snails is on choose k=20000 (more snails - less limitation), k=10000 : strong (less snails - more limitation)
parms$snails$carrying.capacity = case_when(snails == "Absent" ~ 1, #No module
                                           snails == "Mild" ~ 20000,
                                           snails == "Strong" ~ 10000)
#Density-dependent fecundity (DDF)
load("Matched_alphas.RData")
DDF_strength <- "Absent"
parms$parasite$eggs$alpha = case_when(DDF_strength== "Absent" ~ alpha_lin,
                                      DDF_strength== "Mild" ~ alpha_exp,
                                      DDF_strength== "Strong" ~ alpha_strong) #fecundity parameter
parms$parasite$eggs$z = case_when(DDF_strength== "Absent" ~ 0,
                                  DDF_strength== "Mild" ~ 0.0007,
                                  DDF_strength== "Strong" ~ 0.004) #severity of density dependency

#Transmission parameters for tuning endemicity
parms$parasite$zeta = 0.004 #overall exposure rate. 

################
#Set output directory to save results
################
output.dir <- file.path(source.dir, paste("Output")) #/Individual
#dir.create(output.dir) #Add check: this command to be run only if the directory is not existent
#Empty the Output folder (only if needed)
#unlink(file.path(output.dir, "/*")) 

################
#Run the model
################
time.start <- Sys.time()
source("C:\\Users\\Z541213\\Documents\\Project\\Model\\Schisto_model\\04_Model_specification.R")
time.end <- Sys.time()
time.end - time.start

################
#Collating and saving population-level results 
################
#Population-level results
res <- bind_rows(results)
save(res, file = file.path(source.dir, 
                           paste("/Output/Population/Imm", 
                                 imm_strength, "Sn=", snails, "DDF=", DDF_strength, ".RData", sep="")))

################
#Running plotting code
################
source("06_Plotting_single_scenario.R")

################
#Inspecting and saving plots (if desired)
################

#Prevalence timelines
tiff("Plots/Modeling scenarios/IMM_mild/Prevalence_HIGHsetting.tif", width=7, height=6, units = "in", res = 300)
Fig1
dev.off()

#Prevalence of infected snails
tiff("Plots/Prevalence_snail_low.tif", width=7, height=6, units = "in", res = 300)
Fig2
dev.off()

Fig3
Fig4
cerc
mirac
cercVSmirac

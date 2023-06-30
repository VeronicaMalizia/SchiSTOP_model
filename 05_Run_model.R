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

################
#Loadings
################

#Load initial population and conditions
source("01.1_Initial_conditions.R")

#0=male, 1=female
#Checks
hist(cohort$age, breaks = c(0, prob_death$Age_hi[-c(1, nrow(prob_death))]+1),
     main = "Initial age distribution", xlab = "Age (ys)")
table(cohort$sex)

#Load functions
source("01_Handy_functions.R")
#If you want to use Water contact data (duration):
# age_groups <- c(0, 10, 20, 150)
# exposure_rates <- c(0.62, 1, 0.51, 0.51) #Relative Age-specific exposure (minutes/person)
#(frequency)
#exposure_rates <- c(0.75, 1, 0.50, 0.50) #Relative Age-specific exposure (activity/person)

##Load parameters
source("02_Parameters_Smansoni.R")

#Check age-exposure(s) and contribution
# plot(Age_profile_exp(parms$exposure$ICL_derived$ages, 
#                      parms$exposure$ICL_derived$exp, x=0:80, 
#                      method = parms$exposure$ICL_derived$method), 
#      xlim = c(0, 100), ylim = c(0, 1), 
#      type = 's', xlab = "Age", ylab = "Relative exposure")
# lines(Age_profile_exp(parms$exposure$Sow_derived$ages, 
#                       parms$exposure$Sow_derived$exp, x=c(0,80), 
#                       method = parms$exposure$Sow_derived$method), type = 'l',col = "dark green")
# abline(v = c(5, 15), col = 'red')
# #lines(approx(x=c(0, 5, 10, 16, 100), y=c(0.01, 0.61, 1, 0.12, 0), method = "linear"), type = 'l', col='red')
# lines(approx(x=c(0, 10, 200), y=c(1, 1, 1), method = "linear"), col = "brown")

################
#Initializing
################

#Worms are initialized in the cohort with an artificial FOI (1 worm pair per person per month)
init$humans$pop <- init$humans$pop %>%
    mutate(Ind_sus = rgamma(nrow(cohort), shape = parms$parasite$k_w, scale = 1/parms$parasite$k_w),
           complier = as.numeric(rbernoulli(nrow(cohort), 1-parms$mda$fr_excluded)))

################
#Preparing simulation settings
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
#res <- bind_rows(results)
saveRDS(bind_rows(results), file = file.path(pop.output.dir, 
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

#cHECK DEMOGRAPHY
data_all <- list.files(path = ind.output.dir,  # Identify all output CSV files
                       pattern = "^Ind_out_seed", full.names = TRUE) %>% 
  lapply(readRDS) %>%              # Store all files in list
  #lapply(read_csv, show_col_types = F) %>%
  bind_rows %>%
  filter(time == 191)
age_out <- data_all %>%
  group_by(seed, time, age_group, Immunity, Snails, DDF) %>% #sex not for now
  summarise(pop = n())
#Average over seeds
age_out <- age_out %>%
  group_by(time, age_group, Immunity, Snails, DDF) %>% #sex not for now
  summarise(pop = mean(pop))

ggplot(age_out, aes(x = as.factor(age_group), y = pop)) +
  geom_col(aes(colour=DDF, fill=DDF), alpha = 0.7, position = "dodge") +
  facet_grid(Snails ~ Immunity, labeller = labeller(.rows = label_both, .cols = label_both))
            
            
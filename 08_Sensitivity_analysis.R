###############################
# Author: Veronica Malizia
# Date: 29/06/2024
# R version 4.4.1
#
# This script loads simulation outputs and computes predictions for the sensitivity analysis on the MDA drug efficacy
# This parameter is defined in the model as portion of adult worm pairs killed by each MDA round
# Drug efficacy in the main text: 0.86
# Values here explored: 0.8, 0.9
#
# Loaded scripts:
# - 01_Handy_functions.R, list of functions used to define the transmission model
# - 02_Parameters_Smansoni.R, a full list of parameters employed for the simulations
#
# Output:
# - .docx file table with predicted probabilities
################################

rm(list = ls())

library(tidyverse)
library(readxl)
library(ggridges)
library(egg)
library(rstudioapi)
library(scales)
library(tagger)
#to install the package above run: devtools::install_github("eliocamp/tagger")
library(rempsyc)

#Set folder
source.dir <- dirname(getActiveDocumentContext()$path)
setwd(source.dir)

pop.output.dir <- file.path(source.dir, "Output/Population/Rebuttal2")

#Load functions
source("01_Handy_functions.R")

#Load parameters
source("02_Parameters_Smansoni.R")
seeds <- 100
# Adjust accordingly
parms$mda <- list(age.lo = 2, #SAC is 5-15 #all population >= 2ys (WHO)
                  age.hi = 100,
                  start = 150,
                  end = 159,
                  frequency = 1, #annual
                  coverage = 0.75,
                  fr_excluded = 0.05, #systematic non-compliance 
                  efficacy = 0.8)

#Load data of interest: predictions are computed only on those modelling scenarios found successful in reproducing epidemiological patterns for schistosomiasis (see manuscript)
# Successful scenarios: 
# # - Exposure function: based on water contacts ("Sow")
# # - Human-level regulation (Immunity): mild or strong
# # - Snail-level regulation (Snails): mild or strong
# # - Worm-level regulation (DDF): any

res_Sow <- readRDS(file.path(pop.output.dir, "Sow_func_ALLmda_0.9killworm.RDS")) %>%
  filter(Immunity != "Absent" & Snails != "Absent") %>% #Successful scenarios only 
  mutate(Exposure = "Based on water contacts")
res_Sow$Endemicity <- factor(as.factor(res_Sow$Endemicity),
                             levels = c("Low", "Moderate", "High"))

# Reference scenario: 
# # - Exposure function: model-derived ("ICL")
# # - Human-level regulation (Immunity): absent
# # - Snail-level regulation (Snails): absent
# # - Worm-level regulation (DDF): strong
# This is not among successful scenarios, but it is added to computations as reference for the commonly employed assumptions

reference <- readRDS(file.path(pop.output.dir, "ICL_func_ALLmda_0.9killworm.RDS")) %>%
  filter(Immunity == "Absent" & Snails == "Absent" & DDF == "Strong") %>% #Reference scenario for assumptions commonly in use 
  mutate(Exposure = "Model-based")
reference$Endemicity <- factor(as.factor(reference$Endemicity),
                               levels = c("Low", "Moderate", "High"))

#Computing runs that reach EPHP or interruption of transmission:

#Writing the table
options(digits = 1)
pred <- res_Sow %>%
  filter(time==(parms$mda$end+20)*12) %>%
  group_by(Endemicity, Immunity, Snails, DDF, Exposure) %>%
  summarise(ephp_seed = length(which(Heggs_prev<=0.01))*100/seeds,
            eliminated = length(which(eggs_prev_SAC==0))*100/seeds)

ref <- reference %>%
  filter(time==(parms$mda$end+20)*12) %>%
  group_by(Endemicity, Immunity, Snails, DDF, Exposure) %>%
  summarise(ephp_seed = length(which(Heggs_prev<=0.01))*100/seeds,
            eliminated = length(which(eggs_prev_SAC==0))*100/seeds)

#library(flextable)

table <- nice_table(pred, options(digits = 1))
table_ref <- nice_table(ref, options(digits = 1))

print(table, preview ="docx")
print(table_ref, preview ="docx")

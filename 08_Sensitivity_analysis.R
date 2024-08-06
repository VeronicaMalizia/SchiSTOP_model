###############################
# Author: Veronica Malizia
#
# This script load data and computes predictions for the sensitivity analysis on the worm killing rate
# Killing rate in main text: 0.86
# Values explored: 0.8, 0.9
#
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

#Load data of interest
res_Sow <- readRDS(file.path(pop.output.dir, "Sow_func_ALLmda_0.9killworm.RDS")) %>%
  filter(Immunity != "Absent" & Snails != "Absent") %>% #Successful scenarios only 
  mutate(Exposure = "Based on water contacts")
res_Sow$Endemicity <- factor(as.factor(res_Sow$Endemicity),
                             levels = c("Low", "Moderate", "High"))

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

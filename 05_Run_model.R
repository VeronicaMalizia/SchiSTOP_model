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
#library(profvis)
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

##Load parameters
source("02_Parameters_Smansoni.R")

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
res <- bind_rows(results)
saveRDS(bind_rows(results), file = file.path(pop.output.dir, 
                           paste(setting, ".RDS", sep = "")))


################
#Check equilibria
################
faided <- res %>%
  filter(time==(parms$mda$start-1)*12) %>%
  group_by(Immunity, Snails, DDF) %>%
  summarise(faided_seed = length(which(eggs_prev_SAC==0))*100/seeds)
n_faided <- length(which(faided$faided_seed>0))
#Average by seed
if(n_faided>0){
  res <- res[- which(res$time==(parms$mda$start-1)*12 & res$eggs_prev_SAC==0), ]
}

data_avg <- res %>%
  group_by(time, Immunity, Snails, DDF, Endemicity) %>%
  dplyr::summarise(eggs_prev_SAC = mean(eggs_prev_SAC),
            eggs_prev_tot = mean(eggs_prev),
            PHI = mean(Heggs_prev_SAC),
            snail_inf = mean(inf_snail),
            snail_exp = mean(exp_snail),
            snail_prev = mean(inf_snail/(susc_snail+inf_snail+exp_snail)))

#filter(data_avg, Endemicity == "Low") %>%
data_avg %>% 
ggplot(aes(x=time/12, group=interaction(DDF, Endemicity))) +
  # geom_line(data = filter(res, Endemicity == "Low"),
  #           aes(y=eggs_prev_SAC*100,
  #               group = interaction(DDF, seed),
  #               colour = DDF), alpha = 0.01) +
  geom_line(aes(y=eggs_prev_SAC*100, colour = DDF), alpha = 3) +
  geom_hline(yintercept = 30) +
  geom_line(aes(y=PHI*100, colour = DDF), linetype = "longdash") +
  geom_line(aes(y=snail_prev*100, colour = DDF), linetype = "dotted") +
  facet_grid(Snails ~ Immunity, labeller = labeller(.rows = label_both, .cols = label_both)) +
  scale_y_continuous(name = "Prevalence of infection in SAC (%) \n",
                     breaks = seq(0, 100, 10),
                     limits = c(0, 40),
                     expand = c(0, 0)) +
  scale_x_continuous(name = "\n Time [Years]",
                     expand = c(0, 0)) +
  # scale_x_continuous(name = "\n Years since last treatment round",
  #                    breaks = seq(parms$mda$end-10, parms$mda$end+20, 10),
  #                    labels = seq(-10, 20, 10),
  #                    #limits = c(0, 1200),
  #                    expand = c(0, 0)) +
  #coord_cartesian(xlim=c(parms$mda$end-10, parms$mda$end+20)) +
  expand_limits(x = 0,y = 0) +
  theme_bw() +
  theme(legend.position="bottom",
        #axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5),
        plot.margin = margin(5, 10, 0, 10, "pt"),
        panel.spacing = unit(1.5, "lines"))


################
#Running plotting code
################
##Checks when running single scenario
# source("06_Plotting_single_scenario.R")
# 
# data_avg <- res %>%
#   group_by(time, Immunity, Snails, DDF, Endemicity) %>%
#   dplyr::summarise(eggs_prev_SAC = mean(eggs_prev_SAC),
#             eggs_prev_tot = mean(eggs_prev),
#             PHI = mean(Heggs_prev_SAC),
#             snail_inf = mean(inf_snail),
#             snail_exp = mean(exp_snail),
#             snail_prev = mean(inf_snail/(susc_snail+inf_snail+exp_snail)))

# ggplot(data=data_avg, aes(x=time/12, group=interaction(DDF, Endemicity))) +
#   # geom_line(data = filter(res, Endemicity == "Low"), 
#   #           aes(y=eggs_prev_SAC*100,
#   #               group = interaction(DDF, seed),
#   #               colour = DDF), alpha = 0.01) +
#   geom_line(aes(y=eggs_prev_SAC*100, colour = DDF), alpha = 3) +
#   geom_hline(yintercept = 10) +
#   geom_line(aes(y=PHI*100, colour = DDF), linetype = "longdash") +
#   #geom_line(aes(y=snail_prev*100, colour = DDF), linetype = "dotted") +
#   facet_grid(Snails ~ Immunity, labeller = labeller(.rows = label_both, .cols = label_both)) +
#   scale_y_continuous(name = "Prevalence of infection in SAC (%) \n",
#                      breaks = seq(0, 100, 10),
#                      #limits = c(0, 20),
#                      expand = c(0, 0)) +
#   scale_x_continuous(name = "\n Time [Years]",
#                      expand = c(0, 0)) +
#   # scale_x_continuous(name = "\n Years since last treatment round",
#   #                    breaks = seq(parms$mda$end-10, parms$mda$end+20, 10),
#   #                    labels = seq(-10, 20, 10),
#   #                    #limits = c(0, 1200),
#   #                    expand = c(0, 0)) +
#   #coord_cartesian(xlim=c(parms$mda$end-10, parms$mda$end+20)) +
#   expand_limits(x = 0,y = 0) +
#   theme_bw() +
#   theme(legend.position="bottom",
#         #axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5),
#         plot.margin = margin(5, 10, 0, 10, "pt"),
#         panel.spacing = unit(1.5, "lines"))
# 


#Check if demography is at equilibrium
# data_all <- list.files(path = ind.output.dir,  # Identify all output CSV files
#                        pattern = "^Ind_out_seed", full.names = TRUE) %>% 
#   lapply(readRDS) %>%              # Store all files in list
#   #lapply(read_csv, show_col_types = F) %>%
#   bind_rows %>%
#   filter(time == 191)
# age_out <- data_all %>%
#   group_by(seed, time, age_group, Immunity, Snails, DDF) %>% #sex not for now
#   summarise(pop = n())
# #Average over seeds
# age_out <- age_out %>%
#   group_by(time, age_group, Immunity, Snails, DDF) %>% #sex not for now
#   summarise(pop = mean(pop))
# 
# ggplot(age_out, aes(x = as.factor(age_group), y = pop)) +
#   geom_col(aes(colour=DDF, fill=DDF), alpha = 0.7, position = "dodge") +
#   facet_grid(Snails ~ Immunity, labeller = labeller(.rows = label_both, .cols = label_both))
#             
            
# This script prepares data for plotting with 07_PLOTTING_CODE.R
# and computes predictions to reach control targets, for the Results section (Fig 4 and 5)
# It also contains (hereby commented) checks performed along the analysis 

rm(list = ls())

.libPaths(c("C:/Program Files/R/R-4.1.2/library",.libPaths()))
library(tidyverse)
library(readxl)
library(ggridges)
library(patchwork)
library(egg)
library(rstudioapi)
library(scales)
#library(plotly)
library(tagger)
#to install the package above run: devtools::install_github("eliocamp/tagger")
library(rempsyc)

#Set folder
source.dir <- dirname(getActiveDocumentContext()$path)
setwd(source.dir)


#Load functions
source("01_Handy_functions.R")

#Load parameters
source("02_Parameters_Smansoni.R")
source("Setting_simulation_scenario.R") #This is for parameters of MDA and n. of seeds used
######### Setting #########

#Behavior in exposure
# exposure = "Sow" #Choices: "ICL" (model-derived), "Sow" (water contacts)
# setting <- paste(exposure, "func", sep = "")

#######################
#Load population output
#######################
#Load collated results and produce multi-panel plots
pop.output.dir <- file.path(source.dir, "Output/Population")

setting <- paste("ICL", "func_SACmda", sep = "_")
res_ICL <- readRDS(file.path(pop.output.dir, paste(setting, ".RDS", sep = ""))) %>%
  #filter(time/12 > (parms$mda$end - 20) & time/12 < (parms$mda$end + 20)) %>%
  mutate(Exposure = "Model-based function")

setting <- paste("Sow", "func_SACmda", sep = "_")
res_Sow <- readRDS(file.path(pop.output.dir, paste(setting, ".RDS", sep = ""))) %>%
  #filter(time/12 > (parms$mda$end - 20) & time/12 < (parms$mda$end + 20)) %>%
  mutate(Exposure = "Based on water contacts")

res <- bind_rows(res_Sow, res_ICL)
res$Exposure <- factor(as.factor(res$Exposure),
                            levels = c("Model-based function", "Based on water contacts"))
res$Endemicity <- factor(as.factor(res$Endemicity),
                              levels = c("Low", "Moderate", "High"))
rm(res_ICL, res_Sow)
gc()

# Averaging population data
data_avg <- res %>%
  group_by(time, Immunity, Snails, DDF, Endemicity, Exposure) %>%
  summarise(eggs_prev_SAC = mean(eggs_prev_SAC),
            eggs_prev_tot = mean(eggs_prev),
            PHI = mean(Heggs_prev),
            miracidiae = mean(miracidiae),
            cercariae = mean(cercariae),
            snail_inf = mean(inf_snail),
            snail_exp = mean(exp_snail),
            snail_prev = mean(inf_snail/(susc_snail+inf_snail+exp_snail))) 


save(res, data_avg, file = "Population data for Figure2&3_10062024.RData")

#Computing faded-out runs at pre-control (149 ys is soon before MDA)
faided <- res %>%
  filter(time==(parms$mda$start-1)*12) %>%
  group_by(Immunity, Snails, DDF) %>%
  summarise(faided_seed = length(which(eggs_prev_SAC==0))*100/n.un(res$seed))
n_faided <- length(which(faided$faided_seed>0))

#Average by seed
if(n_faided>0){
  res <- res[- which(res$time==149*12 & res$eggs_prev_SAC==0), ]
}

######################
#Load individual data
#####################
#Age-intensity profiles
#"Sow and "ICL" datasets should be loaded separately, changing "exposure" below
exposure = "Sow"
setting <- paste(exposure, "func_SACmda", sep = "_")
ind.output.dir <- file.path(source.dir, paste("Output/Individual/", setting, sep = "")) 
data_all <- list.files(path = ind.output.dir,  # Identify all output CSV files
                       pattern = "^Ind_out_seed", full.names = TRUE) %>% 
  lapply(readRDS) %>%              # Store all files in list
  #lapply(read_csv, show_col_types = F) %>%
  bind_rows %>%
  dplyr::filter(time == parms$mda$start-9) #, # | time == parms$mda$end+1)
                #Snails == "Mild" & DDF == "Mild") 

#Averaging individual data
age_out <- data_all %>% 
  group_by(time, Immunity, Snails, DDF, Endemicity) %>%
  mutate(age_group = cut_number(age, n = 14)) %>%
  group_by(seed, time, age_group, Immunity, Snails, DDF, Endemicity) %>% #sex not for now
  summarise(avg_age_group = mean(age),
            epg = mean(ec), #means over individuals of that age group
            geom_epg = (geom_mean(ec+1)), 
            wp = mean(tot_wp),
            dwp = mean(cum_dwp), 
            rate = mean(rate))

#Average over seeds
data_toplot_Sow <- age_out %>%
  group_by(age_group, time, Immunity, Snails, DDF, Endemicity) %>%
  summarise(avg_age_group = mean(avg_age_group),
            epg_mean = mean(epg),
            epg_lo = ci(epg)[1],
            epg_hi = ci(epg)[2],
            geom_epg_mean = mean(geom_epg),
            geom_epg_lo = ci(geom_epg)[1],
            geom_epg_hi = ci(geom_epg)[2],
            wp_mean = mean(wp),
            wp_lo = ci(wp)[1],
            wp_hi = ci(wp)[2],
            dwp_mean = mean(dwp),
            dwp_lo = ci(dwp)[1],
            dwp_hi = ci(dwp)[2],
            rate_mean = mean(rate),
            rate_lo = ci(rate)[1],
            rate_hi = ci(rate)[2]) %>%
  mutate(Exposure = exposure)

########### ERR check: egg reduction rate after 1 MDA round
epg_pre_Sow <- data_toplot_Sow %>%
  ungroup() %>%
  filter(avg_age_group >= 5 & avg_age_group <= 15) %>%
  group_by(Immunity, Snails, DDF, Endemicity) %>%
  summarise(epg_mean = mean(epg_mean),
            geom_epg_mean = mean(geom_epg_mean))
#Now to have data after 1 round of MDA
data_all <- list.files(path = ind.output.dir,  # Identify all output CSV files
                       pattern = "^Ind_out_seed", full.names = TRUE) %>% 
  lapply(readRDS) %>%              # Store all files in list
  bind_rows %>%
  dplyr::filter(time == parms$mda$start+1,
                age >= 5 & age <= 15) %>%
  dplyr::select(seed, time, Immunity, Snails, DDF, Endemicity, ec) #and repeat averaging above

epg_post_Sow <- data_all %>%
  #filter(avg_age_group >= 5 & avg_age_group <= 15) %>%
  group_by(time, Immunity, Snails, DDF, Endemicity) %>%
  summarise(epg_mean = mean(ec),
            geom_epg_mean = (geom_mean(ec+1)))

ERR_Sow <- (epg_pre_Sow$epg_mean - epg_post_Sow$epg_mean)/epg_pre_Sow$epg_mean
quantile(ERR_Sow, probs = c(0.05, 0.5, 0.95), na.rm = T)

epg_pre_ICL <- data_toplot_ind %>%
  ungroup() %>%
  filter(Exposure == "ICL") %>%
  filter(avg_age_group >= 5 & avg_age_group <= 15) %>%
  group_by(Immunity, Snails, DDF, Endemicity) %>%
  summarise(epg_mean = mean(epg_mean),
            geom_epg_mean = mean(geom_epg_mean))

#Now to have data after 1 round of MDA
exposure = "ICL"
setting <- paste(exposure, "func_SACmda", sep = "_")
ind.output.dir <- file.path(source.dir, paste("Output/Individual/", setting, sep = "")) 
data_all <- list.files(path = ind.output.dir,  # Identify all output CSV files
                       pattern = "^Ind_out_seed", full.names = TRUE) %>% 
  lapply(readRDS) %>%              # Store all files in list
  bind_rows %>%
  dplyr::filter(time == parms$mda$start+1,
                age >= 5 & age <= 15) %>%
  dplyr::select(seed, time, Immunity, Snails, DDF, Endemicity, ec) #and repeat averaging above

epg_post_ICL <- data_all %>%
  #filter(avg_age_group >= 5 & avg_age_group <= 15) %>%
  group_by(time, Immunity, Snails, DDF, Endemicity) %>%
  summarise(epg_mean = mean(ec),
            geom_epg_mean = (geom_mean(ec+1)))

ERR_ICL <- (epg_pre_ICL$epg_mean - epg_post_ICL$epg_mean)/epg_pre_ICL$epg_mean
quantile(ERR_ICL, probs = c(0.05, 0.5, 0.95), na.rm = T)

ERR <- left_join(
  bind_rows(epg_post_ICL %>%
                  rename(epg_mean_post = epg_mean,
                         geom_epg_mean_post = geom_epg_mean) %>%
                  mutate(exposure = "ICL"),
                  
                  epg_post_Sow %>%
                    rename(epg_mean_post = epg_mean,
                           geom_epg_mean_post = geom_epg_mean) %>%
                    mutate(exposure = "Sow")),
  bind_rows(epg_pre_ICL %>%
            rename(epg_mean_pre = epg_mean,
                  geom_epg_mean_pre = geom_epg_mean) %>%
            mutate(exposure = "ICL"),
                             
            epg_pre_Sow %>%
                    rename(epg_mean_pre = epg_mean,
                           geom_epg_mean_pre = geom_epg_mean) %>%
                    mutate(exposure = "Sow")),
  by = c("Immunity", "Snails", "DDF", "Endemicity", "exposure"))

ERR <- ERR %>%
  mutate(ERR = (epg_mean_pre - epg_mean_post)/epg_mean_pre)

quantile(ERR$ERR[ERR$DDF=="Absent"], probs = c(0.05, 0.5, 0.95), na.rm = T)
quantile(ERR$ERR[ERR$DDF=="Mild"], probs = c(0.05, 0.5, 0.95), na.rm = T)
quantile(ERR$ERR[ERR$DDF=="Strong"], probs = c(0.05, 0.5, 0.95), na.rm = T)

##################################################

data_toplot_ind <- bind_rows(data_toplot_Sow, data_toplot_ICL)
save(data_toplot_ind, file = "Individual data for Figure1_10062024.RData")

###############
#Predictions
################
#This code will compute predictions displayed in Fig 4 and 5

#Load data of interest
setting <- paste("Sow", "func_commMDA", sep = "_")
res_Sow <- readRDS(file.path(pop.output.dir, paste(setting, ".RDS", sep = ""))) %>%
  filter(Immunity != "Absent" & Snails != "Absent") %>% #Successful scenarios only 
  mutate(Exposure = "Based on water contacts")
res_Sow$Endemicity <- factor(as.factor(res_Sow$Endemicity),
                             levels = c("Low", "Moderate", "High"))

#Computing runs that reach EPHP or elimination:
#This is to write the table
options(digits = 1)
# res$Exposure <- factor(as.factor(res$Exposure),
#                        levels = c("Model-derived function", "Based on water contacts"),
#                        labels = c("MD", "WC"))
ephp <- res_Sow %>%
  filter(time==(parms$mda$end+50)*12) %>%
  group_by(Snails, Immunity, DDF, Exposure, Endemicity) %>%
  summarise(ephp_seed = length(which(Heggs_prev<=0.01))*100/seeds,
            eliminated = length(which(eggs_prev_SAC==0))*100/seeds)
# ephp <- ephp %>%
#   unite(Snails_DDF, c(Snails,DDF), sep = ", ")
# ephp$Endemicity <- factor(as.factor(ephp$Endemicity),
#                              levels = c("Low", "Moderate", "High"))

library(flextable)

table <- nice_table(ephp, options(digits = 1))
print(table, preview ="docx")

#Check distribution of worms
# #worms <- 
#   ggplot(filter(data_toplot, Exposure == "Based on water contacts")) +   
#     geom_line(aes(x=avg_age_group, y=wp_mean, group = interaction(Immunity, DDF, Snails), 
#                   colour = Snails)) +
#     facet_grid(Endemicity ~ Immunity, labeller = labeller(.rows = label_both, .cols = label_both), 
#                scales = "free_y") + 
#     tag_facets(tag_levels = "a", position = "tr") +
#     scale_y_continuous(name = "Average worm load \n",
#                        expand = expansion(mult = c(0, 0.1), 
#                                           add = c(0, 0))) +
#     scale_x_continuous(name = "\n Age [years]",
#                        breaks = seq(0, 80, 10),
#                        limits = c(0, 70),
#                        expand = c(0, 0)) +
#     expand_limits(x = 0,y = 0) +
#     scale_color_manual(name = "DDF",
#                        values = c(hue_pal()(3)[1], "purple", hue_pal()(3)[3])) +
#     theme_bw() +
#     theme(legend.position="bottom",
#           plot.margin = margin(5, 10, 0, 10, "pt"),
#           panel.spacing = unit(1.5, "lines"),
#           legend.key.width = unit(2, "lines"),
#           strip.text = element_text(size = 13),
#           tagger.panel.tag.text = element_text(size = 14),
#           tagger.panel.tag.background = element_blank()) 

###################
# Check intensities after treatment for successful scenarios only
###################
# Successful scenarios
# res <- readRDS(file.path(pop.output.dir, "Sow_func_successfulscen.RDS"))
# res <- res %>%
#   mutate(Exposure = "Based on water contacts")
# res$Endemicity <- factor(as.factor(res$Endemicity),
#                                  levels = c("Low", "Moderate", "High"))
# 
# # Averaging population data
# data_avg3 <- res %>%
#   filter(time == parms$mda$start*12-36 | time == parms$mda$end*12+12 | time == parms$mda$end*12+60) %>% #pre-control:1y before start, #soon before last round, #2ys after last round
#   group_by(time, Immunity, Snails, DDF, Endemicity, Exposure) %>% 
#   summarise(eggs_prev_SAC = mean(eggs_prev_SAC),
#             eggs_prev_tot = mean(eggs_prev),
#             PHI = mean(Heggs_prev),
#             intensity = mean(avg_intensity_SAC*24)) 
# 
# ggplot(data_avg3, aes(x = intensity, y = eggs_prev_SAC*100, 
#                       group = interaction(time, Immunity, Snails, DDF))) +
#   geom_point(aes(colour = as.factor(time)), size = 3, alpha = 0.8) +
#   facet_wrap(Endemicity ~ ., nrow = 3, scales = "free", strip.position="right")+
#              #labeller = labeller(.rows = label_both, .cols = label_both)) +
#   tag_facets(tag_levels = "a", position = "tr") +
#   scale_y_continuous(name = "Prevalence of infection in school-aged children (%) \n",
#                      expand = expansion(mult = c(0, 0), 
#                                         add = c(0, 10))) + 
#   scale_x_continuous(name = "\n Average intensity of infection in SAC (epg)",
#                      expand = expansion(mult = c(0, 0), 
#                                         add = c(0, 3))) +
#   expand_limits(x = 0,y = 0) +
#   scale_color_manual(name = "Time",
#                      labels = c("Pre-control", "After 1 year", "After 5 years"),
#                      values = c(hue_pal()(3)[1], "purple", hue_pal()(3)[3])) +
#   theme_bw() +
#   theme(legend.position="bottom",
#         plot.margin = margin(5, 10, 0, 10, "pt"),
#         panel.spacing = unit(1, "lines"),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         #legend.key.width = unit(2, "lines"),
#         strip.text = element_text(size = 13),
#         tagger.panel.tag.text = element_text(size = 14),
#         tagger.panel.tag.background = element_blank())


#############
# Check worms in SAC pre- and post-MDA, Sow func, with and without DDF
#############
# data_pre <- data_all %>%
#   filter(Snails == "Mild", Immunity == "Mild", Endemicity == "High") %>%
#   filter(age >= 5 & age <= 15)
# 
# data <- data_pre
# summary(data$tot_wp[which(data$DDF == "Absent")])
# summary(data$tot_wp[which(data$DDF == "Mild")])
# summary(data$tot_wp[which(data$DDF == "Strong")])

#This is at time 141 (soon before MDA)
# > summary(data$tot_wp[which(data$DDF == "Absent")])
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0       2      22      95     109    2151 
# > summary(data$tot_wp[which(data$DDF == "Mild")])
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0       2      22      92     109    1645 
# > summary(data$tot_wp[which(data$DDF == "Strong")])
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0       2      23      92     109    2489 

#This is at time 191 (after MDA)
# > summary(data$tot_wp[which(data$DDF == "Absent")])
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0       2      28     107     131    1956 
# > summary(data$tot_wp[which(data$DDF == "Mild")])
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0       2      26     107     122    1737 
# > summary(data$tot_wp[which(data$DDF == "Strong")])
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0       2      24      98     116    1564 

# bind_rows(mutate(data_pre, Control = "Pre"),
#           mutate(data_post, Control = "Post")) %>%
# ggplot(aes(x = DDF, y = tot_wp)) +
#   #geom_density() +
#   geom_violin(aes(colour = Control)) +
#   #facet_grid(DDF ~ .) +
#   #scale_y_log10() +
#   theme_bw()
#   

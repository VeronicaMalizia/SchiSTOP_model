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
library(rempsyc)

#Set folder
source.dir <- dirname(getActiveDocumentContext()$path)
#"C:/Users/Z541213/Documents/Project/Model/Schisto_model"
setwd(source.dir)

#Load functions
source("01_Handy_functions.R")

#Load parameters
source("02_Parameters_Smansoni.R")
source("Setting_simulation_scenario.R")
######### Setting #########

#Behavior in exposure
exposure = "Sow" #Choices: "ICL" (model-derived), "Sow" (water contacts)
setting <- paste(exposure, "func", sep = "")

#######################
#Load population output
#######################
#Load collated results and produce multi-panel plots
pop.output.dir <- file.path(source.dir, "Output/Population")

setting <- paste("ICL", "func", sep = "_")
res_ICL <- readRDS(file.path(pop.output.dir, paste(setting, ".RDS", sep = ""))) %>%
  #filter(time/12 > (parms$mda$end - 20) & time/12 < (parms$mda$end + 20)) %>%
  mutate(Exposure = "Model-derived function")

setting <- paste("Sow", "func", sep = "_")
res_Sow <- readRDS(file.path(pop.output.dir, paste(setting, ".RDS", sep = ""))) %>%
  #filter(time/12 > (parms$mda$end - 20) & time/12 < (parms$mda$end + 20)) %>%
  mutate(Exposure = "Based on water contacts")

res <- bind_rows(res_Sow, res_ICL)
res$Exposure <- factor(as.factor(res$Exposure),
                            levels = c("Model-derived function", "Based on water contacts"))
res$Endemicity <- factor(as.factor(res$Endemicity),
                              levels = c("Low", "Moderate", "High"))

save(res, data_avg, file = "Population data for Figure2&3.RData")
load("Population data for Figure2&3.RData")

#Computing faded-out runs at pre-control (149 ys is soon before MDA)
faided <- res %>%
  filter(time==149*12) %>%
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
exposure = "ICL"
setting <- paste(exposure, "func", sep = "_")
ind.output.dir <- file.path(source.dir, paste("Output/Individual/All_SACMDA/", setting, sep = "")) 
data_all <- list.files(path = ind.output.dir,  # Identify all output CSV files
                       pattern = "^Ind_out_seed", full.names = TRUE) %>% 
  lapply(readRDS) %>%              # Store all files in list
  #lapply(read_csv, show_col_types = F) %>%
  bind_rows %>%
  dplyr::filter(time == parms$mda$start-9) #, # | time == parms$mda$end+1)
                #Snails == "Mild" & DDF == "Mild") 


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
data_toplot_ICL <- age_out %>%
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

data_toplot_ind <- bind_rows(data_toplot_Sow, data_toplot_ICL)
save(data_toplot_ind, file = "Individual data for Figure1_ALL.RData")


##############
# Age-intensity profiles
##############
data_toplot_ind$Exposure <- factor(as.factor(data_toplot_ind$Exposure),
                                   levels = c("ICL", "Sow"),
                                   labels = c("Model-derived function", "Based on water contacts"))
data_toplot_ind$Endemicity <- factor(as.factor(data_toplot_ind$Endemicity),
                                     levels = c("Low", "Moderate", "High"))

load("Individual data for Figure1_ALL.RData")
data_toplot <- data_toplot_ind %>%
  filter(!(Snails == "Absent" & Exposure == "Model-derived function"))  #No equilibrium
data_toplot$Endemicity <- factor(as.factor(data_toplot$Endemicity),
                                     levels = c("Low", "Moderate", "High"))

# Ietune <- data.frame(Age = c(2,4,6,8,10,12,14,16,18,22,28,38,50,68),
#                      Epg = c(0,10,40,100,120,160,340,360,260,210,120,70,30,50))
Fulford.data <- list.files(path = getwd(),  # Identify all output CSV files
                           pattern = "^Ageprofiles", full.names = TRUE) %>% 
  lapply(read_csv2, show_col_types = F) %>%
  bind_rows 
colnames(Fulford.data) <- c("Age", "Eggs", "Village", "Endemicity")
Fulford.data$Endemicity <- factor(as.factor(Fulford.data$Endemicity),
                                 levels = c("Low", "Moderate", "High"))
eggs <- 
ggplot(data_toplot)+ #group = Immunity
  #geom_pointrange(aes(ymin=epg_lo, ymax=epg_hi, colour = DDF, shape = Immunity)) +
  #geom_line(aes(linetype = Immunity), colour = hue_pal()(3)[2], size = 1) +
  geom_line(aes(x=avg_age_group, y=(geom_epg_mean-1)*24, group = interaction(Immunity, DDF, Snails), colour = Immunity)) +
  geom_line(data = Fulford.data,
            aes(x=Age, y=Eggs, linetype = Village), colour = "black") +
  #labs(title = paste(endem, "endemicity setting", sep = " ")) +
  facet_grid(Endemicity ~ Exposure, labeller = labeller(.rows = label_both, .cols = label_both), 
             scales = "free_y") + 
  tag_facets(tag_levels = "a", position = "tr") +
  scale_y_continuous(name = "Average egg per gram of feaces \n (geometric mean)",
                     #breaks = seq(0, 10000, 200),
                     #limits = c(0.5, 8000), #3000 l-m 8000 high
                     expand = expansion(mult = c(0, 0.1), 
                                        add = c(0, 0))) +
  scale_x_continuous(name = "\n Age [years]",
                     breaks = seq(0, 80, 10),
                     limits = c(0, 70),
                     expand = c(0, 0)) +
  expand_limits(x = 0,y = 0) +
  #guides(shape = "none", colour = "none") +
  scale_color_manual(name = "(Snails, Exposure)",
                     values = c(hue_pal()(3)[1], "purple", hue_pal()(3)[3])) +
  #scale_linetype_manual(values = c("Absent" = "solid", "Mild" = "dotted", "Strong" = "dashed")) +
  theme_bw() +
  theme(legend.position="bottom",
        plot.margin = margin(5, 10, 0, 10, "pt"),
        panel.spacing = unit(1.5, "lines"),
        legend.key.width = unit(2, "lines"),
        strip.text = element_text(size = 13),
        tagger.panel.tag.text = element_text(size = 14),
        tagger.panel.tag.background = element_blank()) 

tiff(paste("Plots/Manuscript/Fig1.tif", sep = ""), 
     compression="lzw", width=12, height=12, units = "in", res = 600)
eggs
dev.off()


##############
# Observed pattern
##############
#Check
# data_avg %>%
#   filter(time/12==parms$mda$start-1 & eggs_prev_SAC < 0.08)
# data_avg2 <- data_avg %>%
#   filter(!(Immunity == "Absent" & Snails == "Absent" & DDF == "Mild" & Exposure == "Based on water contacts")) %>%
#   filter(!(Immunity == "Absent" & Snails == "Absent" & DDF == "Strong" & Exposure == "Based on water contacts"))
  
Fig2 <- 
  #filter(data_avg, Snails == "Mild", Immunity != "Mild", DDF == "Mild") %>%
  ggplot(data_avg2, aes(x=time/12, eggs_prev_SAC*100)) +
  # geom_line(data = filter(res, Snails == "Mild", Immunity != "Mild", DDF == "Mild"),
  #                         aes(group = interaction(DDF, seed, Immunity), colour = DDF), alpha = 0.08) +
  #geom_line(aes(linetype = Immunity, colour = DDF), size = 1) +
  geom_line(aes(group = interaction(Immunity, DDF, Snails), colour = Immunity)) +
  facet_grid(Endemicity ~ Exposure, labeller = labeller(.rows = label_both, .cols = label_both), 
             scales = "free_y") + 
  tag_facets(tag_levels = "a", position = "tr") +
  scale_y_continuous(name = "Prevalence of infection in school-aged children (%) \n",
                     #breaks = seq(0, 100, 20),
                     #limits = c(0, 100),
                     expand = expansion(mult = c(0, 0), 
                                        add = c(0, 6))) +
  scale_x_continuous(name = "\n Round of treatment",
                     breaks = seq(parms$mda$end-10, parms$mda$end, 1),
                     labels = seq(0, 10, 1),
                     #limits = c(0, 1200),
                     expand = c(0, 0)) +
  coord_cartesian(xlim=c(parms$mda$end-10, parms$mda$end+0.5)) +
  #scale_color_manual(values = c("Mild" = hue_pal()(3)[2])) +
  expand_limits(x = 0,y = 0) +
  #scale_linetype_manual(values = c("Absent" = "solid", "Strong" = "dashed")) +
  theme_bw() +
  theme(legend.position="bottom",
        plot.margin = margin(5, 10, 0, 10, "pt"),
        panel.spacing = unit(1.5, "lines"),
        legend.key.width = unit(2, "lines"),
        strip.text = element_text(size = 13),
        tagger.panel.tag.text = element_text(size = 14),
        tagger.panel.tag.background = element_blank())

tiff(paste("Plots/Manuscript/Fig2_ALL.tif", sep = ""), 
     compression = "lzw", width=11, height=11, units = "in", res = 600)
Fig2 
dev.off()


res2 <- res %>%
  #filter(Endemicity == "Moderate", Snails != "Strong", Immunity != "Mild", DDF != "Strong") %>%
  #filter(!(Exposure == "Model-derived function" & Snails == "Absent" & Immunity == "Absent" & DDF == "Mild")) %>%
  #filter(!(Exposure == "Model-derived function" & Snails == "Mild" & Immunity == "Absent" & DDF == "Absent"))
  filter(!(Exposure == "Based on water contacts" & Immunity == "Absent")) %>% #No age-profiles
  filter(!(Snails == "Absent" & Exposure == "Model-derived function"))  #No equilibrium

# Successful scenarios
# res3 <- res %>%
#   filter(Exposure == "Based on water contacts" & Immunity != "Absent" & Snails != "Absent") 
  
# Averaging population data
data_avg2 <- res2 %>%
  group_by(time, Immunity, Snails, DDF, Endemicity, Exposure) %>%
  summarise(eggs_prev_SAC = mean(eggs_prev_SAC),
            eggs_prev_tot = mean(eggs_prev),
            PHI = mean(Heggs_prev),
            miracidiae = mean(miracidiae),
            cercariae = mean(cercariae),
            snail_inf = mean(inf_snail),
            snail_exp = mean(exp_snail),
            snail_prev = mean(inf_snail/(susc_snail+inf_snail+exp_snail))) %>%
  mutate(Snails2 = case_when(Snails == "Absent" ~ "Absent",
                             Snails == "Mild" ~ "Mild / Strong",
                             Snails == "Strong" ~ "Mild / Strong"))
Fig3 <- 
  ggplot(data_avg2, aes(x=time/12, eggs_prev_SAC*100)) +
  # geom_line(data = res,
  #           aes(group = interaction(DDF, seed, Immunity), colour = DDF), alpha = 0.03) +
  #geom_line(aes(linetype = Immunity, colour = DDF), size = 1) +
  geom_line(aes(group = interaction(DDF, Snails, Immunity, Exposure, Endemicity), 
                colour = interaction(Snails2, Exposure))) +
  facet_grid(Endemicity ~ Immunity, labeller = labeller(.rows = label_both, .cols = label_both), 
             scales = "free_y") + 
  tag_facets(tag_levels = "a", position = "tr") +
  scale_y_continuous(name = "Prevalence of infection in school-aged children (%) \n",
                     #breaks = c(10, 30, 60),
                     #limits = c(0, 80),
                     expand = expansion(mult = c(0, 0), 
                                        add = c(0, 5))) + #5
  scale_x_continuous(name = "\n Years since last treatment round",
                     breaks = seq(parms$mda$end-10, parms$mda$end+20, 10),
                     labels = seq(-10, 20, 10),
                     #limits = c(0, 1200),
                     expand = c(0, 0)) +
  coord_cartesian(xlim=c(parms$mda$end-10, parms$mda$end+20)) +
  expand_limits(x = 0,y = 0) +
  #scale_linetype_manual(values = c("Absent" = "solid", "Strong" = "dashed")) +
  scale_color_manual(name = "(Snails, Exposure)",
                     labels = c("(Mild / Strong, Model-derived function)", "(Absent, Based on water contacts)", "(Mild / Strong, Based on water contacts)"),
                     values = c(hue_pal()(3)[1], "purple", hue_pal()(3)[3])) +
  theme_bw() +
  theme(legend.position="bottom",
        plot.margin = margin(5, 10, 0, 10, "pt"),
        panel.spacing = unit(1.5, "lines"),
        legend.key.width = unit(2, "lines"),
        strip.text = element_text(size = 13),
        tagger.panel.tag.text = element_text(size = 14),
        tagger.panel.tag.background = element_blank())

tiff(paste("Plots/Manuscript/Prova.tif", sep = ""), 
     compression = "lzw", width=12, height=10, units = "in", res = 300)
Fig3
dev.off()

###################
# Check intensities after treatment for successful scenarios only
###################
# Successful scenarios
load("Sow_func_successfulscen.RDS")
res <- res %>%
  mutate(Exposure = "Based on water contacts")
res$Endemicity <- factor(as.factor(res$Endemicity),
                                 levels = c("Low", "Moderate", "High"))

# Averaging population data
data_avg3 <- res %>%
  filter(time == parms$mda$end*12-1 | time == parms$mda$end*12+24) %>% #pre-control:1y before start, #soon before last round, #2ys after last round
  group_by(time, Immunity, Snails, DDF, Endemicity, Exposure) %>% #time == parms$mda$start*12-12
  summarise(eggs_prev_SAC = mean(eggs_prev_SAC),
            eggs_prev_tot = mean(eggs_prev),
            PHI = mean(Heggs_prev),
            intensity = mean(avg_geom_intensity_SAC*24)) 

ggplot(data_avg3, aes(x = intensity, y = eggs_prev_SAC*100, 
                      group = interaction(time, Immunity, Snails, DDF))) +
  geom_point(aes(colour = as.factor(time))) +
  facet_grid(Endemicity ~ . , labeller = labeller(.rows = label_both, .cols = label_both), 
             scales = "free_y") +
  tag_facets(tag_levels = "a", position = "tr") +
  scale_y_continuous(name = "Prevalence of infection in school-aged children (%) \n") + 
  scale_x_continuous(name = "\n Average intensity of infection (epg)",
                     expand = c(0, 0)) +
  expand_limits(x = 0,y = 0) +
  theme_bw() 

###############
#Predictions
################
#Figure 4
# setting <- paste("Sow", "func_commMDA", sep = "_")
# res_Sow <- readRDS(file.path(pop.output.dir, paste(setting, ".RDS", sep = ""))) %>%
#   mutate(Exposure = "Based on water contacts")
# res_Sow$Endemicity <- factor(as.factor(res_Sow$Endemicity),
#                              levels = c("Low", "Moderate", "High"))

#Computing runs that reach EPHP or elimination:
options(digits = 1)
res$Exposure <- factor(as.factor(res$Exposure),
                       levels = c("Model-derived function", "Based on water contacts"),
                       labels = c("a.", "b."))
ephp <- res %>%
  filter(time==(parms$mda$end+50)*12) %>%
  group_by(Endemicity, Snails, Immunity, DDF, Exposure) %>%
  summarise(ephp_seed = length(which(Heggs_prev<=0.01))*100/seeds,
            eliminated = length(which(eggs_prev_SAC==0))*100/seeds) 
table <- nice_table(ephp, options(digits = 1))
print(table, preview ="docx")

library(flextable)

#Target eliminated seeds
# el_seeds <- res %>%
#   mutate(eliminated = ifelse((time==(parms$mda$end+50)*12&eggs_prev_SAC==0), 1, 0)) 
# tmp <- el_seeds[which(el_seeds$eliminated==1), c("Snails", "Immunity", "DDF", "seed")] 

#Scenarios of interest for Figure 4
## Imm NO - Snails NO - DDF all
## Imm STRONG - Snails NO - DDF all
## Imm NO - Snails MILD - DDF all
res <- res_Sow %>%
  filter((Immunity == "Absent" & Snails == "Absent") |
           (Immunity == "Strong" & Snails == "Absent") |
           (Immunity == "Absent" & Snails == "Strong"))
#res <- filter(res, !(Immunity == "Absent" & Snails == "Absent" & DDF == "Absent"))

Fig4high <- 
  ggplot(data=filter(res, Endemicity=="Moderate" & DDF!="Mild"), 
         aes(x=time/12, y=eggs_prev_SAC*100, group = seed)) +
  geom_line(aes(colour = DDF), size = 1, alpha = 0.1) +
  facet_grid(Immunity+Snails ~ DDF, labeller = labeller(.rows = label_both, .cols = label_both)) +
    tag_facets(tag_levels = "a", position = "tr") +
    scale_y_continuous(name = "Prevalence of infection in school-aged children (%) \n",
                       #breaks = seq(0, 100, 20),
                       #limits = c(0, 80),
                       expand = expansion(mult = c(0, 0), 
                                          add = c(3, 5))) +
    scale_x_continuous(name = "\n Years since last treatment round",
                       breaks = seq(parms$mda$end-10, parms$mda$end+50, 10),
                       labels = seq(-10, 50, 10),
                       #limits = c(0, 1200),
                       expand = c(0, 0)) +
    coord_cartesian(xlim=c(parms$mda$end-10, parms$mda$end+50)) +
    expand_limits(x = 0,y = 0) +
    scale_color_manual(values = c("Absent" = hue_pal()(3)[1], "Strong" = hue_pal()(3)[3])) +
    guides(color = guide_legend(override.aes = list(alpha = 1))) +
    theme_bw() +
    theme(legend.position="bottom",
          plot.margin = margin(5, 10, 0, 10, "pt"),
          panel.spacing = unit(1.5, "lines"),
          legend.key.width = unit(2, "lines"),
          strip.text = element_text(size = 13),
          tagger.panel.tag.text = element_text(size = 14),
          tagger.panel.tag.background = element_blank())

tiff(paste("Plots/Manuscript/Fig4high.tif", sep = ""), 
     width=12, height=12, units = "in", res = 300)
Fig4high 
dev.off()

#############
# Equilibria (not useful for manuscript I think)
#############

ggplot(data_avg, aes(x=time/12, eggs_prev_SAC*100)) +
  # geom_line(data = res,
  #           aes(group = interaction(DDF, seed, Immunity), colour = DDF), alpha = 0.08) +
  geom_line(aes(colour = interaction(DDF, Immunity, Snails))) +
  # facet_grid( ~ Exposure, labeller = labeller(.rows = label_both, .cols = label_both), 
  #            scales = "free_y") + 
  tag_facets(tag_levels = "a", position = "tr") +
  scale_y_continuous(name = "Prevalence of infection \n in school-aged children (%) \n",
                     breaks = seq(0, 100, 20),
                     limits = c(0, 30),
                     expand = expansion(mult = c(0, 0), 
                                        add = c(0, 6))) +
  scale_x_continuous(name = "\n Time (Years)",
                     #breaks = seq(0, T, 10),
                     #limits = c(0, 1200),
                     expand = c(0, 0)) +
  expand_limits(x = 0,y = 0) +
  #scale_x_break(breaks = c(500, 700)) +
  theme_bw() +
  theme(legend.position="bottom",
        plot.margin = margin(5, 10, 0, 10, "pt"),
        panel.spacing = unit(1.5, "lines"),
        legend.key.width = unit(2, "lines"),
        strip.text = element_text(size = 13),
        tagger.panel.tag.text = element_text(size = 14),
        tagger.panel.tag.background = element_blank())

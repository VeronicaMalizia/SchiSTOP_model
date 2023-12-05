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

setting <- paste("ICL", "func100seeds", sep = "_")
res_ICL <- readRDS(file.path(pop.output.dir, paste(setting, ".RDS", sep = ""))) %>%
  #filter(time/12 > (parms$mda$end - 20) & time/12 < (parms$mda$end + 20)) %>%
  mutate(Exposure = "Model-derived function")

setting <- paste("Sow", "func100seeds", sep = "_")
res_Sow <- readRDS(file.path(pop.output.dir, paste(setting, ".RDS", sep = ""))) %>%
  #filter(time/12 > (parms$mda$end - 20) & time/12 < (parms$mda$end + 20)) %>%
  mutate(Exposure = "Based on water contacts")

res <- bind_rows(res_Sow, res_ICL)
res$Exposure <- factor(as.factor(res$Exposure),
                            levels = c("Model-derived function", "Based on water contacts"))
res$Endemicity <- factor(as.factor(res$Endemicity),
                              levels = c("Low", "Moderate", "High"))
rm(res_ICL, res_Sow)

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


save(res, data_avg, file = "Population data for Figure2&3_100seeds.RData")
load("Population data for Figure2&3.RData")

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
exposure = "Sow"
setting <- paste(exposure, "func", sep = "_")
ind.output.dir <- file.path(source.dir, paste("Output/Individual/All_SACMDA/", setting, sep = "")) 
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

##ICL reference
reference <- data_toplot_ind %>%
  filter(Snails == "Absent" & Immunity == "Absent" & DDF == "Strong" & Exposure == "Model-derived function")

eggs <- 
ggplot(data_toplot)+ #group = Immunity
  #geom_pointrange(aes(ymin=epg_lo, ymax=epg_hi, colour = DDF, shape = Immunity)) +
  #geom_line(aes(linetype = Immunity), colour = hue_pal()(3)[2], size = 1) +
  geom_line(aes(x=avg_age_group, y=(geom_epg_mean-1)*24, 
                group = interaction(Immunity, DDF, Snails), colour = Immunity)) +
  geom_line(data = Fulford.data,
            aes(x=Age, y=Eggs, linetype = Village), colour = "black") +
  geom_line(data = reference,
            aes(x=avg_age_group, y=(geom_epg_mean-1)*24), colour = "grey40", linewidth = 1) +
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
  scale_color_manual(name = "Human-level regulation",
                     values = c(hue_pal()(3)[1], "purple", hue_pal()(3)[3])) +
  #scale_linetype_manual(values = c("dotdashed", "dotted", "dashed")) +
  theme_bw() +
  theme(legend.position="bottom",
        plot.margin = margin(5, 10, 0, 10, "pt"),
        panel.spacing = unit(1.5, "lines"),
        legend.key.width = unit(2, "lines"),
        strip.text = element_text(size = 13),
        tagger.panel.tag.text = element_text(size = 14),
        tagger.panel.tag.background = element_blank()) 

tiff(paste("Plots/Manuscript/Fig1.2.tif", sep = ""), 
     compression="lzw", width=12, height=12, units = "in", res = 600)
eggs
dev.off()

#Check worms
#worms <- 
  ggplot(filter(data_toplot, Exposure == "Based on water contacts")) +   
    geom_line(aes(x=avg_age_group, y=wp_mean, group = interaction(Immunity, DDF, Snails), 
                  colour = Snails)) +
    facet_grid(Endemicity ~ Immunity, labeller = labeller(.rows = label_both, .cols = label_both), 
               scales = "free_y") + 
    tag_facets(tag_levels = "a", position = "tr") +
    scale_y_continuous(name = "Average worm load \n",
                       expand = expansion(mult = c(0, 0.1), 
                                          add = c(0, 0))) +
    scale_x_continuous(name = "\n Age [years]",
                       breaks = seq(0, 80, 10),
                       limits = c(0, 70),
                       expand = c(0, 0)) +
    expand_limits(x = 0,y = 0) +
    scale_color_manual(name = "DDF",
                       values = c(hue_pal()(3)[1], "purple", hue_pal()(3)[3])) +
    theme_bw() +
    theme(legend.position="bottom",
          plot.margin = margin(5, 10, 0, 10, "pt"),
          panel.spacing = unit(1.5, "lines"),
          legend.key.width = unit(2, "lines"),
          strip.text = element_text(size = 13),
          tagger.panel.tag.text = element_text(size = 14),
          tagger.panel.tag.background = element_blank()) 

##############
# Observed pattern
##############
res2 <- res %>%
  #filter(Endemicity == "Moderate", Snails != "Strong", Immunity != "Mild", DDF != "Strong") %>%
  #filter(!(Exposure == "Model-derived function" & Snails == "Absent" & Immunity == "Absent" & DDF == "Mild")) %>%
  #filter(!(Exposure == "Model-derived function" & Snails == "Mild" & Immunity == "Absent" & DDF == "Absent"))
  #filter(!(Exposure == "Based on water contacts" & Immunity == "Absent")) %>% #No age-profiles
  filter(!(Snails == "Absent" & Exposure == "Model-derived function"))  #No equilibrium (DDF strong no equilibrium, but reference)

reference <- res %>%
  filter(Snails == "Absent" & Immunity == "Absent" & DDF == "Strong" & 
           Exposure == "Model-derived function" & Endemicity != "Low") %>%
  group_by(time, Immunity, Snails, DDF, Endemicity, Exposure) %>%
  summarise(eggs_prev_SAC = mean(eggs_prev_SAC),
            eggs_prev_tot = mean(eggs_prev)) %>%
  rename(Humans = Immunity)

# Successful scenarios
# res3 <- res %>%
#   filter(Exposure == "Based on water contacts" & Immunity != "Absent" & Snails != "Absent") 
  
# Averaging population data
data_avg2 <- res2 %>%
  filter(!(Snails == "Absent" & Immunity == "Absent" & DDF == "Absent" & Exposure == "Based on water contacts")) %>% #only if needed
  group_by(time, Immunity, Snails, DDF, Endemicity, Exposure) %>%
  summarise(eggs_prev_SAC = mean(eggs_prev_SAC)) %>%
  mutate(Snails2 = case_when(Snails == "Absent" ~ "Absent",
                             Snails == "Mild" ~ "Mild / Strong",
                             Snails == "Strong" ~ "Mild / Strong")) %>%
  mutate(Immunity2 = case_when(Immunity == "Absent" ~ "Absent / Mild",
                               Immunity == "Mild" ~ "Absent / Mild",
                               Immunity == "Strong" ~ "Strong"))

Fig3 <- 
  ggplot(rename(data_avg2, Humans = Immunity), aes(x=time/12, eggs_prev_SAC*100)) +
  # geom_line(data = res,
  #           aes(group = interaction(DDF, seed, Immunity), colour = DDF), alpha = 0.03) +
  #geom_line(aes(linetype = Immunity, colour = DDF), size = 1) +
  geom_line(aes(group = interaction(DDF, Snails, Humans, Exposure, Endemicity), 
                colour = interaction(Snails2, Exposure))) +
  geom_line(data = reference, aes(x=time/12, eggs_prev_SAC*100), colour = "grey30") +
  facet_grid(Endemicity ~ Humans, labeller = labeller(.rows = label_both, .cols = label_both), 
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
        text = element_text(size = 14),
        tagger.panel.tag.text = element_text(size = 14),
        tagger.panel.tag.background = element_blank())

tiff(paste("Plots/Manuscript/Fig2.tif", sep = ""), 
     compression = "lzw", width=12, height=10, units = "in", res = 300)
Fig3
dev.off()

###################
# Check intensities after treatment for successful scenarios only
###################
# Successful scenarios
res <- readRDS(file.path(pop.output.dir, "Sow_func_successfulscen.RDS"))
res <- res %>%
  mutate(Exposure = "Based on water contacts")
res$Endemicity <- factor(as.factor(res$Endemicity),
                                 levels = c("Low", "Moderate", "High"))

# Averaging population data
data_avg3 <- res %>%
  filter(time == parms$mda$start*12-36 | time == parms$mda$end*12+12 | time == parms$mda$end*12+60) %>% #pre-control:1y before start, #soon before last round, #2ys after last round
  group_by(time, Immunity, Snails, DDF, Endemicity, Exposure) %>% 
  summarise(eggs_prev_SAC = mean(eggs_prev_SAC),
            eggs_prev_tot = mean(eggs_prev),
            PHI = mean(Heggs_prev),
            intensity = mean(avg_intensity_SAC*24)) 

ggplot(data_avg3, aes(x = intensity, y = eggs_prev_SAC*100, 
                      group = interaction(time, Immunity, Snails, DDF))) +
  geom_point(aes(colour = as.factor(time)), size = 3, alpha = 0.8) +
  facet_wrap(Endemicity ~ ., nrow = 3, scales = "free", strip.position="right")+
             #labeller = labeller(.rows = label_both, .cols = label_both)) +
  tag_facets(tag_levels = "a", position = "tr") +
  scale_y_continuous(name = "Prevalence of infection in school-aged children (%) \n",
                     expand = expansion(mult = c(0, 0), 
                                        add = c(0, 10))) + 
  scale_x_continuous(name = "\n Average intensity of infection in SAC (epg)",
                     expand = expansion(mult = c(0, 0), 
                                        add = c(0, 3))) +
  expand_limits(x = 0,y = 0) +
  scale_color_manual(name = "Time",
                     labels = c("Pre-control", "After 1 year", "After 5 years"),
                     values = c(hue_pal()(3)[1], "purple", hue_pal()(3)[3])) +
  theme_bw() +
  theme(legend.position="bottom",
        plot.margin = margin(5, 10, 0, 10, "pt"),
        panel.spacing = unit(1, "lines"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        #legend.key.width = unit(2, "lines"),
        strip.text = element_text(size = 13),
        tagger.panel.tag.text = element_text(size = 14),
        tagger.panel.tag.background = element_blank())

###############
#Predictions
################
#Figure 4
setting <- paste("Sow", "func_commMDA", sep = "_")
res_Sow <- readRDS(file.path(pop.output.dir, paste(setting, ".RDS", sep = ""))) %>%
  filter(Immunity != "Absent" & Snails != "Absent") %>% #Successful scenarios only 
  #filter(Endemicity == "Moderate", Immunity == "Mild", Snails == "Mild", DDF == "Absent") %>%
  mutate(Exposure = "Based on water contacts")
res_Sow$Endemicity <- factor(as.factor(res_Sow$Endemicity),
                             levels = c("Low", "Moderate", "High"))

el_seeds <- res_Sow[which(res_Sow$eggs_prev_SAC==0 & res_Sow$time==(parms$mda$end+50)*12), "seed"]
res_Sow <- res_Sow %>%
  # group_by(Immunity, Endemicity, Snails, DDF) %>%
  # mutate(el_run = ifelse(eggs_prev_SAC==0 & time==(parms$mda$end+50)*12, seed, F)) #%>%
  mutate(eliminated = ifelse(seed %in% el_seeds$seed, "Interruption", "Bounce-back")) 

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
# res <- res_Sow %>%
#   filter((Immunity == "Absent" & Snails == "Absent") |
#            (Immunity == "Strong" & Snails == "Absent") |
#            (Immunity == "Absent" & Snails == "Strong"))
#res <- filter(res, !(Immunity == "Absent" & Snails == "Absent" & DDF == "Absent"))

#Fig4 <- 
  ggplot(data = res_Sow,
         aes(x=time/12, y=eggs_prev_SAC*100, 
         group = interaction(seed, Endemicity, Immunity, Snails, DDF))) +
    #geom_line(aes(colour = eliminated), alpha = 0.3) +
    geom_line(data = filter(res_Sow, eliminated == "Interruption"), colour = hue_pal()(2)[2]) +
    facet_grid(DDF ~ Endemicity, labeller = labeller(.rows = label_both, .cols = label_both)) +
   # tag_facets(tag_levels = "a", position = "tr") +
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
    #coord_cartesian(xlim=c(parms$mda$end-11, parms$mda$end+50)) +
    #coord_cartesian(ylim=c(0, 10)) +
    expand_limits(x = 0,y = 0) +
    scale_color_discrete(name = "") +
    guides(color = guide_legend(override.aes = list(alpha = 1))) +
    theme_bw() +
    theme(legend.position=c(0.8, 0.85),
          plot.margin = margin(5, 10, 0, 10, "pt"),
          panel.spacing = unit(1.5, "lines"),
          legend.key.width = unit(2, "lines"),
          strip.text = element_text(size = 13),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          tagger.panel.tag.text = element_text(size = 14),
          tagger.panel.tag.background = element_blank())

tiff(paste("Plots/Manuscript/Fig4.3.tif", sep = ""), 
     compression = "lzw", width=8, height=6, units = "in", res = 300)
Fig4
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


#############
# Check worms in SAC pre- and post-MDA, Sow func, with and without DDF
#############
data_pre <- data_all %>%
  filter(Snails == "Mild", Immunity == "Mild", Endemicity == "High") %>%
  filter(age >= 5 & age <= 15)

data <- data_pre
summary(data$tot_wp[which(data$DDF == "Absent")])
summary(data$tot_wp[which(data$DDF == "Mild")])
summary(data$tot_wp[which(data$DDF == "Strong")])

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

bind_rows(mutate(data_pre, Control = "Pre"),
          mutate(data_post, Control = "Post")) %>%
ggplot(aes(x = DDF, y = tot_wp)) +
  #geom_density() +
  geom_violin(aes(colour = Control)) +
  #facet_grid(DDF ~ .) +
  #scale_y_log10() +
  theme_bw()
  

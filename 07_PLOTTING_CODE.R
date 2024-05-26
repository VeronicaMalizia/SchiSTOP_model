########################
# Creating Figure 1
#######################
rm(list = ls())

#Loading packages
.libPaths(c("C:/Program Files/R/R-4.1.2/library",.libPaths()))
library(ggplot2)
#library(plotly)
library(tidyverse)
library(patchwork)
library(dplyr)
library(tagger)
library(scales)
library(rstudioapi)

#Set folder
source.dir <- dirname(getActiveDocumentContext()$path)
#"C:/Users/Z541213/Documents/Project/Model/Schisto_model"
setwd(source.dir)

###############################
#Figure 1
##############################

#Loading data
load("Data_for_Fig1.RData")
load("Data_for_Fig1_DDFafterrebuttal.RData")

##### With ggplot and patchwork
A <- ggplot(exposure, aes(x = x, y = y)) +
  geom_step(data = filter(exposure, degree == "Model-based"),
            aes(linetype = degree), linewidth = 1) +
  geom_line(data = filter(exposure, degree == "Water-contacts-based"),
            aes(linetype = degree), linewidth = 1) +
  scale_x_continuous(name = " \n Age (years)",
                     breaks = seq(0, 100, 20),
                     limits = c(0, 100), #3000 l-m 8000 high
                     expand = c(0, 0)) +
  scale_y_continuous(name = "Relative exposure \n",
                     breaks = seq(0, 1, 0.2),
                     limits = c(0, 1),
                     expand = c(0, 0)) +
  expand_limits(x = 0,y = 0) +
  scale_linetype_discrete(name = " ") +
  theme_classic(base_size = 14) +
  theme(legend.position="bottom",
        plot.margin = margin(5, 10, 0, 10, "pt"),
        panel.spacing = unit(1.5, "lines"),
        legend.key.width = unit(2, "lines"),
        strip.text = element_text(size = 13)) 


B <- ggplot(worms_reg, aes(x = wp, y = y)) + 
  geom_line(aes(colour = degree), linewidth = 1) +
  scale_y_continuous(name = "Expected egg counts \n",
                     #breaks = seq(0, 10000, 200),
                     limits = c(0, 300), #3000 l-m 8000 high
                     expand = expansion(mult = c(0, 0.1),
                                        add = c(0, 0))) +
  scale_x_continuous(name = "\n Number of worm pairs",
                     #breaks = seq(0, 80, 10),
                     #limits = c(0, 70),
                     expand = c(0, 0)) +
  expand_limits(x = 0,y = 0) +
  scale_color_manual(name = "Degree of regulation",
                     values = c(hue_pal()(3)[1], "purple", hue_pal()(3)[3])) +
  theme_classic(base_size = 14) +
  theme(legend.position="bottom",
        plot.margin = margin(5, 10, 0, 10, "pt"),
        panel.spacing = unit(1.5, "lines"),
        legend.key.width = unit(2, "lines"),
        strip.text = element_text(size = 13)) 


C <- ggplot(dataf, aes(x = dwp+0.5, y = y)) + 
  geom_line(aes(colour = Immunity), linewidth = 1) +
  scale_y_continuous(name = "Immunity factor \n",
                     breaks = seq(0, 1, 0.2),
                     limits = c(0, 1.1), 
                     expand = c(0, 0)) +
  scale_x_log10(name = "\n Cumulated dead worm pairs + 0.5",
                     #breaks = seq(0, 80, 10),
                     #limits = c(0, 70),
                     ) +
  expand_limits(x = 0,y = 0) +
  scale_color_manual(name = "Degree of regulation",
                     values = c(hue_pal()(3)[1], "purple", hue_pal()(3)[3])) +
  theme_classic(base_size = 14) +
  theme(legend.position="bottom",
        plot.margin = margin(5, 10, 0, 10, "pt"),
        panel.spacing = unit(1.5, "lines"),
        legend.key.width = unit(2, "lines"),
        strip.text = element_text(size = 13))

D <- ggplot(snail_reg, aes(x = x, y = y)) + 
  geom_line(aes(colour = degree), linewidth = 1) +
  scale_y_continuous(name = "Snail birth rate \n",
                     #breaks = seq(0, 10000, 200),
                     limits = c(0, 6000), #3000 l-m 8000 high
                     expand = c(0, 0)) +
  scale_x_continuous(name = "\n Total number of snails",
                     #breaks = seq(0, 80, 10),
                     #limits = c(0, 70),
                     expand = c(0, 0)) +
  expand_limits(x = 0,y = 0) +
  scale_color_manual(name = "Degree of regulation",
                     values = c("purple", hue_pal()(3)[3]),
                     guide = "none") +
  theme_classic(base_size = 14) +
  theme(legend.position="bottom",
        plot.margin = margin(5, 10, 0, 10, "pt"),
        panel.spacing = unit(1.5, "lines"),
        legend.key.width = unit(2, "lines"),
        strip.text = element_text(size = 13)) 

tiff("Fig 1.tif", 
     compression="lzw", width=12, height=10, units = "in", res = 300)
(A + B) / (C + D) + 
  plot_layout(guides = "collect") + plot_annotation(tag_levels = 'A') & theme(legend.position = 'bottom',
                                                                              legend.text = element_text(size=14),
                                                                              plot.tag.position  = c(1, 1))
dev.off()


###############################
#Figure 2
##############################
#Load collated individual-level data
load("Individual data for Figure1_ALL.RData")

#tidy
data_toplot <- data_toplot_ind %>%
  filter(!(Snails == "Absent" & Exposure == "Model-derived function")) %>% #No equilibrium
  mutate(Exposure = ifelse(Exposure == "Model-derived function", "Model-based", "Water-contacts-based"))
data_toplot$Endemicity <- factor(as.factor(data_toplot$Endemicity),
                                 levels = c("Low", "Moderate", "High"))

#Field data
Fulford.data <- list.files(path = getwd(),  # Identify all output CSV files
                           pattern = "^Ageprofiles", full.names = TRUE) %>% 
  lapply(read_csv2, show_col_types = F) %>%
  bind_rows 
colnames(Fulford.data) <- c("Age", "Eggs", "Village", "Endemicity")
Fulford.data$Endemicity <- factor(as.factor(Fulford.data$Endemicity),
                                  levels = c("Low", "Moderate", "High"))

##Common assumptions reference
reference <- data_toplot_ind %>%
  filter(Snails == "Absent" & Immunity == "Absent" & DDF == "Strong" & Exposure == "Model-derived function") %>%
  mutate(Exposure = ifelse(Exposure == "Model-derived function", "Model-based", "Water-contacts-based"))

eggs <- 
  ggplot(data_toplot)+ #group = Immunity
  geom_line(aes(x=avg_age_group, y=(geom_epg_mean-1)*24, 
                group = interaction(Immunity, DDF, Snails), colour = Immunity)) +
  geom_line(data = Fulford.data,
            aes(x=Age, y=Eggs, linetype = Village), colour = "black") +
  geom_line(data = reference,
            aes(x=avg_age_group, y=(geom_epg_mean-1)*24), colour = "grey40", linewidth = 1) +
  facet_grid(Endemicity ~ Exposure, labeller = labeller(.rows = label_both, .cols = label_both), 
             scales = "free_y") + 
  tag_facets(tag_levels = "A", position = "tr") +
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
  theme_bw(base_size = 16) +
  theme(legend.position="bottom",
        plot.margin = margin(5, 10, 0, 10, "pt"),
        panel.spacing = unit(1.5, "lines"),
        legend.key.width = unit(2, "lines"),
        strip.text = element_text(size = 14),
        tagger.panel.tag.text = element_text(size = 14),
        tagger.panel.tag.background = element_blank()) 

tiff(paste("Plots/Manuscript/Fig 2.tif", sep = ""), 
     compression="lzw", width=12, height=12, units = "in", res = 300)
eggs
dev.off()


###############################
#Figure 3
##############################
#Load parameters
#Load functions (in this case it's only needed for expressing time in terms of MDA rounds)
source("01_Handy_functions.R")

#Load parameters
source("02_Parameters_Smansoni.R")
source("Setting_simulation_scenario.R")

#Load collated population-level data
load("Population data for Figure2&3_100seeds.RData")

#Tidy
#Exclude scenarios which do not reproduce pre-control low equilibria 
res2 <- res %>%
  filter(!(Snails == "Absent" & Exposure == "Model-derived function"))  #No equilibrium (DDF strong no equilibrium, but reference)

#Common assumptions reference
reference <- res %>%
  filter(Snails == "Absent" & Immunity == "Absent" & DDF == "Strong" & 
           Exposure == "Model-derived function" & Endemicity != "Low") %>%
  group_by(time, Immunity, Snails, DDF, Endemicity, Exposure) %>%
  summarise(eggs_prev_SAC = mean(eggs_prev_SAC),
            eggs_prev_tot = mean(eggs_prev)) %>%
  rename(Human_level = Immunity) %>%
  mutate(Exposure = ifelse(Exposure == "Model-derived function", "Model-based", "Water-contacts-based"))

# Averaging population data and adjust labels for plotting
data_avg2 <- res2 %>%
  filter(!(Snails == "Absent" & Immunity == "Absent" & DDF == "Absent" & Exposure == "Based on water contacts")) %>% #only if needed
  group_by(time, Immunity, Snails, DDF, Endemicity, Exposure) %>%
  summarise(eggs_prev_SAC = mean(eggs_prev_SAC)) %>%
  mutate(Snails2 = case_when(Snails == "Absent" ~ "Absent",
                             Snails == "Mild" ~ "Mild / Strong",
                             Snails == "Strong" ~ "Mild / Strong")) %>%
  mutate(Immunity2 = case_when(Immunity == "Absent" ~ "Absent / Mild",
                               Immunity == "Mild" ~ "Absent / Mild",
                               Immunity == "Strong" ~ "Strong")) %>%
  mutate(Exposure = ifelse(Exposure == "Model-derived function", "Model-based", "Water-contacts-based"))

#ggplot
Fig3 <- 
  ggplot(rename(data_avg2, Human_level = Immunity), aes(x=time/12, eggs_prev_SAC*100)) +
  geom_line(aes(group = interaction(DDF, Snail, Human_level, Exposure, Endemicity), 
                colour = interaction(Snails2, Exposure))) +
  geom_line(data = reference, aes(x=time/12, eggs_prev_SAC*100), colour = "grey30") +
  facet_grid(Endemicity ~ Human_level, labeller = labeller(.rows = label_both, .cols = label_both), 
             scales = "free_y") + 
  tag_facets(tag_levels = "A", position = "tr") +
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
  scale_color_manual(name = "(Snail-level, Exposure)",
                     labels = c("(Mild / Strong, Model-based)", "(Absent, Water-contacts-based)", "(Mild / Strong, Water-contacts-based)"),
                     values = c(hue_pal()(3)[1], "purple", hue_pal()(3)[3])) +
  theme_bw(base_size = 16) +
  theme(legend.position="bottom",
        plot.margin = margin(5, 10, 0, 10, "pt"),
        panel.spacing = unit(1.5, "lines"),
        legend.key.width = unit(2, "lines"),
        strip.text = element_text(size = 13),
        #text = element_text(size = 14),
        tagger.panel.tag.text = element_text(size = 14),
        tagger.panel.tag.background = element_blank())

tiff(paste("Plots/Manuscript/Fig 3.tif", sep = ""), 
     compression = "lzw", width=13, height=12, units = "in", res = 300)
Fig3
dev.off()







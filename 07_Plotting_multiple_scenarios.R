rm(list = ls())

.libPaths(c("C:/Program Files/R/R-4.1.2/library",.libPaths()))
library(tidyverse)
library(readxl)
library(ggridges)
library(patchwork)
library(egg)
`%!in%` <- Negate(`%in%`)
geom_mean <- function(x){exp(mean(log(x)))}
ci <- function(x){quantile(x, probs=c(0.025, 0.975), na.rm = T)}
################
#SETTING THE SCENARIO

stoch_scenarios <- expand.grid(list(seed = 1:seeds,
                                    DDF_strength = c("Absent", "Mild", "Strong"),
                                    imm_strength = c("Absent", "Mild", "Strong"),
                                    snails = c("Absent", "Mild", "Strong")))

#Load tuned transmission parameters (zetas) for the scenarios above (attention to the order!) 
#Transmission parameters for tuning endemicity
zetas <- read_excel("Zetas.xlsx")
stoch_scenarios <- mutate(stoch_scenarios, zeta = rep(zetas$Zeta, each = seeds))

#Set folder
source.dir <- "C:/Users/Z541213/Documents/Project/Model/Schisto_model"
setwd(source.dir)

################
#Prevalence timelines
#Load population output
#####Load collated results and produce multi-panel plots

load(file.path(source.dir, 
               paste("/Output/Population/Cumulative_pop_results_moderate.RData")))

#Average by seed
data_avg <- res %>%
  group_by(time, Immunity, Snails, DDF) %>%
  summarise(eggs_prev_SAC = mean(eggs_prev_SAC),
            PHI = mean(Heggs_prev),
            snail_prev = mean(inf_snail/(susc_snail+inf_snail+exp_snail)))
#Save collated data
write.csv(data, file = "Population data_collated scenarios.csv")
#Plot prevalence

#Let's create an additional data frame to hold the text annotations:
text <- zetas %>%
  group_by(Snails, Immunity) %>%
  summarise(Snail_prevalence = list(`Snail prevalence`)) 
text$Snail_prevalence[1:3] = NA
my_tag <- paste("Snail prevalence:\n", text$Snail_prevalence)

Fig <- ggplot(data=data_avg, aes(x=time/12)) +
  geom_line(aes(y=eggs_prev_SAC*100, colour = DDF)) +
  labs(title = "Moderate endemicity setting") +
  facet_grid(Snails ~ Immunity, labeller = labeller(.rows = label_both, .cols = label_both)) +
  geom_text(data=text, aes(x=200, y=80, label=my_tag), 
            size=3) +
  scale_y_continuous(name = "Prevalence of infection in SAC (%) \n",
                     breaks = seq(0, 100, 20),
                     limits = c(0, 100),
                     expand = c(0, 0)) +
  scale_x_continuous(name = "\n Time [Years]",
                     #limits = c(0, 1200),
                     expand = c(0, 0)) +
  #coord_cartesian(xlim=c(500, (T*12))) +
  expand_limits(x = 0,y = 0) +
  theme_bw() +
  theme(legend.position="bottom",
        plot.margin = margin(5, 10, 0, 10, "pt"),
        panel.spacing = unit(1, "lines")) 

setwd(source.dir)
tiff("Plots/Modeling_scenarios/PrevalenceSAC.tif", width=10, height=9, units = "in", res = 300)
Fig 
dev.off()
#Zoom on MDA
tiff("Plots/Modeling_scenarios/PrevalenceSAC_MDA.tif", width=10, height=9, units = "in", res = 300)
Fig + 
  # geom_text(data=zetas, aes(x=170, y=80, label=my_tag), 
  #           size=3) +
  coord_cartesian(xlim=c(140, 180)) 
dev.off()

###Heavy intensity prevalence
Fig2 <- ggplot(data=data_avg, aes(x=time/12)) +
  geom_line(aes(y=PHI*100, colour = DDF)) +
  labs(title = "Moderate endemicity setting") +
  facet_grid(Snails ~ Immunity, labeller = labeller(.rows = label_both, .cols = label_both)) +
  geom_text(data=text, aes(x=200, y=80, label=my_tag), 
            size=3) +
  scale_y_continuous(name = "Prevalence of heavy infections (%) \n",
                     breaks = seq(0, 10, 2),
                     limits = c(0, 10),
                     expand = c(0, 0)) +
  scale_x_continuous(name = "\n Time [Years]",
                     #limits = c(0, 1200),
                     expand = c(0, 0)) +
  #coord_cartesian(xlim=c(500, (T*12))) +
  expand_limits(x = 0,y = 0) +
  theme_bw() +
  theme(legend.position="bottom",
        plot.margin = margin(5, 10, 0, 10, "pt"),
        panel.spacing = unit(1, "lines")) 

tiff("Plots/Modeling_scenarios/Prevalence_Heavy.tif", width=10, height=9, units = "in", res = 300)
Fig2 
dev.off()

#PLot prevalence of infected patent snails (only for snails = 'yes')
# tiff("Plots/Modeling_scenarios/I++S+_SnailsPrev.tif", width=7, height=6, units = "in", res = 300)
# ggplot(data_avg, aes(x=time/12, colour = as.factor(Scenario))) +
#   geom_line(aes(y=snail_prev*100)) +
#   labs(title = "Moderate endemicity setting") +
#   guides(col = guide_legend("Scenario")) +
#   scale_y_continuous(name = "Prevalence of infected snails (%)",
#                      breaks = seq(0, 100, 20),
#                      limits = c(0, 100),
#                      expand = c(0, 0)) +
#   scale_x_continuous(name = "Time [Years]",
#                      #limits = c(0, 1200),
#                      expand = c(0, 0)) +
#   expand_limits(x = 0,y = 0) +
#   theme_bw() +
#   theme(legend.position="bottom",
#         plot.margin = margin(5, 10, 0, 10, "pt")) 
# dev.off()

################
#Age-exposure profiles
#Load individual output
#Collate results by seeds
output.dir <- file.path(source.dir, "Output")
# setwd(output.dir)
#Run this below for each of the DDF level
data_all <- list.files(path = output.dir,  # Identify all output CSV files
                       pattern = "*.csv", full.names = TRUE) %>% 
  lapply(read_csv, show_col_types = F) %>%              # Store all files in list
  bind_rows

# Combine data sets into one
#Aggregate by age groups
#Filter one year of interest
#unique(data_all$time) #individual output contains time in years
#years <- c(81, 101, 151, 181) #c(51, 81, 101, 151, 161, 171, 181, 191)
age_out <- data_all %>%
  filter(time == 131) %>%
  #filter(seed %!in% dead.seeds) %>%
  mutate(age_group = case_when(age<=1 ~ "0_1",
                               age>1 & age<=5 ~ "1_5",
                               age>5 & age<=15 ~ "5_15",
                               age>15 & age<=30 ~ "15_30",
                               age>30 & age<=40 ~ "30_40",
                               age>40 ~ "40_90")) %>% 
  group_by(seed, time, age_group, Immunity, Snails, DDF) %>% #sex not for now
  summarise(epg = mean(ec*24), #geom_mean(ec*24+1), #means over individuals of that sex-age group
            wp = mean(tot_wp),
            dwp = mean(cum_dwp), 
            rate = mean(rate)) 

age_out$age_group <- factor(age_out$age_group, levels = c("0_1",
                                                          "1_5",
                                                          "5_15",
                                                          "15_30",
                                                          "30_40",
                                                          "40_90"))

save(age_out, file = file.path(output.dir, "Collated_time131.RData"))


#####Load collated results and produce multi-panel plots
output.dir <- file.path(source.dir, paste("Output/Imm=", imm, "Sn=", snails, sep = ""))
load(file.path(output.dir, "Collated_time131.RData"))

#Average over seeds
data_toplot <- age_out %>%
  #filter(!(Immunity == "Absent" & Snails== "Absent" & DDF== "Absent")) %>%
  group_by(age_group, Immunity, Snails, DDF) %>%
  summarise(epg_mean = mean(epg),
            epg_sd = sd(epg),
            wp_mean = mean(wp),
            wp_sd = sd(wp),
            dwp_mean = mean(dwp),
            dwp_sd = sd(dwp),
            rate_mean = mean(rate),
            rate_sd = sd(rate))

#Mean intensity
eggs <- ggplot(data_toplot, aes(x=age_group, y=epg_mean+1, group = DDF))+
  geom_pointrange(aes(ymin=epg_mean-epg_sd+1, ymax=epg_mean+epg_sd+1, colour = DDF)) +
  geom_line(aes(colour = DDF)) +
  labs(title = "Moderate endemicity setting") +
  facet_grid(Snails ~ Immunity, labeller = labeller(.rows = label_both, .cols = label_both)) +
  scale_y_log10(name = "Observed egg counts (epg) + 1 \n",
                     #breaks = seq(0, 10000, 200),
                limits = c(0.5, 3000),
                     #expand = c(0, 0)
                ) +
  scale_x_discrete(name = "\n Age group [years]") +
  expand_limits(x = 0,y = 0) +
  theme_bw() +
  theme(#axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.position="bottom",
        plot.margin = margin(5, 10, 0, 10, "pt"),
        panel.spacing = unit(1, "lines")) 

tiff("Plots/Modeling_scenarios/Age profiles/Eggs.tif", width=10, height=9, units = "in", res = 300)
eggs
dev.off()

#Parasite establishment rate
rate <- ggplot(data_toplot, aes(x=age_group, y=rate_mean+1, group = DDF))+
  geom_pointrange(aes(ymin=rate_mean-rate_sd+1, ymax=rate_mean+rate_sd+1, colour = DDF)) +
  geom_line(aes(colour = DDF)) +
  labs(title = "Moderate endemicity setting") +
  facet_grid(Snails ~ Immunity, labeller = labeller(.rows = label_both, .cols = label_both)) +
  scale_y_log10(name = "Parasite establishment rate + 1 \n",
                #breaks = seq(0, 10000, 200),
                limits = c(0.5, 10),
                #expand = c(0, 0)
  ) +
  scale_x_discrete(name = "\n Age group [years]") +
  expand_limits(x = 0,y = 0) +
  theme_bw() +
  theme(#axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
    legend.position="bottom",
    plot.margin = margin(5, 10, 0, 10, "pt"),
    panel.spacing = unit(1, "lines")) 

tiff("Plots/Modeling_scenarios/Age profiles/Rate.tif", width=10, height=9, units = "in", res = 300)
rate
dev.off()

#Worms
worms <- ggplot(data_toplot, aes(x=age_group, y=wp_mean+1, group = DDF))+
  geom_pointrange(aes(ymin=wp_mean-wp_sd+1, ymax=wp_mean+wp_sd+1, colour = DDF)) +
  geom_line(aes(colour = DDF)) +
  labs(title = "Moderate endemicity setting") +
  facet_grid(Snails ~ Immunity, labeller = labeller(.rows = label_both, .cols = label_both)) +
  scale_y_log10(name = "Worm load (pairs) + 1 \n",
                #breaks = seq(0, 10000, 200),
                limits = c(0.5, 100),
                #expand = c(0, 0)
  ) +
  scale_x_discrete(name = "\n Age group [years]") +
  expand_limits(x = 0,y = 0) +
  theme_bw() +
  theme(#axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
    legend.position="bottom",
    plot.margin = margin(5, 10, 0, 10, "pt"),
    panel.spacing = unit(1, "lines")) 

tiff("Plots/Modeling_scenarios/Age profiles/Worms.tif", width=10, height=9, units = "in", res = 300)
worms
dev.off()

#Dead worm pairs
dead_worms <- ggplot(data_toplot, aes(x=age_group, y=dwp_mean+1, group = DDF))+
  geom_pointrange(aes(ymin=dwp_mean-dwp_sd+1, ymax=dwp_mean+dwp_sd+1, colour = DDF)) +
  geom_line(aes(colour = DDF)) +
  labs(title = "Moderate endemicity setting") +
  facet_grid(Snails ~ Immunity, labeller = labeller(.rows = label_both, .cols = label_both)) +
  scale_y_log10(name = "Cumulated dead worm (pairs) + 1 \n",
                #breaks = seq(0, 10000, 200),
                limits = c(0.5, 1000),
                #expand = c(0, 0)
  ) +
  scale_x_discrete(name = "\n Age group [years]") +
  expand_limits(x = 0,y = 0) +
  theme_bw() +
  theme(#axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
    legend.position="bottom",
    plot.margin = margin(5, 10, 0, 10, "pt"),
    panel.spacing = unit(1, "lines")) 

tiff("Plots/Modeling_scenarios/Age profiles/Dead_worms.tif", width=10, height=9, units = "in", res = 300)
dead_worms
dev.off()

##Print average worm pairs burden per person
age_out %>% 
  #filter(time==141) %>%
  group_by(time) %>%
  summarise(avg_wp = mean(wp),
            avg_epg = mean(epg))

##############
#Plot demographics
##############
# ##Age distribution by time
#https://r-charts.com/distribution/ggridges/
cohort_plot_age <- cohort %>%
  select(age) %>%
  mutate(ID = 1:nrow(cohort),
         time = 0,
         seed = NA)
bind_rows(cohort_plot_age,
          select(data_all, age, ID, time, seed) %>%
            filter(time %in% c(seq(12, (T*12), 20*12)/12))) %>%
  ggplot(aes(x=round(age))) +
  geom_histogram(aes(y = ..density.., colour=as.factor(seed)), binwidth = 5,
                 fill="white", position = 'dodge') +
  #stat_bin(aes(y=..count.., label=..count..), binwidth = 5, geom="text", vjust=-.5) +
  guides(colour=guide_legend(title="Simulation seed")) +
  scale_color_grey() +
  facet_wrap(~ as.factor(time)) +
  scale_y_continuous(name = "Density",
                     expand = c(0, 0)) +
  scale_x_continuous(name = "Age [years]",
                     expand = c(0, 0)) +
  expand_limits(x = 0,y = 0) +
  theme(legend.position="top")
# 
# 
# data_all %>%
#    filter(time==1, seed==1) %>%
#    ggplot(aes(x=age)) +
#    geom_histogram(binwidth = 5) + 
#    stat_bin(aes(y=..count.., label=..count..), binwidth = 5, geom="text", vjust=-.5) 

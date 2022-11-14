rm(list = ls())

.libPaths(c("C:/Program Files/R/R-4.1.2/library",.libPaths()))
library(tidyverse)
library(readxl)
library(ggridges)
library(patchwork)
`%!in%` <- Negate(`%in%`)
geom_mean <- function(x){exp(mean(log(x)))}
ci <- function(x){quantile(x, probs=c(0.025, 0.975), na.rm = T)}
################
#Set folder
source.dir <- "C:/Users/Z541213/Documents/Project/Model/Schisto_model"
setwd(source.dir)
################
#Prevalence timelines

#Load population output
#SETTING THE SCENARIO
#Anti-reinfection immunity: choose 0=no,	0.00005=mild, 0.0001=strong
imm = 0.00005 #immunity slope parameter
#Snails 
snails = "no" #Choose "yes" or "no"
#If snails is on choose k=xx : mild, k=xx : strong
if(snails == "yes")
  carrying.capacity = 20000
#Density-dependent fecundity
#Choose DDF level: 'no', 'mild', 'strong'      
DDF_strength <- "mild"

load(file.path(source.dir, 
               paste("/Population_output/Res", 
                     "_Imm=", imm, "Sn=", snails, "DDF=", DDF_strength, ".RData", sep="")))
res_worms_0 <- res %>%
  mutate(Scenario = "(Imm++, Snails0, DDF0)")

res_worms_mild <- res %>%
  mutate(Scenario = "(Imm++, Snails0, DDF+)")

res_worms_strong <- res %>%
  mutate(Scenario = "(Imm+, Snails0, DDF++)")

data <- bind_rows(res_worms_0, res_worms_mild, res_worms_strong)
data_avg <- data %>%
  group_by(time, Scenario) %>%
  summarise(eggs_prev_SAC = mean(eggs_prev_SAC))
#Plot prevalence
Fig <- ggplot(data_avg, aes(x=time/12, colour = as.factor(Scenario))) +
  geom_line(aes(y=eggs_prev_SAC*100)) +
  labs(title = "Moderate endemicity setting") +
  guides(col = guide_legend("Scenario")) +
  scale_y_continuous(name = "Prevalence of infection in SAC (%)",
                     breaks = seq(0, 100, 20),
                     limits = c(0, 100),
                     expand = c(0, 0)) +
  scale_x_continuous(name = "Time [Years]",
                     #limits = c(0, 1200),
                     expand = c(0, 0)) +
  #coord_cartesian(xlim=c(500, (T*12))) +
  expand_limits(x = 0,y = 0)

tiff("Plots/Modeling_scenarios/I0S0.tif", width=7, height=6, units = "in", res = 300)
Fig +
  #coord_cartesian(xlim=c(50, 200)) +
  theme_bw()
dev.off()

################
#Age-exposure profiles
#Load individual output
#Collate results by seeds
output.dir <- file.path(source.dir, "/Output/Imm++_Snail0")

#Run this below for each of the DDF level
data_all_DDF <- list.files(path = output.dir,  # Identify all output CSV files
                            pattern = "*_no_.csv", full.names = TRUE) %>% 
  lapply(read_csv, show_col_types = F) %>%              # Store all files in list
  bind_rows                                             # Combine data sets into one
#Aggregate by age groups
#Filter one year of interest
#unique(data_all_DDF0$time) #individual output contains time in years
#years <- c(81, 101, 151, 181) #c(51, 81, 101, 151, 161, 171, 181, 191)
age_out_DDF <- data_all_DDF %>%
  filter(time == 121) %>%
  #filter(seed %!in% dead.seeds) %>%
  mutate(age_group = case_when(age<=1 ~ "0_1",
                               age>1 & age<=5 ~ "1_5",
                               age>5 & age<=15 ~ "5_15",
                               age>15 & age<=30 ~ "15_30",
                               age>30 & age<=40 ~ "30_40",
                               age>40 ~ "40_90")) %>% 
  group_by(seed, time, age_group, sex) %>%
  summarise(epg = mean(ec*24), #geom_mean(ec*24+1), #means over individuals of that sex-age group
            wp = mean(tot_wp),
            dwp = mean(cum_dwp), 
            rate = mean(rate)) %>%
  mutate(Scenario = "(Imm++, Snail0, DDF0)")

age_out_DDF$age_group <- factor(age_out_DDF$age_group, levels = c("0_1",
                                                                  "1_5",
                                                                  "5_15",
                                                                  "15_30",
                                                                  "30_40",
                                                                  "40_90"))
DDF0 <- age_out_DDF
DDFmild <- age_out_DDF
DDFstrong <- age_out_DDF

data <- bind_rows(DDF0, DDFmild, DDFstrong)

box <- ggplot(data)+
  geom_boxplot(aes(x=age_group, y=epg, colour = Scenario), fill = "white", outlier.shape = NA) +
  #facet_wrap(~ time, ncol = 4) +
  labs(title = "Moderate endemicity setting") +
  scale_y_continuous(name = "Observed egg counts (epg) \n",
                     #breaks = seq(0, 10000, 200),
                     #limits = c(0, 1000),
                     expand = c(0, 0)) +
  scale_x_discrete(name = "\n Age group [years]") +
  expand_limits(x = 0,y = 0) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  theme_bw()

tiff("Plots/Modeling scenarios/I++S0_boxplot.tif", width=7, height=6, units = "in", res = 300)
box
dev.off()

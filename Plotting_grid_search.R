# ---
# title: "Explorations of transmission parameters"
# author: "Veronica Malizia"
# date: '25-01-2023'
# output: html_document
# ---
  
# Simulation results

rm(list = ls())

#.libPaths(c("C:/Program Files/R/R-4.1.2/library",.libPaths()))
library(tidyverse)
library(readxl)
library(ggridges)
library(patchwork)
library(egg)
library(rstudioapi)
library(lattice)

`%!in%` <- Negate(`%in%`)
geom_mean <- function(x){exp(mean(log(x)))}
ci <- function(x){quantile(x, probs=c(0.025, 0.975), na.rm = T)}
################
setting <- "Grid_search_refined_panel1_DDFstrong"

#Set folder
source.dir <- dirname(getActiveDocumentContext()$path)
#"C:/Users/Z541213/Documents/Project/Model/Schisto_model"
setwd(source.dir)

#Load population output
#####Load collated results and produce multi-panel plots
grid.output.dir <- file.path(source.dir, "Output/Grid_search")
res <- readRDS(file.path(grid.output.dir, 
               paste(setting, ".RDS", sep = "")))
max.time <- max(res$time)

#Average by seed
data_avg <- res %>%
  filter(time == max.time) %>%
  group_by(time, zeta, worms_aggr, tr_snails, Immunity, Snails, DDF) %>%
  summarise(eggs_prev_SAC = mean(eggs_prev_SAC),
            PHI = mean(Heggs_prev))
            #snail_prev = mean(inf_snail/(susc_snail+inf_snail+exp_snail)))

#Check faded-out simulations:

# data_avg <- data_avg %>%
#   filter(!(time==max.time & eggs_prev_SAC == 0))

# Prevalence heatmap

## Prevalence of any infection in school-aged-children (SAC) 

# Fig <- filter(data_avg, time==max.time) %>% 
#   
#   ggplot(aes(x = as.factor(worms_aggr), y = zeta)) +
#   geom_point(aes(colour = eggs_prev_SAC), size = 6) +
#   #geom_raster(aes(fill = eggs_prev_SAC)) +
#   scale_colour_gradient2(low="yellow", mid="orange", high="red", 
#                        midpoint=0.3, limits=range(data_avg$eggs_prev_SAC)) +
#   theme_classic()

#OR
Fig <- levelplot(eggs_prev_SAC ~ as.factor(tr_snails)*zeta | as.factor(worms_aggr), 
                 data=data_avg[which(data_avg$time==max.time),],
                 at=c(0.1*0:10),
                 xlab="Level of worms aggregation", ylab = "Transmission parameter on humans",
                 main="Heatmap of prevalence in SAC, at the end of simulation.",
                 col.regions = rev(heat.colors(100)))

tiff("Grid search/Panel4_DDFstrong.tif", width=12, height=10, units = "in", res = 300)
Fig 
dev.off()

#Heavy intensity infections
#OR
Fig2 <- levelplot(PHI ~ worms_aggr*zeta, data=data_avg[which(data_avg$time==max.time),],
                  at=c(0.1*0:10),
                  xlab="Level of worms aggregation", ylab = "Transmission parameter on humans",
                  main="Heatmap of heavy prevalence in SAC, at the end of simulation.",
                  col.regions = rev(heat.colors(100)))

tiff("Plots/Grid search/Panel1_DDFstrong_refined_heavy.tif", width=12, height=10, units = "in", res = 300)
Fig2 
dev.off()

##Selection
#Plot prevalence as a function of zeta at the end of simulations
Fig3 <- filter(data_avg, time == max.time) %>%
  ggplot(aes(x=zeta, y=eggs_prev_SAC, colour = interaction(as.factor(worms_aggr),
                                                           as.factor(tr_snails)))) +
  geom_point(size = 3, alpha = 0.7) +
  geom_line() + 
  geom_hline(yintercept = 0.6, linetype = "longdash", size = 1) +
  geom_hline(yintercept = 0.3, linetype = "longdash", size = 1) +
  geom_hline(yintercept = 0.1, linetype = "longdash", size = 1) +
  scale_color_discrete(name = "Level of worms \n aggregation") +
  scale_y_continuous(name = "Prevalence of infection in SAC \n",
                     breaks = seq(0, 1, 0.2),
                     limits = c(0, 1),
                     expand = c(0, 0)) +
  scale_x_continuous(name = "\n Transmission parameter on humans [zeta]",
                     #breaks = seq(0, 4*10^(-4), 0.00005),
                     #limits = c(0, 1200),
                     expand = c(0, 0)) +
  #coord_cartesian(xlim=c(500, (T*12))) +
  expand_limits(x = 0,y = 0) +
  theme_bw()

tiff("Grid search/Panel1_DDFstrong_simualtions.tif", width=12, height=9, units = "in", res = 300)
Fig3
dev.off()

#High (60%)
#I chose:
f <- filter(data_avg, worms_aggr==0.2 #& time == max.time 
            & eggs_prev_SAC >= 0.6 & eggs_prev_SAC <= 0.61)
f$zeta[1]
#zeta = 0.000164 & k_w = 0.2
high <- filter(res, (zeta==f$zeta[1] & worms_aggr == f$worms_aggr[1])) %>%
  mutate(Endemicity = "High")

#Moderate (30%)
f2 <- filter(data_avg, worms_aggr==0.1 & time == max.time 
             & eggs_prev_SAC >= 0.3 & eggs_prev_SAC <= 0.31)
f2
f2$zeta[1]
mod <- filter(res, (zeta==f2$zeta[1] & worms_aggr == f2$worms_aggr[1])) %>%
  mutate(Endemicity = "Moderate")

#Low (10%)
f3 <- filter(data_avg, round(worms_aggr, digits = 2)==0.27 
             & eggs_prev_SAC >= 0.08 & eggs_prev_SAC <= 0.12)
f3
low <- filter(res, (zeta==f3$zeta[1] & worms_aggr == f3$worms_aggr[1])) %>%
  mutate(Endemicity = "Low")

#Check for the equilibrium
Eq <- bind_rows(high, mod, low) %>%
  ggplot(aes(x=time/12, y=eggs_prev_SAC*100, group = interaction(seed, Endemicity))) +
  geom_line(aes(colour = Endemicity), alpha=0.8) +
  labs(title = "Check for the endemic equilibrium") +
  # geom_text(data=text, aes(x=200, y=80, label=my_tag),
  #           size=3) +
  scale_y_continuous(name = "Prevalence of infection in SAC (%) \n",
                     breaks = seq(0, 100, 20),
                     limits = c(0, 100),
                     expand = c(0, 0)) +
  scale_x_continuous(name = "\n Time [Years]",
                     #limits = c(0, 1200),
                     expand = c(0, 0)) +
  expand_limits(x = 0,y = 0) +
  theme_bw() 

tiff("Grid search/Panel1_DDFstrong_equilibrium_onlylow.tif", width=12, height=9, units = "in", res = 300)
Eq
dev.off()


#Check for the equilibrium for LOW endemicity 
#with MULTIPLE worms aggregations
filter(res, zeta==f3$zeta[1] & worms_aggr %in% f3$worms_aggr) %>%
  ggplot(aes(x=time/12, y=eggs_prev_SAC*100, group = seed)) +
  geom_line(aes(colour = worms_aggr), alpha=0.6) +
  labs(title = "Check for the endemic equilibrium") +
  # geom_text(data=text, aes(x=200, y=80, label=my_tag),
  #           size=3) +
  scale_y_continuous(name = "Prevalence of infection in SAC (%) \n",
                     breaks = seq(0, 100, 20),
                     limits = c(0, 100),
                     expand = c(0, 0)) +
  scale_x_continuous(name = "\n Time [Years]",
                     #limits = c(0, 1200),
                     expand = c(0, 0)) +
  #coord_cartesian(xlim=c(500, (T*12))) +
  expand_limits(x = 0,y = 0) +
  theme_bw() 

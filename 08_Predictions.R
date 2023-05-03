rm(list = ls())

.libPaths(c("C:/Program Files/R/R-4.1.2/library",.libPaths()))

library(tidyverse)
library(readxl)
library(ggridges)
library(patchwork)
library(egg)
library(rstudioapi)
library(scales)
library(plotly)

`%!in%` <- Negate(`%in%`)
geom_mean <- function(x){exp(mean(log(x)))}
ci <- function(x){quantile(x, probs=c(0.025, 0.975), na.rm = T)}
################
#SETTING THE SCENARIO
seeds <- 50

endem <- "Low"
setting <- paste(endem, "complete", sep = "_")

#Set folder
source.dir <- dirname(getActiveDocumentContext()$path)
#"C:/Users/Z541213/Documents/Project/Model/Schisto_model"
setwd(source.dir)

#Load parameters
source("02_Parameters_Smansoni.R")

#Load population output
#####Load collated results and produce multi-panel plots
pop.output.dir <- file.path(source.dir, "Output/Population")
res <- readRDS(file.path(pop.output.dir, 
                         paste(setting, ".RDS", sep = "")))


if(endem=="High" | endem=="Moderate"){
  res <- res[- which(res$Immunity=="Absent" & res$Snails=="Absent" & res$DDF=="Absent"),]
}
#Computing faded-out runs at pre-control (149 ys is soon before MDA)
faided <- res %>%
  filter(time==149*12) %>%
  group_by(Immunity, Snails, DDF) %>%
  summarise(faided_seed = length(which(eggs_prev_SAC==0))*100/seeds)
n_faided <- length(which(faided$faided_seed>0))

#Average by seed
if(n_faided>0){
  res <- res[- which(res$time==149*12 & res$eggs_prev_SAC==0), ]
}

data_avg <- res %>%
  group_by(time, Immunity, Snails, DDF) %>%
  summarise(eggs_prev_SAC = mean(eggs_prev_SAC),
            eggs_prev_tot = mean(eggs_prev),
            PHI = mean(Heggs_prev),
            miracidiae = mean(miracidiae),
            cercariae = mean(cercariae),
            snail_inf = mean(inf_snail),
            snail_exp = mean(exp_snail),
            snail_prev = mean(inf_snail/(susc_snail+inf_snail+exp_snail)))

#Computing runs that reach EPHP or elimination:
ephp <- res %>%
  filter(time==(parms$mda$end+40)*12) %>%
  group_by(Snails, Immunity, DDF) %>%
  summarise(ephp_seed = length(which(Heggs_prev<=0.01))*100/seeds,
            eliminated = length(which(eggs_prev_SAC==0))*100/seeds) %>%
  ungroup() %>%
  mutate(x_coord = c(157, 159, rep(c(155, 157, 159), 8)),
         #x_coord = rep(c(155, 157, 159), 7), #low case
         y_coord = 3.5) #10 mod, 25 high

#Predictions
Fig <- ggplot(data=filter(data_avg, Immunity != "Strong"), aes(x=time/12)) +
  geom_line(aes(y=PHI*100, colour = DDF), size = 1) +
  geom_hline(yintercept = 1, linetype = "longdash") +
  geom_line(data = filter(res, Immunity != "Strong"), 
            aes(y = Heggs_prev*100, group = interaction(DDF, seed), colour = DDF), 
            alpha = 0.1) +
  labs(title = paste(endem, "endemicity setting", sep = " ")) +
  facet_grid(Snails ~ Immunity, labeller = labeller(.rows = label_both, .cols = label_both)) +
  # geom_text(data=ephp, aes(x=x_coord, y=y_coord, 
  #                          label=paste(ephp_seed, "%"), colour = DDF), 
  #           size=3, show.legend = FALSE) +
  # geom_text(data=ephp, aes(x=x_coord, y=y_coord-1, 
  #                          label=paste(eliminated, "%"), colour = DDF), 
  #           size=3, show.legend = FALSE) +
  # annotate("text", x = 157, y = ephp$y_coord[1]+0.5, label = "Probability of EPHP:", size =3) +
  # annotate("text", x = 157, y = ephp$y_coord[1]-1+0.5, label = "Probability of elimination:", size =3) +
  scale_y_continuous(name = "Prevalence of haevy intensity \n of infection in SAC (%) \n",
                     breaks = c(1, seq(0, 50, 5)),
                     limits = c(0, 30), #5, 15 mod, 30 high
                     expand = c(0, 0)) +
  scale_x_continuous(name = "\n Time",
                     breaks = c(145, seq(150, (parms$mda$end+40), 10)),
                     #limits = c(0, 1200),
                     labels = c(-15, seq(-10, 40, 10)),
                     expand = c(0, 0)) +
  coord_cartesian(xlim=c(145, (parms$mda$end+40))) +
  expand_limits(x = 0,y = 0) +
  theme_bw() +
  theme(legend.position="bottom",
        plot.margin = margin(5, 10, 0, 10, "pt"),
        panel.spacing = unit(1, "lines"))

tiff(paste("Plots/Predictions/", endem, ".tif", sep = ""), 
     width=12, height=9, units = "in", res = 300)
Fig 
dev.off()

####Main figure of observed pattern
#prev_endPC <- filter(data_avg, time == parms$mda$end*12)

Fig2 <- ggplot(data=data_avg, aes(x=time/12)) +
  geom_line(data = res, aes(y=eggs_prev_SAC*100,
                            group = interaction(DDF, seed),
                            colour = DDF), alpha = 0.01) +
  geom_line(aes(y=eggs_prev_SAC*100, colour = DDF), alpha = 3) +
  #geom_line(aes(y=PHI*100, colour = DDF), linetype = "longdash") +
  # geom_segment(data = prev_endPC,
  #              aes(x = 0, y = eggs_prev_SAC*100, xend = parms$mda$end, yend = eggs_prev_SAC*100, colour = DDF),
  #              linetype = "dashed") +
  labs(title = paste(endem, "endemicity setting", sep = " ")) +
  facet_grid(Snails ~ Immunity, labeller = labeller(.rows = label_both, .cols = label_both)) +
  # geom_text(data=text, aes(x=200, y=80, label=my_tag), 
  #           size=3) +
  scale_y_continuous(name = "Prevalence of infection in SAC (%) \n",
                     breaks = seq(0, 100, 5),
                     limits = c(0, 20),
                     expand = c(0, 0)) +
  scale_x_continuous(name = "\n Time [Years from last round of PC]",
                     breaks = seq(140, parms$mda$end+50, 10),
                     labels = seq(-20, 50, 10),
                     #limits = c(0, 1200),
                     expand = c(0, 0)) +
  coord_cartesian(xlim=c(140, parms$mda$end+50)) +
  expand_limits(x = 0,y = 0) +
  theme_bw() +
  theme(legend.position="bottom",
        #axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5),
        plot.margin = margin(5, 10, 0, 10, "pt"),
        panel.spacing = unit(1.5, "lines")) 

# ggsave(paste("Plots/Modeling_scenarios/", endem, ".tif", sep = ""), 
#        width=12, height=9, units = "in", res = 300)

tiff(paste("Plots/Modeling_scenarios/", endem, ".tif", sep = ""), 
     width=12, height=9, units = "in", res = 300)
Fig2 
dev.off()

#Age-intensity profiles
ind_data <- readRDS(file.path(source.dir, 
                              "Output/Individual/High_complete/Avg_individual_output.RDS"))

#Average over seeds
data_toplot <- age_out %>% #ind_data %>%
  #filter(!(Immunity == "Absent" & Snails== "Absent" & DDF== "Absent")) %>%
  group_by(age_group, Immunity, Snails, DDF) %>%
  summarise(epg_mean = mean(epg),
            epg_lo = ci(epg)[1],
            epg_hi = ci(epg)[2],
            wp_mean = mean(wp),
            wp_lo = ci(wp)[1],
            wp_hi = ci(wp)[2],
            dwp_mean = mean(dwp),
            dwp_lo = ci(dwp)[1],
            dwp_hi = ci(dwp)[2],
            rate_mean = mean(rate),
            rate_lo = ci(rate)[1],
            rate_hi = ci(rate)[2])

eggs <- 
  ggplot(filter(data_toplot, DDF=="Mild"), 
         aes(x=age_group, y=epg_mean+1, group = interaction(DDF, Immunity)))+
  geom_pointrange(aes(ymin=epg_lo+1, ymax=epg_hi+1, colour = DDF, shape = Immunity)) +
  geom_line(data = filter(data_toplot, DDF=="Mild" & Immunity == "Strong"),
            aes(colour = DDF, linetype = Immunity), linetype = "dashed", size = 1) +
  #labs(title = paste(endem, "endemicity setting", sep = " ")) +
  facet_grid( ~ Snails, labeller = labeller(.rows = label_both, .cols = label_both)) +
  scale_y_log10(name = "Observed egg counts (epg) + 1 \n",
                #breaks = seq(0, 10000, 200),
                limits = c(0.5, 8000), #3000 l-m 8000 high
                #expand = c(0, 0)
  ) +
  scale_x_discrete(name = "\n Age group [years]") +
  expand_limits(x = 0,y = 0) +
  guides(linetype = "none") +
  scale_color_manual(values = c("Mild" = hue_pal()(3)[2])) +
  theme_bw() +
  theme(legend.position="bottom",
        plot.margin = margin(5, 10, 0, 10, "pt"),
        panel.spacing = unit(1, "lines")) 

tiff(paste("Plots/Modeling_scenarios/", endem, "age-profile2.tif", sep = ""), width=10, height=6, units = "in", res = 300)
eggs
dev.off()
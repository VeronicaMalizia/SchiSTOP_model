#############################
#Author: Veronica Malizia
#R version: 4.1.2
#
#This script can plot and inspect results for single modelling scenario
#Prevalence plots and diagnostic checks
#############################
#Packages
.libPaths(c("C:/Program Files/R/R-4.1.2/library",.libPaths()))
library(tidyverse)
library(readxl)
library(ggridges)
library(patchwork)
`%!in%` <- Negate(`%in%`)
geom_mean <- function(x){exp(mean(log(x)))}
ci <- function(x){quantile(x, probs=c(0.025, 0.975), na.rm = T)}


#Average results over stochastic simulations
#Check if there are faded out runs (find a better implementation)
# dead.seeds <- unique(res$seed[which
#                               (res$true_prev==0 & res$time!=1)])

avg_res <- res %>%
  #filter(seed %!in% dead.seeds) %>%
  group_by(time, Immunity, Snails, DDF) %>%
  summarise(N = mean(pop_size),
            miracidiae = mean(miracidiae),
            cercariae = mean(cercariae),
            true_prev = mean(true_prev),
            eggs_prev = mean(eggs_prev),
            eggs_prev_SAC = mean(eggs_prev_SAC),
            Heggs_prev = mean(Heggs_prev),
            snail_prev = mean(inf_snail/(susc_snail+inf_snail+exp_snail)))


#Plot prevalence
#Fig1 <- 
  #filter(res, Immunity=="Absent" & DDF=="Absent") %>%
ggplot(res, aes(x=time/12)) +
  geom_line(aes(y=true_prev, group = seed), color = "grey20", alpha = 0.3) +
  geom_line(data=avg_res, aes(y=true_prev, color="True")) +
  geom_line(aes(y=eggs_prev, group = seed), color = "turquoise", alpha = 0.3) +
  geom_line(data=avg_res, aes(y=eggs_prev, color="KK-based")) +
  geom_line(aes(y=eggs_prev_SAC, group = seed), color = "light green", alpha = 0.3) +
  geom_line(data=avg_res, aes(y=eggs_prev_SAC, color="KK-based in SAC")) +
  geom_line(aes(y=Heggs_prev, group = seed), color = "coral", alpha = 0.3) +
  geom_line(data=avg_res, aes(y=Heggs_prev, color="High intensities")) +
  facet_grid(Snails ~ Immunity, labeller = labeller(.rows = label_both, .cols = label_both)) +
  scale_color_manual(name = "Prevalence:",
                     values = c("True" = "black",
                                "KK-based" = "dark blue",
                                "KK-based in SAC" = "dark green",
                                "High intensities" = "dark red")) +
  scale_y_continuous(name = "Prevalence (fraction)",
                     breaks = seq(0, 1, 0.2),
                     limits = c(0, 1),
                     expand = c(0, 0)) +
  scale_x_continuous(name = "Time [Years]",
                     #limits = c(0, 1200),
                     expand = c(0, 0)) +
  #coord_cartesian(xlim=c(500, (T*12))) +
  expand_limits(x = 0,y = 0)

Fig2 <- Fig1 +
  coord_cartesian(xlim=c(130, 150)) 


#Prevalence of infected snails
Fig3 <- ggplot(avg_res) +
  geom_line(aes(x=time, y=snail_prev, colour = DDF)) +
  facet_grid(Snails ~ Immunity, labeller = labeller(.rows = label_both, .cols = label_both)) +
  scale_y_continuous(name = "Prevalence of infected snails (fraction)",
                     breaks = seq(0, 1, 0.2),
                     limits = c(0, 1),
                     expand = c(0, 0)) +
  scale_x_continuous(name = "Time [Months]",
                     #limits = c(0, 1200),
                     expand = c(0, 0)) +
  expand_limits(x = 0,y = 0)

#ODE solutions
Fig4 <- ggplot(res, aes(x=time/12)) +
  geom_line(aes(y=inf_snail, group = seed), color = "red", alpha = 0.3) +
  geom_line(aes(y=susc_snail, group = seed), color = "green", alpha = 0.3) +
  geom_line(aes(y=exp_snail, group = seed), color = "orange", alpha = 0.3) +
  facet_grid(Snails ~ Immunity, labeller = labeller(.rows = label_both, .cols = label_both)) +
  scale_y_continuous(name = "Abundance of snails",
                     expand = c(0, 0)) +
  scale_x_continuous(name = "Time [Months]",
                     #limits = c(0, 1200),
                     expand = c(0, 0)) +
  expand_limits(x = 0,y = 0)

#Plot cercariae and miracidiae in the environmental reservoir
cerc <- ggplot(res) +
  geom_line(aes(x=time, y=cercariae, group = seed), color = "grey20", alpha = 0.3) +
  geom_line(data=avg_res, aes(x=time, y=cercariae), size=1) +
  scale_y_continuous(name = "Environmental reservoir (cercariae)",
                     #breaks = seq(0, 10000, 200),
                     #limits = c(0, 1500),
                     expand = c(0, 0)) +
  scale_x_continuous(name = "Time [years]",
                     limits = c(0, T*12),
                     expand = c(0, 0)) +
  expand_limits(x = 0,y = 0)

mirac <- ggplot(avg_res) +
  #geom_line(aes(x=time, y=miracidiae, group = seed), color = "grey20", alpha = 0.3) +
  geom_point(aes(x=time, y=miracidiae, colour = DDF), size=1) +
  facet_grid(Snails ~ Immunity, labeller = labeller(.rows = label_both, .cols = label_both)) +
  scale_y_continuous(name = "Total output of miracidiae",
                     #breaks = seq(0, 10000, 200),
                     #limits = c(0, 200000),
                     expand = c(0, 0)) +
  scale_x_continuous(name = "Time [years]",
                     #limits = c(500, 2400),
                     expand = c(0, 0)) +
  expand_limits(x = 0,y = 0)

cercVSmirac <- ggplot(avg_res) +
  geom_point(aes(x=miracidiae, y=cercariae, colour = DDF), size=2, alpha=0.5) +
  facet_grid(Snails ~ Immunity, labeller = labeller(.rows = label_both, .cols = label_both)) +
  scale_y_continuous(name = "Total output of cercariae",
                     #breaks = seq(0, 10000, 1000),
                     limits = c(0, 1000000),
                     expand = c(0, 0)) +
  scale_x_continuous(name = "Total output of miracidiae",
                     #limits = c(0, T*12),
                     expand = c(0, 0)) +
  expand_limits(x = 0,y = 0)

##Inspecting worms lifespans
w <- rgamma(100000, shape = 3, rate = 3/57)
plot(density(rgamma(100000, shape = 3, rate = 3/57)), lwd = 2.0, xlab = 'Age of worms [months]', main = "Worms' age distribution", xlim = c(0, 300))
quantile(w, probs = c(0.25, 0.5, 0.75, 0.95, 0.99, 1))
# 25%       50%       75%       95%       99%      100% 
#   34.52039  53.39290  78.43578 127.09506 168.28005 313.55901 
abline(v = c(119, 160), col = c('red', 'brown'), lwd = 2.0)
abline(v = 57, col = c('green'), lty = 5, lwd = 2.0)
legend("topright", 
       legend = c("Mean = 57 months", "95% perc = 10ys", "99% perc = 13ys"), 
       col = c("green", "red", "brown"), 
       pch = 20, 
       bty = "n", 
       pt.cex = 2, 
       cex = 1.2, 
       text.col = "black", 
       horiz = F , 
       inset = c(0.1, 0.1))



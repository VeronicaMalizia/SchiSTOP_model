#############################
#Author: Veronica Malizia
#Date: 27/07/2021
#R version: 3.6.1

#This script runs the model from R source, with human demography (aging/birth/deaths)
#and transmission dynamics with a central reservoir, consisting of snail population dynamcis.
#Code for plotting included.
#Data required: 
# 1. Initial age distribution, after equilibrium is reached in 00_Demography_eq.R
# 2. Death_probabilities from WHO App - South Africa
#############################
rm(list = ls())

.libPaths(c("C:/Program Files/R/R-4.1.2/library",.libPaths()))
library(tidyverse)
library(readxl)
library(ggridges)
library(foreach)
library(doParallel)
library(patchwork)
`%!in%` <- Negate(`%in%`)
geom_mean <- function(x){exp(mean(log(x)))}
ci <- function(x){quantile(x, probs=c(0.025, 0.975), na.rm = T)}
################
#Initial population
#Starting with the initial cohort
################
source.dir <- "C:/Users/Z541213/Documents/Project/Model/Schisto_model"
setwd(source.dir)
#Load age distribution at equilibrium
load("Equilibrium_age_distribution.RData") #the object is called "to.save"
#death_rates <- read.csv("death_rates_Uganda_2019.csv")
prob_death <- read.csv("prob_death_Uganda_2019.csv")

#0=male, 1=female
cohort <- c()
for(i in 1:nrow(to.save)){
  tmp <- tibble(age = rep(to.save$age[i], round(to.save$n[i])))
  cohort <- rbind(cohort, tmp) 
}
cohort <- cohort %>%
  mutate(sex = as.numeric(rbernoulli(nrow(cohort), 0.5)))

hist(cohort$age, main = "Initial age distribution", xlab = "Age (ys)")
table(cohort$sex)

# ggplot(cohort,aes(x=age,group=as.factor(sex),fill= as.factor(sex)))+
#   geom_histogram(position="dodge",binwidth=5)+
#   theme_bw()

#Functions (they can be a separate script)
age_groups <- c(0, 5, 12, 20, 100)
exposure_rates <- c(0.01, 0.61, 1, 0.12, 0.12) #Relative Age-specific exposure rates (activity/person/day) 
#exposure_rates <- c(0.33, 0.44, 0.22, 0) #Relative Age-specific exposure rates (activity/person/day) 
#Data from water contacts computed from Seydou S., De Vlas SJ, et al. (2011)
Age_profile_exp <- function(a){
  approx(x=age_groups, y=exposure_rates, xout=c(a), method = "constant")
}
plot(approx(x=age_groups, y=exposure_rates, method = "constant"), xlim = c(0, 100), ylim = c(0, 1), 
     type = 'l', xlab = "Age", ylab = "Relative exposure rate")
#lines(approx(x=c(0, 5, 10, 16, 100), y=c(0.01, 0.61, 1, 0.12, 0), method = "linear"), type = 'l', col='red')
Age_profile_contr <- function(a){
  approx(x=c(0, 10, 100), y=c(1, 1, 1), xout=c(a), method = "linear")
}
lines(approx(x=c(0, 10, 100), y=c(1, 1, 1), method = "linear"), col = "red")

FOIh <- function(l, zeta, v, a, is){
  Rel_ex <- Age_profile_exp(a)$y * is
  return(l * zeta * v * Rel_ex)
}
hyp_sat <- function(alpha, beta, w){ #Hyperbolic saturating function 
  f <- (alpha*w) / (1 + ((alpha*w) / beta))
  return(f)
}
exp_sat <- function(a, b, c, x){ #Exponential saturating function
  f <- a*(1-exp(-b*x))*(1-exp(-c*x)) 
  return(f)
}
SEI <- function(t, x, parms) {    
  with(as.list(c(parms, x)), {
    N <- S+E+I #is it better to work with densities?
    #Logistic growth
    #For population growing with limited amount of resources
    beta <- beta0*(1-N/k) #Infected snails do not reproduce
    FOIs <- l0*(1-chi^(mir/N))
    
    #Equations
    dS <- beta*(S+E) - (v+FOIs)*S #susceptible
    dE <- FOIs*S - (v+tau)*E #Exposed: snails are invaded, but larvae are not patent yet. Thus, snails do not shed cercariae
    dI <- tau*E - v2*I  #Infected: larvae in the snail are mature and snails shed cercariae
    dC <- lambda*I - m*C #Cercariae (output)
    res <- c(dS, dE, dI, dC)
    list(res)
  })
}

##Parameters
#monthly time step
birth_rate <- 36.5 #is crude annual birth rate Uganda 2019 (same y of available lifetables), 34.8 for Sub-Saharan Africa (per 1000 individuals)
#-1.09 is the net migration rate for 2022 for Uganda (per 1000 individuals)
emig_rate <- 20 #This is calibrated to have constant population
#max.pop <- 1000
k_w <- 0.15 #0.15 Anderson, Turner (2016) #can change for different settings (0.3 Sake) 
v <- 1 #Transmission probability
zeta <- 0.0003 #overall exposure rate. (0.42 water contacts rate per day, Seydou, De Vlas,.. 2011). (changing accordingly to endem. scenario)
ext.foi <- tibble(value = 1, #monthly
                  duration = 3) #years

Tw <- 60 #Average worm's lifespan in host in months (months)(40 m Sake) (5 years for Anderson and May 1985a)
phi1 <- 1-exp(-1/Tw) #(monthly) exponential dying probability of worms within the host

alpha <- 0.14 #expected number of eggs per sample per worm pair (Sake 1996)
gr_stool <- 150 #daily gr of stool produced by each human individual
z <- 0.0007
k_e <- 0.87 #aggregation parameter of egg counts detected (0.1 SCHISTOX; 0.87 Sake1992, but with three months interval and 25gr KK)
#co_rate <- 1 #Average contribution rate (monthly) #to include seasonal patterns
imm = 0.001 #immunity parameter (slope of the decrease in exposure)
mda <- tibble(age.lo = 5,
             age.hi = 15,
             start = 70, #70,
             end = 80, #80,
             frequency = 1, #annual
             coverage = 0.75,
             efficacy = 0.85)

##Parameters for ODEs model for snails (These rates are daily)
max.reproduction.rate = 1 #0.1 d^-1 from Civitello DJ, 2022 #monthly is ~ 1 egg/day 
carrying.capacity = 20000 #arbitrary. To be estimated. #Civitello uses 5 L^-1 (about 30 per m3) 
lifespan = 100 #days, Civitello #Gurarie: about 3 months
lifespan.infected = 30 #days, Gurarie
mortality.rate = 1/lifespan 
#lifespan.reduction=0.8 #arbitrary. Still not enough evidence found.
mortality.rate.infection = 1/lifespan.infected #0.04 d^-1 from Civitello. He works with additional mortality
mu = 30 #1 month: lifespan of larvae within the snail, before shedding cercariae
infection.rate = 1/mu 
rej.prob = 0.5 #probability of rejecting a miracidia, after getting in contact with the snail
# (1-chi)=0.5 for Civitello. OR it is for now computed from a Poisson as P(x=1)=0.8*exp(-0.8) using the infection rate from Anderson & May (1991).
# However, I would consider it arbitrary too and then to be estimated. (Or look for data)
max.invasion = 1/15 #1/d. Rate of sporocyst development in snails, given successful invasion.
cerc.prod.rate = 50 #1/d per infected snail
cerc.mortality = 1 #1/d 

# #Option 1.
# #Population is initialized with a given worms distribution
# #(Initial/external Force of Infection)
# w0= rnbinom(nrow(cohort), size=k_w, mu=6)
# #Compute initial worms' pair
# mw0 <- rbinom(length(w0), w0, 0.5) #male worms
# fw0 <- w0-mw0 #female worms
# cohort <- cohort %>%
#   mutate(wp = pmin(mw0, fw0),
#          Ind_sus = rgamma(nrow(cohort), shape = k_w, scale = 1/k_w))
# 
# ggplot(cohort,aes(x=wp))+
#   geom_histogram(binwidth=5)+ theme_bw()
# hist(cohort$Ind_sus)
# hist(w0)

#Option 2.
# Worms are initialized with an artificial FOI (1 worm pair per person per month)
cohort <- cohort %>%
  mutate(jw1 = 1, 
         Ind_sus = rgamma(nrow(cohort), shape = k_w, scale = 1/k_w))

#Create initial cloud
# #For option 1.
# eggs0 <- alpha*(pmin(mw0, fw0)) #*exp(-z*fw0)
# contributions0 <- eggs0 * Age_profile(cohort$age)$y * cohort$Ind_sus
#For option 2.
eggs0 <- 0
contributions0 <- 0

#Initial cumulative exposure
cum_exp <- sum(Age_profile_exp(cohort$age)$y * cohort$Ind_sus) 
SAC <- which(cohort$age >= 5 & cohort$age <= 15)

#Snail population is initialized
snail.pop=10000
E0=0
I0=0
S0=snail.pop - sum(E0, I0)
C0=0
#If densities: S0=1

T <- 200 #number of years
seeds <- 10

#Run the model
lim_mechanism <- c("worms") #choices: 'worms' (for DDF), 
                             #'snails' (for vect. saturation), 
                             #'humans' (for immunity), 
                             #'no' (none)

time.start <- Sys.time()
source("C:\\Users\\Z541213\\Documents\\Project\\Model\\Schisto_model\\01.a_Model_function.R")
time.end <- Sys.time()
time.end - time.start

#Population-level results
#Collating results (only if seeds>1)
res <- bind_rows(results)

#Average results over stochastic simulations

#Check if there are faded out runs
# dead.seeds <- unique(res$seed[which
#                               (res$true_prev==0 & res$time!=1)])

avg_res <- res %>%
  #filter(seed %!in% dead.seeds) %>%
  group_by(time) %>%
  summarise(N = mean(pop_size),
            miracidiae = mean(miracidiae),
            cercariae = mean(cercariae),
            true_prev = mean(true_prev),
            eggs_prev = mean(eggs_prev),
            eggs_prev_SAC = mean(eggs_prev_SAC),
            Heggs_prev = mean(Heggs_prev),
            snail_prev = mean(inf_snail/tot_snail))


#Plot prevalence
Fig <- ggplot(res) +
  geom_line(aes(x=time, y=true_prev, group = seed), color = "grey20", alpha = 0.3) +
  geom_line(data=avg_res, aes(x=time, y=true_prev, color="True")) +
  geom_line(aes(x=time, y=eggs_prev, group = seed), color = "turquoise", alpha = 0.3) +
  geom_line(data=avg_res, aes(x=time, y=eggs_prev, color="KK-based")) +
  geom_line(aes(x=time, y=eggs_prev_SAC, group = seed), color = "light green", alpha = 0.3) +
  geom_line(data=avg_res, aes(x=time, y=eggs_prev_SAC, color="KK-based in SAC")) +
  geom_line(aes(x=time, y=Heggs_prev, group = seed), color = "coral", alpha = 0.3) +
  geom_line(data=avg_res, aes(x=time, y=Heggs_prev, color="High intensities")) +
  scale_color_manual(name = "Prevalence:",
                     values = c("True" = "black",
                                "KK-based" = "dark blue",
                                "KK-based in SAC" = "dark green",
                                "High intensities" = "dark red")) +
  scale_y_continuous(name = "Prevalence (fraction)",
                     breaks = seq(0, 1, 0.2),
                     limits = c(0, 1),
                     expand = c(0, 0)) +
  scale_x_continuous(name = "Time [Months]",
                     #limits = c(0, 1200),
                     expand = c(0, 0)) +
  #coord_cartesian(xlim=c(500, (T*12))) +
  expand_limits(x = 0,y = 0)

Fig

Fig +
  coord_cartesian(xlim=c(500, 1000)) 

#Saving image
tiff("Plots/Prevalence_humans.tif", width=7, height=6, units = "in", res = 300)
Fig
dev.off()

#PLot prevalence of infected patent snails
tiff("Plots/Prevalence_snail_low.tif", width=7, height=6, units = "in", res = 300)
ggplot(avg_res) +
  geom_line(aes(x=time, y=snail_prev)) +
  scale_y_continuous(name = "Prevalence of infected snails (fraction)",
                     breaks = seq(0, 1, 0.2),
                     limits = c(0, 1),
                     expand = c(0, 0)) +
  scale_x_continuous(name = "Time [Months]",
                     #limits = c(0, 1200),
                     expand = c(0, 0)) +
  expand_limits(x = 0,y = 0)
dev.off()

# ggplot(res) +
#   geom_line(aes(x=time, y=inf_snail, group = seed), color = "purple", alpha = 0.3) +
#   geom_line(aes(x=time, y=tot_snail, group = seed), color = "brown", alpha = 0.3) +
#   scale_y_continuous(name = "Abundance of snails",
#                      expand = c(0, 0)) +
#   scale_x_continuous(name = "Time [Months]",
#                      #limits = c(0, 1200),
#                      expand = c(0, 0)) +
#   expand_limits(x = 0,y = 0)

#Plot cercariae and miracidiae in the environmental reservoir
tiff(file.path(source.dir, "/Plots/Cercariae_moderate.tif"), width=7, height=6, units = "in", res = 300)
ggplot(res) +
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
dev.off()

tiff(file.path(source.dir, "/Plots/Miracidiae_low.tif"), width=7, height=6, units = "in", res = 300)
ggplot(res) +
  geom_line(aes(x=time, y=miracidiae, group = seed), color = "grey20", alpha = 0.3) +
  geom_line(data=avg_res, aes(x=time, y=miracidiae), size=1) +
  scale_y_continuous(name = "Total output of miracidiae",
                     #breaks = seq(0, 10000, 200),
                     #limits = c(0, 200000),
                     expand = c(0, 0)) +
  scale_x_continuous(name = "Time [years]",
                     limits = c(0, T*12),
                     expand = c(0, 0)) +
  expand_limits(x = 0,y = 0)
dev.off()

ggplot(avg_res) +
  geom_point(aes(x=miracidiae, y=cercariae), size=2, alpha=0.5) +
  scale_y_continuous(name = "Total output of cercariae",
                     #breaks = seq(0, 10000, 1000),
                     #limits = c(0, 10000),
                     expand = c(0, 0)) +
  scale_x_continuous(name = "Total output of miracidiae",
                     #limits = c(0, T*12),
                     expand = c(0, 0)) +
  expand_limits(x = 0,y = 0)

#### Individual-level results
#Collate results by seeds
data_all <- list.files(path = file.path(source.dir, "/Output/"),  # Identify all output CSV files
                       pattern = "Ind_out_seed_*", full.names = TRUE) %>% 
            lapply(read_csv, show_col_types = F) %>%              # Store all files in list
            bind_rows                                             # Combine data sets into one
#Aggregate by age groups
#Filter one year of interest
unique(data_all$time) #individual output contains time in years
years <- c(51, 81, 101, 181)
age_out <- data_all %>%
  filter(time %in% years) %>%
  #filter(seed %!in% dead.seeds) %>%
  mutate(age_group = case_when(age<=1 ~ "0_1",
                               age>1 & age<=5 ~ "1_5",
                               age>5 & age<=15 ~ "5_15",
                               age>15 & age<=30 ~ "15_30",
                               age>30 & age<=50 ~ "30_50",
                               age>50 ~ "50_90")) %>% 
  group_by(seed, time, age_group, sex) %>%
  summarise(epg = geom_mean(ec*24+1), #means over individuals of that sex-age group
            wp = mean(wp),
            dwp = mean(cum_dwp), 
            rate = mean(rate))
  
age_out$age_group <- factor(age_out$age_group, levels = c("0_1",
                                                          "1_5",
                                                          "5_15",
                                                          "15_30",
                                                          "30_50",
                                                          "50_90"))
#Eggs by age
eggs <- ggplot(age_out)+
  geom_boxplot(aes(x=age_group, y=epg-1), alpha = 0.3, outlier.shape = NA) +
  facet_wrap(~ time, ncol = 4) +
  scale_y_continuous(name = "Observed egg counts (epg) \n geometric mean",
                     #breaks = seq(0, 10000, 200),
                     #limits = c(0, 200000),
                     expand = c(0, 0)) +
  scale_x_discrete(name = " ") +
  expand_limits(x = 0,y = 0) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

#Intensity by age
worms <- ggplot(age_out)+
  geom_boxplot(aes(x=age_group, y=wp), alpha = 0.3, outlier.shape = NA) +
  facet_wrap(~ time, ncol = 4) +
  scale_y_continuous(name = "Average worm load (pairs) \n",
                     #breaks = seq(0, 10000, 200),
                     #limits = c(0, 200000),
                     expand = c(0, 0)) +
  scale_x_discrete(name = " ") +
  expand_limits(x = 0,y = 0) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

#Establishment rate by age
rate <- ggplot(age_out)+
  geom_boxplot(aes(x=age_group, y=rate), alpha = 0.3, outlier.shape = NA) +
  facet_wrap(~ time, ncol = 4) + #scales = "free_y", 
  scale_y_continuous(name = "Parasite establishment rate \n",
                     #breaks = seq(0, 10000, 200),
                     #limits = c(0, 200000),
                     expand = c(0, 0)) +
  scale_x_discrete(name = "\n Age [years]") +
  expand_limits(x = 0,y = 0) +
  #coord_cartesian(ylim = c(0, 200000)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

tiff(file.path(source.dir, "/Plots/Boxplots_humans.tif"), width=9, height=8, units = "in", res = 300)
eggs / worms / rate + 
  plot_annotation(title = paste("Limiting mechanism on", lim_mechanism, sep = " "))
dev.off()

ggplot(age_out)+
  geom_boxplot(aes(x=age_group, y=dwp), alpha = 0.3, outlier.shape = NA) +
  facet_wrap(~ time, ncol = 4) +
  scale_y_continuous(name = "Cumulative dead worms (pairs) \n",
                     #breaks = seq(0, 10000, 200),
                     #limits = c(0, 200000),
                     expand = c(0, 0)) +
  scale_x_discrete(name = " ") +
  expand_limits(x = 0,y = 0) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

##Print average worm pairs burden per person
age_out %>% 
  filter(time==181) %>%
  group_by(time) %>%
  summarise(avg_wp = mean(wp))

##PLOT demographics

#Plot population size
# tiff(file.path(source.dir, "/Plots/Pop_size_migration.tif"), width=7, height=6, units = "in", res = 300)
# 
# ggplot(res) +
#   geom_line(aes(x=time, y=pop_size, group = seed), color = "grey20", alpha = 0.3) +
#   geom_line(data=avg_res, aes(x=time, y=N), size=1) +
#   scale_y_continuous(name = "Population size (counts)",
#                      breaks = seq(0, 10000, 500),
#                      limits = c(0, 2000),
#                      expand = c(0, 0)) +
#   scale_x_continuous(name = "Time [months]",
#                      #breaks = seq(0, 10000, 500),
#                      limits = c(0, (T*12)+12),
#                      expand = c(0, 0)) +
#   expand_limits(x = 0,y = 0)
# 
# dev.off()
# 
# ##Age distribution by time
# #https://r-charts.com/distribution/ggridges/
# cohort_plot_age <- cohort %>%
#   select(age, sex, wp) %>%
#   mutate(ec = NA,
#          ID = 1:nrow(cohort),
#          time = 0,
#          seed = NA)
# bind_rows(cohort_plot_age,
#   data_all %>%
#   filter(time %in% c(seq(12, (T*12), 20*12)/12))) %>%
#   ggplot(aes(x=round(age))) +
#   geom_histogram(aes(y = ..density.., colour=as.factor(seed)), binwidth = 5, 
#                  fill="white", position = 'dodge') +
#   #stat_bin(aes(y=..count.., label=..count..), binwidth = 5, geom="text", vjust=-.5) +
#   guides(colour=guide_legend(title="Simulation seed")) +
#   scale_color_grey() +
#   facet_wrap(~ as.factor(time)) +
#   scale_y_continuous(name = "Density",
#                      expand = c(0, 0)) +
#   scale_x_continuous(name = "Age [years]",
#                      expand = c(0, 0)) +
#   expand_limits(x = 0,y = 0) +
#   theme(legend.position="top")
# 
# 
# data_all %>%
#    filter(time==1, seed==1) %>%
#    ggplot(aes(x=age)) +
#    geom_histogram(binwidth = 5) + 
#    stat_bin(aes(y=..count.., label=..count..), binwidth = 5, geom="text", vjust=-.5) 

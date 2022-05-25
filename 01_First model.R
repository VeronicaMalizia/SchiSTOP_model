#############################
#Author: Veronica Malizia
#Date: 27/12/2021
#R version: 3.6.1

#This script runs a toy model from R source, with human demography (aging/birth/deaths)
#and basic transmission dynamics with a central reservoir.
#Code for plotting included.
#Data required: 
# 1. Initial cohort population from WORMSIM input file
# 2. Death_probabilities from WHO App - South Africa
#############################
rm(list = ls())

.libPaths(c("C:/Program Files/R/R-4.1.2/library",.libPaths()))
source.dir <- "C:/Users/Z541213/Documents/Project/Model/Schisto_model"
setwd(source.dir)
library(tidyverse)
library(readxl)
library(ggridges)
library(foreach)
library(doParallel)

################
#Initial population
#Starting with the initial cohort
################

pop0 <- read_xlsx("Initial population_Uganda.xlsx")
death_rates <- read.csv("death_rates_Uganda_2019.csv")
prob_death <- read.csv("prob_death_Uganda_2019.csv")
le <- 60 #Life expectancy (Uganda)
sex.ratio <- 94.6 #number of males per 100 females. Uganda, 2014.
N <- 1000 #population size
N_males <- round((sex.ratio/(sex.ratio + 100))*N)
N_females <- N - N_males

cohort <- tibble(age = round(c(runif(round(pop0$N_males[1]*N_males), 0, pop0$Age_limit[1]),
                               runif(round(pop0$N_females[1]*N_females), 0, pop0$Age_limit[1]))),
                 sex = c(rep(0, round(pop0$N_males[1]*N_males)), 
                         rep(1, round(pop0$N_females[1]*N_females))))
#0=male, 1=female
for(i in 2:nrow(pop0)){
  males <- round(pop0$N_males[i]*N_males)
  females <- round(pop0$N_females[i]*N_females)
  tmp <- tibble(age = round(c(runif(males, pop0$Age_limit[i-1]+1, pop0$Age_limit[i]),
                              runif(females, pop0$Age_limit[i-1]+1, pop0$Age_limit[i]))),
                sex = c(rep(0, males), rep(1, females)))
  cohort <- rbind(cohort, tmp) 
}

hist(cohort$age)
table(cohort$sex)

# ggplot(cohort,aes(x=age,group=as.factor(sex),fill= as.factor(sex)))+
#   geom_histogram(position="dodge",binwidth=5)+
#   theme_bw()

##Parameters
T <- 100 #number of years
seeds <- 10
#monthly time step
birth_rate <- 35.7 #35.7 is crude annual birth rate Uganda 2021, 34.8 for Sub-Saharan Africa (per 1000 individuals)
#-1.09 is the net migration rate for 2022 for Uganda (per 1000 individuals)
emig_rate <- 10 #1.09 #This approximates the emigration rate from rural villages (2021) (per 1000 individuals)
#max.pop <- 1000
k_w <- 0.15 #Anderson, Turner (2016)
v <- 1 #Transmission probability
zeta <- 0.9 #overall exposure rate (check the monthly definition)
Tw <- 5 #Average worm's lifespan in host in months (years)(Anderson and May 1985a)
phi1 <- 1-exp(-1/(Tw*12)) #(monthly) dying probability of worms within the host
#phi2 <- exp(-1/(xx*12)) #survival probability of particles in the reservoir
alpha <- 0.28 #expected number of eggs per sample per worm pair (Sake 1996)
z <- 0.0007
k_e <- 0.22 #aggregation parameter of egg counts detected (0.1 SCHISTOX; 0.22 Sake1992)
co_rate <- 1 #Average contribution rate (monthly) #to include seasonal patterns
#expanding_factor <- 1 #Multiplicative factor of miracidia in snails. (i.e. N.of cercariae released)
mda <- tibble(age.lo = 5,
             age.hi = 15,
             start = 100,
             end = 130,
             frequency = 1, #annual
             coverage = 0.75,
             efficacy = 1)

##Parameters for ODEs model for snails (These rates are daily)
max.reproduction.rate = 0.1 #d^-1 from Civitello DJ, 2022 #monthly is ~ 1 egg/day 
carrying.capacity = 1000 #arbitrary. To be estimated. #Civitello uses 5 L^-1 (about 30 per m3) 
lifespan = 100 #days, Civitello #Gurarie: about 3 months
mortality.rate = 1/lifespan 
#lifespan.reduction=0.8 #arbitrary. Still not enough evidence found.
mortality.rate.infection = 0.04 #d^-1 from Civitello #1/((1-lifespan.reduction)*lifespan)
mu = 30 #1 month: lifespan of larvae within the snail, before shedding cercariae
infection.rate = 1/mu 
#c = 2 #contact rate for snails. This should also be a calibrating parameter. 
#For now we assume that snails get in contact with all free-living miracidiae. I think it cannot be >1
chi = 0.5 #probability of a successful invasion for a single miracidia getting in contact with the host
# 0.5 is Civitello. OR it is for now computed from a Poisson as P(x=1)=0.8*exp(-0.8) using the infection rate from Anderson & May (1991).
# However, I would consider it arbitrary too and then to be estimated. (Or look for data)
l0 = 1/15 #1/d. Rate of sporocyst development in snails, given successful invasion.
cerc.prod.rate = 50 #1/d per infected snail
cerc.mortality = 1 #1/d 

#Population is initialized with a given worms distribution
#(Initial/external Force of Infection)
w0= rnbinom(nrow(cohort), size=k_w, mu=6)
#Compute initial worms' pair
mw0 <- rbinom(length(w0), w0, 0.5) #male worms
fw0 <- w0-mw0 #female worms
cohort <- cohort %>%
  mutate(wp = pmin(mw0, fw0),
         Ind_sus = rgamma(nrow(cohort), shape = k_w, scale = 1/k_w))

ggplot(cohort,aes(x=wp))+
  geom_histogram(binwidth=5)+ theme_bw()
hist(cohort$Ind_sus)
hist(w0)

#Functions (they can be a separate script)
age_groups <- c(0, 5, 10, 16, 200)
exposure_rates <- c(0.032, 0.61, 1, 0.06, 0.06) #These should reflect a moderate adult burden setting
Age_profile <- function(a){
  approx(x=age_groups, y=exposure_rates, xout=c(a), method = "constant")
}
plot(approx(x=age_groups, y=exposure_rates, method = "constant"), xlim = c(0, 80), type = 'l', xlab = "Age", ylab = "Relative exposure rate")
lines(approx(x=age_groups, y=c(0.01, 0.61, 1, 0.12, 0.12), method = "constant"), xlim = c(0, 80), type = 'l', col='red')
foi <- function(l, zeta, v, a, is){
  Rel_ex <- Age_profile(a)$y * is
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
    N <- S+E+I #is better to work with densities?
    #Logistic growth
    #For population growing with limited amount of resources
    beta <- beta0*(1-N/k) #Infected snails do not reproduce
    FOIsnails <- FOIs/N
    
    #Equations
    dS <- beta*(S+E) - (v+FOIsnails)*S #susceptible
    dE <- FOIsnails*S - (v+tau)*E #Exposed: snails are invaded, but larvae are not patent yet. Thus, snails do not shed cercariae
    dI <- tau*E - (v+v2)*I  #Infected: larvae in the snail are mature and snails shed cercariae
    dC <- lambda*I - m*C #Cercariae (output)
    res <- c(dS, dE, dI, dC)
    list(res)
  })
}

#Create initial cloud
eggs0 <- alpha*(pmin(mw0, fw0)) #*exp(-z*fw0)
contributions0 <- eggs0 * Age_profile(cohort$age)$y * cohort$Ind_sus
#Initial cumulative exposure
cum_exp <- sum(Age_profile(cohort$age)$y * cohort$Ind_sus) 
SAC <- which(cohort$age >= 5 & cohort$age <= 15)

#Snail population is initialized
snail.pop=10000
E0=0
I0=0
S0=snail.pop - sum(E0, I0)
C0=0
#If densities: S0=1

#Run the model
source("C:\\Users\\Z541213\\Documents\\Project\\Model\\Schisto_model\\01.a_Model_function.R")

#Population-level results
#Collating results (only if seeds>1)
res <- bind_rows(results)

#Average results over simulations
avg_res <- res %>%
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

#Saving image
tiff("Plots/Prevalence_SEIversion.tif", width=7, height=6, units = "in", res = 300)
Fig
dev.off()


#PLot prevalence of infected patent snails
tiff("Plots/Prevalence_snail_SEIversion.tif", width=7, height=6, units = "in", res = 300)
ggplot(res) +
  geom_line(data=avg_res, aes(x=time, y=snail_prev)) +
  scale_y_continuous(name = "Prevalence of infected snails (fraction)",
                     breaks = seq(0, 1, 0.2),
                     limits = c(0, 1),
                     expand = c(0, 0)) +
  scale_x_continuous(name = "Time [Months]",
                     #limits = c(0, 1200),
                     expand = c(0, 0)) +
  expand_limits(x = 0,y = 0)
dev.off()

ggplot(res) +
  geom_line(aes(x=time, y=inf_snail, group = seed), color = "purple", alpha = 0.3) +
  geom_line(aes(x=time, y=tot_snail, group = seed), color = "brown", alpha = 0.3) +
  scale_y_continuous(name = "Abundance of snails",
                     expand = c(0, 0)) +
  scale_x_continuous(name = "Time [Months]",
                     #limits = c(0, 1200),
                     expand = c(0, 0)) +
  expand_limits(x = 0,y = 0)

#Plot population size
tiff(file.path(source.dir, "/Plots/Pop_size_migration.tif"), width=7, height=6, units = "in", res = 300)

ggplot(res) +
  geom_line(aes(x=time, y=pop_size, group = seed), color = "grey20", alpha = 0.3) +
  geom_line(data=avg_res, aes(x=time, y=N), size=1) +
  scale_y_continuous(name = "Population size (counts)",
                     breaks = seq(0, 10000, 1000),
                     #limits = c(0, 1500),
                     expand = c(0, 0)) +
  scale_x_continuous(name = "Time [months]",
                     #breaks = seq(0, 10000, 500),
                     limits = c(0, (T*12)+12),
                     expand = c(0, 0)) +
  expand_limits(x = 0,y = 0)

dev.off()

#Plot cercariae and miracidiae in the environmental reservoir
ggplot(res) +
  geom_line(aes(x=time, y=cercariae, group = seed), color = "grey20", alpha = 0.3) +
  geom_line(data=avg_res, aes(x=time, y=cercariae), size=1) +
  scale_y_continuous(name = "Environmental reservoir (particles)",
                     #breaks = seq(0, 10000, 200),
                     #limits = c(0, 1500),
                     expand = c(0, 0)) +
  scale_x_continuous(name = "Time [years]",
                     limits = c(0, T*12),
                     expand = c(0, 0)) +
  expand_limits(x = 0,y = 0)

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

ggplot(avg_res) +
  geom_line(aes(x=miracidiae, y=cercariae), size=1) +
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
year <- 51
age_out <- data_all %>%
  filter(time==year) %>%
  mutate(age_group = case_when(age<=1 ~ "0_1",
                               age>1 & age<=5 ~ "1_5",
                               age>5 & age<=15 ~ "5_15",
                               age>15 & age<=30 ~ "15_30",
                               age>30 & age<=50 ~ "30_50",
                               age>50 ~ "50_90")) 
  group_by(seed, age_group, sex) %>%
  summarise(ec = mean(ec*24), #means over individuals of that sex-age group
            wp = mean(wp))
  
age_out$age_group <- factor(age_out$age_group, levels = c("0_1",
                                                          "1_5",
                                                          "5_15",
                                                          "15_30",
                                                          "30_50",
                                                          "50_90"))
#Eggs by age
ggplot(age_out)+
  geom_boxplot(aes(x=age_group, y=ec), alpha = 0.3) +
  facet_wrap(~ sex) +
  scale_y_continuous(name = "Observed egg counts (epg)",
                     #breaks = seq(0, 10000, 200),
                     #limits = c(0, 200000),
                     expand = c(0, 0)) +
  scale_x_discrete(name = "Age [years]") +
  expand_limits(x = 0,y = 0)

#Intensity by age
ggplot(age_out)+
  geom_boxplot(aes(x=age_group, y=wp), alpha = 0.3) +
  facet_wrap(~ sex) +
  scale_y_continuous(name = "Average worm load (pairs)",
                     #breaks = seq(0, 10000, 200),
                     #limits = c(0, 200000),
                     expand = c(0, 0)) +
  scale_x_discrete(name = "Age [years]") +
  expand_limits(x = 0,y = 0)

##Age distribution by time
#https://r-charts.com/distribution/ggridges/
cohort_plot_age <- cohort %>%
  select(age, sex, wp) %>%
  mutate(ec = NA,
         ID = 1:nrow(cohort),
         time = 0,
         seed = NA)
bind_rows(cohort_plot_age,
  data_all %>%
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


data_all %>%
   filter(time==1, seed==1) %>%
   ggplot(aes(x=age)) +
   geom_histogram(binwidth = 5) + 
   stat_bin(aes(y=..count.., label=..count..), binwidth = 5, geom="text", vjust=-.5) 

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

library(tidyverse)
library(readxl)
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
emig_rate <- 10 #This approximates the emigration rate from rural villages (2021) (per 1000 individuals)
#max.pop <- 1000
k_w <- 0.15 #Anderson, Turner (2016)
v <- 1 #Transmission probability
zeta <- 0.9 #overall exposure rate
Tw <- 5 #Average worm's lifespan in host in months (years)(Anderson and May 1985a)
phi1 <- 1-exp(-1/(Tw*12)) #(monthly) dying probability of worms within the host
#phi2 <- exp(-1/(xx*12)) #survival probability of particles in the reservoir
alpha <- 0.28 #expected number of eggs per sample per worm pair (Sake 1996)
z <- 0.0007
k_e <- 0.22 #aggregation parameter of egg counts detected (0.1 SCHISTOX; 0.22 Sake1992)
a <- 300 #a, b parameter of the hyperbolic saturating function for cercariae(saturating of available snails)
b <- 3000 #maximum output of cercariae in the environment (it may be difficult to estimate)
# c <- 0.0861
co_rate <- 1 #Average contribution rate (monthly) #to include seasonal patterns
#expanding_factor <- 1 #Multiplicative factor of miracidia in snails. (i.e. N.of cercariae released)
mda <- tibble(age.lo = 5,
             age.hi = 15,
             start = 100,
             end = 130,
             frequency = 1, #annual
             coverage = 0.75,
             efficacy = 1)

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

#Functions
age_groups <- c(0, 5, 10, 16, 200)
exposure_rates <- c(0.032, 0.61, 1, 0.06, 0.06)
Age_profile <- function(a){
  approx(x=age_groups, y=exposure_rates, xout=c(a), method = "constant")
}
plot(approx(x=age_groups, y=exposure_rates, method = "constant"), type = 'l')
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


#Create initial cloud
eggs0 <- alpha*(pmin(mw0, fw0)) #*exp(-z*fw0)
contributions0 <- eggs0 * Age_profile(cohort$age)$y * cohort$Ind_sus
#Initial cumulative exposure
cum_exp <- sum(Age_profile(cohort$age)$y * cohort$Ind_sus) 
SAC <- which(cohort$age >= 5 & cohort$age <= 15)

#Run the model
source("C:\\Users\\Z541213\\Documents\\Project\\Model\\Schisto_model\\01.a_Model_function.R")

#Collating results (only if seeds>1)
res <- bind_rows(results)

#Average results over simulations
avg_res <- res %>%
  group_by(time) %>%
  summarise(N = mean(pop_size),
            miracidiae = mean(miracidiae),
            reservoir = mean(reservoir),
            true_prev = mean(true_prev),
            eggs_prev = mean(eggs_prev),
            eggs_prev_SAC = mean(eggs_prev_SAC),
            Heggs_prev = mean(Heggs_prev))


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
tiff("Prevalence_densdep_V03.tif", width=7, height=6, units = "in", res = 300)
Fig
dev.off()

#Plot population size
ggplot(res) +
  geom_line(aes(x=time, y=pop_size, group = seed), color = "grey20", alpha = 0.3) +
  geom_line(data=avg_res, aes(x=time, y=N), size=1) +
  scale_y_continuous(name = "Population size (counts)",
                     breaks = seq(0, 10000, 500),
                     #limits = c(0, 1500),
                     expand = c(0, 0)) +
  scale_x_continuous(name = "Time [years]",
                     limits = c(0, T*12),
                     expand = c(0, 0)) +
  expand_limits(x = 0,y = 0)

#Plot particles in the environmental reservoir
ggplot(res) +
  geom_line(aes(x=time, y=reservoir, group = seed), color = "grey20", alpha = 0.3) +
  geom_line(data=avg_res, aes(x=time, y=reservoir), size=1) +
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

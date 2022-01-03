#############################
#Author: Veronica Malizia
#Date: 27/10/2021
#R version: 3.6.1

#This will generate a first toy model, without human demography (no age/birth/deaths)
#and basic transmission dynamics with a central reservoir.
#Data: Initial cohort population from WORMSIM input file
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
pop0 <- read_xlsx("Initial population.xlsx")
death_rates <- read.csv("death_rates_Uganda_2019.csv")
prob_death <- read.csv("prob_death_Uganda_2019.csv")
le <- 60 #Life expectancy (Uganda)
#compute death rate from surv_table

cohort <- tibble(age = round(c(runif(pop0$N_males[1], 0, pop0$Age_limit[1]),
                               runif(pop0$N_females[1], 0, pop0$Age_limit[1]))),
                 sex = c(rep(0, pop0$N_males[1]), rep(1, pop0$N_females[1])))
#0=male, 1=female
for(i in 2:nrow(pop0)){
  tmp <- tibble(age = round(c(runif(pop0$N_males[i], pop0$Age_limit[i-1]+1, pop0$Age_limit[i]),
                              runif(pop0$N_females[i], pop0$Age_limit[i-1]+1, pop0$Age_limit[i]))),
                sex = c(rep(0, pop0$N_males[i]), rep(1, pop0$N_females[i])))
  cohort <- rbind(cohort, tmp) 
}

hist(cohort$age)

ggplot(cohort,aes(x=age,group=as.factor(sex),fill= as.factor(sex)))+
  geom_histogram(position="dodge",binwidth=5)+
  theme_bw()

##Parameters
T <- 100 #number of years
#monthly time step
birth_rate <- 37 #Crude birth rate Uganda (per 1000 individuals)
k_w <- 0.24 #Anderson, Turner (2016)
v <- 1 #Transmission probability
zeta <- 0.1 #overall exposure rate
Tw <- 5.7 #Average worm's lifespan in host (years)(Anderson and May 1985a)
phi1 <- exp(-1/(Tw*12)) #(monthly) dyeing probability of worms within the host
#phi2 <- exp(-1/(xx*12)) #survival probability of particles in the reservoir
alpha <- 0.24 #expected number of eggs per sample per worm pair (Sake 1996)
z <- 0.0007
k_e <- 0.22 #aggregation parameter of egg counts detected (0.1 SCHISTOX; 0.22 Sake1992)
# a <- 1.2 #a, b, c, parameter of the exp saturating function for the material the clous (saturating of available snails)
# b <- 0.0213
# c <- 0.0861
co_rate <- 1 #Average contribution rate (monthly) #miracidia infecting snails per month
expanding_factor <- 1 #Multiplicative factor of miracidia in snails. (i.e. N.of cercariae released)
dens_dep <- F
mda <- tibble(age.lo = 5,
             age.hi = 15,
             start = 0, #50
             end = 0, #70
             frequency = 1, #annual
             coverage = 0.75,
             efficacy = 0.80)

#Population is initialized with a given worms distribution
#(Initial/external Force of Infection)
cohort <- cohort %>%
  mutate(w = rnbinom(nrow(cohort), size=k_w, mu=6),
         Ind_sus = rgamma(nrow(cohort), shape = k_w, scale = 1/k_w))

ggplot(cohort,aes(x=w))+
  geom_histogram(binwidth=5)+ theme_bw()
hist(cohort$Ind_sus)

#Functions
age_groups <- c(0, 5, 10, 16, 150)
exposure_rates <- c(0.032, 0.61, 1, 0.06, 0.06)
Age_profile <- function(a){
  approx(x=age_groups, y=exposure_rates, xout=c(a), method = "constant")
}
plot(approx(x=age_groups, y=exposure_rates, method = "constant"), type = 'l')
foi <- function(l, zeta, v, a, is){
  Rel_ex <- Age_profile(a)$y * is
  return(l * zeta * v * Rel_ex)
}
cum_exp <- sum(Age_profile(cohort$age)$y*cohort$Ind_sus) #Cumulative exposure
egg_prod <- function(alpha, beta, w){ #Hyperbolic saturating function for egg production
  f <- (alpha*w) / (1 + ((alpha*w) / beta))
  return(f)
}
exp_sat <- function(a, b, c, x){ #Exponential saturating function
  f <- a*(1-exp(-b*x))*(1-exp(-c*x)) #translating eggs into miracidiae in snails
  return(f)
}

#Compute initial worms' pair
mw0 <- rbinom(length(cohort$w), cohort$w, 0.5) #male worms
fw0 <- cohort$w-mw0 #female worms
#Create initial cloud
eggs <- alpha*(pmin(mw0, fw0))*exp(-z*fw0)
contributions <- eggs * Age_profile(cohort$age)$y * cohort$Ind_sus


seeds <- 5
writeLines(c(""), "Sink.txt") #initiate log file

cluster <- makeCluster(min(parallel::detectCores(logical = FALSE), seeds))
registerDoParallel(cluster)

results <- foreach(k = 1:seeds,
                   .inorder = TRUE,
                   .errorhandling = "remove",
                   #.combine = bind_rows,
                   .packages = c("tidyverse")) %dopar% {
                    #for each seed:

                    #Time step events, for each individual
                    #Initialize
                    pop <- cohort %>%
                      mutate(mw = mw0, #male worms
                             fw = fw0) #female worms
                    N <- nrow(pop)
                    SAC <- which(pop$age >= 5 & pop$age <= 15)
                    
                    true_prev <- length(which(pop$w>0))/N
                    eggs_prev <- c(0)
                    eggs_prev_SAC <- c(0)
                    Heggs_prev <- c(0)
                    
                    #Create initial cloud
                    cloud <- sum(contributions) #External FOI (to start infection)
                    #m_in <- rep(0, 12*T)
                    for(t in 1:(12*T)){
                      
                      sink("Sink.txt", append=TRUE)
                      cat(paste(Sys.time(), ": Starting seed", k, "time step", t, "\n", sep = " "))
                      sink()
                      
                      #Individual events
                      nw <- rep(0, nrow(pop)) #individual exposures
                      ec <- eggs #individual egg counts
                      co <- contributions #individual contributions
                      
                      #Births
                      #for now birth rate does not depend on age-specific female fertility 
                      births <- rpois(1, (birth_rate*(nrow(pop)/1000))/12)
                      nb <- tibble(age = rep(0, births),
                                   sex = as.numeric(rbernoulli(births, 0.5)),
                                   w = rep(0, births),
                                   Ind_sus = rgamma(1, shape = k_w, scale = 1/k_w),
                                   mw = rep(0, births),
                                   fw = rep(0, births)) #new born
                      
                      pop <- bind_rows(pop, nb)
                      
                      #Deaths
                      ag <- as.numeric(cut(pop$age, c(-1, prob_death$Age_hi))) #age groups
                      # pop <- pop %>%
                      #   mutate(ag = ag)
                      dead <- c()
                      for(i in 1:nrow(pop)){
                        if(pop$sex[i]==0){
                          if(rbernoulli(1, p=prob_death[ag[i], "Male"]/12))
                            dead <- c(dead, i)
                        }
                        if(pop$sex[i]==1){
                          if(rbernoulli(1, p=prob_death[ag[i], "Female"]/12))
                            dead <- c(dead, i)
                        }    
                      }
                      
                      if(length(dead)>0)
                        pop <- pop[-dead,]
                      
                      #Update age and population size
                      pop$age <- pop$age + 1/12
                      N <- c(N, nrow(pop))
                      SAC <- which(pop$age >= 5 & pop$age <= 15)
                      
                      for(i in 1:nrow(pop)){
                        #Exposure
                        #The individual assumes new worms nw (FOI)
                        #and only a portion survives from the previous month
                        rate <- foi(l=cloud, zeta, v, a=pop$age[i], is=pop$Ind_sus[i])
                        nw[i] <- rpois(1, rate/cum_exp)
                        #they are assigned sex
                        malesnw <- rbinom(1, nw[i], 0.5)
                        
                        #Worms
                        #Worms' pairs (so mature) produce eggs
                        wp <- min(pop$mw[i], pop$fw[i])
                        eggs <- alpha*wp #*exp(-z*fw[i])
                        #Diagnosis 
                        ec[i] <- rnbinom(1, size=k_e, mu=eggs) 
                        
                        #Individual contribution 
                        #Eggs released that will become miracidia and infect snails
                        #neglect contribution rate for now
                        co[i] <- co_rate * eggs * Age_profile(pop$age[i])$y * pop$Ind_sus[i]
                        
                        #Control (MDA: 75% coverage, annual to SAC, starting to year 50 to 70)
                        if(t %in% c(12*seq(mda$start, mda$end, mda$frequency))){
                          if(pop$age[i] >= mda$age.lo & pop$age[i] <= mda$age.hi){
                            if(rbernoulli(1, mda$coverage)){ #it works in around 75% of the target pop
                              pop$mw[i] <- round((1-mda$efficacy)*pop$mw[i])
                              pop$fw[i] <- round((1-mda$efficacy)*pop$fw[i])
                           }
                         }
                        }
                        
                        #Worms are updated for the next month with:
                        #the newborn, which will be mature the next month
                        #also with survival portion of males and females worms from the previous month
                        pop$mw[i] <- pop$mw[i] - rbinom(1, pop$mw[i], phi1) + malesnw
                        pop$fw[i] <- pop$fw[i] - rbinom(1, pop$fw[i], phi1) + (nw[i] - malesnw)
                      }
                      #Reservoir/cloud
                      #Miracidiae intake at step t from the cloud
                      #m_in[t] <- sum(co)
                      #try first simple cloud without snails
                      #cloud <- m_in[t - 1]*expanding_factor
                      cloud <- sum(co) 
                      # if(t < 12*3) #External foi to start infection (3 years)
                      #   cloud <- cloud + Ext_foi
                      
                      #Summary statistics
                      true_prev[t] <- length(which(pop$mw>0|pop$fw>0))/nrow(pop)
                      eggs_prev[t] <- length(which(ec>0))/nrow(pop)
                      eggs_prev_SAC[t] <- length(which(ec[SAC]>0))/length(SAC)
                      Heggs_prev[t] <- length(which((ec*24)>=400))/nrow(pop)
                    }

                    res <- tibble(time = 1:(12*T),
                                  seed = rep(k, (12*T)),
                                  pop_size = N[-1],
                                  true_prev = true_prev,
                                  eggs_prev = eggs_prev,
                                  eggs_prev_SAC = eggs_prev_SAC,
                                  Heggs_prev = Heggs_prev)
                    #I would add seed number
                   }

stopCluster(cluster)

#Collating results
res <- bind_rows(results)

#Average results over simulations
avg_res <- res %>%
  group_by(time) %>%
  summarise(N = mean(pop_size),
            true_prev = mean(true_prev),
            eggs_prev = mean(eggs_prev),
            eggs_prev_SAC = mean(eggs_prev_SAC),
            Heggs_prev = mean(Heggs_prev))

#Plot population size
ggplot(res) +
  geom_line(aes(x=time, y=pop_size, group = seed), color = "grey20", alpha = 0.3) +
  geom_line(data=avg_res, aes(x=time, y=N), size=1) +
  scale_y_continuous(name = "Population size (counts)",
                     breaks = seq(0, 1000, 200),
                     #limits = c(0, 1000),
                     expand = c(0, 0)) +
  scale_x_continuous(name = "Time [years]",
                     limits = c(0, T*12),
                     expand = c(0, 0)) +
  expand_limits(x = 0,y = 0)

#Plot prevalence
ggplot(res) +
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
                     limits = c(0, 1200),
                     expand = c(0, 0)) +
  expand_limits(x = 0,y = 0)

#Average results
avg_res <- lapply(results, function(x) x$true_prev)
avg_res <- matrix(unlist(avg_res), ncol = 12*T, byrow = TRUE)
#runs <- tibble(True_prev = avg_res)
avg_trueprev <- apply(avg_res, 2, mean)

avg_res <- lapply(results, function(x) x$eggs_prev)
avg_res <- matrix(unlist(avg_res), ncol = 12*T, byrow = TRUE)
#runs <- bind_cols(runs, Eggs_prev = )
avg_eggsprev <- apply(avg_res, 2, mean)

avg_res <- lapply(results, function(x) x$eggs_prev_SAC)
avg_res <- matrix(unlist(avg_res), ncol = 12*T, byrow = TRUE)
avg_eggsprevSAC <- apply(avg_res, 2, mean)

avg_res <- lapply(results, function(x) x$Heggs_prev)
avg_res <- matrix(unlist(avg_res), ncol = 12*T, byrow = TRUE)
avg_Heggsprev <- apply(avg_res, 2, mean)


plot(avg_trueprev, type = 'l', xlab = "Time [Months]", ylim=c(0, 1))
lines(avg_eggsprev, type = 'l', xlab ="Time [Months]", col="dark blue")
lines(avg_eggsprevSAC, type = 'l', xlab ="Time [Months]", col="dark green")
lines(avg_Heggsprev, type = 'l', xlab ="Time [Months]", col="dark red")

# plot(res$true_prev[which(res$time>1400 & res$time<1800)], type = 'l', xlab ="Time [Months]")
# plot(res$eggs_prev[which(res$time>1400 & res$time<1800)], type = 'l', xlab ="Time [Months]")
# plot(res$eggs_prev_SAC[which(res$time>1400 & res$time<1800)], type = 'l',
#       col="red", xlab ="Time [Months]")

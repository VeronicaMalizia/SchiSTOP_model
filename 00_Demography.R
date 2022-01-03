#############################
#Author: Veronica Malizia
#Date: 27/11/2021
#R version: 3.6.1

#This will generate a first toy model for human demography (with age/birth/deaths)
#Data: 
# - Initial cohort population from WORMSIM input file
# - Death_rates_Uganda_2019 specific by age (from: https://apps.who.int/gho/data/view.searo.61730?lang=en)

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

##Parameters
T <- 400 #number of years
birth_rate <- 37 #Crude birth rate Uganda (per 1000 individuals)
#annual time step

seeds <- 1
writeLines(c(""), "Sink_demography.txt") #initiate log file

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
                     pop <- cohort
                     N <- nrow(pop)
                     
                     for(t in 1:T){ #(T*12)
                       
                       sink("Sink_demography.txt", append=TRUE)
                       cat(paste(Sys.time(), ": Starting seed", k, "time step", t, "\n", sep = " "))
                       sink()
                       
                       #Births
                       #for now birth rate does not depend on age-specific female fertility 
                       births <- rpois(1, (birth_rate*(nrow(pop)/1000)))
                       nb <- tibble(age=rep(0, births),
                                    sex=as.numeric(rbernoulli(births, 0.5))) #new born
                       
                       pop <- bind_rows(pop, nb)
                       
                       #Deaths
                       ag <- as.numeric(cut(pop$age, c(-1, prob_death$Age_hi))) #age groups
                       # pop <- pop %>%
                       #   mutate(ag = ag)
                       dead <- c()
                       for(i in 1:nrow(pop)){
                         if(pop$sex[i]==0){
                           if(rbernoulli(1, p=prob_death[ag[i], "Male"]))
                             dead <- c(dead, i)
                         }
                         if(pop$sex[i]==1){
                           if(rbernoulli(1, p=prob_death[ag[i], "Female"]))
                             dead <- c(dead, i)
                         }    
                       }
                       
                       #if(length(dead>0)){
                         pop <- pop[-dead,]
                       
                       #Update age
                       pop$age <- pop$age + 1
                       N <- c(N, nrow(pop))
                     }
                     
                     res <- tibble(time = 0:(T),
                                   seed = rep(k, ((T)+1)),
                                   pop_size = N)
                   }
                      
stopCluster(cluster)

#Collating results
res <- bind_rows(results)

#Average results over simulations
avg_res <- res %>%
  group_by(time) %>%
  summarise(N = mean(pop_size))
            
plot(res$time, res$pop_size, type = 'l')

ggplot(res) +
  geom_line(aes(x=time, y=pop_size, group = seed), color = "grey20", alpha = 0.3) +
  geom_line(data=avg_res, aes(x=time, y=N), size=1) +
  scale_y_continuous(name = "Population size (counts)",
                     breaks = seq(0, 1000, 200),
                     limits = c(0, 1000),
                     expand = c(0, 0)) +
  scale_x_continuous(name = "Time [years]",
                     limits = c(0, T),
                     expand = c(0, 0)) +
  expand_limits(x = 0,y = 0)

hist(cohort$age)
hist(pop$age)

#Saving image
tiff("Demography.tif", width=7, height=6, units = "in", res = 300)
ggplot(res) +
  geom_line(aes(x=time, y=pop_size, group = seed), color = "grey20", alpha = 0.3) +
  geom_line(data=avg_res, aes(x=time, y=N), size=1) +
  scale_y_continuous(name = "Population size (counts)",
                     breaks = seq(0, 1000, 200),
                     limits = c(0, 1000),
                     expand = c(0, 0)) +
  scale_x_continuous(name = "Time [years]",
                     limits = c(0, T),
                     expand = c(0, 0)) +
  expand_limits(x = 0,y = 0)
dev.off()

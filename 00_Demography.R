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
                 sex = c(rep(0, pop0$N_males[1]), rep(1, pop0$N_females[1]))) #Sub-Saharan Africa
#0=male, 1=female
for(i in 2:nrow(pop0)){
  tmp <- tibble(age = round(c(runif(pop0$N_males[i], pop0$Age_limit[i-1]+1, pop0$Age_limit[i]),
                              runif(pop0$N_females[i], pop0$Age_limit[i-1]+1, pop0$Age_limit[i]))),
                sex = c(rep(0, pop0$N_males[i]), rep(1, pop0$N_females[i])))
  cohort <- rbind(cohort, tmp) 
}

hist(cohort$age)

##Parameters
T <- 200 #number of years
ts <- T*12
birth_rate <- 34.8 #37 is crude annual birth rate Uganda, 34.8 for Sub-Saharan Africa (per 1000 individuals)
max.pop <- 700
#annual time step

seeds <- 10
writeLines(c(""), "Sink_demography.txt") #initiate log file
writeLines(c(""), "Find_bug.txt") #initiate log file

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
                     
                     for(t in 1:ts){ 
                       
                       sink("Sink_demography.txt", append=TRUE)
                       cat(paste(Sys.time(), ": Starting seed", k, "time step", t, "\n", sep = " "))
                       sink()
                       
                       #The reaper(acts annually)
                       if(t %in% seq(1, ts, 12) & nrow(pop)>max.pop)
                         pop <- slice_sample(pop, prop = 0.9)
                       
                       #Births
                       #for now birth rate does not depend on age-specific female fertility 
                       births <- rpois(1, (birth_rate*(nrow(pop)/1000))/12)
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
                           if(rbernoulli(1, p=prob_death[ag[i], "Male"]/12))
                             dead <- c(dead, i)
                         }
                         if(pop$sex[i]==1){
                           if(rbernoulli(1, p=prob_death[ag[i], "Female"]/12))
                             dead <- c(dead, i)
                         }    
                       }
                       
                       if(length(dead>0))
                         pop <- pop[-dead,]
                       
                       #Update age
                       pop$age <- pop$age + 1/12
                       N <- c(N, nrow(pop))
                       
                       # sink("Find_bug.txt", append=TRUE)
                       # cat(paste(Sys.time(), ": Time step", t, "NB:",
                       #           births, "deads", length(dead), "Pop size:", nrow(pop), "\n", sep = " "))
                       # sink()
                     }
                     
                     res <- tibble(time = 0:ts,
                                   seed = rep(k, (ts+1)),
                                   pop_size = N)
                   }
                      
stopCluster(cluster)

#Collating results
res <- bind_rows(results)

#Average results over simulations
avg_res <- res %>%
  group_by(time) %>%
  summarise(N = mean(pop_size))
            
ggplot(res) +
  geom_line(aes(x=time, y=pop_size, group = seed), color = "grey20", alpha = 0.3) +
  geom_line(data=avg_res, aes(x=time, y=N), size=1) +
  scale_y_continuous(name = "Population size (counts)",
                     breaks = seq(0, 1000, 200),
                     #limits = c(0, 1000),
                     expand = c(0, 0)) +
  scale_x_continuous(name = "Time [years]",
                     limits = c(0, ts),
                     expand = c(0, 0)) +
  expand_limits(x = 0,y = 0)

hist(cohort$age)
hist(pop$age, breaks = c(0, prob_death$Age_hi[-c(1, nrow(prob_death))]))

#Saving image
tiff("Demography_monthlyts.tif", width=7, height=6, units = "in", res = 300)
ggplot(res) +
  geom_line(aes(x=time, y=pop_size, group = seed), color = "grey20", alpha = 0.3) +
  geom_line(data=avg_res, aes(x=time, y=N), size=1) +
  scale_y_continuous(name = "Population size (counts)",
                     breaks = seq(0, 1000, 200),
                     limits = c(0, 1000),
                     expand = c(0, 0)) +
  scale_x_continuous(name = "Time [months]",
                     limits = c(0, ts),
                     expand = c(0, 0)) +
  expand_limits(x = 0,y = 0)
dev.off()


plot(N, type='l')

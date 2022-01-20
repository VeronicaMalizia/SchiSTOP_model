#############################
#Author: Veronica Malizia
#Date: 20/01/2022
#R version: 3.6.1

#This will generate a first toy model for human demography (with age/birth/deaths)
#To keep constant population, all deaths are replaced with new borns.
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

##Parameters
T <- 500 #number of years
ts <- T*12
birth_rate <- 34.8 #37 is crude annual birth rate Uganda, 34.8 for Sub-Saharan Africa (per 1000 individuals)
#max.pop <- 700
#monthly time step

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
                     pop <- cohort %>%
                       mutate(death_age=150)
                     N <- nrow(pop)
                     
                     for(t in 1:ts){ 
                       
                       sink("Sink_demography.txt", append=TRUE)
                       cat(paste(Sys.time(), ": Starting seed", k, "time step", t, "\n", sep = " "))
                       sink()
                       
                       #The reaper(acts annually)
                       # if(t %in% seq(1, ts, 12) & nrow(pop)>max.pop) #Januaries
                       #   pop <- slice_sample(pop, prop = 0.9)
                       # 
                       # #Births
                       # #for now birth rate does not depend on age-specific female fertility 
                       # births <- rpois(1, (birth_rate*(nrow(pop)/1000))/12)
                       # nb <- tibble(age=rep(0, births),
                       #              sex=as.numeric(rbernoulli(births, 0.5))) #new born
                       # 
                       # pop <- bind_rows(pop, nb)
                       
                       #Deaths
                       
                       if(t %in% seq(1, ts, 12)){ #Januaries
                         ag <- as.numeric(cut(pop$age, c(-1, prob_death$Age_hi))) #age groups
                         # pop <- pop %>%
                         #   mutate(ag = ag)
                         for(i in 1:nrow(pop)){
                           if(pop$sex[i]==0){
                             if(rbernoulli(1, p=prob_death[ag[i], "Male"])){
                               death_month <- runif(1, 1, 12)
                               pop$death_age[i] <- pop$age[i] + death_month/12
                             }
                           }
                           if(pop$sex[i]==1){
                             if(rbernoulli(1, p=prob_death[ag[i], "Female"])){
                               death_month <- runif(1, 1, 12)
                               pop$death_age[i] <- pop$age[i] + death_month/12
                             }
                           }    
                         }
                       }
                       
                       #Replace dead individuals with new borns (reset age and sex)
                       tmp <- which(pop$age >= pop$death_age)
                       if(length(tmp)>0){
                         pop$age[tmp] <- -1/12 #to end up at zero with age updating
                         pop$sex[tmp] <- as.numeric(rbernoulli(length(tmp), 0.5)) #later, update with sex ratio
                         pop$death_age[tmp] <- 150
                       } 
                       
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
                     limits = c(0, 1000),
                     expand = c(0, 0)) +
  scale_x_continuous(name = "Time [years]",
                     limits = c(0, ts),
                     expand = c(0, 0)) +
  expand_limits(x = 0,y = 0)

hist(cohort$age)
#hist(pop$age)
hist(pop$age, breaks = c(0:9*10))

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

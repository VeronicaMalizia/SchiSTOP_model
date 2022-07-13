#############################
#Author: Veronica Malizia
#Date: 20/06/2022
#R version: 3.6.1

#This script initializes human demography (with age/birth/deaths)
#The equilibrium reached by simulations will determine the initial population of the transmission model
#Data: 
# - Age distribution Uganda 2019
# - Death_rates_Uganda_2019 specific by age (from: https://apps.who.int/gho/data/view.searo.61730?lang=en)
#Main use: calibrating birth rate and migration rate to keep population constant. 
#Death probabilities are given to preserve age distribution.
#To save as output:
# - demography parameters
# - final population composition
#
#############################
rm(list = ls())

.libPaths(c("C:/Program Files/R/R-4.1.2/library",.libPaths()))
library(tidyverse)
library(readxl)
library(foreach)
library(doParallel)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source.dir <- "C:/Users/Z541213/Documents/Project/Model/Schisto_model"

################
#Initial plausible Ugandan population
#Starting with the initial cohort
################
age_distr0 <- read_xlsx("Age distribution_Uganda_2019.xlsx")
death_rates <- read.csv("death_rates_Uganda_2019.csv")
prob_death <- read.csv("prob_death_Uganda_2019.csv")

N <- 1000 #population size

#0=male, 1=female
cohort <- c()
for(i in 1:nrow(age_distr0)){
  people <- (age_distr0$Both_sexes[i]/100)*N
  tmp <- tibble(age = round(runif(round(people), age_distr0$Age_lo[i], age_distr0$Age_hi[i])),
                sex = c(as.numeric(rbernoulli(round(people)))))
  cohort <- rbind(cohort, tmp) 
}

hist(cohort$age)
table(cohort$sex)


##Parameters
T <- 100 #number of years
#ts <- T*12
birth_rate <- 36.5 #37.38 is crude annual birth rate Uganda 2019 (same y of available lifetables), 34.8 for Sub-Saharan Africa (per 1000 individuals)
#-1.09 is the net migration rate for 2022 for Uganda (per 1000 individuals)
emig_rate <- 20 #This approximates the emigration rate from rural villages (2021) (per 1000 individuals)
# Annual time step to get equilibrium

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
                     
                     ind_file <- c()
                     for(t in 1:T){ 
                       
                       sink("Sink_demography.txt", append=TRUE)
                       cat(paste(Sys.time(), ": Starting seed", k, "time step", t, "\n", sep = " "))
                       sink()
                       
                       #The reaper(acts annually)
                       # if(t %in% seq(1, ts, 12) & nrow(pop)>max.pop) #Januaries
                       #   pop <- slice_sample(pop, prop = 0.9)
                       
                       #Births (deterministic process)
                       #for now birth rate does not depend on age-specific female fertility 
                       #here, not stochastic. Just a fixed rate.
                       births <- birth_rate*(nrow(pop)/1000)
                       nb <- tibble(age=rep(0, births),
                                    sex=as.numeric(rbernoulli(births, 0.5)),
                                    death_age=150) #new born
                       
                       pop <- bind_rows(pop, nb)
                       
                       #Deaths 
                         ag <- as.numeric(cut(pop$age, c(-1, prob_death$Age_hi))) #age groups
                        
                         for(i in 1:nrow(pop)){
                           if(pop$sex[i]==0){
                             if(rbernoulli(1, p=prob_death[ag[i], "Male"])){
                               #death_month <- runif(1, 1, 12)
                               pop$death_age[i] <- pop$age[i]
                             }
                           }
                           if(pop$sex[i]==1){
                             if(rbernoulli(1, p=prob_death[ag[i], "Female"])){
                               #death_month <- runif(1, 1, 12)
                               pop$death_age[i] <- pop$age[i]
                             }
                           }
                         }
                       
                       tmp <- which(pop$age >= pop$death_age)
                       if(length(tmp)>0)
                         pop <- pop[-tmp,]
                       
                       #Migration
                       #For now we do not account for age-specific emigration
                       lambda <- emig_rate*(nrow(pop)/1000) #net out migration rate (annual)
                       emig_age <- which(pop$age>15&pop$age<55)
                       emigrated <- sample(emig_age, lambda)
                       if(length(emigrated)>0)
                         pop <- pop[-emigrated,]
                       
                       #Update age
                       pop$age <- pop$age + 1
                       N <- c(N, nrow(pop))
                       
                       #Save annual individual output
                       ind_file <- rbind(ind_file,
                                         select(pop, age, sex) %>%
                                         mutate(ID = 1:nrow(pop),
                                                time = t,
                                                seed = k))

                       filename <- paste("Ind_out_seed_", k, ".csv", sep="")
                       write.csv(ind_file, file.path(source.dir, "/Output_demography/", filename), row.names = F)
                       
                       ##ATTENTION!
                       #Need to cancel extra files if running less seeds at the following round
                  
                       # sink("Find_bug.txt", append=TRUE)
                       # cat(paste(Sys.time(), ": Time step", t, "NB:",
                       #           births, "deads", length(dead), "Pop size:", nrow(pop), "\n", sep = " "))
                       # sink()
                     }
                     
                     res <- tibble(time = 0:T,
                                   seed = rep(k, (T+1)),
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
                     breaks = seq(0, 10000, 500),
                     limits = c(0, 2000),
                     expand = c(0, 0)) +
  scale_x_continuous(name = "Time [years]",
                     limits = c(0, T),
                     expand = c(0, 0)) +
  expand_limits(x = 0,y = 0) 

hist(cohort$age)
hist(pop$age, breaks = c(0, prob_death$Age_hi[-c(1, nrow(prob_death))]))

#Saving image
tiff("Demography_monthlyts.tif", width=7, height=6, units = "in", res = 300)
...
dev.off()


plot(N, type='l')

#### Individual-level results
#Collate results by seeds
data_all <- list.files(path = file.path(source.dir, "/Output_demography/"),  # Identify all output CSV files
                       pattern = "Ind_out_seed_*", full.names = TRUE) %>% 
  lapply(read_csv, show_col_types = F) %>%              # Store all files in list
  bind_rows                                             # Combine data sets into one

##Age distribution by time
#https://r-charts.com/distribution/ggridges/
cohort_plot_age <- cohort %>%
  mutate(ID = 1:nrow(cohort),
         time = 0,
         seed = NA)
bind_rows(cohort_plot_age,
          data_all %>%
            filter(time %in% c(20*1:5))) %>%
  ggplot(aes(x=round(age))) +
  geom_histogram(aes(y = ..density.., colour=as.factor(seed)), binwidth = 6, 
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


#Save population at the end of simulation as starting population for the transmission model
#I save it as an age distribution
to.save <- data_all %>%
  filter(time == T) %>%
  group_by(seed, age) %>%
  tally() %>%
  group_by(age) %>% #summarise over seeds
  summarise(n = mean(n))

save(to.save, file = "Equilibrium_age_distribution.RData")


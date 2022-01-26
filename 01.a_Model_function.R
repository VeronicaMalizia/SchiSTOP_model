#############################
#Author: Veronica Malizia
#Date: 27/10/2021
#R version: 3.6.1

#This version will generate a toy model, with human demography (aging/birth/deaths)
#and basic transmission dynamics with a central reservoir.
#Data required: 
# 1. Initial cohort population from WORMSIM input file
# 2. Death_probabilities from WHO App - South Africa
#############################

library(tidyverse)
library(readxl)
library(foreach)
library(doParallel)

writeLines(c(""), "Sink.txt") #initiate log file
writeLines(c(""), "Find_bug.txt") #initiate log file

time.start <- Sys.time()
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
                       mutate(death_age = 150,
                              mw = mw0, #male worms
                              fw = fw0,
                              nw = 0,
                              ec = eggs0,
                              co = contributions0) #female worms
                     N <- nrow(pop)
                     SAC <- which(pop$age >= 5 & pop$age <= 15)
                     
                     true_prev <- c(0)
                     eggs_prev <- c(0)
                     eggs_prev_SAC <- c(0)
                     Heggs_prev <- c(0)
                     
                     #Create initial cloud
                     m_in <- sum(contributions0)
                     cloud <- hyp_sat(alpha=a, beta=b, m_in)
                     for(t in 2:(12*T)){
                       
                       sink("Sink.txt", append=TRUE)
                       cat(paste(Sys.time(), ": Starting seed", k, "time step", t, "\n", sep = " "))
                       sink()
                       
                       # #Beginning of each month, update demography
                       # #The reaper(acts annually)
                       # if(t %in% seq(1, (T*12), 12) & nrow(pop)>max.pop)
                       #   pop <- slice_sample(pop, prop = 0.9)
                       # 
                       # #Births
                       # #for now birth rate does not depend on age-specific female fertility
                       # births <- rpois(1, (birth_rate*(nrow(pop)/1000))/12)
                       # #new born
                       # nb <- tibble(age=rep(0, births),
                       #              sex=as.numeric(rbernoulli(births, 0.5)),
                       #              w = rep(0, births),
                       #              Ind_sus = rgamma(1, shape = k_w, scale = 1/k_w),
                       #              death_age = 150,
                       #              mw = rep(0, births),
                       #              fw = rep(0, births),
                       #              nw = rep(0, births),
                       #              ec = rep(0, births),
                       #              co = rep(0, births))
                       # 
                       # pop <- bind_rows(pop, nb)
                       # 
                       # #Deaths
                       # 
                       # if(t %in% seq(1, (T*12), 12)){ #Januaries
                       #   ag <- as.numeric(cut(pop$age, c(-1, prob_death$Age_hi))) #age groups
                       # 
                       #   for(i in 1:nrow(pop)){
                       #     if(pop$sex[i]==0){
                       #       if(rbernoulli(1, p=prob_death[ag[i], "Male"])){
                       #         death_month <- runif(1, 1, 12)
                       #         pop$death_age[i] <- pop$age[i] + death_month/12
                       #       }
                       #     }
                       #     if(pop$sex[i]==1){
                       #       if(rbernoulli(1, p=prob_death[ag[i], "Female"])){
                       #         death_month <- runif(1, 1, 12)
                       #         pop$death_age[i] <- pop$age[i] + death_month/12
                       #       }
                       #     }
                       #   }
                       # }
                       # 
                       # tmp <- which(pop$age >= pop$death_age)
                       # if(length(tmp)>0)
                       #   pop <- pop[-tmp,]
                       
                       #Replacement of deaths
                       #Deaths
                       if(t %in% seq(1, (T*12), 12)){ #Januaries
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

                       #Replace dead individuals with new borns (reset age, sex, ind quantities)
                       tmp <- which(pop$age >= pop$death_age)
                       if(length(tmp)>0){
                         pop$age[tmp] <- -1/12 #to end up at zero with age updating
                         pop$sex[tmp] <- as.numeric(rbernoulli(length(tmp), 0.5)) #later, update with sex ratio
                         pop$death_age[tmp] <- 150
                         pop$w[tmp] <- 0
                         pop$Ind_sus[tmp] <- rgamma(length(tmp), shape = k_w, scale = 1/k_w)
                         pop$mw[tmp] <- 0
                         pop$fw[tmp] <- 0
                         pop$nw[tmp] <- 0
                         pop$ec[tmp] <- 0
                         pop$co[tmp] <- 0
                       }

                       #Update age, population size, SAC and cumulative exposure
                       pop$age <- pop$age + 1/12
                       N[t] <- nrow(pop)
                       SAC <- which(pop$age >= 5 & pop$age <= 15)
                       
                       cum_exp <- sum(Age_profile(pop$age)$y*pop$Ind_sus) #Cumulative exposure
                       
                       #Parasitology (vectorised)
                       #Exposure
                       #The individual assumes new worms nw (FOI)
                       #and only a portion survives from the previous month
                       #One exposure rate per individual (based on age and ind. sus.)
                       rate <- foi(l=cloud, zeta, v, a=pop$age, is=pop$Ind_sus)
                       pop$nw <- rpois(nrow(pop), rate/cum_exp)
                       #per each human host, worms are assigned with sex
                       malesnw <- rbinom(nrow(pop), pop$nw, 0.5)
                       
                       #Control (MDA: 75% coverage, annual to SAC, starting to year 50 to 70)
                       if(t %in% c(12*seq(mda$start, mda$end, mda$frequency))){
                         for(i in 1:nrow(pop)){
                           if(pop$age[i] >= mda$age.lo & pop$age[i] <= mda$age.hi){
                             if(rbernoulli(1, mda$coverage)){ #it works in around 75% of the target pop
                               pop$mw[i] <- round((1-mda$efficacy)*pop$mw[i])
                               pop$fw[i] <- round((1-mda$efficacy)*pop$fw[i])
                             }
                           }
                         }
                       }
                       
                       #Worms
                       #Worms' pairs (so mature) produce eggs
                       #Shall we consider insemination probability??
                       wp <- pmin(pop$mw, pop$fw)
                       eggs <- alpha*wp #*exp(-z*pop$fw)
                       #Diagnosis 
                       pop$ec <- rnbinom(nrow(pop), size=k_e, mu=eggs) 
                       
                       #Individual contributions 
                       #Eggs released that will become miracidia and infect snails
                       #neglect contribution rate for now (set it to 1, to indicate no seasonal pattern)
                       pop$co <- co_rate * eggs * Age_profile(pop$age)$y * pop$Ind_sus
            
                       
                       #Worms are updated for the next month with:
                       #the newborn, which will be considered mature the next month
                       #also with survival portion of males and females worms from the previous month
                       pop$mw <- pop$mw - rbinom(nrow(pop), pop$mw, phi1) + malesnw
                       pop$fw <- pop$fw - rbinom(nrow(pop), pop$fw, phi1) + (pop$nw - malesnw)
                       
                       #Reservoir/cloud
                       #Miracidiae intake at step t from the cloud
                       m_in[t] <- sum(pop$co)
                       #first simple cloud without explicit snails, but with both larval stages
                       #we assume miracidiae spend 1 month into snails
                       #in the final environmental cloud we have cercariae
                       #we assume particles not infecting humans do not survive from the previous month
                       #Exponential saturation more accurate?
                       cloud <- hyp_sat(alpha=a, beta=b, m_in[t - 1])
                       #cloud <- sum(pop$co) 

                       sink("Find_bug.txt", append=TRUE)
                       cat(paste(Sys.time(), ": Seed:", k, "Time step", t, "cloud", cloud, "\n", sep = " "))
                       sink()
                       
                       #Summary statistics
                       true_prev[t] <- length(which((pop$mw+pop$fw)>0))/nrow(pop)
                       eggs_prev[t] <- length(which(pop$ec>0))/nrow(pop)
                       eggs_prev_SAC[t] <- length(which(pop$ec[SAC]>0))/length(SAC)
                       Heggs_prev[t] <- length(which((pop$ec*24)>=400))/nrow(pop)
                     }
                     
                     res <- tibble(time = 1:(12*T),
                                   seed = rep(k, (12*T)),
                                   pop_size = N,
                                   true_prev = true_prev,
                                   eggs_prev = eggs_prev,
                                   eggs_prev_SAC = eggs_prev_SAC,
                                   Heggs_prev = Heggs_prev)
                   }

stopCluster(cluster)
time.end <- Sys.time()
time.end - time.start

#############################
#Author: Veronica Malizia
#Date: 27/07/2021
#R version: 3.6.1

#This script contains model specification and is called into the main script 01_First\model.R
#
#############################

library(tidyverse)
library(readxl)
library(foreach)
library(doParallel)

writeLines(c(""), "Sink.txt") #initiate log file
writeLines(c(""), "Find_bug.txt") #initiate log file

cluster <- makeCluster(min(parallel::detectCores(logical = FALSE), seeds))
registerDoParallel(cluster)

results <- foreach(k = 1:seeds,
                   .inorder = TRUE,
                   .errorhandling = "pass", #remove or pass
                   .verbose = TRUE,
                   #.combine = bind_rows,
                   .packages = c("tidyverse", "deSolve")) %dopar% {
                     #for each seed:
                     
                     #Time step events, for each individual
                     #Initialize
                     pop <- cohort %>%
                       mutate(death_age = 150,
                              rate = 0,
                              jw1 = 0,
                              jw2 = 0,
                              jw3 = 0,
                              cum_dwp = 0,
                              ec = 0,
                              co = contributions0) #female worms
                     N <- nrow(pop)
                     SAC <- which(pop$age >= 5 & pop$age <= 15)
                     
                     true_prev <- c(0)
                     eggs_prev <- c(0)
                     eggs_prev_SAC <- c(0)
                     Heggs_prev <- c(0)
                     inf_snail <- c(0)
                     tot_snail <- snail.pop
                     
                     #Create initial cloud
                     #Initialize quantities
                     m_in <- sum(contributions0)
                     ## Start values
                     newstart <- c(S0, E0, I0, C0) #Initializing snails
                     cercariae <- 0 #to eventually plot the cercariae (cloud)
                     ind_file <- c()
                     for(t in 2:(12*T)){
                       
                       sink("Sink.txt", append=TRUE)
                       cat(paste(Sys.time(), ": Starting seed", k, "time step", t, "\n", sep = " "))
                       sink()
                       
                       #Beginning of each month, update demography
                       
                       #Births
                       #for now birth rate does not depend on age-specific female fertility
                       births <- birth_rate*(nrow(pop)/1000)/12
                       #new born
                       nb <- tibble(age=rep(0, births),
                                    sex=as.numeric(rbernoulli(births, 0.5)),
                                    wp = rep(0, births),
                                    Ind_sus = rgamma(births, shape = k_w, scale = 1/k_w),
                                    death_age = 150,
                                    rate = rep(0, births),
                                    jw1 = rep(0, births),
                                    jw2 = rep(0, births),
                                    jw3 = rep(0, births),
                                    cum_dwp = rep(0, births),
                                    ec = rep(0, births),
                                    co = rep(0, births))

                       pop <- bind_rows(pop, nb)

                       #Deaths

                       if(t %in% seq(1, (T*12), 12)){ #Januaries
                         ag <- as.numeric(cut(pop$age, c(-1, prob_death$Age_hi))) #age groups

                         for(i in 1:nrow(pop)){
                           if(pop$sex[i]==0){
                             if(rbernoulli(1, p=prob_death[ag[i], "Male"])){
                               death_month <- runif(1, 1, 12)
                               pop$death_age[i] <- pop$age[i] + death_month/12
                             }
                           }
                           if(pop$sex[i]==1){
                             if(rbernoulli(1, p=prob_death[ag[i], "Both_sexes"])){
                               death_month <- runif(1, 1, 12)
                               pop$death_age[i] <- pop$age[i] + death_month/12
                             }
                           }
                         }
                       }

                       tmp <- which(pop$age >= pop$death_age)
                       if(length(tmp)>0)
                         pop <- pop[-tmp,]
                       
                       #Migration
                       #For now we do not account for age-specific emigration
                       lambda <- (emig_rate*(nrow(pop)/1000))/12
                       emig_age <- which(pop$age>5&pop$age<55)
                       emigrated <- sample(emig_age, lambda)
                       if(length(emigrated)>0)
                         pop <- pop[-emigrated,]

                       #Update age, population size, SAC and cumulative exposure
                       pop$age <- pop$age + 1/12
                       N[t] <- nrow(pop)
                       SAC <- which(pop$age >= 5 & pop$age <= 15)
                       
                       cum_exp <- sum(Age_profile_exp(pop$age)$y*pop$Ind_sus) #Cumulative exposure
                       
                       #Parasitology (vectorised)
                       
                       #EXPOSURE and Worms update
                       #The individual assumes new worms nw (FOI)
                       #and only a portion survives from the previous month
                       #One exposure rate per individual (based on age and ind. sus.)
                       pop$rate <- FOIh(l=cercariae[t-1], zeta, v, a=pop$age, is=pop$Ind_sus)/cum_exp
                       if(lim_mechanism == "humans") #Immunity
                         pop$rate <- pop$rate*(1-hyp_sat(alpha = imm, beta = 1, w = pop$cum_dwp))
                       
                       #Worms are updated with:
                       #juveniles
                       #survival portion of males and females worms from the previous month
                       #Juvenile worms at stage 3 do pair and move to stage 4. Who doesn't, do not survive.
                       #per each human host, juvenile worms at stage 3 are assigned with sex
                       malesnw <- rbinom(nrow(pop), pop$jw3, 0.5)
                       new_pairs <- pmin(malesnw, pop$jw3)
                       dying_pairs <- rbinom(nrow(pop), pop$wp, phi1)
                       pop$wp <- pop$wp - dying_pairs + new_pairs
                       pop$cum_dwp <- pop$cum_dwp + dying_pairs 
                       #From here we track worm pairs as units. Not tracked individual male/female worms in this version
                       pop$jw3 <- pop$jw2 - rbinom(nrow(pop), pop$jw2, phi1)
                       pop$jw2 <- pop$jw1 - rbinom(nrow(pop), pop$jw1, phi1)
                       #First stage juveniles(new acquired)
                       pop$jw1 <- rpois(nrow(pop), pop$rate) 
                       if(t < ext.foi$duration*12)
                         pop$jw1 <- pop$jw1 + ext.foi$value
                      
                       #EGGS production (before updating worms)
                       #Worms' pairs (so mature at stage 4) produce eggs
                       #(Shall we consider insemination probability??)
                       #eggs is the daily egg amount passed to the intestine
                       eggs <- eggs_prod*pop$wp
                       
                       #Diagnosis 
                       #mu represents the expected egg load in a stool sample
                       mu <- alpha*pop$wp #*exp(-z*pop$fw) 
                       
                       if(lim_mechanism == "worms"){
                         eggs <- eggs*exp(-z*(pop$wp/2))
                         mu <- mu*exp(-z*(pop$wp/2))
                       }
                          
                       pop$ec <- rnbinom(nrow(pop), size=k_e, mu=mu) 
                       
                       #WORMS
                       #Control (MDA: 75% coverage, annual to SAC, starting to year 50 to 70)
                       killed_worms <- 0
                       if(t %in% c(12*seq(mda$start, mda$end, mda$frequency))){
                         treated <- sample(SAC, mda$coverage*length(SAC))
                         killed_worms <- round(mda$efficacy*pop$wp[treated])
                         pop$wp[treated] <- pop$wp[treated] - killed_worms
                         pop$cum_dwp[treated] <- pop$cum_dwp[treated] + killed_worms
                       }
                       
                       
                       #Individual CONTRIBUTIONS 
                       #Eggs released (daily quantity, since it will enter the SEIC system) that will become miracidiae and will infect snails
                       pop$co <- eggs * Age_profile_contr(pop$age)$y * pop$Ind_sus
                       
                       #RESERVOIR
                       #Miracidiae intake at step t by the reservoir
                       m_in[t] <- sum(pop$co)
                       
                       if(lim_mechanism != "snails")
                         cercariae[t] <- m_in[t-1] # I could add a portion of dying miracidiae
                       
                       if(lim_mechanism == "snails"){
                         #Run snails ODEs model
                         #in the final reservoir we have the output of cercariae
                         #we assume particles not infecting humans do not survive from the previous month
                         
                         #Call parameters and initial conditions
                         #SEIC is at daily time step
                         mirac.input = m_in[t] #chi*miracidiae will be divided by N[t] in the system #Civitello uses a cumulative factor of 0.01
                         parms  <- c(beta0 = max.reproduction.rate, k = carrying.capacity, v = mortality.rate,
                                     mir = mirac.input, l0 = max.invasion, chi = rej.prob,
                                     v2 = mortality.rate.infection, tau = infection.rate,
                                     lambda = cerc.prod.rate, m = cerc.mortality)
                         
                         #Set initial conditions from the last day of previous run/month
                         S0 = newstart[1]
                         E0 = newstart[2]
                         I0 = newstart[3]
                         C0 = newstart[4]
                         ## vector of time steps (30 days)
                         times <- 1:30
                         
                         ## Start values for steady state
                         xstart <- c(S = S0, E = E0, I = I0, C = C0)
                         
                         #Run and solve
                         out <-  lsoda(xstart, times, SEI, parms) 
                         
                         ## Translate the output into a data.frame
                         out2 <- as.data.frame(out)
  
                         #Saving results at t=30, to be saved as input of the next simulation
                         newstart <- c(out2$S[nrow(out2)],
                                       out2$E[nrow(out2)],
                                       out2$I[nrow(out2)],
                                       out2$C[nrow(out2)])
                         
                         #Cercarial production 
                         cercariae[t] <- out2$C[nrow(out2)]
                       }
                       # sink("Find_bug.txt", append=TRUE)
                       # cat(paste(Sys.time(), ": Seed:", k, "Time step", t, "res", cercariae[t], "\n", sep = " "))
                       # sink()
                       
                       #Summary statistics
                       true_prev[t] <- length(which(pop$wp>0))/nrow(pop)
                       eggs_prev[t] <- length(which(pop$ec>0))/nrow(pop)
                       eggs_prev_SAC[t] <- length(which(pop$ec[SAC]>0))/length(SAC)
                       Heggs_prev[t] <- length(which((pop$ec*24)>=400))/nrow(pop)
                       if(lim_mechanism == "snails"){
                       inf_snail[t] <- out2$I[nrow(out2)] 
                       tot_snail[t] <- (out2$S[nrow(out2)]+out2$E[nrow(out2)]+out2$I[nrow(out2)])
                       }
                       
                       #Save annual individual output
                       if(t %in% seq(12, (T*12), 10*12) & t > 500){ #Decembers, every 10 years
                         ind_file <- rbind(ind_file,
                                           select(pop, age, sex, rate, wp, ec, cum_dwp) %>%
                                           mutate(ID = 1:nrow(pop),
                                                  time = t/12,
                                                  seed = k))
                       }
                     }
                     
                     filename <- paste("Ind_out_seed_", k, ".csv", sep="")
                     write.csv(ind_file, 
                               file.path(source.dir, "/Output/", filename),
                               row.names = F)
                     
                     res <- tibble(time = 1:(12*T),
                                   seed = rep(k, (12*T)),
                                   pop_size = N,
                                   miracidiae = m_in/30, #to have a daily output
                                   cercariae = cercariae,
                                   true_prev = true_prev,
                                   eggs_prev = eggs_prev,
                                   eggs_prev_SAC = eggs_prev_SAC,
                                   Heggs_prev = Heggs_prev,
                                   inf_snail = inf_snail,
                                   tot_snail = tot_snail)
                   }

stopCluster(cluster)


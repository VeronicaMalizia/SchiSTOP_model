#############################
#Author: Veronica Malizia
#Date: 27/07/2021
#R version: 3.6.1

#This script contains the model specification 
#It is called and launched into the main script 01_First\model.R
#
#############################

library(tidyverse)
library(readxl)
library(foreach)
library(doParallel)
library(parallelly)

writeLines(c(""), "Sink.txt") #initiate log file
#writeLines(c(""), "Find_bug.txt") #initiate log file

cluster <- makeCluster(min(detectCores(logical = T), nrow(stoch_scenarios)))
#clusterEvalQ(cluster, .libPaths(c("C:\\Users\\NTD-Computing\\Bureaublad\\Veronica\\R-4.1.2\\library",.libPaths())))
#clusterEvalQ(cluster, .libPaths(c("C:/Program Files/R/R-4.1.2/library",.libPaths())))
registerDoParallel(cluster)

results <- foreach(k = 1:nrow(stoch_scenarios),
                   .inorder = TRUE,
                   .errorhandling = "pass", #remove or pass
                   .verbose = TRUE,
                   #.combine = bind_rows, #default is a list
                   .packages = c("tidyverse", "deSolve", "splitstackshape")) %dopar% {
                     
                     #set scenario-specific parameters
                     scen <- stoch_scenarios[k, ]
                     
                     #Load matched alphas for Density-dependent fecundity (DDF) given the endemicity
                     load(paste("Matched_alphas_", scen$endem, ".RData", sep = ""))
                     
                     #Transmission parameters
                     parms$parasite$zeta = scen$zeta #overall exposure rate.
                     parms$parasite$k_w = scen$worms_aggr
                     parms$snails$snail_transmission_rate = scen$tr_snails
                     
                     parms$parasite$ext.foi$value = scen$Ext_foi_value
                     parms$parasite$ext.foi$duration = scen$Ext_foi_duration
                     
                     #Regulating mechanisms
                     parms$immunity$imm = case_when(scen$imm_strength== "Absent" ~ 0,
                                                    scen$imm_strength== "Mild" ~ 0.0005,
                                                    scen$imm_strength== "Strong" ~ 0.002) #immunity slope parameter
                     parms$snails$carrying.capacity = case_when(scen$snails == "Absent" ~ 1, #No module
                                                                scen$snails == "Mild" ~ 20000,
                                                                scen$snails == "Strong" ~ 10000)
                     parms$parasite$eggs$alpha = case_when(scen$DDF_strength== "Absent" ~ alpha_lin,
                                                           scen$DDF_strength== "Mild" ~ alpha_mild,
                                                           scen$DDF_strength== "Strong" ~ alpha_strong) #fecundity parameter
                     parms$parasite$eggs$z = case_when(scen$DDF_strength== "Absent" ~ 0,
                                                       scen$DDF_strength== "Mild" ~ 0.00022,
                                                       scen$DDF_strength== "Strong" ~ 0.0007) #severity of density dependency
                     
                     
                     sink("Sink.txt", append=TRUE)
                     cat(paste(Sys.time(), ": Scenario nr.", k, "\n", sep = " "))
                     sink()
                     
                     #for each seed:
                     
                     #profvis({   
                     #Time step events, for each individual
                     #Initialize population
                     pop <- init$humans$pop 
                     N <- nrow(pop)
                     ever_lived <- N
                     SAC <- init$humans$SAC
                     
                     true_prev <- c(0)
                     eggs_prev <- c(0)
                     eggs_prev_SAC <- c(0)
                     Heggs_prev <- c(0)
                     Heggs_prev_SAC <- c(0)
                     avg_geom_intensity <- c(0)
                     avg_geom_intensity_SAC <- c(0)
                     avg_intensity <- c(0)
                     avg_intensity_SAC <- c(0)
                     inf_snail <- init$snails$I0
                     susc_snail <- init$snails$S0
                     exp_snail <- init$snails$E0
                     
                     #Create initial cloud
                     #Initialize quantities
                     m_in <- sum(init$environment$contributions0)
                     ## Start values
                     newstart <- c(init$snails$S0, init$snails$E0, init$snails$I0, init$snails$C0) #Initializing snails
                     cercariae <- 0 #to eventually plot the cercariae (cloud)
                     ind_file <- c()
                     for(t in 2:(12*T)){
                       
                       ###########################################
                       #Demography
                       ###########################################
                       
                       #Beginning of each month, update demography
                       ########## AGING
                       pop$age = pop$age + 1/12
                       
                       ########## BIRTHS
                       # (births and migration: deterministic processed -> fixed rate)
                       #for now birth rate does not depend on age-specific female fertility
                       births <- round(parms$demography$birth_rate*(nrow(pop)/1000)/12)
                       #new born
                       
                       pop <- bind_rows(pop, 
                                        data.frame(age=0,
                                                   sex=as.numeric(rbernoulli(births, 0.5)),
                                                   jw1 = 0,
                                                   Ind_sus = rgamma(births, shape = parms$parasite$k_w, scale = 1/parms$parasite$k_w),
                                                   complier = as.numeric(rbernoulli(births, 1-parms$mda$fr_excluded)),
                                                   ID = c((ever_lived+1):(ever_lived+births)),
                                                   death_age = 100,
                                                   age_group = 1,
                                                   rate = 0,
                                                   wp1 = 0,
                                                   wp2 = 0,
                                                   wp3 = 0,
                                                   jw2 = 0,
                                                   jw3 = 0,
                                                   cum_dwp = 0,
                                                   mu = 0,
                                                   ec = 0))
                       ever_lived <- ever_lived+births
                       
                       ########## DEATHS
                       if(t %in% seq(1, (T*12), 12)){ #Januaries
                         ag <- as.numeric(cut(pop$age, c(-1, prob_death$Age_hi))) #age groups
                         pop <- mutate(pop, age_group = ag)
                         
                         deaths <- as.tibble(table(ag)) %>%
                           rename(Ag = ag) %>%
                           mutate(n.deaths = round(prob_death[Ag, "Both_sexes"]*n)) %>%
                           filter(n.deaths > 0)
                         names(deaths$n.deaths) <- deaths$Ag
                         
                         # if(length(which(deaths$n.alive==0))>0) #we already remove the age group where nobody is left (if any)
                         #   pop <- filter(pop, age_group != deaths$Ag[deaths$n.alive==0])
                         
                         deads <- stratified(pop[pop$age_group %in% deaths$Ag, ], "age_group", 
                                             deaths$n.deaths, keep.rownames = T) %>%
                           mutate(death_age = age + runif(n(), 1, 12)/12) #starts from 1/12, but it reaches 12/12, that is the following January
                         pop[deads$rn,"death_age"] <- deads$death_age #we assign death age within that year
                       }
                       
                       pop <- filter(pop, age < death_age)
                       
                       ########## MIGRATION
                       #For now we do not account for age-specific emigration
                       lambda <- round(parms$demography$emig_rate*(nrow(pop)/1000)/12)
                       emigrated <- sample(which(pop$age>5&pop$age<55), lambda)
                       if(length(emigrated)>0)
                         pop <- pop[-emigrated,]
                       
                       ########## UPDATE population size, SAC and cumulative exposure
                       N[t] <- nrow(pop)
                       SAC <- which(pop$age >= 5 & pop$age <= 15)
                       
                       #Relative exposure
                       if(exposure=="ICL")
                         rel_exp <- Age_profile_exp(parms$exposure$ICL_derived$ages,
                                                     parms$exposure$ICL_derived$exp,
                                                     pop$age,
                                                     parms$exposure$ICL_derived$method)
                       if(exposure=="Sow")
                         rel_exp <- Age_profile_exp(parms$exposure$Sow_derived$ages,
                                                     parms$exposure$Sow_derived$exp,
                                                     pop$age,
                                                     parms$exposure$Sow_derived$method)
                       
                       cum_exp <- sum(rel_exp$y*pop$Ind_sus) #Cumulative exposure
                       
                       ###########################################
                       #Parasitology (vectorised)
                       ###########################################
                       
                       ########### 1. CONTROL (MDA: 75% coverage, annual to SAC, starting to year 50 to 70)
                       # If MDA applies, it is scheduled at the beginning of the month
                       killed_worms1 <- 0
                       killed_worms2 <- 0
                       killed_worms3 <- 0
                       if(t %in% c(12*seq(parms$mda$start, parms$mda$end, parms$mda$frequency))){
                         target <- which(pop$age >= parms$mda$age.lo & pop$age <= parms$mda$age.hi)
                         compliers <- which(pop$age >= parms$mda$age.lo & pop$age <= parms$mda$age.hi & pop$complier==1)
                         fr_compliers <- sum(pop$complier[target])/length(target)
                         real_coverage <- parms$mda$coverage/fr_compliers
                         treated <- sample(compliers, real_coverage*length(compliers)) #index of individuals
                         n_treated <- length(treated)
                         #Killing of worms in the three age baskets:
                         killed_worms1 <- rbinom(n_treated, size = pop$wp1[treated], prob = parms$mda$efficacy)
                         killed_worms2 <- rbinom(n_treated, size = pop$wp2[treated], prob = parms$mda$efficacy)
                         killed_worms3 <- rbinom(n_treated, size = pop$wp3[treated], prob = parms$mda$efficacy)
                         pop$wp1[treated] = pop$wp1[treated] - killed_worms1
                         pop$wp2[treated] = pop$wp2[treated] - killed_worms2
                         pop$wp3[treated] = pop$wp3[treated] - killed_worms3
                         pop$cum_dwp[treated] = pop$cum_dwp[treated] + 
                           killed_worms1 + killed_worms2 + killed_worms3
                       }
                       
                       ########### 2. EGGS production (before updating worms)
                       #Worms' pairs (so mature at stage 4) produce eggs
                       #(Shall we consider insemination probability??)
                       
                       #'mu' represents the expected egg load in a stool sample (41.7mg)
                       #'eggs' is the daily egg amount passed to the environment
                       # 24 is the conversion factor from egg counts and epg, because stoolxday is expressed in grams
                       Tot_wp <- pop$wp1+pop$wp2+pop$wp3
                       
                       pop$mu <- parms$parasite$eggs$alpha*Tot_wp*expon_reduction(parms$parasite$eggs$z, w=Tot_wp) 
                       
                       eggs <- round(pop$mu*24*parms$parasite$eggs$gr_stool) #daily quantity, since particles in the environment are short-lived
                       
                       ########## 3. DIAGNOSIS 
                       pop$ec = rnbinom(nrow(pop), size=parms$parasite$eggs$k_e, mu=pop$mu)  
                       
                       ########## 4. CONTRIBUTION to the environment
                       #Individual contributions
                       #Eggs are passed to the intestine and released in the environment
                       # TIP: I can add a delay in eggs' excretion
                       
                       #We can skip contribution, we assume here that contribution is not dependent on age
                       #pop$co <- round(eggs * Age_profile_contr(pop$age)$y)
                       
                       ########## 5. LIFE IN THE ENVIRONMENTAL RESERVOIR
                       #Miracidiae intake at step t by the reservoir
                       m_in[t] <- sum(eggs)
                       
                       if(scen$snails == "Absent")
                         cercariae[t] <- m_in[t-1] 
                       
                       if(scen$snails != "Absent"){
                         #Run snails ODEs model
                         #in the final reservoir we have the output of cercariae
                         #we assume particles not infecting humans do not survive from the previous month
                         
                         #Call parameters and initial conditions
                         #SEIC runs at daily time step
                         mirac.input = m_in[t-1] #chi*miracidiae will be divided by N[t] in the system #Civitello uses a unique factor of 0.01
                         parms.ODE  <- c(beta0 = parms$snails$max.reproduction.rate, 
                                         k = parms$snails$carrying.capacity, 
                                         v = parms$snails$mortality.rate,
                                         mir = mirac.input, 
                                         b = parms$snails$snail_transmission_rate,
                                         v2 = parms$snails$mortality.rate.infection, 
                                         tau = parms$snails$infection.rate,
                                         lambda = parms$snails$cerc.prod.rate, 
                                         m = parms$snails$cerc.mortality)
                         
                         #Set initial conditions from the last day of previous run/month
                         ## Start values for steady state
                         xstart <- c(S = newstart[1], E = newstart[2], I = newstart[3], C = newstart[4])
                         
                         ## vector of time steps (30 days)
                         times <- 1:30
                         
                         #Run and solve
                         out <-  lsoda(xstart, times, SEI, parms.ODE, atol = 1e-4, rtol = 1e-4) #, verbose = T
                         
                         ## Translate the output into a data.frame
                         out2 <- as.data.frame(out)
                         
                         #Saving results at t=30, to be saved as input of the next simulation
                         newstart <- c(out2$S[nrow(out2)],
                                       out2$E[nrow(out2)],
                                       out2$I[nrow(out2)],
                                       out2$C[nrow(out2)])
                         if(out2$E[nrow(out2)]<1e-10)
                           newstart[2] <- 0
                         if(out2$I[nrow(out2)]<1e-10)
                           newstart[3] <- 0
                         
                         #Cercarial production 
                         cercariae[t] <- out2$C[nrow(out2)]
                         
                         if(out2$C[nrow(out2)]<1e-10){
                           newstart[4] <- 0
                           cercariae[t] <- 0
                         }
                         
                       }
                       
                       ########## 6. EXPOSURE OF HUMANS 
                       #Each individual assumes new worms (FOIh) and immunity applies
                       #One exposure rate per individual (based on age and ind. sus.)
                       pop$rate = FOIh(l=cercariae[t], zeta=parms$parasite$zeta, rel_exp=rel_exp, is=pop$Ind_sus)*expon_reduction(parms$immunity$imm, pop$cum_dwp)/cum_exp
                       #pop$rate <- pop$rate*logistic(k=imm, w0=w0_imm, w = pop$cum_dwp)
                       
                       ########## 7. SAVE OUTPUT
                       # sink("Find_bug.txt", append=TRUE)
                       # cat(paste(Sys.time(), ": Seed:", k, "Time step", t, "res", cercariae[t], "\n", sep = " "))
                       # sink()
                       
                       #Summary statistics
                       tot_worms <- pop$jw1+ pop$jw2 + pop$jw3 + pop$wp1 + pop$wp2 + pop$wp3
                       true_prev[t] <- length(which(tot_worms>0))/nrow(pop)
                       eggs_prev[t] <- length(which(pop$ec>0))/nrow(pop)
                       Heggs_prev[t] <- length(which((pop$ec*24)>=400))/nrow(pop)
                       eggs_prev_SAC[t] <- length(which(pop$ec[SAC]>0))/length(SAC)
                       Heggs_prev_SAC[t] <- length(which((pop$ec[SAC]*24)>=400))/length(SAC)
                       avg_geom_intensity[t] <- geom_mean(pop$ec+1)-1
                       avg_geom_intensity_SAC[t] <- geom_mean(pop$ec[SAC]+1)-1
                       avg_intensity[t] <- mean(pop$ec)
                       avg_intensity_SAC[t] <- mean(pop$ec[SAC])
                       if(scen$snails != "Absent"){
                         inf_snail[t] <- out2$I[nrow(out2)] 
                         susc_snail[t] <- out2$S[nrow(out2)]
                         exp_snail[t] <- out2$E[nrow(out2)]
                       }
                       
                       #Save annual individual output
                       if(write.output == TRUE){
                         if(t %in% seq(12, (T*12), fr*12) & t > 500){ #Decembers, every 10 years
                           ind_file <- rbind(ind_file,
                                             select(pop, ID, age, age_group, rate, wp1, wp2, wp3, ec, mu, cum_dwp) %>%
                                               mutate(tot_wp = wp1+wp2+wp3,
                                                      time = t/12, #years
                                                      seed = scen$seed,
                                                      Endemicity = scen$endem,
                                                      zeta = scen$zeta,
                                                      worms_aggr = scen$worms_aggr,
                                                      tr_snails = scen$tr_snails,
                                                      Immunity = scen$imm_strength,
                                                      Snails = scen$snails,
                                                      DDF = scen$DDF_strength))
                         }
                       }
                       
                       ########## 8. WORMS UPDATE (for next month)
                       
                       ###Juveniles worms
                       #Juvenile worms at stage 3 do pair and move to stage 4. Who doesn't, do not survive.
                       #per each human host, juvenile worms at stage 3 are assigned with sex
                       malesnw <- rbinom(nrow(pop), pop$jw3, 0.5)
                       new_pairs <- pmin(malesnw, pop$jw3) #From here on we track worm pairs as units. 
                       
                       pop$jw3 = pop$jw2 
                       pop$jw2 = pop$jw1 
                       
                       #First stage juveniles(new acquired)
                       pop$jw1 = rpois(nrow(pop), pop$rate) 
                       if(t <= parms$parasite$ext.foi$duration*12)
                         pop$jw1 = pop$jw1 + round(parms$parasite$ext.foi$value*pop$Ind_sus)
                       
                       ### Adults worm pairs
                       # Aging of worms: Erlang distributed lifespans
                       # Three baskets: wp1, wp2, wp3, three aging groups
                       aging_1 <- rbinom(nrow(pop), pop$wp1, parms$parasite$phi)
                       aging_2 <- rbinom(nrow(pop), pop$wp2, parms$parasite$phi)
                       aging_3 <- rbinom(nrow(pop), pop$wp3, parms$parasite$phi) #dying_pairs
                       
                       pop$wp1 = pop$wp1 + new_pairs - aging_1
                       pop$wp2 = pop$wp2 + aging_1 - aging_2
                       pop$wp3 = pop$wp3 + aging_2 - aging_3
                       
                       # Cumulated death worm pairs
                       pop$cum_dwp = pop$cum_dwp + aging_3
                     }
                     
                     if(write.output==TRUE){
                       filename <- paste("Ind_out_seed_", scen$seed,
                                         "_", scen$endem,
                                         "_Imm=", scen$imm_strength,
                                         "Sn=", scen$snails,
                                         "DDF=", scen$DDF_strength, ".RDS", sep="")
                       #Write
                       saveRDS(ind_file,
                                 file.path(ind.output.dir, filename))
                     }
                     
                     res <- tibble(time = 1:(12*T),
                                   seed = rep(scen$seed, (12*T)),
                                   Endemicity = rep(scen$endem, (12*T)),
                                   zeta = rep(scen$zeta, (12*T)),
                                   worms_aggr = rep(scen$worms_aggr, (12*T)),
                                   tr_snails = rep(scen$tr_snails, (12*T)),
                                   Immunity = rep(scen$imm_strength, (12*T)),
                                   Snails = rep(scen$snails, (12*T)),
                                   DDF= rep(scen$DDF_strength, (12*T)),
                                   pop_size = N,
                                   miracidiae = m_in, 
                                   cercariae = cercariae,
                                   true_prev = true_prev,
                                   eggs_prev = eggs_prev,
                                   eggs_prev_SAC = eggs_prev_SAC,
                                   Heggs_prev = Heggs_prev,
                                   Heggs_prev_SAC = Heggs_prev_SAC,
                                   avg_geom_intensity = avg_geom_intensity,
                                   avg_geom_intensity_SAC = avg_geom_intensity_SAC,
                                   avg_intensity = avg_intensity,
                                   avg_intensity_SAC = avg_intensity_SAC,
                                   inf_snail = inf_snail,
                                   susc_snail = susc_snail,
                                   exp_snail = exp_snail)
                     #})
                   }

stopCluster(cluster)

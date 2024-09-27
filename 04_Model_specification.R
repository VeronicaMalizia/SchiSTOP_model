#############################
#Author: Veronica Malizia
#Date: 27/07/2021
#R version: 3.6.1

# This is the main script containing the model's formulation and specifications
# It contains the function defining the ABM, with integrated ODEs module for snails dynamics
# This script is called and launched from the "console script" 05_Run_model.R
# It runs in parallel simulations using "doParallel" and "foreach" packages
#
# The functions called within the script are stored in 01_Handy_functions.R
#############################

library(tidyverse)
library(readxl)
library(foreach)
#library(doParallel)
library(doFuture)
library(future)
#library(parallelly)

#writeLines(c(""), "Sink.txt") #initiate log file

#This is for doParallel
#num_cores <- as.numeric(Sys.getenv("SLURM_CPUS_ON_NODE", "1"))
#cluster <- makeCluster(min(num_cores, nrow(stoch_scenarios)))
# line below only if library folder is different than default
#clusterEvalQ(cluster, .libPaths(c("C:/Program Files/R/R-4.1.2/library",.libPaths())))
#registerDoParallel(cluster)

#This is for doFuture
#registerDoFuture()
plan(multicore, workers = 16)

results <- foreach(k = 1:nrow(stoch_scenarios),
                   .inorder = TRUE,
                   .errorhandling = "pass", #remove or pass
                   .verbose = TRUE,
                   .combine = 'list', #default is a list
                   .options.future = list(seed = TRUE, packages = c("tidyverse", "deSolve", "splitstackshape"))) %dofuture% {
                     
                     #set scenario-specific parameters
                     scen <- stoch_scenarios[k, ]

                     #Load matched alphas for Density-dependent fecundity (DDF) given the endemicity
                     load(paste("Matched_alphas_", scen$endem, ".RData", sep = ""))
                     
                     #Transmission parameters
                     parmsk <- parms
                     parmsk$parasite$zeta = scen$zeta #overall exposure rate.
                     parmsk$parasite$k_w = scen$worms_aggr
                     parmsk$snails$snail_transmission_rate = scen$tr_snails
                     
                     parmsk$parasite$ext.foi$value = scen$Ext_foi_value
                     parmsk$parasite$ext.foi$duration = scen$Ext_foi_duration
                     
                     #Regulating mechanisms
                     parmsk$immunity$imm = case_when(scen$imm_strength== "Absent" ~ 0,
                                                    scen$imm_strength== "Mild" ~ 0.0005,
                                                    scen$imm_strength== "Strong" ~ 0.002) #immunity slope parameter
                     parmsk$snails$carrying.capacity = case_when(scen$snails == "Absent" ~ 1, #No ODEs module is called
                                                                scen$snails == "Mild" ~ 20000,
                                                                scen$snails == "Strong" ~ 10000)
                     parmsk$parasite$eggs$alpha = case_when(scen$DDF_strength== "Absent" ~ alpha_lin,
                                                           scen$DDF_strength== "Mild" ~ alpha_mild,
                                                           scen$DDF_strength== "Strong" ~ alpha_strong) #fecundity parameter
                     parmsk$parasite$eggs$z = case_when(scen$DDF_strength== "Absent" ~ 0,
                                                       scen$DDF_strength== "Mild" ~ 0.00022,
                                                       scen$DDF_strength== "Strong" ~ 0.0007) #severity of density dependency
 
                    #  log_file <- file.path(source.dir, paste("/Sink_worker_", k, ".txt", sep = ""))
                    #  writeLines(c(""), log_file) #initiate log file

                    #  sink(log_file, append=TRUE)
                    #  cat(paste(Sys.time(), ": Scenario nr.", k, "\n", sep = " "))
                    #  flush.console()
                    #  sink()

                     #for each seed:
                     
                     #profvis({   
                     #Time step events, for each individual
                     #Initialize population
                     initk <- init
                     pop <- initk$humans$pop 
                     N <- nrow(pop)
                     ever_lived <- N
                     SAC <- initk$humans$SAC
                     
                     true_prev <- c(0)
                     eggs_prev <- c(0)
                     eggs_prev_SAC <- c(0)
                     Heggs_prev <- c(0)
                     Heggs_prev_SAC <- c(0)
                     avg_geom_intensity <- c(0)
                     avg_geom_intensity_SAC <- c(0)
                     avg_intensity <- c(0)
                     avg_intensity_SAC <- c(0)
                     inf_snail <- initk$snails$I0
                     susc_snail <- initk$snails$S0
                     exp_snail <- initk$snails$E0
                     
                     #Create initial cloud
                     #Initialize quantities
                     m_in <- sum(initk$environment$contributions0)
                     ## Start values
                     newstart <- c(initk$snails$S0, initk$snails$E0, initk$snails$I0, initk$snails$C0) #Initializing snails
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
                       births <- round(parmsk$demography$birth_rate*(nrow(pop)/1000)/12)
                       #new born
                       
                       pop <- bind_rows(pop, 
                                        data.frame(age=0,
                                                   sex=as.numeric(rbernoulli(births, 0.5)),
                                                   jw1 = 0,
                                                   Ind_sus = rgamma(births, shape = parmsk$parasite$k_w, scale = 1/parmsk$parasite$k_w),
                                                   complier = as.numeric(rbernoulli(births, 1-parmsk$mda$fr_excluded)),
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
                       #For now we do not account for age-specific emigration. We only select an eligible age group (5-55)
                       lambda <- round(parmsk$demography$emig_rate*(nrow(pop)/1000)/12)
                       emigrated <- sample(which(pop$age>5&pop$age<55), lambda)
                       if(length(emigrated)>0)
                         pop <- pop[-emigrated,]
                       
                       ########## UPDATE population size, SAC vector, and cumulative exposure
                       N[t] <- nrow(pop)
                       SAC <- which(pop$age >= 5 & pop$age <= 15)
                       
                       #Relative exposure
                       if(exposure=="ICL")
                         rel_exp <- Age_profile_exp(parmsk$exposure$ICL_derived$ages,
                                                     parmsk$exposure$ICL_derived$exp,
                                                     pop$age,
                                                     parmsk$exposure$ICL_derived$method)
                       if(exposure=="Sow")
                         rel_exp <- Age_profile_exp(parmsk$exposure$Sow_derived$ages,
                                                     parmsk$exposure$Sow_derived$exp,
                                                     pop$age,
                                                     parmsk$exposure$Sow_derived$method)
                       
                       cum_exp <- sum(rel_exp$y*pop$Ind_sus) #Cumulative exposure
                       
                       ###########################################
                       #Parasitology (vectorised)
                       ###########################################
                       
                       ########### 1. CONTROL (MDA: 75% coverage, annual to target population)
                       # If MDA occurs, it is scheduled at the beginning of the month
                       killed_worms1 <- 0
                       killed_worms2 <- 0
                       killed_worms3 <- 0
                       if(t %in% c(12*seq(parmsk$mda$start, parmsk$mda$end, parmsk$mda$frequency))){
                         target <- which(pop$age >= parmsk$mda$age.lo & pop$age <= parmsk$mda$age.hi)
                         compliers <- which(pop$age >= parmsk$mda$age.lo & pop$age <= parmsk$mda$age.hi & pop$complier==1)
                         fr_compliers <- sum(pop$complier[target])/length(target)
                         real_coverage <- parmsk$mda$coverage/fr_compliers
                         treated <- sample(compliers, real_coverage*length(compliers)) #index of individuals
                         n_treated <- length(treated)
                         #Killing of worms in the three age baskets:
                         killed_worms1 <- rbinom(n_treated, size = pop$wp1[treated], prob = parmsk$mda$efficacy)
                         killed_worms2 <- rbinom(n_treated, size = pop$wp2[treated], prob = parmsk$mda$efficacy)
                         killed_worms3 <- rbinom(n_treated, size = pop$wp3[treated], prob = parmsk$mda$efficacy)
                         pop$wp1[treated] = pop$wp1[treated] - killed_worms1
                         pop$wp2[treated] = pop$wp2[treated] - killed_worms2
                         pop$wp3[treated] = pop$wp3[treated] - killed_worms3
                         pop$cum_dwp[treated] = pop$cum_dwp[treated] + 
                           killed_worms1 + killed_worms2 + killed_worms3
                       }
                       
                       ########### 2. EGGS production (before updating worms)
                       #Worms' pairs (so mature at stage 4) produce eggs
                       #Can be improved with introduction of insemination probability
                       
                       #'mu' represents the expected egg load in a stool sample (41.7mg)
                       #'eggs' is the daily egg amount passed to the environment
                       # 24 is the Kato-Katz conversion factor from egg counts to epg, because (stool x day) is expressed in grams
                       Tot_wp <- pop$wp1+pop$wp2+pop$wp3
                       
                       pop$mu <- parmsk$parasite$eggs$alpha*Tot_wp*expon_reduction(parmsk$parasite$eggs$z, w=Tot_wp) 
                       
                       eggs <- round(pop$mu*24*parmsk$parasite$eggs$gr_stool) #daily quantity, since particles in the environment are short-lived
                       
                       ########## 3. DIAGNOSIS 
                       # 'ec' represents the egg counts detected with a simulated Kato-Katz test. Random draw from a Negative Binomial distribution
                       pop$ec = rnbinom(nrow(pop), size=parmsk$parasite$eggs$k_e, mu=pop$mu)  
                       
                       ########## 4. CONTRIBUTION to the environment
                       #Individual contributions
                       #Eggs are passed to the intestine and released in the environment
                       # TIP: this part can be improved with a delay in eggs' excretion
                       
                       #This is commented in this version of SchiSTOP, as we assume that contribution is not dependent on age
                       #but the feature can be turned on by the user
                       #pop$co <- round(eggs * Age_profile_contr(pop$age)$y)
                       
                       ########## 5. LIFE IN THE ENVIRONMENTAL RESERVOIR
                       #Miracidiae intake at step t by the reservoir
                       m_in[t] <- sum(eggs)
                       
                       if(scen$snails == "Absent")
                         cercariae[t] <- m_in[t-1] 
                       
                       if(scen$snails != "Absent"){
                         #Run snails ODEs model
                         #in the final reservoir we have the total output of cercariae
                         #we assume particles not infecting humans do not survive from the previous month
                         
                         #Call parameters and initial conditions
                         #SEI-system runs at daily time step
                         mirac.input = m_in[t-1] #miracidial input
                         
                         #for details about the parameters see 02_Parameters_Smansoni.R
                         parms.ODE  <- c(beta0 = parmsk$snails$max.reproduction.rate, 
                                         k = parmsk$snails$carrying.capacity, 
                                         v = parmsk$snails$mortality.rate,
                                         mir = mirac.input, 
                                         b = parmsk$snails$snail_transmission_rate,
                                         v2 = parmsk$snails$mortality.rate.infection, 
                                         tau = parmsk$snails$infection.rate,
                                         lambda = parmsk$snails$cerc.prod.rate, 
                                         m = parmsk$snails$cerc.mortality)
                         
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
                       #Each individual acquires new worms (FOIh) and anti-reinfection immunity here applies
                       #One exposure rate per individual (based on age and individual susceptibility)
                       pop$rate = FOIh(l=cercariae[t], zeta=parmsk$parasite$zeta, rel_exp=rel_exp, is=pop$Ind_sus)*expon_reduction(parmsk$immunity$imm, pop$cum_dwp)/cum_exp
                       #pop$rate <- pop$rate*logistic(k=imm, w0=w0_imm, w = pop$cum_dwp)
                       
                       ########## 7. SAVE OUTPUT
 
                       #Summary statistics for population-level output 
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
                       
                       #Save annual individual-level output
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
                       
                       ########## 8. WORMS UPDATE (for next monthly time step)
                       
                       ###Juveniles worms
                       #Juvenile worms at stage 3 do pair and move to stage 4. The ones who don't, do not survive.
                       #per each human host, juvenile worms at stage 3 are assigned with sex
                       malesnw <- rbinom(nrow(pop), pop$jw3, 0.5)
                       new_pairs <- pmin(malesnw, pop$jw3) #From here on we track worm pairs as infective units. 
                       
                       pop$jw3 = pop$jw2 
                       pop$jw2 = pop$jw1 
                       
                       #First stage juveniles (new acquired)
                       pop$jw1 = rpois(nrow(pop), pop$rate) 
                       if(t <= parmsk$parasite$ext.foi$duration*12)
                         pop$jw1 = pop$jw1 + round(parmsk$parasite$ext.foi$value*pop$Ind_sus)
                       
                       ### Adults worm pairs
                       # Aging of worms: Erlang distributed lifespans
                       # Three baskets: wp1, wp2, wp3, three aging groups
                       aging_1 <- rbinom(nrow(pop), pop$wp1, parmsk$parasite$phi)
                       aging_2 <- rbinom(nrow(pop), pop$wp2, parmsk$parasite$phi)
                       aging_3 <- rbinom(nrow(pop), pop$wp3, parmsk$parasite$phi) #these are the dying_pairs
                       
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
                    #prova = 2
                    #return(prova) #print di base funziona
                    #return(res) 
                    #write.csv(res, file = file.path(pop.output.dir, paste("/Res_", k, ".RDS", sep = "")))
                     #}) #This come from profvis
                    }
#print(results)
#stopCluster(cluster)
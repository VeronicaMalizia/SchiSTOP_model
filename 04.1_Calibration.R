#############################
#Author: Veronica Malizia
#R version: 4.1.2

#Calibration code for the three transmission parameters:
## transmission on humans (zeta)
## aggregation of worms
## transmission on snails

#Figures to match: 
## Moderate endemicity: 30% prevalence in SAC, 6% HI prevalence in SAC, 6% snail prevalence
## Low endemicity: 10% prevalence in SAC, ~ 0% HI prevalence in SAC, 6% snail prevalence
## High endemicity: 60% prevalence in SAC, 20% HI prevalence in SAC, 6% snail prevalence
#############################
rm(list = ls())

#Loading packages
.libPaths(c("C:/Program Files/R/R-4.1.2/library",.libPaths()))
library(tidyverse)
library(readxl)
library(rstudioapi)
library(beepr)
library(EasyABC)
library(deSolve)
library(splitstackshape)

# Working directory
source.dir <- dirname(getActiveDocumentContext()$path)
setwd(source.dir)

# 'use_seed' must be turned to TRUE.
#n=10

###########
#Priors
###########
#theta = c(zeta, worms_aggr, tr_snails)

theta_prior = list(c("normal", -7, 0.5), # zeta (> 0) #-7.5
                   c("exponential", 1/0.15)) # worms_aggr (in (0,1)) or beta(che può essere skewed)
                   #c("lognormal", -23, 1)) # tr_snails (> 0, very small) 
#kw exponential 
theta_prior

#######################
# OUTPUT
#######################
#Select the three prevalence indicators at the end of simulation

#Load functions
source("01_Handy_functions.R")
##Load parameters
source("02_Parameters_Smansoni.R")
#Load initial conditions
source("01.1_Initial_conditions.R")

# Input (user choices)
#Through this file the user can set the desired modelling scenario
source("Setting_simulation_scenario.R")
#set scenario-specific parameters
scen <- stoch_scenarios
#parms$parasite$zeta = scen$zeta #overall exposure rate.  
parms$immunity$imm = case_when(scen$imm_strength== "Absent" ~ 0,
                               scen$imm_strength== "Mild" ~ 0.0005,
                               scen$imm_strength== "Strong" ~ 0.005) #immunity slope parameter
parms$snails$carrying.capacity = case_when(scen$snails == "Absent" ~ 1, #No module
                                           scen$snails == "Mild" ~ 20000,
                                           scen$snails == "Strong" ~ 10000)
parms$parasite$eggs$alpha = case_when(scen$DDF_strength== "Absent" ~ alpha_lin,
                                      scen$DDF_strength== "Mild" ~ alpha_mild,
                                      scen$DDF_strength== "Strong" ~ alpha_strong) #fecundity parameter
parms$parasite$eggs$z = case_when(scen$DDF_strength== "Absent" ~ 0,
                                  scen$DDF_strength== "Mild" ~ 0.0005,
                                  scen$DDF_strength== "Strong" ~ 0.0007) #severity of density dependency
# parms$parasite$k_w = scen$worms_aggr
# parms$snails$snail_transmission_rate = scen$tr_snails

#theta = c(-5, 0.25, 10^(-10))

#OUTPUT (model) FUNCTION
schisto_model <- function(theta){
  
  # parms$parasite$zeta = exp(theta[1])
  # parms$parasite$k_w = theta[2]
  # if(stoch_scenarios$snails != "Absent")
  #   parms$snails$snail_transmission_rate = theta[3]
  
  #Run model
  # Will give results as object 'results'
  #Step 1
  pop <- init$humans$pop %>%
    mutate(Ind_sus = rgamma(nrow(cohort), shape = theta[2], scale = 1/theta[2]))
  N <- nrow(pop)
  ever_lived <- N
  SAC <- init$humans$SAC
  
  eggs_prev_SAC <- c(0)
  Heggs_prev_SAC <- c(0)
  inf_snail <- init$snails$I0
  susc_snail <- init$snails$S0
  exp_snail <- init$snails$E0
  snail_prev <- inf_snail/(susc_snail+exp_snail+inf_snail) 
  
  #Create initial cloud
  m_in <- sum(init$environment$contributions0)
  ## Start values
  newstart <- c(init$snails$S0, init$snails$E0, init$snails$I0, init$snails$C0) 
  cercariae <- 0 #to eventually plot the cercariae (cloud)
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
    
    if(births>0){
      pop <- bind_rows(pop, 
                       data.frame(age=0,
                                  sex=as.numeric(rbernoulli(births, 0.5)),
                                  jw1 = 0,
                                  Ind_sus = rgamma(births, shape = theta[2], scale = 1/theta[2]),
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
                                  ec = 0))
      ever_lived <- ever_lived+births
    }
    
    ########## DEATHS
    if(t %in% seq(1, (T*12), 12)){ #Januaries
      ag <- as.numeric(cut(pop$age, c(-1, prob_death$Age_hi))) #update age groups
      pop <- mutate(pop, age_group = ag)
      
      deaths <- as_tibble(table(ag)) %>%
        rename(Ag = ag) %>%
        mutate(n.deaths = round(prob_death[Ag, "Both_sexes"]*n)) %>%
        filter(n.deaths > 0)
      names(deaths$n.deaths) <- deaths$Ag
      
      deads <- stratified(pop[pop$age_group %in% deaths$Ag, ], "age_group", 
                          deaths$n.deaths, keep.rownames = T) %>%
        mutate(death_age = age + runif(n(), 1, 12)/12) #starts from 1/12, but it reaches 12/12, that is the following January
      pop[deads$rn,"death_age"] <- deads$death_age #we assign death age within that year
    }
    
    pop <- filter(pop, age < death_age)
    
    ########## MIGRATION
    #For now we do not account for age-specific emigration
    lambda <- round(parms$demography$emig_rate*(nrow(pop)/1000)/12)
    emig_age <- which(pop$age>5&pop$age<55)
    emigrated <- sample(emig_age, lambda)
    if(length(emigrated)>0)
      pop <- pop[-emigrated,]
    
    ########## UPDATE population size, SAC and cumulative exposure
    N[t] <- nrow(pop)
    SAC <- which(pop$age >= 5 & pop$age <= 15)
    
    cum_exp <- sum(Age_profile_exp(pop$age)$y*pop$Ind_sus) #Cumulative exposure
    
    ###########################################
    #Parasitology (vectorised)
    ###########################################
    
    ########### 1. CONTROL (No for grid search)
    
    ########### 2. EGGS production (before updating worms)
    #Worms' pairs (so mature at stage 4) produce eggs
    #(Shall we consider insemination probability??)
    
    #'mu' represents the expected egg load in a stool sample (41.7mg)
    #'eggs' is the daily egg amount passed to the environment
    # 24 is the conversion factor from egg counts and epg
    Tot_wp <- pop$wp1+pop$wp2+pop$wp3
    
    mu <- parms$parasite$eggs$alpha*Tot_wp*expon_reduction(parms$parasite$eggs$z, w=Tot_wp/2) 
    
    eggs <- round(mu*24*parms$parasite$eggs$gr_stool) #daily quantity, since particles in the environment are short-lived
    
    ########## 3. DIAGNOSIS 
    pop$ec = rnbinom(nrow(pop), size=parms$parasite$eggs$k_e, mu=mu)  
    
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
                      b = theta[3],
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
      ### SMT to replace NAs with zeros
      
      #Cercarial production 
      cercariae[t] <- out2$C[nrow(out2)]
    }
    
    ########## 6. EXPOSURE OF HUMANS 
    #Each individual assumes new worms (FOIh) and immunity applies
    #One exposure rate per individual (based on age and ind. sus.)
    pop$rate = FOIh(l=cercariae[t], exp(theta[1]), parms$parasite$v, a=pop$age, is=pop$Ind_sus)*expon_reduction(parms$immunity$imm, pop$cum_dwp)/cum_exp
    
    ########## 7. SAVE OUTPUT
    # sink("Find_bug.txt", append=TRUE)
    # cat(paste(Sys.time(), ": Seed:", k, "Time step", t, "res", cercariae[t], "\n", sep = " "))
    # sink()
    
    #Summary statistics
    # eggs_prev[t] <- length(which(pop$ec>0))/nrow(pop)
    eggs_prev_SAC[t] <- length(which(pop$ec[SAC]>0))/length(SAC)
    Heggs_prev_SAC[t] <- length(which((pop$ec[SAC]*24)>=400))/length(SAC)
    if(scen$snails != "Absent"){
      inf_snail <- out2$I[nrow(out2)] 
      susc_snail <- out2$S[nrow(out2)]
      exp_snail <- out2$E[nrow(out2)]
      snail_prev[t] <- inf_snail/(susc_snail+exp_snail+inf_snail) 
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
    if(t < parms$parasite$ext.foi$duration*12)
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
  
  res <- c(mean(eggs_prev_SAC[720:1200]), #any prev in SAC
           mean(Heggs_prev_SAC[720:1200]), #high intensity prev in SAC
           mean(snail_prev[720:1200])) #infection prevalence in snails
  
  return(res[1:2])
  # Save if needed
  # saveRDS(res, file = file.path(pop.output.dir, 
  #                               paste(setting, ".RDS", sep = "")))
  #then read with readRDS
}

#Check if it works
schisto_model()

###########################
#Desired output of vector 3
###########################
## Moderate endemicity: 30% prevalence in SAC, 6% HI prevalence in SAC, 6% snail prevalence
## Low endemicity: 10% prevalence in SAC, ~ 0% HI prevalence in SAC, 6% snail prevalence
## High endemicity: 60% prevalence in SAC, 20% HI prevalence in SAC, 6% snail prevalence

high_stat_obs = c(0.6, 0.2) #c(0.6, 0.2, 0.06)
#Weigths for summary statistics
weigths = c(0.7, 0.30)

###### Specific parameters #######
## to perform the Beaumont et al. (2009)'s method:
#Ex. desired max distance 0.01 (trying to back-calculate target tolerance from formula)
# sum(c(0.01, 0.01)^2/c(0.3, 0.2)^2) #0.0036
tolerance=3 #0.01 


ABC_Delmoral<-ABC_sequential(method="Delmoral", model=schisto_model, prior=theta_prior, 
                             #n_cluster = parallel::detectCores(logical = FALSE),
                             #dist_weights = weigths,
                             nb_simul=10, M=5, use_seed=FALSE, verbose = FALSE, progress_bar = TRUE,
                             summary_stat_target=high_stat_obs, tolerance_target=tolerance)
ABC_Delmoral




sd(ABC_Delmoral$stats[,1])
sd(ABC_Delmoral$stats[,2])

method="Beaumont"
model=schisto_model 
prior=theta_prior #inside_prior = FALSE,
nb_simul=10 
use_seed=TRUE 
verbose = TRUE 
progress_bar = TRUE
summary_stat_target=high_stat_obs 
tolerance_tab=tolerance


toy_model3 <- function(theta){
  #set.seed(1234)
  shape = theta[2]
  scale = 1/theta[2]
  theta[1] + rgamma(1, shape = shape, scale = scale)
}

schisto_model_cluster <- function(theta){
  
  library(parallel)
  library(tidyverse)
  library(deSolve)
  library(splitstackshape)
  library(EasyABC)
  library(readxl)
  
  #source.dir <- dirname(getActiveDocumentContext()$path)
  setwd("C:/Users/MAGE-NTD-07/Desktop/Federica/Schisto_model")
  
  #Load functions
  `%!in%` <- Negate(`%in%`)
  geom_mean <- function(x){exp(mean(log(x)))}
  ci <- function(x){quantile(x, probs=c(0.025, 0.975), na.rm = T)}
  
  #Functions (they can be a separate script)
  age_groups <- c(0, 5, 10, 16, 200)
  exposure_rates <- c(0.032, 0.61, 1, 0.06, 0.06) #Relative Age-specific exposure rates (activity/person/day) 
  #exposure_rates <- c(0.33, 0.44, 0.22, 0) #Relative Age-specific exposure rates (activity/person/day) 
  #Data from water contacts computed from Seydou S., De Vlas SJ, et al. (2011)
  Age_profile_exp <- function(a){
    approx(x=age_groups, y=exposure_rates, xout=c(a), method = "constant")
  }
  Age_profile_contr <- function(a){
    approx(x=c(0, 10, 100), y=c(1, 1, 1), xout=c(a), method = "linear")
  }
  FOIh <- function(l, zeta, v, a, is){
    Rel_ex <- Age_profile_exp(a)$y * is
    return(l * zeta * v * Rel_ex)
  }
  hyp_sat <- function(alpha, beta, w){ #Hyperbolic saturating function 
    f <- (alpha*w) / (1 + ((alpha*w) / beta))
    return(f)
  }
  expon_reduction <- function(alpha, w){
    f <- exp(-alpha*w)
    return(f)
  }

  #ODE system diefinition for snails' module
  SEI <- function(t, x, parms) {    
    with(as.list(c(parms, x)), {
      N_s <- S+E+I #is it better to work with densities?
      #Logistic growth
      #For population growing with limited amount of resources
      beta <- beta0*(1-N_s/k) #Infected snails do not reproduce
      FOIs <- b*mir 
      #l0*(1-chi^(mir/N)) #Gurarie - non linear FOIs
      
      #Equations
      dS <- beta*(S+E) - (v+FOIs)*S #susceptible
      dE <- FOIs*S - (v+tau)*E #Exposed: snails are invaded, but larvae are not patent yet. Thus, snails do not shed cercariae
      dI <- tau*E - v2*I  #Infected: larvae in the snail are mature and snails shed cercariae
      dC <- lambda*I - m*C #Cercariae (output)
      res <- c(dS, dE, dI, dC)
      list(res)
    })
  }
  ##Load parameters
  #The time step of the ABM is monthly
  parms <- list(#Demography
    demography = list(birth_rate = 36.5,
                      emig_rate = 18.6),
    #This is crude annual birth rate Uganda 2019 (same y of available life tables), 34.8 for Sub-Saharan Africa (per 1000 individuals)
    #-1.09 is the net migration rate for 2022 for Uganda (per 1000 individuals)
    #The emigration rate is calibrated to have constant population
    
    parasite = list(k_w = 0.15, #0.15 Anderson, Turner (2016) #can change for different settings (0.3 Sake) 
                    v = 1, #Transmission probability
                    zeta = 0.004, #overall exposure rate. (0.42 water contacts rate per day per individual, Seydou, De Vlas,.. 2011). (changing accordingly to endem. scenario)
                    ext.foi = list(value = 1, #monthly
                                   duration = 3), #years
                    Tw = 60, #Average worm's lifespan in host in months (months)(40 m Sake) (5 years for Anderson and May 1985a)
                    pp = 3, #Pre-patent period (months)
                    eggs = list(alpha = 0.14, #expected number of eggs per sample per worm pair (Sake 1996)
                                max = 100, #eggs plateau of hyp saturation in DDF
                                gr_stool = 150, #daily gr of stool produced by each human individual
                                z = 0.0007, #severity of density dependent fecundity
                                k_e = 0.87)), #aggregation parameter of egg counts detected (0.1 SCHISTOX; 0.87 Sake1992, but with three months interval and 25gr KK)
    #co_rate <- 1 #Average contribution rate (monthly) #to include seasonal patterns
    mda = list(age.lo = 5,
               age.hi = 15,
               start = 150, #70,
               end = 160, #80,
               frequency = 1, #annual
               coverage = 0.75,
               efficacy = 0.86), #Turner (2017)
    
    immunity = list(imm = 0), #immunity parameter
    
    ##Parameters for ODEs model for snails (These rates are daily rates)
    snails = list(max.reproduction.rate = 1, #0.1 d^-1 from Civitello DJ, 2022 #monthly is ~ 1 egg/day 
                  carrying.capacity = 20000, #arbitrary. To be estimated. #Civitello uses 5 L^-1 (about 30 per m3) 
                  mortality.rate = 1/100, #1/days of lifespan, Civitello #Gurarie: about 3 months
                  mortality.rate.infection = 1/30, #0.04 1/lifespan.infected, from Civitello. He works with additional mortality
                  infection.rate = 1/30, #1/lifespan of larvae within the snail, before shedding cercariae
                  snail_transmission_rate = 0.00005, #a combined version of exposure rate and probability of success. invasion
                  #rej.prob = 0.5 #probability of rejecting a miracidia, after getting in contact with the snail
                  # (1-chi)=0.5 for Civitello. OR it is for now computed from a Poisson as P(x=1)=0.8*exp(-0.8) using the infection rate from Anderson & May (1991).
                  cerc.prod.rate = 50, #1/d per infected snail
                  cerc.mortality = 1) #1/d 
  )
  
  parms$parasite$phi = 1-exp(-3/(parms$parasite$Tw - parms$parasite$pp)) #(monthly) proportion of adult worm pairs aging from age basket i to i+1, assuming 3 baskets 
  
  #Load initial conditions
  ## Load initial population 
  load("Equilibrium_age_distribution.RData") 
  prob_death <- read.csv("prob_death_Uganda_2019.csv")
  #Initial human cohort
  # cohort  <-  cohort %>%
  #   mutate(jw1 = 1, 
  #          Ind_sus = rgamma(nrow(cohort), shape = 0.15, scale = 1/0.15))
  
  ## Initial conditions
  init <- list(# Human cohort
    humans = list(pop = cbind(cohort, #%>%
                              jw1 = 1,
                              ID = c(1:nrow(cohort)),
                              death_age = 100,
                              age_group = as.numeric(cut(cohort$age, c(-1, prob_death$Age_hi))),
                              rate = 0,
                              wp1 = 0,
                              wp2 = 0,
                              wp3 = 0,
                              jw2 = 0,
                              jw3 = 0,
                              cum_dwp = 0,
                              ec = 0),
                  #Initial cumulative exposure
                  #cum_exp = sum(Age_profile_exp(cohort$age)$y * cohort$Ind_sus),
                  #Index of SAC in the cohort
                  SAC = which(cohort$age >= 5 & cohort$age <= 15)),
    #Environment
    environment = list(eggs0 = 0,
                       contributions0 = 0),
    #Snail population 
    snails = list(snail.pop=10000,
                  E0=0,
                  I0=0,
                  C0=0))
  init$snails$S0 = init$snails$snail.pop - sum(init$snails$E0 + init$snails$I0) 
  
  
  # Input (user choices)
  #Through this file the user can set the desired modelling scenario
  ################
  #Simulation settings
  ################
  T <- 100 #number of years simulated
  seeds <- 10
  fr <- 10 #frequency for printing to file the individual output [years]
  write.output <- FALSE #disable individual output for grid search/calibration (saving time)
  endem <- "high"
  
  ################
  #SETTING THE MODELLING SCENARIO: grid search
  ################
  #Combinations of modelling scenarios and stochastic seed
  
  #No mda
  parms$mda$start <- 0
  parms$mda$end <- 0
  
  #Grid search of endemic parameters
  stoch_scenarios <- expand.grid(list(#seed = 1:seeds,
    DDF_strength = c("Strong"), #Multiple combinations are possible
    imm_strength = c("Strong"), #Multiple combinations are possible
    snails = c("Absent"))) #Multiple combinations are possible
  ## The three below are the varying parameters
  # Need to assign values along the calibration
  # worms_aggr = 0.15,
  # zeta = 0.15,
  # tr_snails = 10^(-10)))
  
  
  #Load matched alphas for Density-dependent fecundity (DDF) given the endemicity
  load(paste("Matched_alphas_", endem, ".RData", sep = ""))
  
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
  
  # parms$parasite$zeta = exp(theta[2])
  # parms$parasite$k_w = theta[3]
  # if(stoch_scenarios$snails != "Absent")
  #   parms$snails$snail_transmission_rate = theta[4]
  
  #Stochastic seed
  set.seed(theta[1])
  
  #Run model
  # Will give results as object 'results'
  #Step 1
  pop <- init$humans$pop %>%
    mutate(Ind_sus = rgamma(nrow(cohort), shape = theta[3], scale = 1/theta[3]))
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
                                  Ind_sus = rgamma(births, shape = theta[3], scale = 1/theta[3]),
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
                      b = theta[4],
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
    pop$rate = FOIh(l=cercariae[t], exp(theta[2]), parms$parasite$v, a=pop$age, is=pop$Ind_sus)*expon_reduction(parms$immunity$imm, pop$cum_dwp)/cum_exp
    
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
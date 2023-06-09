T <- 300 #number of years simulated
seeds <- 20
fr <- 10 #frequency for printing to file the individual output [years]
write.output <- TRUE #disable individual output for grid search (saving time)

################
#SETTING THE MODELLING SCENARIO: limiting mechanism
################
#For each limiting mechanism, choose level: 'No', 'Mild', 'Strong'
#Combinations of modelling scenarios and stochastic seed

parms$mda <- list(age.lo = 5, #SAC is 5-15 #all population >= 2ys (WHO)
                 age.hi = 15,
                 start = 150,
                 end = 160,
                 frequency = 1, #annual
                 coverage = 0.75,
                 fr_excluded = 0.05, #systematic non-compliance 
                 efficacy = 0.86)
#Endemicity
endem <- "Moderate"

#Behavior in exposure
exposure = "Sow" #Choices: "ICL" (model-derived), "Sow" (water contacts)

stoch_scenarios <- expand.grid(list(seed = 1:seeds,
                                    DDF_strength = c("Absent", "Mild", "Strong"),
                                    imm_strength = c("Absent", "Mild", "Strong"),
                                    snails = c("Absent", "Mild", "Strong")))

#Load tuned transmission parameters (zetas) for the scenarios above (attention to the order!) 
#Transmission parameters for tuning endemicity

######
#IF LOW ENDEM
if(endem=="Low"){
  stoch_scenarios <- stoch_scenarios[(seeds*6+1):nrow(stoch_scenarios),] #61 or 301
  parms$parasite$ext.foi$value = 0.1 #0.1
  parms$parasite$ext.foi$duration = 0.25 #1/12 #years
}
######
#IF MODERATE ENDEM
if(endem=="Moderate"){
  parms$parasite$ext.foi$value = 1
  parms$parasite$ext.foi$duration = 2 #years
}
######
#######
#IF HIGH ENDEM
if(endem=="High"){
  parms$parasite$ext.foi$value = 5 #0.1
  parms$parasite$ext.foi$duration = 2
}
######

zetas <- read_excel("Zetas_Sow_func.xlsx") %>%
  filter(Endemicity == endem) # & Snails == "Strong" & Immunity == "Mild" & DDF != "Absent")
#Remove scenario with "all absent" because not at eq.
stoch_scenarios <- mutate(stoch_scenarios, 
                          zeta = rep(zetas$Zeta_grid_search, each = seeds),
                          worms_aggr = rep(zetas$Kw, each = seeds),
                          tr_snails = rep(zetas$`Transmission on snails`, each = seeds)) %>%
  filter(!(DDF_strength == "Absent" & snails == "Absent" & imm_strength == "Absent"))

#Load matched alphas for Density-dependent fecundity (DDF) given the endemicity
load(paste("Matched_alphas_", endem, ".RData", sep = ""))

#Specific parameters to be changed
# parms$snails$snail_transmission_rate = 5e-10
#parms$parasite$k_w = 0.1

setting <- paste(endem, "retuning_Sowfunc", sep = "_")

################
#Set output directory to save results
################
#This will be the directory where the individual output is automatically saved throughout the simulations
if(write.output == TRUE){
  ind.output.dir <- file.path(source.dir, paste("Output/Individual/", setting, sep = "")) 
  if(!file.exists(ind.output.dir)){
    dir.create(ind.output.dir) #Add check: this command to be run only if the directory is not existent
  }
  # if(file.exists(ind.output.dir)){
  #   #Empty the Output folder (only if needed)
  #   unlink(file.path(ind.output.dir, "/*")) 
  # }
}


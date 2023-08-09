T <- 300 #number of years simulated
seeds <- 30
fr <- 10 #frequency for printing to file the individual output [years]
write.output <- TRUE #disable individual output for grid search (saving time)

################
#SETTING THE MODELLING SCENARIO: limiting mechanism
################
#For each limiting mechanism, choose level: 'No', 'Mild', 'Strong'
#Combinations of modelling scenarios and stochastic seed

parms$mda <- list(age.lo = 2, #SAC is 5-15 #all population >= 2ys (WHO)
                 age.hi = 100,
                 start = 150,
                 end = 159,
                 frequency = 1, #annual
                 coverage = 0.75,
                 fr_excluded = 0.05, #systematic non-compliance 
                 efficacy = 0.86)
#Endemicity
#endem <- "High"

#Behavior in exposure
exposure = "Sow" #Choices: "ICL" (model-derived), "Sow" (water contacts)

stoch_scenarios <- expand.grid(list(seed = 1:seeds,
                                    DDF_strength = c("Absent", "Mild", "Strong"),
                                    imm_strength = c("Absent", "Mild", "Strong"),
                                    snails = c("Absent", "Mild", "Strong"),
                                    endem = c("Moderate", "High", "Low") #same order of the Zetas file
))

#Load tuned transmission parameters (zetas) for the scenarios above (attention to the order!) 
#Transmission parameters for tuning endemicity
zetas <- read_excel(paste("Zetas_", exposure, "_func.xlsx", sep = "")) 
#Remove scenario not at equilibrium
stoch_scenarios <- mutate(stoch_scenarios, 
                          zeta = rep(zetas$Zeta_grid_search, each = seeds),
                          worms_aggr = rep(zetas$Kw, each = seeds),
                          tr_snails = rep(zetas$`Transmission on snails`, each = seeds),
                          equilibrium = rep(zetas$Equilibrium, each = seeds)) %>%
  filter(equilibrium==TRUE) %>%
  mutate(Ext_foi_value = case_when(endem == "Low" ~ 0.1,
                                   endem == "Moderate" ~ 1,
                                   endem == "High" ~ 5),
         Ext_foi_duration = case_when(endem == "Low" ~ 0.5,
                                      endem == "Moderate" ~ 2,
                                      endem == "High" ~ 2))

setting <- paste(exposure, "func_commMDA", sep = "_")

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


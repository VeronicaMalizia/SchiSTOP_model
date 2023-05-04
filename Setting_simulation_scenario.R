################
#Simulation settings
################
T <- 100 #number of years simulated
seeds <- 10
fr <- 10 #frequency for printing to file the individual output [years]
write.output <- FALSE #disable individual output for grid search (saving time)

################
#SETTING THE MODELLING SCENARIO: grid search
################
#For each limiting mechanism, choose level: 'No', 'Mild', 'Strong'
#Combinations of modelling scenarios and stochastic seed

#No mda
parms$mda$start <- 0
parms$mda$end <- 0

#Grid search of endemic parameters
stoch_scenarios <- expand.grid(list(seed = 1:seeds,
                                    DDF_strength = c("Absent", "Strong"),
                                    imm_strength = c("Strong"),
                                    snails = c("Absent"),
                                    # zeta = exp(runif(50, -9.9, -9.5)),
                                    worms_aggr = c(0.05, seq(0.1, 0.5, 0.1)),
                                    zeta = seq(10^(-5), 0.5, 5*10^(-3)),
                                    #worms_aggr = c(0.1, 0.25),
                                    tr_snails = c(10^(-10), 2*10^(-10), 5*10^(-10), 10^(-9), 5*10^(-9))))



#transmission_snails <- seq(0.000001, 0.0001, 0.00001)

#Load tuned transmission parameters (zetas) for the scenarios above (attention to the order!) 
#Transmission parameters for tuning endemicity
# No need in case of grid search
# zetas <- read_excel("Zetas.xlsx")
# stoch_scenarios <- mutate(stoch_scenarios, zeta = rep(zetas$Zeta, each = seeds))

#Load matched alphas for Density-dependent fecundity (DDF) given the endemicity
load("Matched_alphas_moderate.RData")

#Heterogeneity of worms
#parms$parasite$k_w = 0.3

################
#Set output directory to save results
################
setting <- "Grid_search_Panel3_watercontact"
#"Moderate_setting_k=030_Age-int"

#This will be the directory where the individual output is automatically saved throughout the simulations
if(write.output == TRUE){
  ind.output.dir <- file.path(source.dir, paste("Output/Individual/", setting, sep = "")) 
  dir.create(ind.output.dir) #Add check: this command to be run only if the directory is not existent
  #Empty the Output folder (only if needed)
  #unlink(file.path(output.dir, "/*")) 
}
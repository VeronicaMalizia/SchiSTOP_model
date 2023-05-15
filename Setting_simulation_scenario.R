################
#Simulation settings
################
T <- 100 #number of years simulated
seeds <- 10
fr <- 10 #frequency for printing to file the individual output [years]
write.output <- FALSE #disable individual output for grid search/calibration (saving time)
endem <- "High"

################
#SETTING THE MODELLING SCENARIO: grid search
################
#Combinations of modelling scenarios and stochastic seed

#No mda
parms$mda$start <- 0
parms$mda$end <- 0

#Grid search of endemic parameters
stoch_scenarios <- expand.grid(list(seed = 1:seeds,
                                    DDF_strength = c("Strong"), #Multiple combinations are possible
                                    imm_strength = c("Strong"), #Multiple combinations are possible
                                    snails = c("Mild"), #Multiple combinations are possible
                                    ## The three below are the varying parameters
                                    # Need to assign values along the calibration
                                    worms_aggr = 0.15,
                                    zeta = 0.15,
                                    tr_snails = 10^(-10)))


#Load matched alphas for Density-dependent fecundity (DDF) given the endemicity
load(paste("Matched_alphas_", endem, ".RData", sep = ""))

################
#Set output directory to save results
################
#The name will be given to saved file
setting <- "XXX"

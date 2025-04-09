#############################
#Author: Veronica Malizia
#Date: 01/07/2022
#R version: 4.1.2
#
# Customize simulation scenario-specific parameters. 
# This will prepare the setting to run simulations.
# Species of interest: Schistosoma mansoni
#############################

T <- 300 #number of simulated years #300 for paper
seeds <- 100 #stochastic seeds #100 for paper
fr <- 10 #frequency for printing to file the individual output [years] - this is to save memory
write.output <- FALSE #disable individual output if not needed (saving time and memory)

################
#SETTING THE MODELLING SCENARIO: combinations of regulating mechanisms and age-exposure function
################
#For each regulating mechanism, choose level: 'No', 'Mild', 'Strong'
#Create combinations of modelling scenarios and stochastic seed

# Set mda target population as desired
# Two strategies for the manuscript: #SAC (5-15ys) #community-wide >= 2ys (WHO)
parms$mda <- list(age.lo = 5, 
                 age.hi = 15,
                 start = 150,
                 end = 159,
                 frequency = 1, #annual
                 coverage = 0.75,
                 fr_excluded = 0.05, #systematic non-compliance 
                 efficacy = 0.86)

#Age-specific behavior in exposure
# Two choices:
# 1. model-derived from Toor2018 and Turner2017 
# 2. based on water-contacts here estimated, based on Sow2011 and Fulford1996
exposure = "Sow" #Choices: "ICL" (model-derived), "Sow" (water contacts)

stoch_scenarios <- expand.grid(list(seed = 1:seeds, ## For simplicity, DO NOT change the order here if you do not want to simulate certain scenarios, but apply this change at line 55
                                    DDF_strength = c("Absent", "Mild", "Strong"),
                                    imm_strength = c("Absent", "Mild", "Strong"),
                                    snails = c("Absent", "Mild", "Strong"),
                                    endem = c("Moderate", "High", "Low") #MUST be the same order as in the Zetas file
))

#Load tuned transmission parameters for all the scenarios above (order important!) 
zetas <- read_excel(paste("Zetas_", exposure, "_func.xlsx", sep = "")) 
#Remove those scenarios which do not reach the endemic equilibrium at pre-control
#  this is optional to save memory and computing time
stoch_scenarios <- mutate(stoch_scenarios, 
                          zeta = rep(zetas$Zeta_grid_search, each = seeds),
                          worms_aggr = rep(zetas$Kw, each = seeds),
                          tr_snails = rep(zetas$`Transmission on snails`, each = seeds),
                          equilibrium = rep(zetas$Equilibrium, each = seeds)) %>%
## Filter HERE the scenarios you want to simulate
  filter(equilibrium==TRUE) %>%
  mutate(Ext_foi_value = case_when(endem == "Low" ~ 0.5,
                                   endem == "Moderate" ~ 1,
                                   endem == "High" ~ 5),
         Ext_foi_duration = case_when(endem == "Low" ~ 0.5,
                                      endem == "Moderate" ~ 2,
                                      endem == "High" ~ 2))

# Define the title of the simulation setting. This will give the name to all the output files.
setting <- paste(exposure, "func_SACmda", sep = "_") 

################
#Set output directory to save individual results
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


T <- 200 #number of years simulated
seeds <- 10
fr <- 10 #frequency for printing to file the individual output [years]
write.output <- TRUE #disable individual output for grid search (saving time)

################
#SETTING THE MODELLING SCENARIO: limiting mechanism
################
#For each limiting mechanism, choose level: 'No', 'Mild', 'Strong'
#Combinations of modelling scenarios and stochastic seed

parms$mda <- list(age.lo = 5, #SAC is 5-15 #all population >= 2ys (WHO)
                 age.hi = 15,
                 start = 150, #70,
                 end = 160, #80,
                 frequency = 1, #annual
                 coverage = 0.75,
                 efficacy = 0.86)

endem <- "High"
stoch_scenarios <- expand.grid(list(seed = 1:seeds,
                                    DDF_strength = c("Absent", "Mild", "Strong"),
                                    imm_strength = c("Absent", "Mild", "Strong"),
                                    snails = c("Absent", "Mild", "Strong")))

#Load tuned transmission parameters (zetas) for the scenarios above (attention to the order!) 
#Transmission parameters for tuning endemicity

######
#IF LOW ENDEM
if(endem=="Low"){
  stoch_scenarios <- stoch_scenarios[301:nrow(stoch_scenarios),]
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

stoch_scenarios <- filter(stoch_scenarios, #DDF_strength != "Absent" &
                            snails == "Mild" & imm_strength == "Strong")
zetas <- read_excel("Zetas_new.xlsx") %>%
  filter(Endemicity == endem & Snails == "Mild" & Immunity == "Strong") # & DDF != "Absent")
stoch_scenarios <- mutate(stoch_scenarios, 
                          zeta = rep(zetas$Zeta_grid_search, each = seeds),
                          worms_aggr = rep(zetas$Kw, each = seeds),
                          tr_snails = rep(zetas$`Transmission on snails`, each = seeds))

#Load matched alphas for Density-dependent fecundity (DDF) given the endemicity
load(paste("Matched_alphas_", endem, ".RData", sep = ""))

#Specific parameters to be changed
# parms$snails$snail_transmission_rate = 5e-10
#parms$parasite$k_w = 0.1

setting <- paste(endem, "watercontacts", sep = "_")


#############################
#Author: Veronica Malizia
#Date: 20/06/2022
#R version: 4.1.2
#
#Complete list of parameters to prepare and run the simulations
#Species of interest: Schistosoma mansoni
#############################

#The time step of the ABM is monthly
parms <- list(#Demography
              demography = list(birth_rate = 36.5,
                                emig_rate = 18.6),
              #This is crude annual birth rate Uganda 2019 (same y of available life tables), 34.8 for Sub-Saharan Africa (per 1000 individuals)
              #-1.09 is the net migration rate for 2022 for Uganda (per 1000 individuals)
              #The emigration rate is calibrated in script 00_Demography.R to have constant population
              
              parasite = list(k_w = 0.15, #0.15 Anderson, Turner (2016) #can change for different settings (0.3 de Vlas (1992)) 
                              zeta = 0.004, #overall exposure rate. This is a varying parameter. (0.42 water contacts rate per day per individual according to Seydou, De Vlas,.. 2011).
                               ext.foi = list(value = 1, #monthly
                                              duration = 2), #years
                               Tw = 60, #Average worm's lifespan in host in months (months)(40mo de Vlas (1996)) (5 years for Anderson and May 1985a)
                               pp = 3, #Pre-patent period (months)
                               eggs = list(alpha = 0.14, #expected number of eggs per sample per worm pair (de Vlas 1996)
                                            max = 100, #eggs plateau if using hyperbolic saturating function for DDF
                                            gr_stool = 150, #daily gr of stool produced by each human individual
                                            z = 0.0007, #severity of density-dependent fecundity
                                            n.slides = 1, #number of slides per stool sample
                                            k_e = 0.87)), #aggregation parameter of egg counts detected (0.1 SCHISTOX; 0.87 de Vlas1992, but with three months interval and 25gr KK)
                                            #co_rate <- 1 #Average contribution rate (monthly) #can be improved to include seasonal patterns
              exposure = list(#model-derived from Toor2018 and Turner2017
                              ICL_derived = list(ages = c(0, 5, 10, 16, 200),
                                                 exp = c(0.032, 0.61, 1, 0.06, 0.06),
                                                 method = "constant"), #Relative Age-specific exposure rates (activity/person/day))
                              #based on water-contacts here estimated, based on Sow2011 and Fulford1996
                              Sow_derived = list(ages = c(0, 5, 15, 40, 100),
                                                 exp = c(0, 0.62, 1, 0.51, 0.51),
                                                 method = "linear")),
              mda = list(age.lo = 5,
                        age.hi = 15,
                        start = 150, 
                        end = 160, 
                        frequency = 1, #annual
                        coverage = 0.75,
                        fr_excluded = 0.05, #systematic non-compliance                       
                        efficacy = 0.86), #Turner (2017)
              
              immunity = list(imm = 0), #immunity parameter. Varying across models.
              
              ##Parameters for ODEs model for snails (These rates are daily rates)
              snails = list(max.reproduction.rate = 1, #0.1 d^-1 from Civitello DJ, 2022 #monthly is ~ 1 egg/day 
                            carrying.capacity = 20000, #arbitrary. Varying across models. #Civitello et al. use 5 L^-1 (about 30 per m3) 
                            mortality.rate = 1/100, #1/days of lifespan, Civitello and Gurarie: about 3 months
                            mortality.rate.infection = 1/30, #0.04 1/lifespan.infected, from Civitello. They work with additional mortality
                            infection.rate = 1/30, #1/lifespan of larvae within the snail, before shedding cercariae
                            snail_transmission_rate = 1e-09, #a combined version of exposure rate and probability of successful invasion
                            #rej.prob = 0.5 #probability of rejecting a miracidia, after getting in contact with the snail. Here not in use.
                            # (1-chi)=0.5 for Civitello. OR it is for now computed from a Poisson as P(x=1)=0.8*exp(-0.8) using the infection rate from Anderson & May (1991).
                            cerc.prod.rate = 50, #1/d per infected snail
                            cerc.mortality = 1) #1/d 
)

parms$parasite$phi = 1-exp(-3/(parms$parasite$Tw - parms$parasite$pp)) #(monthly) proportion of adult worm pairs aging from age basket i to i+1, assuming 3 baskets for the Erlang distribution. 

# age_groups <- c(0, 5, 10, 16, 200)
# exposure_rates <- c(0.032, 0.61, 1, 0.06, 0.06) #Relative Age-specific exposure rates (activity/person/day) 
# Moderate exposure in adults

#exposure_rates <- c(0.33, 0.44, 0.22, 0) #Relative Age-specific exposure rates (activity/person/day) 
#Data from water contacts computed from Seydou S., De Vlas SJ, et al. (2011)

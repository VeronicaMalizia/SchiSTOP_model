#############################
#Author: Veronica Malizia
#Date: 10/03/2022
#R version: 3.6.1

#Snails are infected with a FOIs: human-to-snail
#This modules implements the SEI model for freshwater snail population
#The model provides as outcome the FOIh: snail-to-human

#Daily timesteps for snail population dynamics

#This script is explorative for snails dynamics. In the main Schisto model it is already included in line
#############################
rm(list = ls())

library(deSolve)
SEI <- function(t, x, parms) {    
    with(as.list(c(parms, x)), {
        N <- S+E+I #is better to work with densities?
        #Logistic growth
        #For population growing with limited amount of resources
        beta <- beta0*(1-N/k) #Infected snails do not reproduce

        #Equations
        dS <- beta*(S+E) - (v+FOIs/N)*S #susceptible
        dE <- (FOIs/N)*S - (v+tau)*E #Exposed: snails are invaded, but larvae are not patent yet. Thus, snails do not shed cercariae
        dI <- tau*E - (v+v2)*I  #Infected: larvae in the snail are mature and snails shed cercariae
        res <- c(dS, dE, dI)
        list(res)
    })
}


## Parameters
#ENV = 500 #L, total volume of water
max.reproduction.rate = 0.1 #d^-1 from Civitello DJ, 2022 #monthly is ~ 1 egg/day 
carrying.capacity = 5000 #arbitrary. To be estimated. #Civitello uses 5 L^-1 (about 30 per m3) 
lifespan = 100 #days, Civitello #Gurarie: about 3 months
mortality.rate = 1/lifespan 
#lifespan.reduction=0.8 #arbitrary. Still not enough evidence found.
mortality.rate.infection = 0.04 #d^-1 from Civitello #1/((1-lifespan.reduction)*lifespan)
mu = 30 #1 month: lifespan of larvae within the snail, before shedding cercariae
infection.rate = 1/mu 
#c = 2 #contact rate for snails. This should also be a calibrating parameter. 
#For now we assume that snails get in contact with all free-living miracidiae. I think it cannot be >1
chi = 0.5 #probability of a successful invasion for a single miracidia getting in contact with the host
# 0.5 is Civitello. OR it is for now computed from a Poisson as P(x=1)=0.8*exp(-0.8) using the infection rate from Anderson & May (1991).
# However, I would consider it arbitrary too and then to be estimated. (Or look for data)
miracidiae = 1000
l0 = 1/15 #1/d. Rate of sporocyst development in snails, given successful invasion.
cerc.prod.rate = 50 #1/d per infected snail
cerc.mortality = 1 #1/d

FOIs= chi*miracidiae #c*chi*miracidiae #Civitello uses 0.01
parms  <- c(beta0 = max.reproduction.rate, k = carrying.capacity, v = mortality.rate,
            FOIs = FOIs, v2 = mortality.rate.infection, tau = infection.rate)

## vector of time steps
#I have to stick to monthly time step as in the main module.
ndays=180
times <- 1:ndays

## initial conditions
#M0 = 100
pop.size=10000
E0=0
I0=0
S0=pop.size - sum(E0, I0)
#If densities: S0=1

## Start values for steady state
xstart <- c(S = S0, E = E0, I = I0)

## Solving
out <-  lsoda(xstart, times, SEI, parms) #tra gli argomenti ha la mia funzione che definisce il sistema differenziale

#Solving iteratively for different carring capacities
# ks <- c(1000, 5000, 10000)
# for(i in 1:3){
#   parms  <- c(beta0 = max.reproduction.rate, k = ks[i], v = mortality.rate,
#               FOIs = FOIs, v2 = mortality.rate.infection, tau = infection.rate)
#   assign(paste("out_", ks[i], sep = ""),
#   as.data.frame(lsoda(xstart, times, SEI, parms)))
# }


## Translate the output into a data.frame
out2 <- as.data.frame(out)

## Plotting
yname <- paste("Density of snails", expression(L^(-1)), sep=" ")
TOT <- out2$S + out2$E + out2$I
plot(out2$time, out2$S, col='blue', ylim=c(0, max(out[,-1])), type='l', xlab = "Time in days",
     ylab = yname) #Susceptibles
lines(out2$time, out2$E, col='dark green') #Exposed (infected but not shedding larvae)
lines(out2$time, out2$I, col='red') #Infected (larvae are mature and snails excrete cercariae)
#str(out2)
legend("topright", 
       legend = c("Scusceptibles", "Exposed", "Infected"), 
       col = c("blue", "dark green", "red"), 
       pch = 20, 
       bty = "n", 
       pt.cex = 2, 
       cex = 1.2, 
       text.col = "black", 
       horiz = F , 
       inset = c(0.1, 0.1))

plot(out2$time, TOT, type='l', 
     xlab = "Time in days", ylab = "Total snail population size") 

#Plotting prevalence
plot(out2$time, out2$I/TOT, ylim=c(0, 1), type='l', 
     xlab = "Time in days", ylab = "Prevalence of infection in snails") 
lines(out2$time, (out2$I+out2$E)/out2$N, ylim=c(0, 1), type='l', col='violet',
     xlab = "Time in days", ylab = "Prevalence of infection in snails")
legend("topright", 
       legend = c("I/N", "(E+I)/N"), 
       col = c("black", "violet"), 
       pch = 20, 
       bty = "n", 
       pt.cex = 2, 
       cex = 1.2, 
       text.col = "black", 
       horiz = F , 
       inset = c(0.1, 0.1))


# plot(out_10000$time, out_10000$I, col='purple', type='l', xlab = "Time in days",
#      ylab = yname) #Susceptibles
# lines(out_1000$time, out_1000$I, col='blue', type='l', xlab = "Time in days",
#       ylab = yname) #Susceptibles
# lines(out_5000$time, out_5000$I, col='red', type='l', xlab = "Time in days",
#       ylab = yname) #Susceptibles

#############################
#Author: Veronica Malizia
#Date: 10/03/2022
#R version: 3.6.1

#Snails are infecte with a FOIs: human-to-snail
#This modules implements the SEI model for freshwater snail population
#The model provides as outcome the FOIh: snail-to-human

#Daily timesteps for snail population dynamics
#############################
rm(list = ls())

library(deSolve)
SEI <- function(t, x, parms) {    
    with(as.list(c(parms, x)), {
        N <- S+E+I #is better to work with densities?
        #Logistic growth
        #For population growing with limited amount of resources
        beta <- beta0*(1-N/k) #Infected snails do not reproduce
        dS <- beta*(S+E) - (v+FOIs)*S #susceptible
        dE <- FOIs*S - (v+tau)*E #Exposed: snails are invaded, but larvae are not patent yet. Thus, snails do not shed cercariae
        dI <- tau*E - (v+v2)*I  #Infected: larvae in the snail are mature and snails shed cercariae
        res <- c(dS, dE, dI)
        list(res)
    })
}


## Parameters
ENV = 500 #L, total volume of water
max.reproduction.rate = 0.1 #d^-1 from Civitello DJ, 2022 #monthly is ~ 1 egg/day 
carrying.capacity = 5 #Civitello, L^-1 per m3 #arbitrary. To be estimated.
lifespan = 100 #days, Civitello #Another source: 15-18 months in natural conditions (sure??)
mortality.rate = 1/lifespan 
#lifespan.reduction=0.8 #arbitrary. Still not enough evidence found.
mortality.rate.infection = 0.04 #d^-1 from Civitello #1/((1-lifespan.reduction)*lifespan)
mu = 30 #1 month: lifespan of larvae within the snail, before shedding cercariae
infection.rate = 1/mu 
c = 2 #contact rate for snails. This should also be a calibrating parameter. 
#For now we assume that snails get in contact with all free-living miracidiae. I think it cannot be >1
chi = 0.5 #probability of a successful invasion for a single miracidia getting in contact with the host
# It is for now computed from a Poisson as P(x=1)=0.8*exp(-0.8) using the infection rate from Anderson & May (1991).
# However, I would consider it arbitrary too and then to be estimated. (Or look for data)
miracidiae = 1 #if we work with densities, we will need to define miracidiae with 1/L
inf.contact.rate = c*chi
FOIs= 0.01 #c*chi*miracidiae
parms  <- c(beta0 = max.reproduction.rate, k = carrying.capacity, v = mortality.rate,
            FOIs = FOIs, v2 = mortality.rate.infection, tau = infection.rate)

## vector of timesteps
#I have to stick to monthly timestep as in the main module.
nmonths=36 #3 years
ndays=180
times <- seq(0, nmonths, length = 30*nmonths)
times <- 1:ndays

## initial conditions
#M0 = 100
pop.size=1
E0=0
I0=0
S0=pop.size - sum(E0, I0)
#If densities: S0=1

## Start values for steady state
xstart <- c(S = S0, E = E0, I = I0)

## Solving
out <-  lsoda(xstart, times, SEI, parms) #tra gli argomenti ha la mia funzione che definisce il sistema differenziale

## Translate the output into a data.frame
out2 <- as.data.frame(out)

## Plotting
yname <- paste("Density of snails", expression(L^(-1)), sep=" ")
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

out2$N = out2$S+out2$E+out2$I
plot(out2$time, out2$I/out2$N, ylim=c(0, 1), type='l', 
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

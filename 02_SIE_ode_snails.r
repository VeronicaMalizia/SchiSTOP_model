#############################
#Author: Veronica Malizia
#Date: 10/03/2022
#R version: 3.6.1

#Snails are infecte with a FOIs: human-to-snail
#This modules implements the SEI model for freshwater snail population
#The model provides as outcome the FOIh: snail-to-human
#############################

library(deSolve)
SIE <- function(t, x, parms) {    
    with(as.list(c(parms, x)), {
        N <- S+E+I #is better to work with densities?
        #Logistic growth
        #For population growing with limited amount of resources
        beta <- beta0*(1-N/k)*(S+E) #Infected snails do not reproduce
        dS <- beta - (v+FOIs)*S #susceptible
        dE <- FOIs*S - (v+tau)*E #Exposed: snails are invaded, but larvae are not patent yet. Thus, snails do not shed cercariae
        dI <- tau*E - (v+v2)*I  #Infected: larvae in the snail are mature and snails shed cercariae
        res <- c(dS, dE, dI)
        list(res)
    })
}


## Parameters
max.reproduction.rate=30 #monthly # ~ 1 egg/day 
carrying.capacity=1000 #arbitrary. To be estimated.
lifespan=16 #15-18 months in natural conditions
mortality.rate=0.3 #1/lifespan 
lifespan.reduction=0.8 #arbitrary. Still not enough evidence found.
mortality.rate.infection=1.2 #1/((1-lifespan.reduction)*lifespan)
mu=1 #1 month: lifespan of larvae within the snail, before shedding cercariae
infection.rate=1/mu 
c=0.3 #exposure rate for snails. This should also be a calibrating parameter. 
#For now we assume that snails get in contact with all free-living miracidiae. I think it cannot be >1
chi=0.4 #probability of a successful invasion for a single miracidia getting in contact with the host
# It is for now computed from a Poisson as P(x=1)=0.8*exp(-0.8) using the infection rate from Anderson & May (1991).
# However, I would consider it arbitrary too and then to be estimated. (Or look for data)
miracidiae=176
FOIs=c*chi*miracidiae
parms  <- c(beta0 = max.reproduction.rate, k = carrying.capacity, v = mortality.rate,
            FOIs = FOIs, v2 = mortality.rate.infection, tau = infection.rate)

## vector of timesteps
#I have to stick to monthly timestep as in the main module.
nmonths=36 #3 years
times <- seq(0, nmonths, length = 30*nmonths)

## initial conditions
pop.size=500
E0=1
I0=1
S0=pop.size - 2

## Start values for steady state
xstart <- c(S = S0, E = E0, I = I0)

## Solving
out <-  lsoda(xstart, times, SIE, parms) #tra gli argomenti ha la mia funzione che definisce il sistema differenziale

## Translate the output into a data.frame
out2 <- as.data.frame(out)

## Plotting
plot(out2$time, out2$S, col='blue', ylim=c(0, max(out[,-1])), type='l') #Susceptibles
lines(out2$time, out2$E, col='dark green') #Exposed (infected but not shedding larvae)
lines(out2$time, out2$I, col='red') #Infected (larvae are mature and snails excrete cercariae)
str(out2)

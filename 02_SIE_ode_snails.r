#############################
#Author: Veronica Malizia
#Date: 10/03/2022
#R version: 3.6.1

#Snails are infected with a FOIs: human-to-snail
#This modules implements the SEI model for freshwater snail population
#The model provides as outcome the FOIh: snail-to-human

#Daily time steps for snail population dynamics

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
        FOIs <- l0*(1-chi^(mir/N))
        
        #Equations
        dS <- beta*(S+E) - (v+FOIs)*S #susceptible
        dE <- FOIs*S - (v+tau)*E #Exposed: snails are invaded, but larvae are not patent yet. Thus, snails do not shed cercariae
        dI <- tau*E - v2*I  #Infected: larvae in the snail are mature and snails shed cercariae
        dC <- lambda*I - m*C #Cercariae (output)
        res <- c(dS, dE, dI, dC)
        list(res)
    })
}

#Checks:
## 1. N is not supposed go asymptotically move towards k, because respectively to a normal logistic growth function (beta*(1-N/k)*N)
##    here only S + E contribute to the reproduction. Moreover, I have an extra mortality factor.
##    The population dynamics are different than a normal logistic growth population.
##
## 2. Is (and if so why) the snail infection prevalence affected by the choice of k? (see plotting code at the end of the document)

## Parameters
#ENV = 500 #L, total volume of water
max.reproduction.rate = 1 #0.1 d^-1 from Civitello DJ, 2022 #monthly is ~ 1 egg/day 
carrying.capacity = 20000 #arbitrary. To be estimated. #Civitello uses 5 L^-1 (about 30 per m3) 
lifespan = 100 #days, Civitello #Gurarie: about 3 months
lifespan.infected = 30 #days, Gurarie
mortality.rate = 1/lifespan 
#lifespan.reduction=0.8 #arbitrary. Still not enough evidence found.
mortality.rate.infection = 1/lifespan.infected #0.04 d^-1 from Civitello. He works with additional mortality
mu = 30 #1 month: lifespan of larvae within the snail, before shedding cercariae
infection.rate = 1/mu 
rej.prob = 0.5 #probability of rejecting a miracidia, after getting in contact with the snail
# (1-chi)=0.5 for Civitello. OR it is for now computed from a Poisson as P(x=1)=0.8*exp(-0.8) using the infection rate from Anderson & May (1991).
# However, I would consider it arbitrary too and then to be estimated. (Or look for data)
max.invasion = 1/15 #1/d. Rate of sporocyst development in snails, given successful invasion.
cerc.prod.rate = 50 #1/d per infected snail
cerc.mortality = 1 #1/d 

mirac.input = 30000 #chi*miracidiae will be divided by N[t] in the system #Civitello uses a cumulative factor of 0.01
parms  <- c(beta0 = max.reproduction.rate, k = carrying.capacity, v = mortality.rate,
            mir = mirac.input, l0 = max.invasion, chi = rej.prob, 
            v2 = mortality.rate.infection, tau = infection.rate,
            lambda = cerc.prod.rate, m = cerc.mortality)

## vector of time steps
#I have to stick to monthly time step as in the main module.
ndays=100
times <- 1:ndays

## initial conditions
#M0 = 100
snail.pop=10000
E0=0
I0=0
S0=snail.pop - sum(E0, I0)
C0=0
#If densities: S0=1

## Start values for steady state
xstart <- c(S = S0, E = E0, I = I0, C = C0)

## Solving single run
# solve ODEs
out <-  lsoda(xstart, times, SEI, parms) #tra gli argomenti ha la mia funzione che definisce il sistema differenziale
# Translate the output into a data.frame
out2 <- as.data.frame(out)

#Solving iteratively for different carring capacities
ks <- c(1000, 5000, 10000, 20000)
for(i in 1:length(ks)){
  parms  <- c(beta0 = max.reproduction.rate, k = ks[i], v = mortality.rate,
              FOIs = FOIs, v2 = mortality.rate.infection, tau = infection.rate,
              lambda = cerc.prod.rate, m = cerc.mortality)
  assign(paste("out_", ks[i], sep = ""),
  as.data.frame(lsoda(xstart, times, SEI, parms)))
}

## Plotting
TOT <- out2$S + out2$E + out2$I
yname <- "Abundances of snails"
plot(out2$time, out2$S, col='blue', ylim=c(0, max(out[,2:4])), type='l', xlab = "Time in days",
     ylab = yname) #Susceptibles
lines(out2$time, out2$E, col='dark green') #Exposed (infected but not shedding larvae)
lines(out2$time, out2$I, col='red') #Infected (larvae are mature and snails excrete cercariae)
#str(out2)
legend("topright", 
       legend = c("Susceptibles", "Exposed", "Infected"), 
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
plot(out2$time, out2$C, type='l', 
     xlab = "Time in days", ylab = "Cercarial output") 

#Plotting prevalence
plot(out2$time, out2$I/TOT, ylim=c(0, 1), type='l', 
     xlab = "Time in days", ylab = "Prevalence of infection in snails") 
lines(out2$time, (out2$I+out2$E)/TOT, ylim=c(0, 1), type='l', col='violet',
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


#Plotting if we run for different values of k
#Infected
plot(out_20000$time, out_20000$I, col='dark green', type='l', xlab = "Time in days",
     ylab = "Infected snail population") 
lines(out_10000$time, out_10000$I, col='purple', type='l', xlab = "Time in days",
     ylab = "Infected snail population") 
lines(out_1000$time, out_1000$I, col='blue', type='l', xlab = "Time in days",
      ylab = "Infected snail population")
lines(out_5000$time, out_5000$I, col='red', type='l', xlab = "Time in days",
      ylab = "Infected snail population") 
legend("topright", 
       legend = c("k=20000", "k=10000", "k=5000", "k=1000"), 
       col = c("dark green", "purple", "red", "blue"), 
       pch = 20, 
       bty = "n", 
       pt.cex = 2, 
       cex = 1.2, 
       text.col = "black", 
       horiz = F , 
       inset = c(0.1, 0.1))

#Total
plot(out_20000$time, (out_20000$S + out_20000$E + out_20000$I), col='dark green', type='l', xlab = "Time in days",
     ylab = "Total snail population", ylim = c(0, 20000)) 
lines(out_10000$time, (out_10000$S + out_10000$E + out_10000$I), col='purple', type='l', xlab = "Time in days",
     ylab = "Total snail population") 
lines(out_1000$time, (out_1000$S + out_1000$E + out_1000$I), col='blue', type='l', xlab = "Time in days",
      ylab = "Total snail population") #Susceptibles
lines(out_5000$time, (out_5000$S + out_5000$E + out_5000$I), col='red', type='l', xlab = "Time in days",
      ylab = "Total snail population") #Susceptibles
legend("topright", 
       legend = c("k=20000", "k=10000", "k=5000", "k=1000"), 
       col = c("dark green", "purple", "red", "blue"), 
       pch = 20, 
       bty = "n", 
       pt.cex = 2, 
       cex = 1.2, 
       text.col = "black", 
       horiz = F , 
       inset = c(0.1, 0.1))


#Prevalence
plot(out_20000$time, out_20000$I/(out_20000$S + out_20000$E + out_20000$I), col='dark green', type='l', xlab = "Time in days",
     ylab = "Prevalence of infected snails", ylim = c(0, 1)) 
lines(out_10000$time, out_10000$I/(out_10000$S + out_10000$E + out_10000$I), col='purple', type='l', xlab = "Time in days",
     ylab = "Prevalence of infected snails", ylim = c(0, 1)) 
lines(out_1000$time, out_1000$I/(out_1000$S + out_1000$E + out_1000$I), col='blue', type='l', xlab = "Time in days",
      ylab = "Prevalence of infected snails") 
lines(out_5000$time, out_5000$I/(out_5000$S + out_5000$E + out_5000$I), col='red', type='l', xlab = "Time in days",
      ylab = "Prevalence of infected snails") 
legend("topright", 
       legend = c("k=20000", "k=10000", "k=5000", "k=1000"), 
       col = c("dark green", "purple", "red", "blue"), 
       pch = 20, 
       bty = "n", 
       pt.cex = 2, 
       cex = 1.2, 
       text.col = "black", 
       horiz = F , 
       inset = c(0.1, 0.1))

#Cercariae
plot(out_20000$time, out_20000$C, col='dark green', type='l', xlab = "Time in days",
     ylab = "Cercarial output") 
lines(out_10000$time, out_10000$C, col='purple', type='l', xlab = "Time in days",
      ylab = "Cercarial output") 
lines(out_1000$time, out_1000$C, col='blue', type='l', xlab = "Time in days",
      ylab = "Cercarial output") 
lines(out_5000$time, out_5000$C, col='red', type='l', xlab = "Time in days",
      ylab = "Cercarial output") 
legend("topright", 
       legend = c("k=20000", "k=10000", "k=5000", "k=1000"), 
       col = c("dark green", "purple", "red", "blue"), 
       pch = 20, 
       bty = "n", 
       pt.cex = 2, 
       cex = 1.2, 
       text.col = "black", 
       horiz = F , 
       inset = c(0.1, 0.1))


final_equil <- data.frame(k = ks,
                          S = c(out_1000[180, "S"], out_5000[180, "S"], out_10000[180, "S"], out_20000[180, "S"]),
                          E = c(out_1000[180, "E"], out_5000[180, "E"], out_10000[180, "E"], out_20000[180, "E"]),
                          I = c(out_1000[180, "I"], out_5000[180, "I"], out_10000[180, "I"], out_20000[180, "I"]),
                          C = c(out_1000[180, "C"], out_5000[180, "C"], out_10000[180, "C"], out_20000[180, "C"]))

ggplot(final_equil, aes(x = k)) +
  geom_line(aes(y=S), colour = "green", size=2) +
  geom_line(aes(y=E), colour = "orange", size=2) +
  geom_line(aes(y=I), colour = "red", size=2) +
  geom_line(aes(y=E+I+S), size=2) +
  scale_x_continuous(name = "Carrying capacity (k)",
                     breaks = ks) +
  scale_y_continuous(name = "Snail population size")

ggplot(final_equil, aes(x=k)) + 
  geom_line(aes(y=I/(S+E+I)), size = 2) +
  scale_x_continuous(name = "Carrying capacity (k)",
                     breaks = ks) +
  scale_y_continuous(name = "Prevalence of infection in snails (fraction)",
                     limits = c(0, 1))

ggplot(final_equil, aes(x=k)) + 
  geom_line(aes(y=C), size = 2) +
  scale_x_continuous(name = "Carrying capacity (k)",
                     breaks = ks) +
  scale_y_continuous(name = "Cercarial production")
ggplot(final_equil, aes(x=k)) + 
  geom_line(aes(y=I), size = 2) +
  scale_x_continuous(name = "Carrying capacity (k)",
                     breaks = ks) +
  scale_y_continuous(name = "Shedding snails")
                          
                          
                          
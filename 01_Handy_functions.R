#############################
#Author: Veronica Malizia
#R version: 4.1.2
#
#Handy functions employed throughout the model
#############################
`%!in%` <- Negate(`%in%`)
geom_mean <- function(x){exp(mean(log(x)))}
ci <- function(x){quantile(x, probs=c(0.025, 0.975), na.rm = T)}

#Functions (they can be a separate script)
age_groups <- c(0, 5, 10, 16, 200)
exposure_rates <- c(0.032, 0.61, 1, 0.06, 0.06) #Relative Age-specific exposure rates (activity/person/day) 
#exposure_rates <- c(0.33, 0.44, 0.22, 0) #Relative Age-specific exposure rates (activity/person/day) 
#Data from water contacts computed from Seydou S., De Vlas SJ, et al. (2011)
Age_profile_exp <- function(a){
  approx(x=age_groups, y=exposure_rates, xout=c(a), method = "constant")
}
Age_profile_contr <- function(a){
  approx(x=c(0, 10, 100), y=c(1, 1, 1), xout=c(a), method = "linear")
}
FOIh <- function(l, zeta, v, a, is){
  Rel_ex <- Age_profile_exp(a)$y * is
  return(l * zeta * v * Rel_ex)
}
hyp_sat <- function(alpha, beta, w){ #Hyperbolic saturating function 
  f <- (alpha*w) / (1 + ((alpha*w) / beta))
  return(f)
}
logistic <- function(k, w0, w){ #Logistic reduction function 
  f <- 1 / (1 + exp(k*(w-w0)))
  return(f)
}
expon_reduction <- function(alpha, w){
  f <- exp(-alpha*w)
  return(f)
}
#plot(logistic(k=imm, w0=w0_imm, c(0:8000)), ylim = c(0, 1), ylab = "Reduction", xlab = "Number of worm pairs")

#ODE system diefinition for snails' module
SEI <- function(t, x, parms) {    
  with(as.list(c(parms, x)), {
    N_s <- S+E+I #is it better to work with densities?
    #Logistic growth
    #For population growing with limited amount of resources
    beta <- beta0*(1-N_s/k) #Infected snails do not reproduce
    FOIs <- b*mir #/N_s
    #l0*(1-chi^(mir/N)) #Gurarie - non linear FOIs
    
    #Equations
    dS <- beta*(S+E) - (v+FOIs)*S #susceptible
    dE <- FOIs*S - (v+tau)*E #Exposed: snails are invaded, but larvae are not patent yet. Thus, snails do not shed cercariae
    dI <- tau*E - v2*I  #Infected: larvae in the snail are mature and snails shed cercariae
    dC <- lambda*I - m*C #Cercariae (output)
    res <- c(dS, dE, dI, dC)
    list(res)
  })
}

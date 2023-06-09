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
Age_profile_exp <- function(age_groups, exposure_rates, a, method){
  approx(x=age_groups, y=exposure_rates, xout=c(a), method = method)
}
Age_profile_contr <- function(a){
  approx(x=c(0, 10, 100), y=c(1, 1, 1), xout=c(a), method = "linear")
}
FOIh <- function(l, zeta, rel_exp, is){
  Rel_ex <- rel_exp$y * is
  return(l * zeta * Rel_ex)
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

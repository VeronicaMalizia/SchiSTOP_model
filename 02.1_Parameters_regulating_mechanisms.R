#############################
#Author: Veronica Malizia
#Date: 13/06/2024
#R version: 4.1.2
#
# The script explores functions and parameters for the implementation of three regulating mechanisms and one age-exposure function
# Species of interest: Schistosoma mansoni

# Output:
# - "Matched_alphas_low.RData", tuned alphas (i.e. expected number of eggs per sample per worm pair) to pass on the model in the setting of low endemicity 
# - "Matched_alphas_mod.RData", tuned alphas to pass on the model in the setting of moderate endemicity
# - "Matched_alphas_high.RData", tuned alphas to pass on the model in the setting of high endemicity
# - "Data_for_Fig1.RData", data to create Fig1 of the manuscript

#############################

rm(list = ls())

#Loading packages
# line below only if library folder is different than default
#.libPaths(c("C:/Program Files/R/R-4.1.2/library",.libPaths()))
library(ggplot2)
library(tidyverse)
library(patchwork)
library(dplyr)

#############################
# Worm-level regulation: density-dependent fecundity
# Parametrise function
#############################

# Final choice: exponential saturating function
exp_sat <- function(alpha, z, x){
  f <- alpha*x*exp(-z*x)
  return(f)
}

# Set a range of worm pairs
wp <- c(0:2000)

# Function's parameters:
# - z: severity of density-dependency
# - alpha: expected number of eggs per sample per worm pair
z1 <- 0.0007 #strong degree of regulation - from SCHISTOX
z2 <- 0.00022 #mild degree of regulation
alpha_strong <- 0.14 #strong degree of regulation - from SCHISTOX
M <- 150 #Average worm pair burden (50 in low endem, 150 in moderate endem, 500 in high endem)

# Match the function's initial growth to parametrise alpha for mild and absent degree of regulation
y = exp_sat(alpha_strong, z1, M) 
alpha_lin = y/M #absent degree of regulation
alpha_mild = y/(M*exp(-z2*M)) #mild degree of regulation

#Create dataset for ggplot
worms_reg <- bind_rows(data.frame(wp = wp,
                              y = exp_sat(alpha_lin, 0, wp),
                              degree = "Absent"),
                   data.frame(wp = wp,
                              y = exp_sat(alpha_mild, z2, wp),
                              degree = "Mild"),
                   data.frame(wp = wp,
                              y = exp_sat(alpha_strong, z1, wp),
                              degree = "Strong"))

#Basic plot for inspection
plot(wp, exp_sat(alpha_lin, 0, wp), type = 'l', col = 'red', #log = 'xy', 
     axes = F, lwd = 2, ylim = c(0,300),
     xlab = "Number of worm pairs", ylab = "Expected egg counts")
axis(side = 1, at = seq(0,2000,500), pos = 0)
axis(side = 2, at = seq(0,300,50), pos = 0)
lines(wp, exp_sat(alpha_mild, z2, wp), type = 'l', col = "purple", lwd = 2)
lines(wp, exp_sat(alpha_strong, z1, wp), type = 'l', col = "blue", lwd = 2)
#abline(v=M, col="dark grey", lwd = 2, lty=2)
legend("bottomright", 
       legend = c("Absent", "Mild", "Strong"), 
       col = c("red", "purple", "blue"), 
       lty = c(1, 1),
       bty = "n", 
       pt.cex = 2, 
       cex = 0.8, 
       text.col = "black", 
       horiz = F , 
       inset = c(0.1, 0.1))
mtext(c(paste("M=", M)), side=1, line=1, at=c(M), col="grey30", adj = 0)
text(1800, 280, c("(A)"), adj = 0)

# Save alphas for absent, mild, and strong degree of regulation 
# Repeat for low, moderate, and high endemicity by varying the value of M
save(alpha_mild, alpha_lin, alpha_strong, file = "Matched_alphas_high.RData")

################################################
# Human-level regulation: immunity function
###############################################

# Define a range of dead worm pairs
dwp <- c(0:10000) 

# Final choice: immunity factor defined as exponential reduction on the exposure
expon_reduction <- function(alpha, w){
  f <- exp(-alpha*w)
}
  
# immunity coefficient for absent, mild, strong degree of regulation
imm = c(0, 0.0005, 0.002)

# Data for ggplot
data <- c()
for(i in 1:length(imm)){
  data <- bind_rows(data,
                        bind_cols(dwp = dwp, imm_coefficient = imm[i], 
                                  y = expon_reduction(imm[i], dwp),
                                  func = "exponential"))
}

data$imm_coefficient <- as.factor(data$imm_coefficient)

data <- data %>%
  mutate(Immunity = case_when(imm_coefficient == imm[1] ~ "Absent",
                              imm_coefficient == imm[2] ~ "Mild",
                              imm_coefficient == imm[3] ~ "Strong"))
dataf <- as.data.frame(data)

# Basic plot for inspection
plot(dataf[which(dataf$Immunity=="Absent"), "dwp"]+0.5, 
     dataf[which(dataf$Immunity=="Absent"), "y"], 
     type = 'l', col = 'red', log = 'x', 
     axes = F, lwd = 2, ylim = c(0,1),
     xlab = "Cumulated dead worm pairs + 0.5", ylab = "Immunity factor")
axis(side = 1, pos = 0)
axis(side = 2, at = seq(0,1,0.2), labels = seq(0,1,0.2), pos = 0.5)
lines(dataf[which(dataf$Immunity=="Mild"), "dwp"]+0.5, 
      dataf[which(dataf$Immunity=="Mild"), "y"], 
      type = 'l', col = "purple", lwd = 2, log = 'x')
lines(dataf[which(dataf$Immunity=="Strong"), "dwp"]+0.5, 
      dataf[which(dataf$Immunity=="Strong"), "y"], 
      type = 'l', col = "blue", lwd = 2, log = 'x')
#abline(v=M, col="dark grey", lwd = 2, lty=2)
legend("topright", 
       legend = c("Absent", "Mild", "Strong"), 
       col = c("red", "purple", "blue"), 
       lty = c(1, 1),
       bty = "n", 
       pt.cex = 2, 
       cex = 0.8, 
       text.col = "black", 
       horiz = F , 
       inset = c(0.1, 0.1))
text(6000, 0.8, c("(B)"), adj = 0)


####################################
# Snail-level regulation: carrying capacity on snails' population growth
####################################

# Define function of population logistic growth
beta0 = 1 #max.reproduction.rate
beta <- function(beta0, K, x, i){
  return(beta0*(1-(x/K))*x)
}

# Define range of total number of snails
N = c(0:20000) 

#K = carrying capacity, x = infected snails
i = 0.06 #portion of infected

#Create dataset for ggplot
snail_reg <- bind_rows(data.frame(x = N,
                                  y = beta(beta0, 20000, x=N),
                                  degree = "Mild"),
                       data.frame(x = N,
                                  y = beta(beta0, 10000, x=N),
                                  degree = "Strong"))

#Basic plot for inspection
plot(N, beta(beta0, 20000, x=N),
     type = 'l', col = 'purple',  
     axes = F, lwd = 2, ylim = c(0,5000),
     xlab = "Total number of snails", ylab = "Snail birth rate")
axis(side = 1, pos = 0)
axis(side = 2, pos = 0)
lines(N, beta(beta0, 10000, x=N),
      type = 'l', col = "blue", lwd = 2)
legend("topright", 
       legend = c("Mild", "Strong"), 
       col = c("purple", "blue"), 
       lty = c(1, 1),
       bty = "n", 
       pt.cex = 2, 
       cex = 0.8, 
       text.col = "black", 
       horiz = F , 
       inset = c(0.1, 0.1))
text(18000, 4800, c("(C)"), adj = 0)

##############################  
# Age-exposure to infection: two scenarios
############################## 

# Define function: 
Age_profile_exp <- function(x, y, a){
  approx(x=x, y=y, xout=c(a), method = "constant")
}

# Model-derived function 
age_groups <- c(0, 5, 10, 16, 100)
exposure_rates <- c(0.032, 0.61, 1, 0.06, 0.06) #Relative Age-specific exposure rates as in [Turner 2017] (medium adult burden)
# Based on water contacts
# age categories and relative Age-specific exposure rates computed and extracted from [Sow 2011] and [Fulford 1996]
age <- c(0, 5, 15, 40, 100)
exp <- c(0, 0.62, 1, 0.51, 0.51) 

# Data for ggplot
exposure <- bind_rows(data.frame(x = age_groups,
                                     y = exposure_rates,
                                     degree = "Model-based"),
                          data.frame(x = age,
                                     y = exp,
                                     degree = "Water-contacts-based"))

# Basic plot for inspection
plot(x=age_groups, y=exposure_rates, 
     #main = "Age-exposure functions",
     axes = F, xlim = c(0, 100), ylim = c(0, 1),
     type = 's', col = "dark orange", lwd = 2,
     xlab = "Age", ylab = "Relative exposure") 
axis(side = 1, at = seq(0,100,20), pos = 0)
axis(side = 2, at = seq(0,1,0.2), pos = 0)
lines(x=age, y=exp, 
      type = 'l', col = "turquoise", lwd = 2)
#text(x = 60, y = 0.2, "AUC = 57", cex = 1, col = "dimgray")
legend("topright", 
       legend = c("Extracted from water contact data", 
                  "Available from literature (model-derived)"), 
       col = c("turquoise", "darkblue"), 
       #pch = c(1:2), 
       lty = c(1, 1),
       bty = "n", 
       pt.cex = 2, 
       cex = 0.8, 
       text.col = "black", 
       horiz = F)

################################
# Save data to generate Fig 1 of the manuscript
# This is done in script 07_PLOTTING_CODE.R 
################################
save(exposure, dataf, worms_reg, snail_reg, file = "Data_for_Fig1.RData")

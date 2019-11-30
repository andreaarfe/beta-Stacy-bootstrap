
# Script to run the simulation study.

# Prepares the R session
rm(list=ls())
library(survival)
source('./Programs/functions.R')
set.seed(467853)

# Loads the data
load('./Data/pbc_data.Rdata')

# beta-Stacy process prior
# Centering measure: exponential distribution with median survival of 10 years
# Concentration parameter: c=1
c_bs <- function(x){ rep(1, length(x)) } # prior precision function
lambda <- log(2)/(10) # rate parameter of the centering exponential distribution
F_bs <- function(x){ pexp(x, rate=lambda) } # cenetering distribution function
f_bs <- function(x){ dexp(x, rate=lambda) } # density of the centering distribution function

# numer of samples to generate
NREPS <- 10000 

# Simulates random survival functions using the BSB approximation of the 
# beta-Stacy posterior;
# generates the same number of replications NREPS defined above
# Uses different values of the BSB parameter m
# Placebo arm
#bsb_0_samp_1 <- vector(mode = "list", length = NREPS)
bsb_0_samp_10 <- vector(mode = "list", length = NREPS)
bsb_0_samp_100 <- vector(mode = "list", length = NREPS)
bsb_0_samp_1000 <- vector(mode = "list", length = NREPS)
for(i in 1:NREPS){
  #svMisc::progress(i, max.value = NREPS)
  #bsb_0_samp_1[[i]] <- bsb(1, time[I0], event[I0], c_bs, F_bs, f_bs)
  bsb_0_samp_10[[i]] <- bsb(10, time[I0], event[I0], c_bs, F_bs, f_bs)
  bsb_0_samp_100[[i]] <- bsb(100, time[I0], event[I0], c_bs, F_bs, f_bs)
  bsb_0_samp_1000[[i]] <- bsb(1000, time[I0], event[I0], c_bs, F_bs, f_bs)
}
# Experimental arm
#bsb_1_samp_1 <- vector(mode = "list", length = NREPS)
bsb_1_samp_10 <- vector(mode = "list", length = NREPS)
bsb_1_samp_100 <- vector(mode = "list", length = NREPS)
bsb_1_samp_1000 <- vector(mode = "list", length = NREPS)
for(i in 1:NREPS){
  #svMisc::progress(i, max.value = NREPS)
  #bsb_1_samp_1[[i]] <- bsb(1, time[I1], event[I1], c_bs, F_bs, f_bs)
  bsb_1_samp_10[[i]] <- bsb(10, time[I1], event[I1], c_bs, F_bs, f_bs)
  bsb_1_samp_100[[i]] <- bsb(100, time[I1], event[I1], c_bs, F_bs, f_bs)
  bsb_1_samp_1000[[i]] <- bsb(1000, time[I1], event[I1], c_bs, F_bs, f_bs)
}

save.image(file="./Results/simulations.Rdata")









  

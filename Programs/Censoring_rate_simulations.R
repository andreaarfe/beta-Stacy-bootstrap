
# Prepares the R session
rm(list=ls())
library(survival)
library(tidyverse)
library(ggfortify)
library(gridExtra)
library(compiler)
source('./Programs/functions.R')
set.seed(467853)
enableJIT(3)

# Loads the data
load('./Data/pbc_data.Rdata')

# beta-Stacy process prior
# Centering measure: exponential distribution with median survival of 10 years
# Concentration parameter: c=1
c_bs <- function(x){ rep(1, length(x)) } # prior precision function
lambda <- log(2)/(10) # rate parameter of the centering exponential distribution
F_bs <- function(x){ pexp(x, rate=lambda) } # cenetering distribution function
f_bs <- function(x){ dexp(x, rate=lambda) } # density of the centering distribution function

# Event rate in the placebo arm
lambda <- sum(event[I0])/sum(time[I0])

# numer of samples to generate
NREPS <- 10000 

# Generate non-censored data
d <- vector(mode = 'list', length = NREPS)
for(i in 1:NREPS){
  d[[i]] <- rexp(length(I0), rate = lambda)
}

# Proportions of censored observations
pcens <- c(0.00, 0.25, 0.50, 0.75)

# Simulations
sims <- vector(mode = 'list', length = length(pcens))
grid <- seq(0,10,length.out = 5000) # Grid for GvdVa algorithm
bsb_10   <- vector(mode = 'list', length = NREPS)
bsb_100  <- vector(mode = 'list', length = NREPS)
bsb_1000 <- vector(mode = 'list', length = NREPS)
GvdVa    <- vector(mode = 'list', length = NREPS)
# start <- Sys.time()
for(j in 1:length(pcens)){
  # Generates posterior samples
  print(j)
  for(i in 1:NREPS){
    delta <- rbinom(length(I0), 1, prob = 1 - pcens[j]) # 0 = censored, # 1 = event observed
    bsb_10[[i]]   <- bsb(10,   d[[i]], delta, c_bs, F_bs, f_bs)
    bsb_100[[i]]  <- bsb(100,  d[[i]], delta, c_bs, F_bs, f_bs)
    bsb_1000[[i]] <- bsb(1000, d[[i]], delta, c_bs, F_bs, f_bs)
    GvdVa[[i]]    <- betastacy_sim(1, c_bs, F_bs, f_bs, grid, d[[i]], delta)
  }
  sims[[j]] <- list(pcens    = pcens[j],
                    bsb_10   = bsb_10, 
                    bsb_100  = bsb_100,
                    bsb_1000 = bsb_1000,
                    GvdVa    = GvdVa)
}
# end <- Sys.time()
# end - start

saveRDS(sims, './Results/sims_censoring_final.rds')


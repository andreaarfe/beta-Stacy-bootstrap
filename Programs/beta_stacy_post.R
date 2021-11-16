
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
F_bs <- function(x){ pexp(x, rate=lambda) } # centering distribution function
f_bs <- function(x){ dexp(x, rate=lambda) } # density of the centering distribution function

# numer of samples to generate
NREPS <- 10000 

# Simulates random survival functions over the interval 0-12 years
# from the beta-Stacy posterior law; 
# uses discretization method with step = 1 day
x <- seq(0,10,length.out = 5000) # discretization grid
# Placebo arm
bstacy_0 <- betastacy_sim(NREPS, c_bs, F_bs, f_bs, x, time[I0], event[I0])
# Experimental arm
bstacy_1 <- betastacy_sim(NREPS, c_bs, F_bs, f_bs, x, time[I1], event[I1])

save.image(file = './Results/beta_stacy_post.Rdata')


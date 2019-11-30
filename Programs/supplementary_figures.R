
# Supplementary figures

# Prepares the R session
rm(list=ls())
library(survival)
library(pracma)
library(ggplot2)
library(ggfortify)
library(gridExtra)
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

#### Supplementary Figure S1

# Plot of Kaplan-Meier curves and posterior means
grid <- seq(0,12,0.01)
Fstar_0 <- F_bs_post(grid, time[I0], event[I0], c_bs, F_bs, f_bs)
Fstar_1 <- F_bs_post(grid, time[I1], event[I1], c_bs, F_bs, f_bs)
KM_0 <- survfit(Surv(time[I0],event[I0])~1)
KM_1 <- survfit(Surv(time[I1],event[I1])~1)
KM_0_values <- stepfun(KM_0$time, 1-c(1, KM_0$surv))(grid) 
KM_1_values <- stepfun(KM_1$time, 1-c(1, KM_1$surv))(grid) 
# Kolmogorv-smirnov statistics
#max(abs(Fstar_0-KM_0_values))
#max(abs(Fstar_1-KM_1_values))

pdf(file = './Results/Supplementary_Fig_1.pdf',
    width = 9, height = 5)
par(mfrow=c(1,2))
plot(grid, Fstar_0, type='s', ylim=c(0,0.8), lwd=4, col='blue',
     xlab='Time (in years)', ylab='Cumulative probability',
     main='a) Placebo arm')
points(grid, KM_0_values, type='s', lwd=4, col='darkorange', lty=2)
legend(x='topleft', lwd=c(3,3), lty=c(1,2), col=c('blue', 'darkorange'),
       legend=c('Posterior mean', 'Kaplan-Meier estimate'))
plot(grid, Fstar_1, type='s', ylim=c(0,0.8), lwd=4, col='blue',
     xlab='Time (in years)', ylab='Cumulative probability',
     main='b) D-penicillammine')
points(grid, KM_1_values, type='s', lwd=4, col='darkorange', lty=2)
legend(x='topleft', lwd=c(3,3), lty=c(1,2), col=c('blue', 'darkorange'),
       legend=c('Posterior mean', 'Kaplan-Meier estimate'))
par(mfrow=c(1,1))
dev.off() 

#### Supplementary Figure S2

# Sample of NREPS paths from the beta-Stacy prior
NREPS <- 100 # numer of samples to generate
x <- seq(0,12,1/365.25) # discretization grid
Fsamp <- betastacy_sim(NREPS, c_bs, F_bs, f_bs, x, time = NULL)

# Plot of sample paths and prior mean
pdf(file = './Results/Supplementary_Fig_2.pdf')
plot(Fsamp$grid,Fsamp$cumul_probs[1,], ylim=c(0,1), col='white',
     xlab='Time (in years)', ylab='Cumulative probability')
for(i in 1:NREPS){
  points(Fsamp$grid,Fsamp$cumul_probs[i,], type='s', col='black')
}
points(x, F_bs(x), col='red', type='l', lwd=3)
dev.off()


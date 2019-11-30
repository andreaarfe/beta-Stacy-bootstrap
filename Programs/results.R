
# Script to compute results for the paper

# Prepares the R session
rm(list=ls())
library(survival)
library(pracma)
library(ggplot2)
library(ggfortify)
library(gridExtra)
source('./Programs/functions.R')
load("./Results/simulations.Rdata")
load('./Data/pbc_data.Rdata')
load('./Results/beta_stacy_post.Rdata')
set.seed(467853)

# Numer of deaths per arm
sum(event[I0]) # placebo arm
sum(event[I1]) # experimental arm

# Kolmogorov-Smirnov distance between posterior means and Kapla-Meier curves
grid <- seq(0,12,0.01)
Fstar_0 <- F_bs_post(grid, time[I0], event[I0], c_bs, F_bs, f_bs)
Fstar_1 <- F_bs_post(grid, time[I1], event[I1], c_bs, F_bs, f_bs)
KM_0 <- survfit(Surv(time[I0],event[I0])~1)
KM_1 <- survfit(Surv(time[I1],event[I1])~1)
D_0 <- max(abs(stepfun(KM_0$time, 1-c(1, KM_0$surv))(grid) - Fstar_0))
D_1 <- max(abs(stepfun(KM_1$time, 1-c(1, KM_1$surv))(grid) - Fstar_1))
D_0 # placebo arm
D_1 # experimental arm

# cut-point for RMST computations
tau <- 10

# RMST estimates in the placebo group based on the beta-Stacy process
x <- bstacy_0$grid
delta <- diff(c(x[x<tau], tau))
bstacy_rmst_0 <- apply(1-bstacy_0$cumul_probs,1, 
                       function(s) sum(s[1:length(delta)]*delta))

# RMST estimates in the placebo group based on the BSB
f_rmst <- function(g, tau) sum(ifelse(g$atoms<=tau,g$atoms,tau)*g$weights)
#bsb_rmst_0_1 <- sapply(bsb_0_samp_1, f_rmst)
bsb_rmst_0_10 <- sapply(bsb_0_samp_10, f_rmst, tau)
bsb_rmst_0_100 <- sapply(bsb_0_samp_100, f_rmst, tau)
bsb_rmst_0_1000 <- sapply(bsb_0_samp_1000, f_rmst, tau)

# Survival probability estimates in the placebo group based on the beta-Stacy process
bstacy_surv_0 <- apply(1-bstacy_0$cumul_probs,1, 
                       function(s) stepfun(x,c(0,s))(tau))

# Survival probability estimates in the placebo group based on the BSB
f_surv <- function(g, tau) 1-sum(g$weights[g$atoms<=tau])
#bsb_surv_0_1 <- sapply(bsb_0_samp_1, f_surv)
bsb_surv_0_10 <- sapply(bsb_0_samp_10, f_surv, tau)
bsb_surv_0_100 <- sapply(bsb_0_samp_100, f_surv, tau)
bsb_surv_0_1000 <- sapply(bsb_0_samp_1000, f_surv, tau)

# Expected survival times in the two arms
bsb_exp_0_10 <- sapply(bsb_0_samp_10, f_rmst, Inf)
bsb_exp_0_100 <- sapply(bsb_0_samp_100, f_rmst, Inf)
bsb_exp_0_1000 <- sapply(bsb_0_samp_1000, f_rmst, Inf)
bsb_exp_1_10 <- sapply(bsb_1_samp_10, f_rmst, Inf)
bsb_exp_1_100 <- sapply(bsb_1_samp_100, f_rmst, Inf)
bsb_exp_1_1000 <- sapply(bsb_1_samp_1000, f_rmst, Inf)

# Difference in expected survival times (D-penicilamine - Placebo)
diff_10 <- bsb_exp_1_10 - bsb_exp_0_10
diff_100 <- bsb_exp_1_100 - bsb_exp_0_100
diff_1000 <- bsb_exp_1_1000 - bsb_exp_0_1000

# Two-sample Kolmogorv-Smirnov statistics
# Survival probabilities
ks.test(bstacy_surv_0,bsb_surv_0_10)
ks.test(bstacy_surv_0,bsb_surv_0_100)
ks.test(bstacy_surv_0,bsb_surv_0_1000)
# RMSTS
ks.test(bstacy_rmst_0,bsb_rmst_0_10)
ks.test(bstacy_rmst_0,bsb_rmst_0_100)
ks.test(bstacy_rmst_0,bsb_rmst_0_1000)
# Mean survival time
ks.test(diff_1000,diff_100)

# Save the results needed to prepare Figure 2
outfile <- c('bstacy_rmst_0', 'bsb_rmst_0_10', 'bsb_rmst_0_100', 'bsb_rmst_0_1000',
             'bstacy_surv_0', 'bsb_surv_0_10', 'bsb_surv_0_100', 'bsb_surv_0_1000',
             'bsb_exp_0_10', 'bsb_exp_0_100', 'bsb_exp_0_1000', 
             'bsb_exp_1_10', 'bsb_exp_1_100', 'bsb_exp_1_1000',
             'diff_10', 'diff_100', 'diff_1000')
save(list=outfile, file='./Results/results.Rdata')



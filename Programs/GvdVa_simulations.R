# Prepares the R session
rm(list=ls())
library(survival)
library(tidyverse)
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

# numer of samples to generate
NREPS <- 10000 

# Simulates random survival functions over the interval 0-12 years
# from the beta-Stacy posterior law; 
# uses discretization method with step = 1 day

# Placebo arm
bstacy_5 <- betastacy_sim(NREPS, c_bs, F_bs, f_bs, 
                           seq(0,10,length.out = 5), 
                          time[I0], event[I0])

bstacy_50 <- betastacy_sim(NREPS, c_bs, F_bs, f_bs, 
                           seq(0,10,length.out = 50), 
                           time[I0], event[I0])

bstacy_500 <- betastacy_sim(NREPS, c_bs, F_bs, f_bs, 
                          seq(0,10,length.out = 500),
                          time[I0], event[I0])

bstacy_5000 <- betastacy_sim(NREPS, c_bs, F_bs, f_bs, 
                          seq(0,10,length.out = 5000), 
                          time[I0], event[I0])


# cut-point for RMST computations
tau <- 10

# RMST estimates in the placebo group based on the beta-Stacy process
fn <- function(bstacy){
  x <- bstacy$grid
  delta <- diff(c(x[x<tau], tau))
  apply(1-bstacy$cumul_probs,1, 
      function(s) sum(s[1:length(delta)]*delta))
}
bstacy_rmst_5     <- fn(bstacy_5)
bstacy_rmst_50    <- fn(bstacy_50)
bstacy_rmst_500   <- fn(bstacy_500)
bstacy_rmst_5000  <- fn(bstacy_5000)

# boxplot(bstacy_rmst_5, ylim=c(5,9))
# boxplot(bstacy_rmst_50, ylim=c(5,9))
# boxplot(bstacy_rmst_500, ylim=c(5,9))
# boxplot(bstacy_rmst_5000, ylim=c(5,9))


# 10-years survival probability
bstacy_surv_5     <- apply(1-bstacy_5$cumul_probs,1, 
                        function(s) stepfun(bstacy_5$grid,c(0,s))(tau))
bstacy_surv_50    <- apply(1-bstacy_50$cumul_probs,1, 
                           function(s) stepfun(bstacy_50$grid,c(0,s))(tau))
bstacy_surv_500   <- apply(1-bstacy_500$cumul_probs,1, 
                           function(s) stepfun(bstacy_500$grid,c(0,s))(tau))
bstacy_surv_5000  <- apply(1-bstacy_5000$cumul_probs,1, 
                           function(s) stepfun(bstacy_5000$grid,c(0,s))(tau))

# boxplot(bstacy_surv_5, ylim=c(0,1))
# boxplot(bstacy_surv_50, ylim=c(0,1))
# boxplot(bstacy_surv_500, ylim=c(0,1))
# boxplot(bstacy_surv_5000, ylim=c(0,1))


dset <- data.frame(N = rep(c(5, 50, 500, 5000), each = NREPS),
                   samps = c(bstacy_rmst_5, bstacy_rmst_50,
                             bstacy_rmst_500, bstacy_rmst_5000))
p <- dset %>% 
  ggplot(aes(x = factor(N), y = samps, fill = factor(N))) +
  geom_violin(show.legend = FALSE, alpha=0.6) + theme_minimal() +
  theme(plot.title = element_text(face = 'bold')) +
  xlab('N') +
  ylab('10-years RMST (in years)') +
  scale_fill_brewer(palette="Blues") +
  geom_boxplot(width=0.1, show.legend = FALSE,
               outlier.shape = NA, alpha=0) +
  labs(title='b) Restricted mean survival time',
       subtitle = '(Placebo arm)') +
  ylim(3,10)
  
dset2 <- data.frame(N = rep(c(5, 50, 500, 5000), each = NREPS),
                   samps = c(bstacy_surv_5, bstacy_surv_50,
                             bstacy_surv_500, bstacy_surv_5000))  
q <- dset2 %>% 
  ggplot(aes(x = factor(N), y = samps, fill = factor(N))) +
  geom_violin(show.legend = FALSE, alpha=0.6) + theme_minimal() +
  theme(plot.title = element_text(face = 'bold')) +
  xlab('N') +
  ylab('10-years survival probability') +
  scale_fill_brewer(palette="Blues") +
  geom_boxplot(width=0.1, show.legend = FALSE,
               outlier.shape = NA, alpha=0) +
  labs(title='a) Survival probability',
       subtitle = '(Placebo arm)') +
  ylim(0,1)

out <- grid.arrange(q, p, nrow = 1)
ggsave('./Results/GvdVa.pdf', out,
       width = 21, height = 15/2, units = 'cm')

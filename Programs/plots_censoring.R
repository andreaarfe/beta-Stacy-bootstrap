
# Prepares the R session
rm(list=ls())
library(survival)
library(tidyverse)
library(ggfortify)
library(gridExtra)
#library(compiler)
source('./Programs/functions.R')
set.seed(467853)
#enableJIT(3)

# Loads the simulation results
sims <- readRDS('./Results/sims_censoring_final.rds')

# Proportions of censored observations
pcens <- c(0.00, 0.25, 0.50, 0.75)

# cut-point for RMST computations
tau <- 10

### Results: RMST

rmst_sims <- vector(mode = 'list', length = length(pcens))
for(j in 1:length(pcens)){
  d <- matrix(0, nrow = length(sims[[j]]$GvdVa), ncol = 5)
  colnames(d) <- c('pcens', 'GvdVa', 'BSB-10', 
                   'BSB-100', 'BSB-1000')
  d[, 1] <- pcens[j]
  
  # RMST estimates in the placebo group based on the beta-Stacy process
  x <- sims[[j]]$GvdVa[[1]]$grid
  for(i in 1:length(sims[[j]]$GvdVa)){
    cumul_probs <- sims[[j]]$GvdVa[[i]]$cumul_probs
    delta <- diff(c(x[x<tau], tau))
    d[i, 'GvdVa'] <- apply(1-cumul_probs,1, 
                           function(s) sum(s[1:length(delta)]*delta))
  }
  
  # RMST estimates in the placebo group based on the BSB
  f_rmst <- function(g, tau) sum(ifelse(g$atoms<=tau,g$atoms,tau)*g$weights)
  #bsb_rmst_0_1 <- sapply(bsb_0_samp_1, f_rmst)
  d[, 'BSB-10'] <- sapply(sims[[j]]$bsb_10, f_rmst, tau)
  d[, 'BSB-100'] <- sapply(sims[[j]]$bsb_100, f_rmst, tau)
  d[, 'BSB-1000'] <- sapply(sims[[j]]$bsb_1000, f_rmst, tau)
  rmst_sims[[j]] <- as.data.frame(d)
}
rmst_sims <- do.call(rbind, rmst_sims)
rmst_sims$pcens <- factor(rmst_sims$pcens,
                          labels = paste0('Censoring probability: ', round(pcens, 2)))


p1 <- rmst_sims %>% 
  pivot_longer(cols = -pcens) %>% 
  mutate(name = factor(name, 
                       levels = c('BSB-10', 'BSB-100', 
                                  'BSB-1000', 'GvdVa'),
                       labels = c('BSB-10', 'BSB-100', 
                                  'BSB-1000', expression('GvdVa \n (reference)'))
                       )) %>%  
  ggplot(aes(x = name, y = value, fill=name)) +
  geom_violin(show.legend = FALSE, alpha=0.6) +
  facet_wrap(vars(pcens)) +
  theme_bw() +
  scale_fill_brewer(palette="Blues") +
  theme(plot.title = element_text(face = 'bold')) +
  xlab('Method') + 
  ylab('10-years RMST (in years)') +
  ggtitle('Simulation results: RMST') +
  geom_boxplot(width=0.1,show.legend = FALSE,
               outlier.shape = NA, alpha=0)

ggsave('./Results/Censoring_rmst.pdf', plot = p1,
       width = 15, height = 15,  unit = 'cm')
 
# Two-sample Kolmogorv-Smirnov statistics
# RMSTS
ks.test(rmst_sims %>% 
          filter(pcens == 'Censoring probability: 0') %>% 
          select('BSB-1000') %>% 
          unlist(),
        rmst_sims %>% 
          filter(pcens == 'Censoring probability: 0') %>% 
          select('GvdVa') %>% 
          unlist())
ks.test(rmst_sims %>% 
          filter(pcens == 'Censoring probability: 0.25') %>% 
          select('BSB-1000') %>% 
          unlist(),
        rmst_sims %>% 
          filter(pcens == 'Censoring probability: 0.25') %>% 
          select('GvdVa') %>% 
          unlist())
ks.test(rmst_sims %>% 
          filter(pcens == 'Censoring probability: 0.5') %>% 
          select('BSB-1000') %>% 
          unlist(),
        rmst_sims %>% 
          filter(pcens == 'Censoring probability: 0.5') %>% 
          select('GvdVa') %>% 
          unlist())
ks.test(rmst_sims %>% 
          filter(pcens == 'Censoring probability: 0.75') %>% 
          select('BSB-1000') %>% 
          unlist(),
        rmst_sims %>% 
          filter(pcens == 'Censoring probability: 0.75') %>% 
          select('GvdVa') %>% 
          unlist())

### Results: Survival probability

surv_sims <- vector(mode = 'list', length = length(pcens))
for(j in 1:length(pcens)){
  d2 <- matrix(0, nrow = length(sims[[j]]$GvdVa), ncol = 5)
  colnames(d2) <- c('pcens', 'GvdVa', 'BSB-10', 
                   'BSB-100', 'BSB-1000')
  d2[, 1] <- pcens[j]
  
  # Survival probability estimates in the placebo group based on the beta-Stacy process
  x <- sims[[j]]$GvdVa[[1]]$grid
  for(i in 1:length(sims[[j]]$GvdVa)){
    cumul_probs <- sims[[j]]$GvdVa[[i]]$cumul_probs
    d2[i, 'GvdVa'] <- apply(1-cumul_probs,1, 
                           function(s) stepfun(x,c(0,s))(tau))
  }
  
  # Survival probability estimates in the placebo group based on the BSB
  f_surv <- function(g, tau) 1-sum(g$weights[g$atoms<=tau])
  #bsb_surv_0_1 <- sapply(bsb_0_samp_1, f_surv)
  d2[, 'BSB-10']   <- sapply(sims[[j]]$bsb_10  , f_surv, tau)
  d2[, 'BSB-100']  <- sapply(sims[[j]]$bsb_100 , f_surv, tau)
  d2[, 'BSB-1000'] <- sapply(sims[[j]]$bsb_1000, f_surv, tau)
  surv_sims[[j]] <- as.data.frame(d2)
}
surv_sims <- do.call(rbind, surv_sims)
surv_sims$pcens <- factor(rmst_sims$pcens,
                          labels = paste0('Censoring probability: ', round(pcens, 2)))

p2 <- surv_sims %>% 
  pivot_longer(cols = -pcens) %>% 
  mutate(name = factor(name, 
                       levels = c('BSB-10', 'BSB-100', 
                                  'BSB-1000', 'GvdVa'),
                       labels = c('BSB-10', 'BSB-100', 
                                  'BSB-1000', expression('GvdVa \n (reference)'))
  )) %>%  
  ggplot(aes(x = name, y = value, fill=name)) +
  geom_violin(show.legend = FALSE, alpha=0.6) +
  facet_wrap(vars(pcens)) +
  theme_bw() +
  scale_fill_brewer(palette="Blues") +
  theme(plot.title = element_text(face = 'bold')) +
  xlab('Method') + 
  ylab('10-years survival probability') +
  ggtitle('Simulation results: survival probability') +
  geom_boxplot(width=0.1,show.legend = FALSE,
               outlier.shape = NA, alpha=0)

ggsave('./Results/Censoring_surv.pdf', plot = p2,
       width = 15, height = 15,  unit = 'cm')

# Two-sample Kolmogorv-Smirnov statistics
# survival probability
ks.test(surv_sims %>% 
          filter(pcens == 'Censoring probability: 0') %>% 
          select('BSB-1000') %>% 
          unlist(),
        surv_sims %>% 
          filter(pcens == 'Censoring probability: 0') %>% 
          select('GvdVa') %>% 
          unlist())
ks.test(surv_sims %>% 
          filter(pcens == 'Censoring probability: 0.25') %>% 
          select('BSB-1000') %>% 
          unlist(),
        surv_sims %>% 
          filter(pcens == 'Censoring probability: 0.25') %>% 
          select('GvdVa') %>% 
          unlist())
ks.test(surv_sims %>% 
          filter(pcens == 'Censoring probability: 0.5') %>% 
          select('BSB-1000') %>% 
          unlist(),
        surv_sims %>% 
          filter(pcens == 'Censoring probability: 0.5') %>% 
          select('GvdVa') %>% 
          unlist())
ks.test(surv_sims %>% 
          filter(pcens == 'Censoring probability: 0.75') %>% 
          select('BSB-1000') %>% 
          unlist(),
        surv_sims %>% 
          filter(pcens == 'Censoring probability: 0.75') %>% 
          select('GvdVa') %>% 
          unlist())





# prepares the data for analyses
# Prepares the R session
rm(list=ls())
library(survival)

# Prepares the data
dset <- pbc[1:312,] # selects only the randomized patients
time <- dset$time/365.25 # follow-up times in years
event <- ifelse(dset$status==2,1,0) # 0 = censored, 1 = death 
arm <- dset$trt
I0 <- which(arm==2) # placebo patients
I1 <- which(arm==1) # D-penicillammine patients

save(list=ls(), file='./Data/pbc_data.Rdata')
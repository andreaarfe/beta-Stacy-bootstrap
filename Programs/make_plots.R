
# Script prepare the figures for the paper

# Prepares the R session
rm(list=ls())
library(survival)
library(pracma)
library(ggplot2)
library(ggfortify)
library(gridExtra)
source('./Programs/functions.R')
load("./Data/pbc_data.Rdata")
load("./Results/results.Rdata")
set.seed(467853)

# Kaplan-Meier curve
S <- Surv(time, event)
arm_f <- factor(arm, labels = c('D-penicillammine','Placebo'))
KM <- survfit(S~arm_f)
p1 <- autoplot(KM, conf.int = FALSE) + 
  theme_minimal() + theme(plot.title = element_text(face = 'bold')) +
  theme(legend.position = c(0.3, 0.3),
        legend.background = element_rect(fill="white",
                                         size=0.1, linetype="solid", 
                                         colour ="black")) +
  xlab('Time since randomization (in years)') +
  ylab('Survival probability') +
  ylim(0,1) +
  labs(title='a) Kaplan-Meier curves', subtitle = '(Marks show censored follow-up times)') +
  scale_color_discrete(name = 'Arm') 

# Violin plot - comparison of the posterior distributions of the RMST in 
# the control arm
d_rmst <- data.frame(surv=c(bstacy_rmst_0,
                       bsb_rmst_0_10,
                       bsb_rmst_0_100,
                       bsb_rmst_0_1000),
                Method=rep(c('Beta-Stacy', 'BSB-10', 'BSB-100', 'BSB-1000'), 
                           each=length(bstacy_rmst_0)))
d_rmst$Method <- factor(d_rmst$Method, levels = c('BSB-10', 'BSB-100', 'BSB-1000', 'Beta-Stacy'))
p2 <- ggplot(d_rmst, aes(y=surv, x=Method, fill=Method)) +
  geom_violin(show.legend = FALSE, alpha=0.6) +
  geom_boxplot(width=0.1, show.legend = FALSE,
               outlier.shape = NA, alpha=0) +
  #stat_summary(fun.y=mean, geom="point", size=2, shape=23, color='red',
  #             show.legend = FALSE) +
  scale_fill_brewer(palette="Blues") + 
  theme_minimal() + 
  theme(plot.title = element_text(face = 'bold')) +
  xlab('Method') +
  ylab('10-years RMST (in years)') + 
  labs(title='c) Restricted mean survival time', subtitle='(Placebo arm)')

# Violin plot - comparison of the posterior distributions of the survival 
# probability in the control arm
d_surv <- data.frame(surv=c(bstacy_surv_0,
                            bsb_surv_0_10,
                            bsb_surv_0_100,
                            bsb_surv_0_1000),
                     method=rep(c('Beta-Stacy', 'BSB-10', 'BSB-100', 'BSB-1000'), 
                                each=length(bstacy_surv_0)))
d_surv$method <- factor(d_surv$method, levels = c('BSB-10', 'BSB-100', 'BSB-1000', 'Beta-Stacy'))
p3 <- ggplot(d_surv, aes(y=surv, x=method, fill=method)) + 
  geom_violin(show.legend = FALSE, alpha=0.6) +
  geom_boxplot(width=0.1,show.legend = FALSE,
               outlier.shape = NA, alpha=0) +
  #stat_summary(fun.y=mean, geom="point", size=2, shape=23, color='red',
  #             show.legend = FALSE) +
  scale_fill_brewer(palette="Blues") + theme_minimal() +
  theme(plot.title = element_text(face = 'bold')) +
  xlab('Method') + ylim(0,1) +
  ylab('10-years survival probability') +
  labs(title='b) Survival probability', subtitle =  '(Placebo arm)')

# Histograms - Difference in mean survival times 
d_diff <- data.frame(diff=c(diff_10,
                            diff_100,
                            diff_1000),
                     Method=rep(c('BSB-10', 'BSB-100', 'BSB-1000'), 
                                each=length(bstacy_surv_0)))
d_diff$Method <- factor(d_diff$Method, levels = c('BSB-10', 'BSB-100', 'BSB-1000'))
p4 <- ggplot(d_diff, aes(x=Method, y=diff,
                   #y=..density.., 
                   fill=Method)) + 
  geom_violin(show.legend = FALSE, alpha=0.6) +
  geom_boxplot(width=0.1,show.legend = FALSE,
               outlier.shape = NA, alpha=0) +
  #geom_density(alpha = 0.6) + 
  #stat_summary(fun.y=mean, geom="point", size=2, shape=23, color='red',
  #             show.legend = FALSE) +
  scale_fill_brewer(palette="Blues") + theme_minimal() +
  theme(plot.title = element_text(face = 'bold')) +
  xlab('Methods') + 
  ylab('Difference in mean survival time') +
  labs(title='d) Difference in mean survival times',
       subtitle = '(D-penicillammine vs. placebo)')

# Final plot
p_final <- grid.arrange(p1,p3,p2,p4,nrow=2, ncol=2)
ggsave('./Results/Figure.pdf', p_final,
       width = 21, height = 15, units = 'cm')

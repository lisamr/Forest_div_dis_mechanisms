#Community competency vs diversity
rm(list=ls())
library(tidyverse)
library(bayesplot)
library(rstan)

#functions----------------------------------------
zscore <- function(x) (x-mean(x))/(2*sd(x))
unzscore <- function(z, x)  z * 2*sd(x) + mean(x)


#load data + model--------------------------------
plot_level <- read_csv('data/plot_level_data.csv') #plot level data
CC_stan <- stan_model('Stan_models/CommunityCompetence.stan')


#prep data----------------------------------------
df <- plot_level %>% 
  mutate(ccomp = ccompsqrt^2,
         ccomplog = log(ccomp)) %>% 
  select(BSPlotNumber, ccomp, Richness, ccomplog, Richness_z, ForestAllianceType)
  
#data list for the model
newdat <- expand_grid(Richness_z = unique(df$Richness_z),
                     ForestAllianceType = c(0,1))
df_list <- list(
  N = nrow(df),
  K = 2,
  y = df$ccomplog,
  X = cbind(df$Richness_z, df$ForestAllianceType),
  Nsim = nrow(newdat),
  Xsim = as.matrix(newdat)
)
str(df_list)

#model it!----------------------------------------

#checking out priors
N=10000
a0 <- rnorm(N, 6, 2)
xrange <- seq(-2, 2, length.out = 100)
b <- rnorm(N, 0, 1)
plot(NULL, xlim=c(-2,2), ylim=c(-2,15))
for(i in 1:100) lines(xrange, a0[i] + b[i]*xrange, col=col.alpha('black'))
#sanity check--plot min and max observed community competency
abline(h = c(min(df$ccomplog), max(df$ccomplog)), col='red')

#fit model
set.seed(2020)
fit <- sampling(CC_stan, data = df_list, chains = 4, iter = 2000, seed = 2020)

#quick look at estimates
rethinking::precis(fit, depth = 2, pars = c('a0', 'B', 'sigma'), prob = .9)
#       mean   sd    5%   95% n_eff Rhat4
#a0     5.91 0.08  5.78  6.05  2681     1 #global intercept
#B[1]  -0.32 0.12 -0.52 -0.12  3362     1 #richness effect
#B[2]   0.45 0.13  0.24  0.65  2699     1 #forest type
#sigma  0.72 0.04  0.65  0.79  3775     1 #SD



#posterior predictive check-----------------------
post <- extract.samples(fit)
ppc_dens_overlay(df$ccomplog, post$y_rep[1:100,])+
  labs(x = 'log(community competency)') +
  theme(legend.position = 'none')

ggsave('figures/PPC_community_compentency.pdf', units = 'in', device = 'pdf', width = 3, height = 2)




#export posterior predictions---------------------

#put into a dataframe
sim_median <- apply(post$mu_sim, 2, median)
sim_CI <- apply(post$mu_sim, 2, HPDI, .9) #90% CI
sim_df <- newdat %>% 
  mutate(CCmedian = sim_median, 
         lower = sim_CI[1,], 
         upper = sim_CI[2,],
         Richness = unzscore(Richness_z, df$Richness))

ggplot(sim_df, aes(Richness, CCmedian, group = ForestAllianceType, color = ForestAllianceType)) +
  geom_line() +
  geom_jitter(data = df, aes(Richness, ccomp), height = 0, width = .1)


write_csv(sim_df, 'outputs/CC_postpredictions.csv')

#Gaussian process disease models. Not kfold, so can't contrast model perfomance of plot-level binomial mods.

#just running the models with richness + key host densities since those were consistently the best performing models. compared posteriors with non-gp models--nearly identical posteriors.


rm(list=ls())
library(rethinking)
library(tidyverse)
library(bayesplot)
library(loo)


#load data------------------------------------------------
plot_level <- read_csv('data/plot_level_data.csv') #data on 151 plots
tree_level <- read_csv('data/tree_level_data_HS.csv') #tree level data of highly susceptible species
dmat <- read_rds('data/dmat.RDS')



#create data lists----------------------------------------

#function for the binomial models
plot_vars <- c("ForestAllianceType", "SampleYear", "ppt_z", "psi_z", "forest.200m_z")
make_datlist <- function(predictors, all=T){
  plot_level_mat <- as.matrix(select(plot_level, plot_vars, predictors))
  if(all == T){
    I <- as.integer(plot_level$I_all)
    n <- as.integer(plot_level$n_all)
  }else{
    I <- as.integer(plot_level$I_HS)
    n <- as.integer(plot_level$n_HS)
  }
  datlist <- list(
    N=nrow(plot_level_mat), #number of plots
    K=ncol(plot_level_mat), #number of predictors
    X=plot_level_mat, #covariate matrix
    I=I, #response (number of infections)
    n=n, #total host plants
    D = dmat#distance matrix
  ) 
  return(datlist)
}

#function for Bernoulli models
make_datlist_Bern <- function(predictors=NULL){
  if(is.null(predictors)){
    DF <- plot_level %>% select(plot_vars) %>% as.matrix()
  }else{
    DF <- plot_level %>% select(plot_vars, predictors) %>% as.matrix()
  }
  datlist <- list(
    N = nrow(tree_level), #number of plots
    S = as.integer(max(tree_level$spID)), #number of species
    J = nrow(DF), #number of observations (individuals)
    Kp = ncol(DF), #number of plot-level predictors
    Xp = DF, #plot-level covariate matrix
    I = as.integer(tree_level$I), #infected or not
    SpID = as.integer(tree_level$spID), #species index
    PlotID = as.integer(tree_level$plotID), #plot index
    XS = as.matrix(cbind(1, plot_level$Richness_z[tree_level$plotID])), #model matrix, varies by species
    BA = tree_level$basal_area_z, #basal area of individual
    D = dmat
  )
  return(datlist)
}


# 1. Plot-level models, all hosts
datalists_plotlevel_all <- make_datlist(c('Richness_z', 'BA_UMCAsqrt_z', 'BA_LIDEsqrt_z'), all = T)


# 3. individual-level models, highly symptomatic hosts
datalists_indlevel <- make_datlist_Bern(c('BA_UMCAsqrt_z', 'BA_LIDEsqrt_z'))






#models---------------------------------------------------

#gaussian process beta-binomial model
stan_betaBinom <- stan_model(file = 'Stan_models/beta_binomialGP.stan')

#gaussian process Bernoulli model with varying diversity slopes and individual-level basal area
stan_bern <- stan_model(file = 'Stan_models/bernoulli_varyingslopeGPnc.stan')




#fit models-----------------------------------------------


#Run the models (just doing ones with richness + key host densities since they performed best without the GP)

#plot-level all plots
set.seed(2021)
post_plotlevel <- sampling(stan_betaBinom, datalists_plotlevel_all, iter = 2000, chains = 4, cores = 4, seed = 2021)


#individual-level 
#CAUTION: TAKES A LOT OF TIME AND MEMORY. BUDGET A FEW HOURS!!!
set.seed(2021)
post_indlevel <- sampling(stan_bern, datalists_indlevel, iter = 2000, chains = 4, cores = 4, control = list(adapt_delta = .99), seed = 2021)


#export models--------------------------------------------
saveRDS(post_plotlevel_allGP, '../../../../Box/Stan_model_outputs/Big_Sur/post_plotlevel_all_GP.RDS')
saveRDS(post_indlevel, '../../../../Box/Stan_model_outputs/Big_Sur/post_indlevel_GP.RDS')


#read in models-----------------------------------
#I ran these previously when it was set up in a list. 
post_plotlevel <- read_rds('../../../../Box/Stan_model_outputs/Big_Sur/post_plotlevel_all_GP.RDS')
post_indlevel <- read_rds('../../../../Box/Stan_model_outputs/Big_Sur/post_indlevel_GP.RDS')



#compare posterior with non-GP models-----------------------

#plot level
post_plotlevel_nonGP <- read_rds('../../../../Box/Stan_model_outputs/Big_Sur/post_plotlevel_all.RDS')
precis(post_plotlevel[[1]], pars = c(c('a0', 'beta', 'theta', 'eta_sq', 'rho_sq')), depth = 2, prob = .9)
precis(post_plotlevel_nonGP[[2]], pars = c(c('a0', 'beta', 'theta')), depth = 2, prob = .9)


#individual-level
precis(post_indlevel[[1]], pars = c(c('bBA', 'betap', 'betaS', 'zS', 'eta_sq', 'rho_sq', 'sig_sq')), depth = 3, prob = .9)
#need to add zS to mean to get the values you report in the paper
samp_indlevel <- extract.samples(post_indlevel[[1]])
apply(samp_indlevel$betaS[,2] + samp_indlevel$zS[,2,], 2, HPDI, .9)
#[,1]      [,2]      [,3]      [,4]
#|0.9 -1.7020691 -1.155407 -0.583679 0.1998709
#0.9|  0.3483102  1.122726  1.940522 2.3254298
post_indlevel_nonGP <- read_rds('../../../../Box/Stan_model_outputs/Big_Sur/post_indlevel.RDS')
precis(post_indlevel_nonGP[[2]], pars = c(c('bBA', 'betap', 'betaS', 'zS')), depth = 3, prob = .9)
samp_indlevel <- extract.samples(post_indlevel_nonGP[[2]])
apply(samp_indlevel$betaS[,2] + samp_indlevel$zS[,2,], 2, HPDI, .9)
#[,1]      [,2]       [,3]      [,4]
#|0.9 -1.7510532 -1.087464 -0.6691851 0.1361247
#0.9|  0.2661887  1.098156  1.9489762 2.3431677



#contrast models with loo-------------------------

#plot-level models would need to be compared wih k-fold. My computer doesn't have enough memory to do htat. the posteriors are nearly identical, so it would be moot.
loo_plotlevel <- lapply(list(post_plotlevel[[1]], post_plotlevel_nonGP[[2]]), loo)
loo_compare(loo_plotlevel)
#> loo_compare(loo_plotlevel_all)
#elpd_diff se_diff
#model1  0.0       0.0   (GP) #9 obs with pareto k > .7. 
#model2 -1.5       0.8   (non-GP)

#individual-level
loo_indlevel <- lapply(list(post_indlevel[[1]], post_indlevel_nonGP[[2]]), loo)
loo_compare(loo_indlevel)
#> loo_compare(loo_indlevel)
#elpd_diff se_diff
#model1  0.0       0.0  (GP) 
#model2 -0.4       0.4  (Non-GP)







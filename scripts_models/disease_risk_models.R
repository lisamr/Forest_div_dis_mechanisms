#Disease risk models assessed 1) at the plot-level for all host species, 2) at the plot-level for highly symptomatic species, and 3) at the individual-level for highly symptomatic species. 
#For each of the analyses, 3 competing models are contrasted. Each one includes various plot-level variables, in addition to: i) richness, ii) richness + tanoak BA + bay laurel BA, and iii) richness + community competency.

rm(list=ls())
library(rethinking)
library(tidyverse)
library(bayesplot)
library(loo)
library(flextable)
library(wesanderson)
library(cowplot)
theme_set(theme_classic())#set ggplot theme

#load data------------------------------------------------
plot_level <- read_csv('data/plot_level_data.csv')
tree_level <- read_csv('data/tree_level_data.csv')



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
    n=n) #total host plants
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
    BA = tree_level$basal_area_z #basal area of individual
  )
  return(datlist)
}

#additional variables used in the contrasting models
contrasting_models <- list(c('Richness_z'), c('Richness_z', 'BA_UMCAsqrt_z', 'BA_LIDEsqrt_z'), c('Richness_z', 'ccompsqrt_z'))
contrasting_models_bern <- list(NULL, c('BA_UMCAsqrt_z', 'BA_LIDEsqrt_z'), 'ccompsqrt_z')

# 1. Plot-level models, all hosts
datalists_plotlevel_all <- lapply(contrasting_models, function(x) make_datlist(x, all = T))

# 2. Plot-level models, highly symptomatic hosts only
datalists_plotlevel_HS <- lapply(contrasting_models, function(x) make_datlist(x, all = F))

# 3. individual-level models, highly symptomatic hosts
datalists_indlevel <- lapply(contrasting_models_bern, make_datlist_Bern)



#models---------------------------------------------------

#beta-binomial model
stan_betaBinom <- stan_model(file = 'Stan_models/beta_binomial.stan')

#Bernoulli model with varying diversity slopes and individual-level basal area
stan_bern <- stan_model(file = 'Stan_models/bernoulli_varyingslopeBA_nc.stan')




#fit models-----------------------------------------------

ffit <- function(Dlist, Stanmod, iter = 2000, chains = 4, ...){
  set.seed(2021)
  fit <- sampling(Stanmod, Dlist, iter = iter, chains = chains, cores = chains, seed = 2021, ...)
  return(fit)
}

#Run the models
post_plotlevel_all <- lapply(datalists_plotlevel_all, function(x) ffit(x, stan_betaBinom))
post_plotlevel_HS <- lapply(datalists_plotlevel_HS, function(x) ffit(x, stan_betaBinom))
post_indlevel <- lapply(datalists_indlevel, function(x) ffit(x, stan_bern, control = list(adapt_delta = .99)))


#Contrast model performance
loo_plotlevel_all <- lapply(post_plotlevel_all, loo)
loo_plotlevel_HS <- lapply(post_plotlevel_HS, loo)
loo_indlevel <- lapply(post_indlevel, loo)
loo_compare(loo_plotlevel_all)
loo_compare(loo_plotlevel_HS)
loo_compare(loo_indlevel)

#> loo_compare(loo_plotlevel_all)
#elpd_diff se_diff
#model2   0.0       0.0  
#model3 -13.9       3.6  
#model1 -20.1       5.8  
#> loo_compare(loo_plotlevel_HS)
#elpd_diff se_diff
#model2   0.0       0.0  
#model3 -16.8       4.7  
#model1 -23.7       6.7 
#> loo_compare(loo_indlevel)
#elpd_diff se_diff
#model2  0.0       0.0   
#model3 -3.0       2.2   
#model1 -4.7       2.5 


#plot goodness of fit for best performing model
post <- extract.samples(post_plotlevel_all[[2]])
ppc1 <- ppc_dens_overlay(datalists_plotlevel_all[[2]]$I, post$y_rep[1:100,]) + labs(title = 'plot-level, \nall hosts')
post <- extract.samples(post_plotlevel_HS[[2]])
ppc2 <- ppc_dens_overlay(datalists_plotlevel_HS[[2]]$I, post$y_rep[1:100,]) + labs(title = 'plot-level, \nhighly symtomatic hosts')
post <- extract.samples(post_indlevel[[2]])
ppc3 <- ppc_dens_overlay(datalists_indlevel[[2]]$I, post$y_rep[1:100,]) + labs(title = 'individual-level, \nhighly symtomatic hosts')
PPC <- cowplot::plot_grid(ppc1, ppc2, ppc3, 
                   nrow = 1)

#quick look at parameters
precis(post_plotlevel_all[[1]], pars = c(c('a0', 'beta', 'theta')), depth = 2, prob = .9)
precis(post_indlevel[[3]], pars = c(c('bBA', 'betap', 'betaS', 'zS')), depth = 3, prob = .9)


#export models--------------------------------------------
saveRDS(post_plotlevel_all, '../../../Box/Stan_model_outputs/Big_Sur/post_plotlevel_all.RDS')
saveRDS(post_plotlevel_HS, '../../../Box/Stan_model_outputs/Big_Sur/post_plotlevel_HS.RDS')
saveRDS(post_indlevel, '../../../Box/Stan_model_outputs/Big_Sur/post_indlevel2.RDS')

pdf('figures/PPC_diseaserisk.pdf', width = 6.5, height = 3.5)
PPC
dev.off()


#NOT USED!

#Disease risk models assessed 1a) infection prevalence at the plot-level for all host species, 1b) infection prevalence at the plot-level for highly symptomatic species, 2a) absolute density of infected hosts, 2b) absolute density of bay laurels only, and 3) infection risk at the individual-level for highly symptomatic species. 
#For each of the analyses, 3 competing models are contrasted. Each one includes various plot-level variables, in addition to: i) richness, ii) richness + tanoak BA + bay laurel BA, and iii) richness + community competency.

rm(list=ls())
library(rethinking)
library(tidyverse)
library(bayesplot)
library(loo)
library(flextable)
library(wesanderson)
library(cowplot)


#load data------------------------------------------------
plot_level <- read_csv('data/plot_level_data.csv') #data on 151 plots
tree_level <- read_csv('data/tree_level_data_HS.csv') #tree level data of highly susceptible species

#for model 2B-D, add a column to plot level data for...
# - number of infected bay laurels 
# - basal area of infected plants
# - basal area of infected bay laurel
newcols <- tree_level %>%
  group_by(BSPlotNumber) %>% 
  summarise(I_UMCA = sum(I[Species == 'UMCA']),
            BA_I = sum(basal_area[I == 1]),
            BA_I.UMCA = sum(basal_area[I == 1 & Species == 'UMCA']))

plot_level <- left_join(plot_level, newcols)


#create data lists----------------------------------------


#function for the negative binomial and gamma models (density of infected plants)
make_datlist_density <- function(predictors, response, df = plot_level){
  plot_level_mat <- as.matrix(select(df, plot_vars, predictors))
  datlist <- list(
    N=nrow(plot_level_mat), #number of plots
    K=ncol(plot_level_mat), #number of predictors
    y = response, #response (number of infections)
    X=plot_level_mat #covariate matrix
    ) 
  return(datlist)
}


#additional variables used in the contrasting models
contrasting_models <- list(c('Richness_z'), c('Richness_z', 'BA_UMCAsqrt_z', 'BA_LIDEsqrt_z'), c('Richness_z', 'ccompsqrt_z'))
contrasting_models_bern <- list(NULL, c('BA_UMCAsqrt_z', 'BA_LIDEsqrt_z'), 'ccompsqrt_z')


# 2A. Total number of infected individuals
datalists_number_I <- lapply(contrasting_models, function(x) make_datlist_density(x, plot_level$I_all))
# 2B. Total number of infected UMCA individuals
tmpdf <- plot_level %>% 
  filter(BA_I.UMCA != 0)
datalists_number_IUMCA <- lapply(contrasting_models, function(x) make_datlist_density(x, tmpdf$I_UMCA, tmpdf))
# 2C. Basal area of infected plants
datalists_BA_I <- lapply(contrasting_models, function(x) make_datlist_density(x, plot_level$BA_I))
# 2D. Basal area of infected UMCA plants. need to remove zeros (gamma dist != 0)
datalists_BA_I_UMCA <- lapply(contrasting_models, function(x) make_datlist_density(x, tmpdf$BA_I.UMCA, tmpdf))



#models---------------------------------------------------


#Negative binomial model #make sure to comment out Nsim and Xsim. I'll proabably just make another script without simulating data.
stan_NB <- stan_model(file = 'Stan_models/nindividuals_NB.stan')

#gamma model#make sure to comment out Nsim and Xsim. 
stan_gamma <- stan_model(file = 'Stan_models/basal_area_gamma.stan')


#fit models-----------------------------------------------

ffit <- function(Dlist, Stanmod, iter = 2000, chains = 4, ...){
  set.seed(2021)
  fit <- sampling(Stanmod, Dlist, iter = iter, chains = chains, cores = chains, seed = 2021, ...)
  return(fit)
}


#run the new models 
#NB models... undefined values:  log_lik[1], etc.
post_number_I <- lapply(datalists_number_I, function(x) ffit(x, stan_NB))#loo issues
post_number_IUMCA <- lapply(datalists_number_IUMCA, function(x) ffit(x, stan_NB))#restricted dataset including only 96 plots (where BA.I.UMCA > 0)
post_BA_I <- lapply(datalists_BA_I, function(x) ffit(x, stan_gamma))#loo issues
post_BA_IUMCA <- lapply(datalists_BA_I_UMCA, function(x) ffit(x, stan_gamma)) #restricted dataset including only 96 plots (where BA.I.UMCA > 0)



#new contrasts
loo_number_I <- lapply(post_number_I, loo)
loo_compare(loo_number_I)
loo_number_I <- lapply(post_number_IUMCA, loo)
loo_compare(loo_number_I)
loo_BA_I <- lapply(post_BA_I, loo)
loo_compare(loo_BA_I)
loo_BA_I <- lapply(post_BA_IUMCA, loo)
loo_compare(loo_BA_I)

precis(post_number_I[[3]], pars = c(c('a0', 'B', 'phi')), depth = 2, prob = .9)
precis(post_number_IUMCA[[1]], pars = c(c('a0', 'B', 'phi')), depth = 2, prob = .9)
precis(post_BA_I[[1]], pars = c(c('a0', 'B', 'va')), depth = 2, prob = .9)
precis(post_BA_IUMCA[[1]], pars = c(c('a0', 'B', 'va')), depth = 2, prob = .9)

post <- extract.samples(post_number_I[[2]])
ppc_dens_overlay(datalists_number_I[[2]]$y, post$y_rep[1:100,]) 
post <- extract.samples(post_number_IUMCA[[2]])
ppc_dens_overlay(datalists_number_IUMCA[[2]]$y, post$y_rep[1:100,]) 
post <- extract.samples(post_BA_I[[1]])
ppc_dens_overlay(datalists_BA_I[[1]]$y, post$y_rep[1:100,]) 

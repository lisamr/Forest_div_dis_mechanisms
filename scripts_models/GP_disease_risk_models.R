#Gaussian process disease models, bernoulli only. beta-binomial models don't converge with an extra GP intercept, unless the priors are inappropriately tight.

#just running the models with richness + key host densities since those were consistently the best performing models. compared posteriors with non-gp models--nearly identical posteriors and no difference in ELPD. 


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


# individual-level models, highly symptomatic hosts
datalists_indlevel <- make_datlist_Bern(c('BA_UMCAsqrt_z', 'BA_LIDEsqrt_z'))



#viz GP priors------------------------------------

#neighborhood distance cutoff in Haas et al. 2011 was 5km
N <- 10000
eta <- rnorm(N, 0, 1.5)
rho <- rnorm(N, 0, 3)
x <- seq(0, 10, by = .1)

tmp <- as.data.frame(sapply(1:N, function(i) (eta[i])^2 * exp(-.5/((rho[i])^2) * x^2))) %>% 
  mutate(x = x)
mediandf <- data.frame(x, median = apply(tmp[,-ncol(tmp)], 1, median))
tmpdf <- tmp %>% 
  select(V1:V100, x) %>% 
  pivot_longer(cols = contains("V")) 
ggplot(tmpdf, aes(x, value)) +
  geom_line(alpha = .2, aes(group = name)) +
  geom_line(data = mediandf, aes(x, median), color = 'red3', lwd = 1) #+coord_cartesian(xlim = c(0, 5))




#models---------------------------------------------------

#gaussian process Bernoulli model with varying diversity slopes and individual-level basal area
stan_bern <- stan_model(file = 'Stan_models/bernoulli_varyingslopeGPnc.stan')




#fit models-----------------------------------------------


#Run the models (just doing ones with richness + key host densities since they performed best without the GP)


#individual-level 
#CAUTION: TAKES A LOT OF TIME AND MEMORY. BUDGET A FEW HOURS!!!
set.seed(2021)
post_indlevel <- sampling(stan_bern, datalists_indlevel, iter = 2000, chains = 4, cores = 4, control = list(adapt_delta = .99), seed = 2021)


#export models--------------------------------------------
saveRDS(post_indlevel, '../../../Box/Stan_model_outputs/Big_Sur/post_indlevel_GP.RDS')


#read in models-----------------------------------
post_indlevel <- read_rds('../../../Box/Stan_model_outputs/Big_Sur/post_indlevel_GP.RDS')



#compare posterior with non-GP models-----------------------

#individual-level M2 (GP)
precis(post_indlevel, pars = c(c('bBA', 'betap', 'betaS', 'zS', 'eta_sq', 'rho_sq', 'sig_sq')), depth = 3, prob = .9)
#need to add zS to mean to get the values you report in the paper
samp_indlevel <- extract.samples(post_indlevel)
apply(samp_indlevel$betaS[,2] + samp_indlevel$zS[,2,], 2, HPDI, .9)
#[,1]      [,2]       [,3]      [,4]
#|0.9 -1.7948924 -1.159082 -0.6247144 0.1393561
#0.9|  0.2823888  1.169797  1.9877729 2.3907762

#individual-level M2 (non-GP)
post_indlevel_nonGP <- read_rds('../../../Box/Stan_model_outputs/Big_Sur/post_indlevel.RDS')
precis(post_indlevel_nonGP[[2]], pars = c(c('bBA', 'betap', 'betaS', 'zS')), depth = 3, prob = .9)
samp_indlevel <- extract.samples(post_indlevel_nonGP[[2]])
apply(samp_indlevel$betaS[,2] + samp_indlevel$zS[,2,], 2, HPDI, .9)
#[,1]      [,2]       [,3]      [,4]
#|0.9 -1.7510532 -1.087464 -0.6691851 0.1361247
#0.9|  0.2661887  1.098156  1.9489762 2.3431677



#contrast models with loo-------------------------

loo_indlevel <- lapply(list(post_indlevel, post_indlevel_nonGP[[2]]), loo) #all pareto k values < 0.7
loo_compare(loo_indlevel)
#> loo_compare(loo_indlevel)
#elpd_diff se_diff
#model2  0.0       0.0  (non-GP) 
#model1 -0.3       0.4  (GP)







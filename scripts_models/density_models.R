#models for assessing various measurements of density against richness in the 2 forest types. All models include predictors for plot richness and forest type. Response variables below:
# 1. Basal area of all species (gamma)
# 2. Basal area of tanoak and bay laurel, run seperately (hurdle model--bernoulli, gamma). 
# 3. Number of individuals from highly susceptible hosts (negative binomial)
# 4. Number of individuals from minimally susceptible hosts (negative binomial)

library(rethinking)
library(rstan)
library(tidyverse)
library(bayesplot)

rm(list=ls())
zscore <- function(x) (x-mean(x))/(2*sd(x))
unzscore <- function(z, x)  z *2*sd(x) + mean(x)



#load Stan models-------------------------------------------------

stan_gamma <- stan_model(file = 'Stan_models/basal_area_gamma.stan') 
stan_hurdle_singlesp <- stan_model(file = 'Stan_models/basal_area_gamma_hurd_singlespecies.stan') #for sinlge species 
stan_NB <- stan_model(file = 'Stan_models/nindividuals_NB.stan')




#Wrangle data-----------------------------------------------------

#load data
plot_level <- read_csv('data/plot_level_data.csv') #plot level data
tree_level <- read_csv('data/tree_level_data_all.csv') #tree level data on all species


#Calculate various forms of density
# 1. BA of all species
BA_all <- tree_level %>% 
  group_by(BSPlotNumber) %>% 
  summarise(BA_all = sum(basal_area)) %>% 
  left_join(plot_level %>% select(BSPlotNumber, ForestAllianceType, Richness_z))

# 2. BA of bay and tanoak, seperate dataframe contained in a list
top2 <- c("UMCA", "LIDE")
BA_spp_separated <- expand_grid(BSPlotNumber = plot_level$BSPlotNumber, 
  Species = top2) %>% 
  mutate(plotID = as.integer(as.factor(BSPlotNumber)),
         spID = as.integer(as.factor(Species))) %>% 
  left_join(tree_level %>% 
              filter(Species %in% top2) %>% 
              group_by(BSPlotNumber, Species) %>% 
              summarise(BA = sum(basal_area))) %>% 
  left_join(plot_level %>% select(BSPlotNumber, ForestAllianceType, Richness_z)) %>% 
  replace_na(list(BA = 0)) %>% 
  split(as.factor(.$Species))

# 3, 4. number of individuals from highly and mininally susceptible species
# already included in plot_level dataframe: plot_level$n_HS, plot_level$n_MS




#prep datalists----------------------------------------------------

# Datalists for GLM models (gamma, NB, hurdle model with single species)
#create new dataframe to predict posteriors to
richness_key <- plot_level %>% 
  distinct(Richness, Richness_z) %>% 
  arrange(Richness)
Xsim <- data.frame(
  Forest = c(rep(0, 8), rep(1, 6)),
  Richness_z = c(richness_key$Richness_z[1:8], richness_key$Richness_z[1:6]))


make_datalists_GLM <- function(response, covariates){
  list(
    N = nrow(covariates),
    K = ncol(covariates),
    y = response,
    X = covariates,
    Nsim = nrow(Xsim),
    Xsim = as.matrix(Xsim)
  )
}

covariates <-  cbind(plot_level$ForestAllianceType, plot_level$Richness_z) #the same for all of the following models
datalist_basal_area_all <- make_datalists_GLM(BA_all$BA_all, covariates) # 1. BA of all species
datalist_top2 <- lapply(BA_spp_separated, function(x) make_datalists_GLM(x$BA, covariates)) # 2. hurdle model run for each species separately
datalist_HS <- make_datalists_GLM(as.integer(plot_level$n_HS), covariates) # 3. Number of highly susceptible species
datalist_MS <- make_datalists_GLM(as.integer(plot_level$n_MS), covariates) # 4. Number of highly susceptible species




#run models------------------------------------------------------------

#fit models
fit_basal_area_all <- sampling(stan_gamma, datalist_basal_area_all, iter = 2000, chains = 4, cores = 4)
fit_top2 <- lapply(datalist_top2, function(x) sampling(stan_hurdle_singlesp, x, iter = 2000, chains = 4, cores = 4))
fit_HS <- sampling(stan_NB, datalist_HS, iter = 2000, chains = 4, cores = 4)
fit_MS <- sampling(stan_NB, datalist_MS, iter = 2000, chains = 4, cores = 4)

#posterior predictive checks
ppc <- function(datalist, fit, title, legend = F){
  post <- extract.samples(fit)
  p <- ppc_dens_overlay(datalist$y, post$y_rep[1:50,]) +
    labs(title = title) + 
    theme(title = element_text(size = 8))
  if(legend == F) p <- p + theme(legend.position = 'none')
  return(p)
}
ppc_int <- function(datalist, fit, title, legend = F){
  post <- extract.samples(fit)
  p <- ppc_intervals(datalist$y, post$y_rep) +
    labs(title = title, x = 'Plot index') + 
    theme(title = element_text(size = 8))
  if(legend == F) p <- p + theme(legend.position = 'none')
  return(p)
}

p1 = ppc(datalist_basal_area_all, fit_basal_area_all, 'Total basal area of all species')
p2 = ppc(datalist_HS, fit_HS, 'Number of individuals from \nhighly susceptible species')
p3 = ppc(datalist_MS, fit_MS, 'Number of individuals from \nminimally susceptible species')
ppc_sepspp = list(NULL)
spp_names <- c('Tanoak', 'Bay laurel')
for(i in 1:length(datalist_top2)){
  ppc_sepspp[[i]] <- ppc_int(datalist_top2[[i]], fit_top2[[i]], spp_names[i])
}
ppc_multiplot1 <- cowplot::plot_grid(plotlist = ppc_sepspp)
ppc_multiplot2 <- cowplot::plot_grid(p1, p2, p3, nrow = 1)


#Export---------------------------------------------------------------

#save models
saveRDS(fit_basal_area_all, '../../../Box/Stan_model_outputs/Big_Sur/fit_basal_area_all.RDS')
saveRDS(fit_top2, '../../../Box/Stan_model_outputs/Big_Sur/fit_hurdle_top2spp.RDS')
saveRDS(fit_HS, '../../../Box/Stan_model_outputs/Big_Sur/fit_nindividuals_HS.RDS')
saveRDS(fit_MS, '../../../Box/Stan_model_outputs/Big_Sur/fit_nindividuals_MS.RDS')


#save posterior predictive checks
pdf('figures/PPC_density_models1.pdf', width = 4.5, height = 2)
ppc_multiplot1
dev.off()
pdf('figures/PPC_density_models2.pdf', width = 6.5, height = 2)
ppc_multiplot2
dev.off()






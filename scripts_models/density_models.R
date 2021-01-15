#models for assessing various measurements of density against richness in the 2 forest types. All models include predictors for plot richness and forest type. Response variables below:
# 1. Basal area of all species (gamma)
# 2. Basal area of top 6 most common species (hurdle model--bernoulli, gamma). Needed to run redwood seperately beucase it's large size was causing problems. 
# 3. Number of individuals from highly susceptible hosts (negative binomial)
# 4. Number of individuals from minimally susceptible hosts (negative binomial)

library(rethinking)
library(rstan)
library(tidyverse)
library(cowplot)
library(wesanderson)
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

# 2. BA of each of the top 6 common species, seperate dataframe contained in a list
top6 <- c("UMCA", "LIDE", "QUAG", "QUPA", "ARME", 'SESE')
BA_spp_separated <- expand_grid(BSPlotNumber = plot_level$BSPlotNumber, 
  Species = top6) %>% 
  mutate(plotID = as.integer(as.factor(BSPlotNumber)),
         spID = as.integer(as.factor(Species))) %>% 
  left_join(tree_level %>% 
              filter(Species %in% top6) %>% 
              group_by(BSPlotNumber, Species) %>% 
              summarise(BA = sum(basal_area))) %>% 
  left_join(plot_level %>% select(BSPlotNumber, ForestAllianceType, Richness)) %>% 
  mutate(Richness_z = zscore(Richness)) %>% 
  replace_na(list(BA = 0)) %>% 
  split(as.factor(.$Species))

# 3, 4. number of individuals from highly and mininally susceptible species
# already included in plot_level dataframe: plot_level$n_HS, plot_level$n_MS




#prep datalists----------------------------------------------------

# Datalists for GLM models (gamma, NB, hurdle model with single species)
Xsim <- expand_grid( #create new dataframe to predict posteriors to
  Forest = c(0,1),
  Richness_z = seq(min(plot_level$Richness_z), max(plot_level$Richness_z), length.out = 10)
)
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
datalist_gamma <- make_datalists_GLM(BA_all$BA_all, covariates) # 1. BA of all species
datalist_hurd_sepspp <- lapply(BA_spp_separated, function(x) make_datalists_GLM(x$BA, covariates)) # 2. hurdle model run for each species separately
datalist_NB_HS <- make_datalists_GLM(as.integer(plot_level$n_HS), covariates) # 3. Number of highly susceptible species
datalist_NB_MS <- make_datalists_GLM(as.integer(plot_level$n_MS), covariates) # 4. Number of highly susceptible species





#run models------------------------------------------------------------

#fit models
fit_gamma <- sampling(stan_gamma, datalist_gamma, iter = 2000, chains = 4, cores = 4)
fit_hurd_sepspp <- lapply(datalist_hurd_sepspp, function(x) sampling(stan_hurdle_singlesp, x, iter = 2000, chains = 4, cores = 4))
fit_NB_HS <- sampling(stan_NB, datalist_NB_HS, iter = 2000, chains = 4, cores = 4)
fit_NB_MS <- sampling(stan_NB, datalist_NB_MS, iter = 2000, chains = 4, cores = 4)

#posterior predictive checks
ppc <- function(datalist, fit, title, legend = F){
  post <- extract.samples(fit)
  p <- ppc_dens_overlay(datalist$y, post$y_rep[1:50,]) +
    labs(title = title) + 
    theme(title = element_text(size = 8))
  if(legend == F) p <- p + theme(legend.position = 'none')
  return(p)
}

p1 = ppc(datalist_gamma, fit_gamma, 'Total basal area of all species')
p2 = ppc(datalist_NB_HS, fit_NB_HS, 'Number of individuals from \nhighly susceptible species')
p3 = ppc(datalist_NB_MS, fit_NB_MS, 'Number of individuals from \nminimally susceptible species')
ppc_sepspp = list(NULL)
for(i in 1:length(datalist_hurd_sepspp)){
  ppc_sepspp[[i]] <- ppc(datalist_hurd_sepspp[[i]], fit_hurd_sepspp[[i]], names(datalist_hurd_sepspp)[i])
}
ppc_multiplot1 <- cowplot::plot_grid(plotlist = ppc_sepspp)
ppc_multiplot2 <- cowplot::plot_grid(p1, p2, p3, nrow = 1)


#Export---------------------------------------------------------------

#save models
saveRDS(fit_gamma, '../../../Box/Stan_model_outputs/Big_Sur/fit_basal_area_all.RDS')
saveRDS(fit_hurd_sepspp, '../../../Box/Stan_model_outputs/Big_Sur/fit_hurdle_top6spp.RDS')
saveRDS(fit_NB_HS, '../../../Box/Stan_model_outputs/Big_Sur/fit_nindividuals_HS.RDS')
saveRDS(fit_NB_MS, '../../../Box/Stan_model_outputs/Big_Sur/fit_nindividuals_MS.RDS')

#save posterior predictive checks
pdf('figures/PPC_density_models1.pdf', width = 6.5, height = 3)
ppc_multiplot1
dev.off()
pdf('figures/PPC_density_models2.pdf', width = 6.5, height = 2)
ppc_multiplot2
dev.off()



#unused model------------------------------------------------------------
#I tried to run the hurdle model with many species together (only 5, excludes redwood), but it wasn't running well (inefficient, not amazing fit). I include it here just in case I'll need it in the future.

#data to be used
top5 <- c('UMCA', 'LIDE', 'QUAG', 'QUPA', 'ARME')
BA_top5 <- expand_grid(BSPlotNumber = plot_level$BSPlotNumber, Species = top5) %>% 
  mutate(plotID = as.integer(as.factor(BSPlotNumber)),
         spID = as.integer(as.factor(Species))) %>% 
  left_join(tree_level %>% 
              filter(Species %in% top5) %>% 
              group_by(BSPlotNumber, Species) %>% 
              summarise(BA = sum(basal_area))) %>% 
  left_join(plot_level %>% select(BSPlotNumber, ForestAllianceType, Richness)) %>% 
  mutate(Richness_z = zscore(Richness))
BA_top5$BA[is.na(BA_top5$BA)] <- 0

#data list
cov_top5 <- cbind(1, BA_top5$ForestAllianceType, BA_top5$Richness_z) #intercept + covariates
Xsim2 <- expand_grid( #data for posterior predictions 
  a0 = 1,
  Forest = c(0,1),
  Richness_z = seq(min(BA_top5$Richness_z), max(BA_top5$Richness_z), length.out = 10)
) %>% 
  slice(rep(row_number(), max(BA_top5$spID)))
datalist_hurdletop5 <- list(
  N = nrow(cov_top5),
  K = ncol(cov_top5),
  J = max(BA_top5$plotID),
  S = max(BA_top5$spID),
  y = BA_top5$BA,
  X = cov_top5,
  spID = BA_top5$spID,
  plotID = BA_top5$plotID,
  Nsim = nrow(Xsim2),
  Xsim = as.matrix(Xsim2),
  spID_sim = rep(1:max(BA_top5$spID), each = nrow(Xsim2)/max(BA_top5$spID))
)

#run model
stan_hurdle_manyspp <- stan_model(file = 'Stan_models/basal_area_gamma_hurd_species_studT.stan') #for top 5 species all together
fit_hurdletop5 <- sampling(stan_hurdle_manyspp, datalist_hurdletop5, iter = 2000, chains = 4, cores = 4, control = list(adapt_delta = .99, max_treedepth = 15), seed = 2021)
#No divergent transitions or other warnings

#check out
precis(fit_hurdletop5, depth =3, pars = c('B_bar', 'B','Bz_bar', 'Bz', 'sigma_p', 'sigma_pz', 'sd_j', 'va'), prob = .95)
pairs(fit_hurdletop5, pars = 'B_bar')
traceplot(fit_hurdletop5, pars = c('B_bar', 'Bz_bar', 'sigma_p', 'sigma_pz', 'sd_j', 'va'))
post_species <- extract.samples(fit_hurdletop5)
ppc_dens_overlay(datalist_hurdletop5$y, post_species$y_rep[1:50,]) + coord_cartesian(xlim = c(0,2))
ppc_intervals_grouped(datalist_hurdletop5$y, post_species$y_rep, group = datalist_hurdletop5$spID)

saveRDS(fit_hurdletop5, '../../../Box/Stan_model_outputs/Big_Sur/fit_hurdletop5.RDS')








#AGGREGATED SPECIES GROUPS-------

#groups = susceptibles, all hosts, all species. Run models seperately for each group. Include forest type as an intercept. 

#Basal area: gamma likelihood
#No. of Invididuals: Negative binomial likelihood


#load stan model
stan_BA_GAM <- stan_model(file = 'scripts_2020/Stan_models/Density_vs_diversity/basal_area_gamma.stan')

#create new dataframe to predict posteriors to
newdf <- expand_grid(
  Forest = c(0,1),
  Richness_z = seq(-1.8, 3.1, length.out = 10)
)
newdf_I <- expand_grid( #includes interaction term
  Forest = c(0,1),
  Richness_z = seq(-1.8, 3.1, length.out = 10)
) %>% 
  mutate(R_F = Richness_z*Forest)

#prep data for aggregated species.
prepdata_AGG <- function(df, predictors, density = 'BA'){
  cov <- df %>% #covariate matrix
    select(predictors) %>% 
    as.matrix() 
  dlist <- list( 
    N = nrow(df),
    K = ncol(cov),
    y = pull(df, density),
    X = cov,
    Nsim = nrow(newdf),
    Xsim = newdf
  )  
  return(dlist)
}

d_spp <- prepdata_AGG(BA %>% filter(Species == 'all_species'), c('ForestAllianceType', 'Richness_z'))
d_sus <- prepdata_AGG(BA %>% filter(Species == 'susceptibles'),  c('ForestAllianceType', 'Richness_z'))
d_hosts <- prepdata_AGG(BA %>% filter(Species == 'all_hosts'), c('ForestAllianceType', 'Richness_z'))


#run model
fit_spp <- sampling(stan_BA_GAM, data = d_spp, chains = 4, cores = 4, iter=2000)
fit_sus <- sampling(stan_BA_GAM, data = d_sus, chains = 4, cores = 4, iter=2000)
fit_hosts <- sampling(stan_BA_GAM, data = d_hosts, chains = 4, cores = 4, iter=2000)


precis(fit_spp, prob = .95, depth = 2, pars = c('a0', 'B', 'va'))
precis(fit_sus, prob = .95, depth = 2, pars = c('a0', 'B', 'va'))
precis(fit_hosts, prob = .95, depth = 2, pars = c('a0', 'B', 'va'))


#goodness of fits (might want to allow variance to vary by forest type?)
ppc <- function(datalist, fit){
  post <- extract.samples(fit)
  p1 <- ppc_dens_overlay(datalist$y, post$y_rep[1:50,]) 
  p2 <- ppc_intervals_grouped(datalist$y, post$y_rep, group = datalist$X[,'ForestAllianceType'])
  return(cowplot::plot_grid(p1, p2, nrow = 2))
}
ppc(d_spp, fit_spp)
ppc(d_sus, fit_sus)
ppc(d_hosts, fit_hosts)


#posterior predictive plots
f_PP <- function(post, newdat = newdf){ 
  pmu <- apply(post, 2, mean)
  pCI <- apply(post, 2, HPDI, .95)
  newdat %>% 
    mutate(mean = pmu, lower = pCI[1,], upper = pCI[2,],
           ForestAllianceType = as.factor(Forest), 
           Richness = unzscore(Richness_z, covariates$Richness))
}

post_spp <- extract.samples(fit_spp)
post_sus <- extract.samples(fit_sus)
post_hosts <- extract.samples(fit_hosts)

post_sppdf <- f_PP(post_spp$mu_sim)
post_susdf <- f_PP(post_sus$mu_sim)
post_hostsdf <- f_PP(post_hosts$mu_sim)

plot_BA <- function(post_df, dat, Title, legend.loc = 'none'){ 
  ggplot(post_df, aes(Richness, mean, group = ForestAllianceType)) +
    geom_line(aes(color = ForestAllianceType)) +
    geom_ribbon(aes(ymin = lower, ymax = upper, fill = ForestAllianceType), alpha = .2) +
    geom_jitter(data = dat, 
                aes(Richness, BA, color = as.factor(ForestAllianceType)), 
                alpha = .7, width = .15, height = 0) +
    scale_color_manual(values = c("#46ACC8", "#DD8D29"), labels = c('Mixed evergreen','Redwood')) +
    scale_fill_manual(values = c("#46ACC8", "#DD8D29"), labels = c('Mixed evergreen','Redwood')) +
    labs(y = expression(paste("Basal area (", m^2, ')')),
         title = Title, 
         color = 'Forest type', fill = 'Forest type') +
    theme(legend.position = legend.loc, legend.justification='center', legend.direction='vertical', plot.title = element_text(hjust = 0.5)) 
}

p_all <- plot_BA(post_sppdf, filter(BA, Species == 'all_species'), "All species", legend.loc = 'right')
p_all
plot_BA(post_susdf, filter(BA, Species == 'susceptibles'), 'Susceptible species')
plot_BA(post_hostsdf, filter(BA, Species == 'all_hosts'), 'All host species')



#Do for No. of individuals------

#load model
stan_NB <- stan_model(file = 'scripts_2020/Stan_models/Density_vs_diversity/density_NB.stan')

#prep data
inds$Species %>% unique
d_hosts2 <- prepdata_AGG(inds %>% filter(Species == 'all_hosts'), c('ForestAllianceType', 'Richness_z'), 'ninds')
d_sus2 <- prepdata_AGG(inds %>% filter(Species == 'susceptibles'), c('ForestAllianceType', 'Richness_z'), 'ninds')
d_oth2 <- prepdata_AGG(inds %>% filter(Species == 'other_hosts'), c('ForestAllianceType', 'Richness_z'), 'ninds')

#run models
fit_hosts2 <- sampling(stan_NB, data = d_hosts2, chains = 4, cores = 4, iter=2000)
fit_sus2 <- sampling(stan_NB, data = d_sus2, chains = 4, cores = 4, iter=2000)
fit_oth2 <- sampling(stan_NB, data = d_oth2, chains = 4, cores = 4, iter=2000)

#check out
precis(fit_hosts2, depth = 2, pars =c('a0', 'B', 'phi'), prob = .95)
precis(fit_sus2, depth = 2, pars =c('a0', 'B', 'phi'), prob = .95)
precis(fit_oth2, depth = 2, pars =c('a0', 'B', 'phi'), prob = .95)
ppc(d_sus, fit_sus)
ppc(d_hosts, fit_hosts)
ppc(d_oth2, fit_oth2)

#plot posterior predictions
plot_ninds <- function(post_df, dat, Title, legend.loc = 'none'){ 
  ggplot(post_df, aes(Richness, mean, group = ForestAllianceType)) +
    geom_line(aes(color = ForestAllianceType)) +
    geom_ribbon(aes(ymin = lower, ymax = upper, fill = ForestAllianceType), alpha = .2) +
    geom_jitter(data = dat, 
                aes(Richness, ninds, color = as.factor(ForestAllianceType)), 
                alpha = .7, width = 0, height = 0) +
    scale_color_manual(values = c("#46ACC8", "#DD8D29"), labels = c('Mixed evergreen','Redwood')) +
    scale_fill_manual(values = c("#46ACC8", "#DD8D29"), labels = c('Mixed evergreen','Redwood')) +
    labs(y = 'No. of individuals',
         title = Title, 
         color = 'Forest type', fill = 'Forest type') +
    theme(legend.position = legend.loc, legend.justification='center', legend.direction='horizontal', plot.title = element_text(hjust = 0.5)) 
}

post_sus2 <- extract.samples(fit_sus2)
post_oth2 <- extract.samples(fit_oth2)
post_hosts2 <- extract.samples(fit_hosts2)

#plot it!
p_nsus <- plot_ninds(f_PP(post_sus2$mu_sim), filter(inds, Species == 'susceptibles'), 'Susceptible hosts', legend.loc = 'bottom')
p_noth <- plot_ninds(f_PP(post_oth2$mu_sim), filter(inds, Species == 'other_hosts'), 'Effectively non-susceptible hosts') 
p_nall <- plot_ninds(f_PP(post_hosts2$mu_sim), filter(inds, Species == 'all_hosts'), 'All hosts') 






#Hurdle model approach------

#load stan model (bernoulli, gamma processes)
stan_BA_GAM.HURD <- stan_model(file = 'scripts_2020/Stan_models/Density_vs_diversity/basal_area_gamma_hurd.stan')

#prep data
prepdatahurd <- function(df, predictors, newDF = newdf){
  cov <- df %>% #covariate matrix. using same ones for binomial and gamma process
    select(predictors) %>% 
    as.matrix() 
  dlist <- list( 
    N = nrow(df),
    K = ncol(cov),
    y = df$BA,
    X = cov, 
    Nsim = nrow(newdf),
    Xsim = newDF
  )  
  return(dlist)
}

d_UMCA <- prepdatahurd(BA %>% filter(Species == 'UMCA'),  c('ForestAllianceType', 'Richness_z', 'R_F'), newdf_I)
d_LIDE <- prepdatahurd(BA %>% filter(Species == 'LIDE'),  c('ForestAllianceType', 'Richness_z', 'R_F'), newdf_I)
d_SESE <- prepdatahurd(BA %>% filter(Species == 'SESE'),  c('ForestAllianceType', 'Richness_z'))
d_QUAG <- prepdatahurd(BA %>% filter(Species == 'QUAG'),  c('ForestAllianceType', 'Richness_z'))
d_QUPA <- prepdatahurd(BA %>% filter(Species == 'QUPA'),  c('ForestAllianceType', 'Richness_z'))
d_ARME <- prepdatahurd(BA %>% filter(Species == 'ARME'),  c('ForestAllianceType', 'Richness_z'))

#fit models
fit_UMCA <- sampling(stan_BA_GAM.HURD, data = d_UMCA, chains = 4, cores = 4, iter=2000)
fit_LIDE <- sampling(stan_BA_GAM.HURD, data = d_LIDE, chains = 4, cores = 4, iter=2000)
fit_SESE <- sampling(stan_BA_GAM.HURD, data = d_SESE, chains = 4, cores = 4, iter=2000)
fit_QUAG <- sampling(stan_BA_GAM.HURD, data = d_QUAG, chains = 4, cores = 4, iter=2000)
fit_QUPA <- sampling(stan_BA_GAM.HURD, data = d_QUPA, chains = 4, cores = 4, iter=2000)
fit_ARME <- sampling(stan_BA_GAM.HURD, data = d_ARME, chains = 4, cores = 4, iter=2000)

#check out
precis(fit_UMCA, prob = .95, depth = 2, pars = c('a0', 'theta', 'B', 'Bz', 'va'))
precis(fit_LIDE, prob = .95, depth = 2, pars = c('a0', 'theta', 'B', 'Bz', 'va'))
precis(fit_SESE, prob = .95, depth = 2, pars = c('a0', 'theta', 'B', 'Bz', 'va'))
precis(fit_QUAG, prob = .95, depth = 2, pars = c('a0', 'theta', 'B', 'Bz', 'va'))
precis(fit_QUPA, prob = .95, depth = 2, pars = c('a0', 'theta', 'B', 'Bz', 'va'))
precis(fit_ARME, prob = .95, depth = 2, pars = c('a0', 'theta', 'B', 'Bz', 'va'))

#goodness of fits
ppc(d_UMCA, fit_UMCA)
ppc(d_LIDE, fit_LIDE)
ppc(d_SESE, fit_SESE) #definitely not as good
ppc(d_QUAG, fit_QUAG)
ppc(d_QUPA, fit_QUPA)
ppc(d_ARME, fit_ARME)

#posterior predictive plots
post_UMCA <- extract.samples(fit_UMCA)
post_LIDE <- extract.samples(fit_LIDE)
post_SESE <- extract.samples(fit_SESE)
post_QUAG <- extract.samples(fit_QUAG)
post_QUPA <- extract.samples(fit_QUPA)
post_ARME <- extract.samples(fit_ARME)

#check out richness effects with hpdi for umca, lide
HPDI(post_UMCA$B[,2], .95) #NS
HPDI(post_UMCA$Bz[,2], .95) #neg
HPDI(post_LIDE$B[,2], .95) #NS
HPDI(post_LIDE$Bz[,2], .95) #NS

#First, estimated basal area (when present)
postdf_UMCA <- f_PP(post_UMCA$mu_sim)
postdf_LIDE <- f_PP(post_LIDE$mu_sim)
postdf_SESE <- f_PP(post_SESE$mu_sim)
postdf_QUAG <- f_PP(post_QUAG$mu_sim)
postdf_QUPA <- f_PP(post_QUPA$mu_sim)
postdf_ARME <- f_PP(post_ARME$mu_sim)

p1a = plot_BA(postdf_UMCA, filter(BA, Species == 'UMCA', BA > 0), 'Bay laurel') + scale_y_continuous(limits = c(0,3.5))
p2a = plot_BA(postdf_LIDE, filter(BA, Species == 'LIDE', BA > 0), 'Tanoak') + scale_y_continuous(limits = c(0,3.5))
p3a = plot_BA(postdf_SESE, filter(BA, Species == 'SESE', BA > 0), 'Redwood')
p4a = plot_BA(postdf_QUAG, filter(BA, Species == 'QUAG', BA > 0), 'CLO')
p5a = plot_BA(postdf_QUPA, filter(BA, Species == 'QUPA', BA > 0), 'shreves')
p6a = plot_BA(postdf_ARME, filter(BA, Species == 'ARME', BA > 0), 'madrone')

#Next, prob of occurence
postdf_UMCA2 <- f_PP(post_UMCA$p_sim)
postdf_LIDE2 <- f_PP(post_LIDE$p_sim)
postdf_SESE2 <- f_PP(post_SESE$p_sim)
postdf_QUAG2 <- f_PP(post_QUAG$p_sim)
postdf_QUPA2 <- f_PP(post_QUPA$p_sim)
postdf_ARME2 <- f_PP(post_ARME$p_sim)


plot_BA_Occ <- function(post, Title, legend.loc = 'none'){
  ggplot(post, aes(Richness, mean, group = Forest)) +
    geom_line(aes(color = as.factor(Forest)))+
    geom_ribbon(aes(ymin = lower, ymax = upper, fill = as.factor(Forest)), alpha = .2) +
    scale_color_manual(values = c("#46ACC8", "#DD8D29"), labels = c('Mixed evergreen','Redwood')) +
    scale_fill_manual(values = c("#46ACC8", "#DD8D29"), labels = c('Mixed evergreen','Redwood'))+
    labs(y = 'P(occurrence)',
         title = Title) +
    scale_y_continuous(limits = c(0,1)) +
    theme(legend.position = legend.loc, legend.title = element_blank(), legend.justification='center', legend.direction='horizontal', plot.title = element_text(hjust = 0.5)) 
}

p1b = plot_BA_Occ(postdf_UMCA2, 'Bay laurel')
p2b = plot_BA_Occ(postdf_LIDE2, 'Tanoak')
p3b = plot_BA_Occ(postdf_SESE2, 'Redwood')
p4b = plot_BA_Occ(postdf_QUAG2, 'CLO')
p5b = plot_BA_Occ(postdf_QUPA2, 'shreves')
p6b = plot_BA_Occ(postdf_ARME2, 'madrone')

plot_grid(plotlist = list(p1a, p1b, p2a, p2b, p3a, p3b, p4a, p4b, p5a, p5b, p6a, p6b),
          nrow = 3)
plot_grid(plotlist = list(p1a, p2a, p3a, p4a, p5a, p6a),
          nrow = 3)






#run one mega model for all species------

#there will be a random intercept for species and slopes varying by species in both the gamma and bernoulli models. Also a random effect for plots. Even with a MV-studentT, can't include redwood--too much of an outlier and screws with the fits for all the other species.

#load stan model (bernoulli, lognormal processes, species random effects)
stan_BA_Species <- stan_model(file = 'scripts_2020/Stan_models/Density_vs_diversity/basal_area_gamma_hurd_species.stan')
stan_BA_Species_studT <- stan_model(file = 'scripts_2020/Stan_models/Density_vs_diversity/basal_area_gamma_hurd_species_studT.stan')

#prep data 
tmp <- BA %>% 
  filter(Species %in% c('UMCA', 'LIDE', 'QUAG', 'QUPA')) %>% 
  select(ForestAllianceType, Richness_z, Richness, BA, BSPlotNumber, Species) %>% 
  droplevels() %>% 
  mutate(plotID = as.integer(as.factor(BSPlotNumber)),
         spID = as.integer(as.factor(Species)))
tmpcov <- tmp %>% 
  mutate(a0 = 1L) %>% 
  select(a0, ForestAllianceType, Richness_z) %>% 
  as.matrix()
newdf2 <- expand_grid(
  a0 = 1,
  Forest = c(0,1),
  Richness_z = seq(-1.8, 3.1, length.out = 10)
) %>% 
  slice(rep(row_number(), 5))
d_species <- list(
  N = nrow(tmpcov),
  K = ncol(tmpcov),
  J = max(tmp$plotID),
  S = max(tmp$spID),
  y = tmp$BA,
  X = tmpcov, 
  spID = tmp$spID,
  plotID = tmp$plotID,
  Nsim = nrow(newdf2),
  Xsim = newdf2,
  spID_sim = rep(1:max(tmp$spID), each = nrow(newdf2)/max(tmp$spID))
)
str(d_species)

#fit models
fit_species <- sampling(stan_BA_Species, data = d_species, chains = 4, cores = 4, iter=2000, control=list(adapt_delta=0.9))
fit_species_studT <- sampling(stan_BA_Species_studT, data = d_species, chains = 4, cores = 4, iter=2000, control=list(adapt_delta=0.9))


#check out
#fit <- fit_species_studT
fit <- fit_species
sp_names <- tmp %>% distinct(Species, spID)
#precis(fit, depth =3, pars = c('B_bar', 'B','Bz_bar', 'Bz', 'sigma_p', 'sigma_pz', 'sd_j', 'va', 'nu', 'nuz'), prob = .95)
precis(fit, depth =3, pars = c('B_bar', 'B','Bz_bar', 'Bz', 'sigma_p', 'sigma_pz', 'sd_j', 'va'), prob = .95)
traceplot(fit, pars = c('B_bar', 'Bz_bar', 'sigma_p', 'sigma_pz', 'sd_j', 'va'))
post_species <- extract.samples(fit)
ppc_dens_overlay(d_species$y, post_species$y_rep[1:50,]) + coord_cartesian(xlim = c(0,2))
ppc_intervals_grouped(d_species$y, post_species$y_rep, group = d_species$spID) 


#HPDI estimates of richness effects on occurence and BA
apply(post_species$B[,,3], 2, HPDI, .95)
apply(post_species$Bz[,,3], 2, HPDI, .95)
apply(post_species$B_bar, 2, HPDI, .95)
apply(post_species$Bz_bar, 2, HPDI, .95)

#posterior predictive plots
#first basal area, when present
postdf_species <- f_PP(post_species$mu_sim, newdf2) %>% 
  mutate(spID = d_species$spID_sim) %>% 
  left_join(sp_names) %>% 
  select(Forest:Richness, Species) %>% 
  bind_rows(postdf_SESE %>% mutate(Species = 'SESE'))

psp1 <- ggplot(postdf_species, aes(Richness, mean, group = ForestAllianceType)) +
  geom_line(aes(color = ForestAllianceType)) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = ForestAllianceType), alpha = .2) +
  geom_jitter(data =BA %>% filter(Species %in% c('UMCA', 'LIDE', 'QUAG', 'QUPA', 'SESE'), BA>0), 
              aes(Richness, BA, color = as.factor(ForestAllianceType)), 
              alpha = .8, width = .15, height = 0) +
  scale_shape_manual(values = c(19, 4)) +
  facet_wrap(~Species,nrow = 1, scales = 'free_y') +
  scale_color_manual(values = c("#46ACC8", "#DD8D29"), labels = c('Mixed evergreen','Redwood')) +
  scale_fill_manual(values = c("#46ACC8", "#DD8D29"), labels = c('Mixed evergreen','Redwood')) +
  labs(y = expression(paste("Basal area (", m^2, ')'))) +
  theme(legend.position = 'none') 

#next occurrence
postdf_species2 <- f_PP(post_species$pz_sim, newdf2) %>% 
  mutate(spID = d_species$spID_sim) %>% 
  left_join(sp_names) %>% 
  select(Forest:Richness, Species) %>% 
  bind_rows(postdf_SESE2 %>% mutate(Species = 'SESE'))
psp2 <- ggplot(postdf_species2, aes(Richness, mean, group = ForestAllianceType)) +
  geom_line(aes(color = ForestAllianceType)) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = ForestAllianceType), alpha = .2) +
  #geom_jitter(data = tmp %>% mutate(PA = as.integer(BA>0)), aes(Richness, PA, color = as.factor(ForestAllianceType)), height = 0,  width = .15, alpha = .4) +
  facet_grid(cols = vars(Species), scales = 'free_y') +
  scale_color_manual(values = c("#46ACC8", "#DD8D29"), labels = c('Mixed evergreen','Redwood')) +
  scale_fill_manual(values = c("#46ACC8", "#DD8D29"), labels = c('Mixed evergreen','Redwood')) +
  scale_y_continuous(limits = c(0,1)) +
  labs(y = "P(occurrence)") +
  theme(legend.position = 'none') 

p_5spp <- plot_grid(psp1, psp2, nrow = 2)
#



#replot with just bay laurel and tanoak. might want to include the results from this model isntead of the first hurdle model.
p_baytoak1 <- ggplot(postdf_species %>% filter(Species %in% c("UMCA", 'LIDE')) , aes(Richness, mean, group = ForestAllianceType)) +
  geom_line(aes(color = ForestAllianceType)) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = ForestAllianceType), alpha = .2) +
  geom_jitter(data =BA %>% filter(Species %in% c('UMCA', 'LIDE'), BA>0), 
              aes(Richness, BA, color = as.factor(ForestAllianceType)), 
              alpha = .8, width = .15, height = 0) +
  facet_wrap(~Species, labeller=labeller(Species = c('LIDE'='Tanoak', 'UMCA' = 'Bay laurel')) ) +
  scale_color_manual(values = c("#46ACC8", "#DD8D29"), labels = c('Mixed evergreen','Redwood')) +
  scale_fill_manual(values = c("#46ACC8", "#DD8D29"), labels = c('Mixed evergreen','Redwood')) +
  labs(y = expression(paste("Basal area (", m^2, ')'))) +
  theme(legend.position = 'none', strip.background = element_blank(), strip.text = element_text(size=10)) 

p_baytoak2 <- ggplot(postdf_species2 %>% filter(Species %in% c("UMCA", 'LIDE')), aes(Richness, mean, group = ForestAllianceType)) +
  geom_line(aes(color = ForestAllianceType)) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = ForestAllianceType), alpha = .2) +
  facet_wrap(~Species, labeller=labeller(Species = c('LIDE'='Tanoak', 'UMCA' = 'Bay laurel'))) +
  scale_color_manual(values = c("#46ACC8", "#DD8D29"), labels = c('Mixed evergreen','Redwood')) +
  scale_fill_manual(values = c("#46ACC8", "#DD8D29"), labels = c('Mixed evergreen','Redwood')) +
  scale_y_continuous(limits = c(0,1)) +
  labs(y = "P(occurrence)") +
  theme(legend.position = 'none', strip.background = element_blank(), strip.text = element_text(size=10)) 


#Figures for paper------

#put plot together for main paper
#BA of bay, tanoak, all species
leg <- get_legend(p_all)
pg1 <- plot_grid(p1a, p2a,p_all + theme(legend.position = 'none'), nrow = 1)
pg2 <- plot_grid(p1b, p2b, leg, nrow=1)
pg_both <- plot_grid(pg1, pg2, nrow = 2)

ggsave(plot = pg_both, filename = 'output/figures/final_density/baytoakall.png', width = 173, height = 110, units = 'mm', dpi = 600)

pg3 <- plot_grid(p_baytoak1, p_baytoak2, nrow = 2)
pg4 <- plot_grid(p_all + theme(axis.title.y = element_blank(), legend.position = 'none', title = element_text(size=8)), leg, nrow = 2, rel_heights = c(1, .9), scale = .9)
pg_both2 <- plot_grid(pg3, pg4, nrow = 1, rel_widths = c(1.7,1))
ggsave(plot = pg_both2, filename = 'output/figures/final_density/baytoakall2.png', width = 110, height = 80, units = 'mm', dpi = 600)
ggsave(plot = pg_both2, filename = 'output/figures/final_density/baytoakall2.pdf', width = 110, height = 80, units = 'mm', dpi = 600)

#Supplemental
#top 5 species (SESE model run alone), no. of all hosts
ggsave(plot = p_5spp, filename = 'output/figures/final_density/species_BA.png', width = 173, height = 80, units = 'mm', dpi = 600)
leg2 <- get_legend(p_nsus)
p_ninds <- plot_grid(p_nsus + theme(legend.position = 'none'), p_noth, leg2, nrow = 3, rel_heights = c(1,1,.2))
p_ninds
ggsave(p_ninds, filename = 'output/figures/final_density/ninds_groups.png', width = 82, height = 120, units = 'mm', dpi = 600)













#lognormal with a transformation?-----

#I'd like to model the overall change in density and avoid the hurdle model as much as possible. maybe adding a small constant can get rid of the zeros?

#ultimately, I feel like the constant is too arbitrary and can lead to different results. I think I'm going with the hurdle model :/
stan_logN <- stan_model(file = 'scripts_2020/Stan_models/Density_vs_diversity/basal_area_gamma.stan')


#data to predict to
newdf <- expand_grid(
  Forest = c(0,1),
  Richness = seq(-2,3,length.out = 10),
)

prepdatagamma <- function(df, predictors){
  cov <- df %>% #covariate matrix
    select(predictors) %>% 
    as.matrix() 
  dlist <- list( 
    N = nrow(df),
    Nsim = nrow(newdf),
    K = ncol(cov),
    y = df$BA,
    X = cov,
    simdf = as.matrix(newdf)
  )  
  return(dlist)
}



#constant to add
C <- .01
#proportion of plots with BA <= constant
1 - sum(d_UMCA$y >= C) / sum(d_UMCA$y > 0) 
1 - sum(d_LIDE$y >= C) / sum(d_LIDE$y > 0) 

#create data lists
d_UMCA2 <- prepdatagamma(BA %>% filter(Species == 'Bay laurel'), c('ForestAllianceType', 'Richness'))
d_UMCA2$y <- d_UMCA2$y + C
d_LIDE2 <- prepdatagamma(BA %>% filter(Species == 'Tanoak'), c('ForestAllianceType', 'Richness'))
d_LIDE2$y <- d_LIDE2$y + C

#fit models
fitlogN_UMCA <- sampling(stan_logN, data = d_UMCA2, chains = 1, cores = 1, iter=2000)
precis(fitlogN_UMCA, prob = .95, depth = 2, pars = c('a0', 'B',  'va'))
ppc(d_UMCA2, fitlogN_UMCA)

fitlogN_LIDE <- sampling(stan_logN, data = d_LIDE2, chains = 1, cores = 1, iter=2000)
precis(fitlogN_LIDE, prob = .95, depth = 2, pars = c('a0', 'B',  'va'))
ppc(d_UMCA2, fitlogN_LIDE)














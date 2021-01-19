
#unused model

#I tried to run the hurdle model with many species together (all 6, including redwood), but it wasn't running well (inefficient, not amazing fit). I include it here just in case I'll need it in the future.
library(rethinking)
library(rstan)
library(tidyverse)
library(bayesplot)
theme_set(theme_classic())

rm(list=ls())
zscore <- function(x) (x-mean(x))/(2*sd(x))
unzscore <- function(z, x)  z *2*sd(x) + mean(x)
#posterior predictive plots
f_PP <- function(post, newdat = Xsim){ 
  pmu <- apply(post, 2, mean)
  pCI <- apply(post, 2, HPDI, .9)
  newdat %>% 
    mutate(mean = pmu, lower = pCI[1,], upper = pCI[2,])
}
#for summarizing arrays of posteriors. the tidybayes function was behaving well.
medianHPDI <- function(post_matrix, width = .9){ 
  y <- apply(post_matrix, 2, median)
  yCI <- t(apply(post_matrix, 2, HPDI, width))
  cbind(y, yCI)
}



#load data---------
plot_level <- read_csv('data/plot_level_data.csv') #plot level data
tree_level <- read_csv('data/tree_level_data_all.csv') #tree level data on all species


#data wrangling-------
#data to be used
top6 <- c('UMCA', 'LIDE', 'QUAG', 'QUPA', 'ARME', 'SESE')
BA_top6 <- expand_grid(BSPlotNumber = plot_level$BSPlotNumber, Species = top6) %>% 
  mutate(plotID = as.integer(as.factor(BSPlotNumber)),
         spID = as.integer(as.factor(Species))) %>% 
  left_join(tree_level %>% 
              filter(Species %in% top6) %>% 
              group_by(BSPlotNumber, Species) %>% 
              summarise(BA = sum(basal_area))) %>% 
  left_join(plot_level %>% select(BSPlotNumber, ForestAllianceType, Richness)) %>% 
  mutate(Richness_z = zscore(Richness))
BA_top6$BA[is.na(BA_top6$BA)] <- 0
species_names <- BA_top6 %>% 
  distinct(Species, spID) %>% 
  arrange(spID) %>% pull(Species)




#data list-------
cov_top6 <- cbind(1, BA_top6$ForestAllianceType, BA_top6$Richness_z) #intercept + covariates
Xsim2 <- expand_grid( #data for posterior predictions 
  a0 = 1,
  Forest = c(0,1),
  Richness_z = seq(min(BA_top6$Richness_z), max(BA_top6$Richness_z), length.out = 10)
) %>% 
  slice(rep(row_number(), max(BA_top6$spID)))
datalist_hurdletop6 <- list(
  N = nrow(cov_top6),
  K = ncol(cov_top6),
  J = max(BA_top6$plotID),
  S = max(BA_top6$spID),
  y = BA_top6$BA,
  X = cov_top6,
  spID = BA_top6$spID,
  plotID = BA_top6$plotID,
  Nsim = nrow(Xsim2),
  Xsim = as.matrix(Xsim2),
  spID_sim = rep(1:max(BA_top6$spID), each = nrow(Xsim2)/max(BA_top6$spID))
)

#run model----------
stan_hurdle_manyspp <- stan_model(file = 'Stan_models/basal_area_gamma_hurd_species_studT.stan') #for top 6 species all together
fit_hurdletop6 <- sampling(stan_hurdle_manyspp, datalist_hurdletop6, iter = 2000, chains = 4, cores = 4, control = list(adapt_delta = .99, max_treedepth = 15), seed = 2021)
#No divergent transitions or other warnings, but chains sampled at really different rates and ppc shows not amazing fit still

saveRDS(fit_hurdletop6, '../../../Box/Stan_model_outputs/Big_Sur/fit_hurdleglmm_top6.RDS')


#check out----------
fit_hurdletop6 <- read_rds('../../../Box/Stan_model_outputs/Big_Sur/fit_hurdleglmm_top6.RDS')

precis(fit_hurdletop6, depth =3, pars = c('B_bar', 'B','Bz_bar', 'Bz', 'sigma_p', 'sigma_pz', 'sd_j', 'va'), prob = .9)

traceplot(fit_hurdletop6, pars = c('B_bar', 'Bz_bar', 'sigma_p', 'sigma_pz', 'sd_j', 'va'))
post_species <- extract.samples(fit_hurdletop6)
ppc_dens_overlay(datalist_hurdletop6$y, post_species$y_rep[1:50,]) + coord_cartesian(xlim = c(0,2))
ppc_intervals_grouped(datalist_hurdletop6$y, post_species$y_rep, group = species_names[datalist_hurdletop6$spID])
ggsave('figures/PPC_hurdleGLMM_top6.pdf')

#look at dens overlay for each species
ppc_dens_group <- function(i) {
  use <- datalist_hurdletop6$spID == i
  ppc_dens_overlay(datalist_hurdletop6$y[use], post_species$y_rep[1:50,use])
}
plots <- lapply(1:5, ppc_dens_group)
plot_grid(plotlist = plots) 


#look at posteriors--------


#basal area of top 6 species, hurdle glmm
#species in order: arme, lide, quag, qupa, umca
post_top6 <- extract.samples(fit_hurdletop6)
str(post_top6) 

#average richness effect on basal area
medianHPDI(cbind(post_top6$B_bar[,3]))
#y      |0.9      0.9|
#[1,] -0.1442722 -0.3155603 0.04520396

#average richness effect on occurrence
medianHPDI(cbind(post_top6$Bz_bar[,3]))
#y      |0.9      0.9|
#[1,] -1.163769 -1.858278 -0.346282

#species-specific richness effects on BA (difference from mean)
medianHPDI(post_top6$B[,,3])
#             y       |0.9         0.9|
#[1,] -0.09229648 -0.2732276  0.136832828
#[2,] -0.13228434 -0.3148570  0.026282211
#[3,] -0.19575401 -0.5144508  0.055641073
#[4,] -0.14863110 -0.3741549  0.046787453
#[5,] -0.14760822 -0.2842807 -0.005827378

#species-specific richness effects on occurence (difference from mean)
medianHPDI(post_top6$Bz[,,3])
#             y       |0.9         0.9|
#[1,] -1.9053401 -2.639215 -1.1301344
#[2,] -0.7815364 -1.425846 -0.1217705
#[3,] -0.7368778 -1.366579 -0.0907069
#[4,] -1.3738479 -1.986729 -0.7772336
#[5,] -1.9591398 -2.859596 -1.0868373

#sigma parameters for intercept and 2 covariates
medianHPDI(post_top6$sigma_p) 
#             y        |0.9      0.9|
#[1,] 0.2852823 0.092726537 0.5619644
#[2,] 0.4622918 0.177122634 0.8887390
#[3,] 0.1118801 0.002441554 0.3525422



#density of top 5 species using hurdle GLMM
#solid or dashed lines
CI_BA <- t(apply(post_top6$B[,,3] + post_top6$B_bar[,3], 2, HPDI, .9))
crosses_zero_BA <- data.frame(linetypeBA = unname(ifelse(CI_BA[,1] * CI_BA[,2] > 0, 1, 2)), spID = 1:nrow(CI_BA)) 
CI_Occ <- t(apply(post_top6$Bz[,,3] + post_top6$Bz_bar[,3], 2, HPDI, .9))
crosses_zero_Occ <- data.frame(linetypeOcc = unname(ifelse(CI_Occ[,1] * CI_Occ[,2] > 0, 1, 2)), spID = 1:nrow(CI_Occ)) 

#dataframe used for posterior predictions
Xsim2 <- expand_grid( #data for posterior predictions 
  a0 = 1,
  ForestAllianceType = c(0,1),
  Richness_z = seq(min(BA_top6$Richness_z), max(BA_top6$Richness_z), length.out = 10)
) %>% 
  slice(rep(row_number(), max(BA_top6$spID))) %>% 
  mutate(Richness = unzscore(Richness_z, BA_top6$Richness),
         spID = datalist_hurdletop6$spID_sim) %>% 
  left_join(crosses_zero_BA) %>% 
  left_join(crosses_zero_Occ) %>% 
  left_join(BA_top6 %>% distinct(Species, spID))





#Basal area
post_df <- f_PP(post_top6$mu_sim, Xsim2)
p1 <- ggplot(post_df, aes(Richness, mean, group = ForestAllianceType)) +
  geom_line(aes(color = as.factor(ForestAllianceType)), lty = post_df$linetypeBA) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = as.factor(ForestAllianceType)), alpha = .2) +
  geom_jitter(data = BA_top6 %>% filter(BA>0), aes(Richness, BA, color = as.factor(ForestAllianceType)), alpha = .7, width = .15, height = 0) +
  facet_wrap(~Species, scales = 'free_y', nrow = 1)+
  scale_color_manual(values = c("#46ACC8", "#DD8D29")) +
  scale_fill_manual(values = c("#46ACC8", "#DD8D29")) +
  theme(legend.position = 'none')

#occurence
post_df2 <- f_PP(post_top6$pz_sim, Xsim2)
p2 <- ggplot(post_df2, aes(Richness, mean, group = ForestAllianceType)) +
  geom_line(aes(color = as.factor(ForestAllianceType)), lty = post_df$linetypeOcc) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = as.factor(ForestAllianceType) ), alpha = .2) +
  facet_grid(~Species)  +
  scale_color_manual(values = c("#46ACC8", "#DD8D29")) +
  scale_fill_manual(values = c("#46ACC8", "#DD8D29")) +
  theme(legend.position = 'none')

cowplot::plot_grid(p1, p2, nrow = 2)
ggsave('figures/hurdleGLMM_top6_unused.pdf', width = 6, height = 3)


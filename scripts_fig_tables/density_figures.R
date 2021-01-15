#Figures for density models
library(tidyverse)
library(wesanderson)
rm(list = ls())

theme_set(theme_classic())#set ggplot theme
unzscore <- function(z, x)  z *2*sd(x) + mean(x)


#import model fits
fit_basal_area_all <- readRDS('../../../Box/Stan_model_outputs/Big_Sur/fit_basal_area_all.RDS')
fit_top6 <- readRDS('../../../Box/Stan_model_outputs/Big_Sur/fit_hurdle_top6spp.RDS')
fit_HS <- readRDS('../../../Box/Stan_model_outputs/Big_Sur/fit_nindividuals_HS.RDS')
fit_MS <- readRDS('../../../Box/Stan_model_outputs/Big_Sur/fit_nindividuals_MS.RDS')

#import data
plot_level <- read_csv('data/plot_level_data.csv')
tree_level <- read_csv('data/tree_level_data_all.csv')





#Data wrangling----------------------------------------------------------------

#dataframe used for posterior predictions
Xsim <- expand_grid( 
  ForestAllianceType = as.factor(c(0,1)),
  Richness_z = seq(min(plot_level$Richness_z), max(plot_level$Richness_z), length.out = 10)
) %>% 
  mutate(Richness = unzscore(Richness_z, plot_level$Richness))

#Data for plotting
# 1. BA of all species
BA_all <- tree_level %>% 
  group_by(BSPlotNumber) %>% 
  summarise(BA = sum(basal_area)) %>% 
  left_join(plot_level %>% select(BSPlotNumber, ForestAllianceType, Richness, Richness_z))

# 2. BA and occurence of each of the top 6 common species, seperate dataframe contained in a list
top6 <- c("UMCA", "LIDE", "QUAG", "QUPA", "ARME", 'SESE')
spp_separated <- expand_grid(BSPlotNumber = plot_level$BSPlotNumber, 
                                Species = top6) %>% 
  mutate(plotID = as.integer(as.factor(BSPlotNumber)),
         spID = as.integer(as.factor(Species))) %>% 
  left_join(tree_level %>% 
              filter(Species %in% top6) %>% 
              group_by(BSPlotNumber, Species) %>% 
              summarise(BA = sum(basal_area))) %>% 
  left_join(plot_level %>% select(BSPlotNumber, ForestAllianceType, Richness, Richness_z)) %>% 
  replace_na(list(BA = 0))

BA_spp_separated <- spp_separated %>% 
  filter(BA > 0) %>% 
  split(as.factor(.$Species))
occurence_spp_separated <- spp_separated %>% 
  mutate(Present = as.numeric(BA > 0))%>% 
  split(as.factor(.$Species))





#get posteriors + summaries------------------------------------------------------

#basal area of all species
post_all <- extract.samples(fit_basal_area_all) 
tidybayes::median_hdci(post_all$B[,2], .width = .9)
#y       ymin      ymax .width .point .interval
#1 -0.02501688 -0.1618321 0.1067124    0.9 median      hdci

#basal area of top 6 species
post_spp <- lapply(fit_top6, extract.samples) 
sapply(post_spp, function(x) tidybayes::median_hdci(x$B[,2], .width = .9)) %>% t
#y           ymin       ymax        .width .point   .interval
#ARME 0.08486771  -0.1976215 0.3581919   0.9    "median" "hdci"   
#LIDE -0.08482223 -0.2917723 0.1351496   0.9    "median" "hdci"   
#QUAG -0.4623458  -0.831266  -0.09885835 0.9    "median" "hdci"   
#QUPA -0.1821664  -0.5287347 0.1858088   0.9    "median" "hdci"   
#SESE -0.06469372 -0.4776241 0.343455    0.9    "median" "hdci"   
#UMCA -0.1524417  -0.3136554 0.0263978   0.9    "median" "hdci"   

#probability of occurrence for top 6 species
sapply(post_spp, function(x) tidybayes::median_hdci(x$Bz[,2], .width = .9)) %>% t
#y          ymin       ymax       .width .point   .interval
#ARME -1.848206  -2.576201  -1.142565  0.9    "median" "hdci"   
#LIDE -0.5257579 -1.126284  0.04436544 0.9    "median" "hdci"   
#QUAG -0.5155101 -1.117727  0.07879479 0.9    "median" "hdci"   
#QUPA -1.286342  -1.894791  -0.7171985 0.9    "median" "hdci"   
#SESE -0.1720305 -0.9321409 0.5661623  0.9    "median" "hdci"   
#UMCA -1.846928  -2.634078  -1.102452  0.9    "median" "hdci"   


#Highly susceptible individuals
post_HS <- extract.samples(fit_HS) 
tidybayes::median_hdci(post_HS$B[,2], .width = .9)
#y       ymin      ymax .width .point .interval
#1 0.1345237 -0.1027169 0.3515275    0.9 median      hdci


#Minimally susceptible individuals
post_MS <- extract.samples(fit_MS) 
tidybayes::median_hdci(post_MS$B[,2], .width = .9)
#y     ymin     ymax .width .point .interval
#1 1.344762 1.054569 1.679182    0.9 median      hdci







#Figure for Basal area or No. individuals--------------------------------------------

#posterior predictive plots
f_PP <- function(post, newdat = Xsim){ 
  pmu <- apply(post, 2, mean)
  pCI <- apply(post, 2, HPDI, .95)
  newdat %>% 
    mutate(mean = pmu, lower = pCI[1,], upper = pCI[2,])
}

#does the richness slope cross zero?
cross_zero <- function(posterior, occurrence = NULL){
  if(!is.null(occurrence)) CI <- HPDI(posterior$Bz[,2], .9)
  else CI <- HPDI(posterior$B[,2], .9)
  crosses_zero <- unname(ifelse(CI[1] * CI[2] > 0, 1, 2)) 
  return(crosses_zero)
}

#plot figures
plot_figures <- function(posterior, dat, response, Title, y.axis, legend.loc = 'none', occurrence = NULL){ 
  if(!is.null(occurrence)) {
    post_df <- f_PP(posterior$p_sim) } 
  else{post_df <- f_PP(posterior$mu_sim)} 
  response <- enquo(response)
  
  p <- ggplot(post_df, aes(Richness, mean, group = ForestAllianceType)) +
    geom_line(aes(color = ForestAllianceType), lty = cross_zero(posterior, occurrence)) +
    geom_ribbon(aes(ymin = lower, ymax = upper, fill = ForestAllianceType), alpha = .2) +
    scale_color_manual(values = c("#46ACC8", "#DD8D29"), labels = c('Mixed evergreen','Redwood')) +
    scale_fill_manual(values = c("#46ACC8", "#DD8D29"), labels = c('Mixed evergreen','Redwood')) +
    labs(y = y.axis,
         title = Title, 
         color = 'Forest type', fill = 'Forest type') +
    theme(legend.position = legend.loc, legend.justification='center', plot.title = element_text(hjust = 0.5)) 
  
  if(is.null(occurrence)){
    p <- p + geom_jitter(data = dat, 
                  aes(Richness, !!response, color = as.factor(ForestAllianceType)),
                  alpha = .7, width = .15, height = 0)
  } 
  return(p)
}

#all species basal area
basal_area_y <- expression(paste("Basal area (", m^2, ')'))
plot_BA_all <- plot_figures(post_all, BA_all, BA, 'All species', basal_area_y, legend.loc = 'right')

#top 6 species basal area
post_spp <- lapply(fit_top6, extract.samples)
spp_names <- c('Madrone', 'Tanoak', 'Coast live oak', 'Shreve oak', 'Redwood', 'Bay laurel')
BA_spp_plots <- list(NULL)
for(i in 1:length(spp_names)){
  BA_spp_plots[[i]] <- plot_figures(post_spp[[i]], BA_spp_separated[[i]], BA, spp_names[i], basal_area_y)
}
plot_grid(plotlist = BA_spp_plots) 

#top 6 species occurrence probability
Occurrence_spp_plots <- list(NULL)
for(i in 1:length(spp_names)){
  Occurrence_spp_plots[[i]] <- plot_figures(post_spp[[i]], Title = spp_names[i], y.axis = 'P(Occurrence)', occurrence = T)
}
plot_grid(plotlist = Occurrence_spp_plots) 


#Number of Highly and Minimally susceptible species
pHS <- plot_figures(post_HS, plot_level, n_HS, 'Highly susceptible species', 'No. of individuals')
pMS <- plot_figures(post_MS, plot_level, n_MS, 'Minimally susceptible species', 'No. of individuals', legend.loc = 'bottom')






#Export figures for paper--------

#Basal area + occurence figure
pgrid1 <- plot_grid(
  BA_spp_plots[[2]]+ xlab(""),
  BA_spp_plots[[6]] + xlab("") + theme(axis.title.y = element_blank()),
  plot_BA_all + xlab("") + theme(legend.position = 'none', axis.title.y = element_blank()),
  nrow = 1
)
pgrid2 <- plot_grid(
  Occurrence_spp_plots[[2]]+ xlab(""), 
  Occurrence_spp_plots[[6]]+ theme(axis.title.y = element_blank()),
  get_legend(plot_BA_all),
  nrow = 1
)
Fig1 <- plot_grid(pgrid1, pgrid2, nrow = 2)


#number of individuals
Fig2 <- plot_grid(pHS + xlab(''), 
          pMS + theme(legend.position = 'none'),
          get_legend(pMS),
          rel_heights = c(1,1,.15),
          nrow =3)


#supplemental figure of all of the top 6 species 
FigS2 <- plot_grid(plotlist = c(BA_spp_plots, Occurrence_spp_plots), nrow = 2)


#save figures
ggsave('figures/density_BA_Occurence.pdf', Fig1, width = 6.5, height = 3.5, device = 'pdf')
ggsave('figures/density_ninds.pdf', Fig2, width = 3, height = 4, device = 'pdf')
ggsave('figures/density_top6.pdf', FigS2, width = 10, height = 4, device = 'pdf')
 
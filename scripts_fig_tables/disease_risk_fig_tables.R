#creating supplementary tables with posterior values for the various disease risk models. Also figures of key plot-level and species-level variables. 

rm(list=ls())
library(rethinking)
library(tidyverse)
library(bayesplot)
library(flextable)
library(wesanderson)
library(cowplot)
theme_set(theme_classic())#set ggplot theme

#import model results 
post_plotlevel_all <- readRDS('../../../Box/Stan_model_outputs/Big_Sur/post_plotlevel_all.RDS')
post_plotlevel_HS <- readRDS('../../../Box/Stan_model_outputs/Big_Sur/post_plotlevel_HS.RDS')
post_indlevel <- readRDS('../../../Box/Stan_model_outputs/Big_Sur/post_indlevel.RDS')



#Table of Posterior estimates---------------------------

#functions for summarizing binomial models
summarize_post_binom <- function(posterior, modelname, predictors){
  x <- c('Forest type**', 'Year*', 'Precipitation', 'Potential solar induction', 'Host vegetative coverage', predictors)
  postdf <- with(posterior, cbind(a0, beta = beta, theta) %>% data.frame) %>% rename_at(vars(starts_with('V')), ~x) 
  summary_post <- postdf %>% pivot_longer(everything()) %>% group_by(name) %>% tidybayes::median_hdci(.width = .90) %>% mutate_at(vars(value:.upper), round, 2) %>% select(variable = name, median = value, lower = .lower, upper = .upper) %>% mutate(model = modelname) 
  return(summary_post)
}

make_table <- function(summary_list, Title, parnames){
  #put into a table
  tab_binom <- bind_rows(summary_list) %>%
    mutate('median (90% CI)' = paste0(median, ' \n(', lower, ', ', upper, ')')) %>% 
    pivot_wider(id_cols = variable, values_from = c('median (90% CI)'), names_from = model) %>% 
    replace(is.na(.), '–') 
  tab_binom <- tab_binom[match(parnames, tab_binom$variable),]
  
  #make it pretty for a Word table
  flextab_binom <- tab_binom %>% 
    flextable() %>% 
    align(part = "all") %>% # left align
    set_caption(caption = Title) %>% 
    add_footer_lines(values = c("*Estimate for plots sampled in 2007.", "**Estimate for redwood forests.") ) %>% 
    font(fontname = "Times new roman", part = "all") %>% 
    fontsize(size = 11, part = "body") %>% 
    bold(part = 'header') 
  flextab_binom <- width(flextab_binom, width = 1.0)
  return(flextab_binom)
}

#define some stuff for making a nice table
model_names <- list('Richness only', 'Richness + density of key hosts', 'Richness + community competency')
model_variables <- list('Richness', c('Richness', 'Bay laurel basal area', 'Tanoak basal area'), c('Richness', 'Community competency'))
parnames_binom <- c("a0",  "Year*", "Forest type**", "Host vegetative coverage","Precipitation",  "Potential solar induction", "Richness", "Bay laurel basal area", "Tanoak basal area",  "Community competency", "theta")


#Beta-binomial models for all hosts
draws_plotlevel_all <- lapply(post_plotlevel_all, extract.samples)
summary_plotlevel_all <- list(NULL)
for(i in 1:3){
  summary_plotlevel_all[[i]] <- summarize_post_binom(draws_plotlevel_all[[i]], model_names[[i]], model_variables[[i]])
}
flextab_binom_all <- make_table(summary_plotlevel_all, 'Posterior estimates (median log-odds, 90% HDPI): Infection prevalence aggregated among all hosts', parnames_binom)

#Beta-binomial models for highly symptomatic hosts
draws_plotlevel_HS <- lapply(post_plotlevel_HS, extract.samples)
summary_plotlevel_HS <- list(NULL)
for(i in 1:3){
  summary_plotlevel_HS[[i]] <- summarize_post_binom(draws_plotlevel_HS[[i]], model_names[[i]], model_variables[[i]])
}
flextab_binom_HS <- make_table(summary_plotlevel_HS, 'Posterior estimates (median log-odds, 90% HDPI): Infection prevalence aggregated among highly susceptible hosts', parnames_binom)

#Bernoulli models
summarize_post_bern <- function(posterior, modelname, predictors){
  if(is.null(predictors) ){
    x <- c('Tanoak intercept', 'Coast live oak intercept', 'Shreve oak intercept', 'Bay laurel intercept', 'Tanoak-specific richness', 'Coast live oak-specific richness', 'Shreve oak-specific richness', 'Bay laurel-specific richness', 'Basal area of individual',  'Forest type**', 'Year*', 'Precipitation', 'Potential solar induction', 'Host vegetative coverage', 'Mean effect of richness', 'Species SD', 'Richness slope SD', 'Plot SD')
  }else{
    x <- c('Tanoak intercept', 'Coast live oak intercept', 'Shreve oak intercept', 'Bay laurel intercept', 'Tanoak-specific richness', 'Coast live oak-specific richness', 'Shreve oak-specific richness', 'Bay laurel-specific richness', 'Basal area of individual',  'Forest type**', 'Year*', 'Precipitation', 'Potential solar induction', 'Host vegetative coverage', predictors, 'Mean effect of richness', 'Species SD', 'Richness slope SD', 'Plot SD')
  }
  sp_ints <- posterior$betaS[,1]  + posterior$zS[,1,]
  sp_rich <- posterior$betaS[,2]  + posterior$zS[,2,]
  postdf <- data.frame(cbind(sp_ints, sp_rich, posterior$bBA, posterior$betap, posterior$betaS[,2], posterior$sigmaS, posterior$sdPlot))
  #postdf <- with(posterior, cbind(aSp, bd, bBA, betap, bbar, sigma_pars, sdPlot) %>% data.frame)
  colnames(postdf) <- x 
  
  #summarize posteriors
  summary_postdf <- postdf %>% 
    pivot_longer(everything()) %>% group_by(name) %>% tidybayes::median_hdci(.width = .9) %>% mutate_at(vars(value:.upper), round, 2) %>% select(variable = name, median = value, lower = .lower, upper = .upper) %>% mutate(model = modelname) 
  return(summary_postdf)
}

#Bernoulli models for highly symptomatic hosts
draws_indlevel <- lapply(post_indlevel, extract.samples)
model_variables_bern <- list(NULL, c('Bay laurel basal area', 'Tanoak basal area'), c('Community competency'))
summary_indlevel <- list(NULL)
for(i in 1:3){
  summary_indlevel[[i]] <- summarize_post_bern(draws_indlevel[[i]], model_names[[i]], model_variables_bern[[i]])
}
parnames_bern <- c('Tanoak intercept', 'Coast live oak intercept', 'Shreve oak intercept', 'Bay laurel intercept','Basal area of individual', 'Forest type**', 'Year*', 'Precipitation', 'Potential solar induction', 'Host vegetative coverage', 'Mean effect of richness', 'Tanoak-specific richness', 'Coast live oak-specific richness', 'Shreve oak-specific richness', 'Bay laurel-specific richness', 'Bay laurel basal area', 'Tanoak basal area', 'Community competency', 'Species SD', 'Richness slope SD', 'Plot SD')
flextab_bern <- make_table(summary_indlevel, 'Posterior estimates (median log-odds, 90% HDPI): Individual-level infection risk for susceptible species', parnames_bern)


save_as_docx(flextab_binom_all, path = 'tables/diseaserisk_S1.docx')
save_as_docx(flextab_binom_HS, path = 'tables/diseaserisk_S2.docx')
save_as_docx(flextab_bern, path = 'tables/diseaserisk_S3.docx')



# Coefplot of Posterior estimates----


#Plot-level coefs----

#plot ODDS
plotlevel_coefs_ODDS <- function(summary_list, Title){
  plotdf <- bind_rows(summary_list) %>% 
    filter(variable %in% c('Richness', "Bay laurel basal area", 'Tanoak basal area', 'Community competency')) %>% 
    mutate_at(vars(median:upper), exp) %>% 
    mutate(significant = ifelse(lower<1 & upper>1, F, T)) %>% 
    mutate_at(vars(model, variable), as.factor) %>% 
    mutate(variable = fct_relevel(variable, 'Richness', "Bay laurel basal area", 'Tanoak basal area', 'Community competency'),
           variable = fct_recode(variable, "Bay laurel \nbasal area" = "Bay laurel basal area", 'Tanoak \nbasal area' = 'Tanoak basal area', 'Community \ncompetency' = 'Community competency'),
           model = fct_relevel(model, 'Richness only', "Richness + density of key hosts"))
  
  ggplot(plotdf, aes(median, fct_rev(variable), color = model)) +
    geom_vline(xintercept = 1, lty = 1, lwd = .3, color = grey(.7)) +
    geom_pointrange(aes(xmin = lower, xmax = upper, 
                        lty = significant==F,
                        pch = significant==F), 
                    position = position_dodge2(.55, reverse = T), 
                    fatten = 2.5) +
    scale_shape_manual(guide = 'none', values = c(19, 1)) +
    scale_linetype_manual(guide = 'none', values = c('solid', '11')) +
    scale_color_manual(values = wes_palettes$FantasticFox1[c(1,2,3)], guide = guide_legend(reverse = F)) +
    labs(color = 'Model', x = 'Odds ratio', y = 'Coefficient', title = Title) +
    theme(legend.position = 'bottom', legend.justification='left', legend.direction='horizontal') 
}

p1 <- plotlevel_coefs_ODDS(summary_plotlevel_all, 'Plot infection prevalence \n(all hosts)')
p2 <- plotlevel_coefs_ODDS(summary_plotlevel_HS, 'Plot infection prevalence \n(commonly symptomatic hosts)')
for(i in 1:3){
  summary_indlevel[[i]]$variable <- recode(summary_indlevel[[i]]$variable,  'Mean effect of richness' = "Richness" )
}
p3 <- plotlevel_coefs_ODDS(summary_indlevel, 'Individual infection risk \n(commonly symptomatic hosts)')




#species-level coefs-----

#species intercept and diversity slopes
plotdf_slope_ODDS <- bind_rows(summary_indlevel) %>% 
  filter(grepl('specific', variable)) %>%  
  mutate_at(vars(median:upper), exp) %>% 
  mutate(significant = ifelse(lower<1 & upper>1, F, T)) %>% 
  mutate_at(vars(model), as.factor) %>% 
  mutate(model = fct_relevel(model, 'Richness only', "Richness + density of key hosts"),
         species = as.factor(gsub("-.*$", '', variable)))

plotdf_int <- bind_rows(summary_indlevel) %>% 
  filter(grepl('intercept', variable)) %>% 
  mutate_at(vars(model), as.factor) %>% 
  mutate_at(vars(median:upper), inv_logit) %>% 
  mutate(model = fct_relevel(model, 'Richness only', "Richness + density of key hosts"),
         species = as.factor(gsub(" intercept.*$", '', variable)))


p4 <- plotdf_int %>% 
  ggplot(., aes(median, fct_rev(species), color = model)) +
  geom_pointrange(aes(xmin = lower, xmax = upper), position = position_dodge2(.55, reverse = T), fatten = 2.5, pch = 19) +
  scale_color_manual(values = wes_palettes$FantasticFox1[c(1,2,3)], guide = guide_legend(reverse = F)) +
  labs(color = 'Model', x = 'Probability of infection', y = 'Species', title = 'Species-specific \ninfection rates') +
  theme(legend.position = 'none') +
  scale_x_continuous(limits = c(0,1))

 
p5_ODDS <- plotdf_slope_ODDS %>% 
  ggplot(., aes(median, fct_rev(species), group = model, color = model)) +
  geom_vline(xintercept = 1, lty = 1, lwd = .3, color = grey(.7)) +
  geom_pointrange(aes(xmin = lower, xmax = upper, 
                      lty = significant==F,
                      pch = significant==F),
                      position = position_dodge2(.55, reverse = T), 
                      fatten = 2.5) +
  scale_linetype_manual(values = c('solid', '11')) +
  scale_shape_manual(values = c(19, 1)) +
  scale_color_manual(values = wes_palettes$FantasticFox1[c(1,2,3)], guide = guide_legend(reverse = T)) +
  labs(color = 'Model', x = 'Odds ratio',  title = 'Species-specific \nrichness effects') +
  theme(legend.position = 'none', legend.justification='left', legend.direction='horizontal', 
        axis.title.y = element_blank(), axis.text.y = element_blank())
p5_ODDS



#Put plots together-----


#PLOT LEVEL EFFECTS

# arrange the plots in a single row
prow <- plot_grid(
  p1 + 
    theme(legend.position="none", plot.title = element_text(hjust = 0.5, size = 10)) +
    xlab(''),
  p2 + 
    theme(legend.position="none", 
          plot.title = element_text(hjust = 0.5, size = 10),
          axis.title.y = element_blank(), axis.text.y = element_blank()),
  p3 + 
    theme(legend.position="none", 
          plot.title = element_text(hjust = 0.5, size = 10),
          axis.title.y = element_blank(), axis.text.y = element_blank())+
    xlab(''),
  labels = c("A", "B", "C"),
  hjust = c(-6.5, -.5, -.5), vjust = 2,
  nrow = 1, rel_widths = c(1,.7,.7), scale = .9
)

# add the legend 
legend_h <- get_legend(p1 + theme(legend.justification="center", legend.text = element_text(size = 8), legend.title = element_text(size = 9))) # extract legend 
plot_finished <- plot_grid(prow, legend_h, nrow = 2, rel_heights = c(1, .1))
plot_finished

ggsave2(plot_finished, width = 180, height = 110, units = 'mm', dpi = 600, filename = 'figures/disease_risk_plotlevel.pdf')


#SPECIES-LEVEL EFFECTS 
# arrange the plots in a single row
prow2 <- plot_grid(
  p4 + 
    theme(legend.position="none", plot.title = element_text(hjust = 0.5, size = 10)),
  p5_ODDS + 
    theme(legend.position="none", 
          plot.title = element_text(hjust = 0.5, size = 10),
          axis.title.y = element_blank(), axis.text.y = element_blank()),
  labels = c("A", "B"),
  hjust = c(-6.5, -.5), vjust = 2,
  nrow = 1, rel_widths = c(1,.8), scale = .9
)
plot_finished2 <- plot_grid(prow2, legend_h, nrow = 2, rel_heights = c(1, .1))
plot_finished2

ggsave2(plot_finished2, width = 140, height = 110, units = 'mm', dpi = 600, filename = 'figures/disease_risk_specieslevel.pdf')


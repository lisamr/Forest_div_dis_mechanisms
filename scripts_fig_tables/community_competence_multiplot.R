#create multi-plot figure for showing how competency scales with diversity
#1. nestedness
#2. sporulation (competence)
#3. contribuiton to community compeence
#4. community competence

rm(list = ls())

pal <- wesanderson::wes_palettes$FantasticFox1[c(1,2,3,5)]
library(tidyverse)
library(cowplot)
library(ggrepel)

theme_set(theme_classic())#set ggplot theme

#NESTEDNESS---------------------------------------
PA_mat <- read.table('data/PA_matrix.csv')
plotdf <- read_csv('data/plot_level_data.csv')
plotdf$ForestAllianceType <- as.factor(plotdf$ForestAllianceType + 1)

forests <- plotdf %>% 
  distinct(BSPlotNumber, ForestAllianceType) 

#order the matrix to be packed as most as possible
Y <- names(sort(colSums(PA_mat), T))
X <- names(sort(rowSums(PA_mat), T))
PA_mat2 <- PA_mat[X, Y]

#get forest ID for each plot
FID <- left_join(
  data.frame(BSPlotNumber = colnames(PA_mat2), 
             plotID = 1:ncol(PA_mat2)), 
  forests) 
#rename the columns so they order properly in the plot
colnames(PA_mat2) <- 1:ncol(PA_mat2)

#plot with ggplot
PAdf <- as.data.frame(PA_mat2) %>%
  mutate(species = rownames(.), 
         speciesID = 1:nrow(.)) %>% 
  pivot_longer(-c(species, speciesID), names_to = 'plotID') %>% 
  mutate(plotID = as.integer(plotID)) %>% 
  left_join(FID) %>% 
  mutate(value2 = ifelse(value==1, ForestAllianceType, 0))
print(PAdf, width = Inf)

P1 <- ggplot(PAdf, aes(plotID, rev(speciesID))) +
  geom_tile(aes(fill = as.factor(value2))) +
  scale_x_continuous(expand = c(.05, 0)) +
  scale_fill_manual(values = c('white', pal[3], pal[1])) +
  labs(x='Plot \noccurrence', y='') + 
  theme(legend.position = 'none') +
  theme(axis.text.x = element_text(color = 'white'), axis.ticks.x = element_line(color=NA), panel.border = element_rect(colour = "black", fill = NA))
P1

#latin names
rev(unique(PAdf$species))
latin <- c('S. mexicanus', 'R. ilicifolia', 'P. menziesii', 'Morella cerifera', 'Lathyrus vestitus', 'P. racemosa', 'C. papillosus', 'Ceanothus sp.', 'Q. lobata', 'Q. kelloggii', 'Arctostaphylos sp.', 'A. californica', 'R. californica', 'P. ponderosa', 'L. hispidula', 'C. oliganthus', 'H. arbutifolia', 'T. diversilobum', 'A. macrophyllum', 'Q. chrysolepis', 'A. menziesii', 'Q. parvula', 'Q. agrifolia', 'S. sempervirens', 'N. densiflorus', 'U. californica')



#Sporulation--------------------------------------

#Sporulation values come from the Competency publication estimates (Rosenthal et al. in press, Plant Disease)
spores <- read_csv('data/sporangia_estimates.csv')

#fill in unknowns with 70.7 (median of known estimates)
not_measured <- X[!X %in% spores$species] 
spores2 <- bind_rows(spores,
                     data.frame(species = not_measured,
                                mu = 70.7, SD = NA) )

#get dataframe of species with sporulation and rank
x <- 1:length(X)
names(x) <- (X)
spores2 <- spores2 %>% 
  mutate(SE = SD / sqrt(32), #SE=SD/sqrt(n), n=32 leaves
         rank = recode(species, !!!x),
         unmeasured = ifelse(species %in% not_measured, T, F))

#plot
P2 <- ggplot(spores2, aes(mu, fct_relevel(species, rev(X)) )) +
  geom_col(aes(fill = unmeasured)) +
  geom_linerange(aes(xmin = mu - SE, xmax = mu + SE)) +
  scale_fill_manual(values = c(pal[4], 'grey80'), labels = c('measured', 'unmeasured')) +
  labs( y = '', x = "Competence\n(Spores per cm^2)") +
  scale_x_continuous(breaks = c(0, 300, 600, 900), limits = c(0,900)) +
  theme(legend.position = 'none')




#Contribuion to community competence----
#p_i = totalBA_i*sporulation_i for species i
Ind <- read_csv('data/tree_level_data_all.csv')

Ind2 <- Ind %>% 
  select(BSPlotNumber, Species, basal_area) %>%
  group_by(Species, BSPlotNumber) %>% 
  summarise(BA = sum(basal_area)) %>% 
  left_join(plotdf %>% select(BSPlotNumber, ForestAllianceType)) %>% left_join(spores2 %>% select(Species = species, spores = mu)) %>% 
  mutate(p = BA*spores) %>% 
  group_by(Species, ForestAllianceType) %>% 
  summarise(mu = mean(p), 
            SE = sd(p) / sqrt(length(p))) %>% 
  mutate(rank = recode(Species, !!!x),
         unmeasured = ifelse(Species %in% not_measured, T, F))
  
#plot it
P3 <- ggplot(Ind2, aes(mu, fct_relevel(Species, rev(X)), group = ForestAllianceType )) +
  geom_col(aes(fill = as.factor(ForestAllianceType)),width = .7, position = position_dodge(width = .8)) +
  geom_linerange(aes(xmin = mu - SE, xmax = mu + SE), position = position_dodge(width = .8), lwd = .25) +
  scale_fill_manual(values = c(pal[3], pal[1])) +
  labs(x = 'Community competence \ncontribution x 10^1', y = '') +
  #xlab(expression(italic(p[i])))+
  scale_x_continuous(breaks = c(0, 150, 300, 450), labels = c(0, 15, 30, 45))+
  theme(legend.position = 'none')





#COMMUNITY COMPETENCY----
CC <- read_csv('outputs/CC_postpredictions.csv')

#Community competency: units = CC x 10^-2
CCsmaller <- CC %>% 
  mutate(ForestAllianceType = as.factor(ForestAllianceType + 1)) %>% 
  mutate_at(vars(CCmedian:upper), function(x) x/100)
plotdf$CC <- (plotdf$ccompsqrt^2)/100

#filter out predictions above 6 for redwood forests
CCsmaller <- CCsmaller %>% 
  filter(!(Richness > 6 & ForestAllianceType == 2))
plot_CC_Rich <- ggplot(CCsmaller, aes(Richness, CCmedian, group = ForestAllianceType)) +
  geom_line(aes(color = ForestAllianceType)) +
  geom_jitter(data = plotdf, aes(Richness, CC, color = ForestAllianceType), alpha = .5,width = .15, height = 0 ) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = ForestAllianceType), alpha = .2) +#, fill = '#2171b5'
  geom_line(aes(color = ForestAllianceType)) + #'#2171b5'
  labs(x = 'Richness') +
  ylab(expression(paste("Community \ncompetence"%*% 10^2))) +
  theme(axis.title.y = element_text(vjust=-1, size = 8.5),
        axis.title.x = element_text(size = 8.5), 
        legend.position = 'none') +
  scale_fill_manual(values = c(pal[3], pal[1])) +
  scale_color_manual(values = c(pal[3], pal[1])) 
plot_CC_Rich








#create legend plot-----
mockdf <- data.frame(Competency = c('Measured', 'Unmeasured'), Forest=c('Mixed evergreen', 'Redwood'), value = 1)
colorsC <- c(pal[4], 'grey80')
colorsF <- c(pal[3], pal[1])
legC <- ggplot(mockdf, aes(value, Competency, fill = Competency)) +
  geom_col() +
  scale_fill_manual(values = colorsC) +
  theme(legend.justification = 0, legend.text = element_text(size = 7), legend.title = element_text(size = 8), legend.key.size = unit(3, "mm"))
legF <- ggplot(mockdf, aes(value, Forest, fill = Forest)) +
  geom_col() +
  scale_fill_manual(values = colorsF) +
  labs(fill = 'Forest type')+
  theme(legend.justification = 0, legend.text = element_text(size = 7), legend.title = element_text(size = 8), legend.key.size = unit(3, "mm"))


#PUT PLOTS TOGETHER----

#with CC vs. richness AND shannon

grid1 <- plot_grid(
  P1 + 
    theme(axis.title.x = element_text(size = 8.5), axis.text.y = element_text(size = 7.5, face = 'italic')) +
    # scale_y_continuous(expand = c(.005,.005), breaks = 1:26, labels = rev(unique(PAdf$species))),
    scale_y_continuous(expand = c(.005,.005), breaks = 1:26, labels = (latin)),
  P2 + theme(axis.text.y = element_blank(), axis.title.x = element_text(size = 8.5)), 
  P3 + theme(axis.text.y = element_blank(), axis.title.x = element_text(size = 8.5)),
  nrow = 1, rel_widths = c(1, .5, .5),
  labels = LETTERS[1:3], vjust = .3, hjust = c(-2, -.5, -.5))
#grid1


#without the CC vs. shannon plot
blankplot <- ggplot(NULL) + theme_void()
legends <- plot_grid(get_legend(legF), get_legend(legC), ncol = 1)
legends2 <- plot_grid(blankplot, legends, nrow = 1, rel_widths = c(.4, 1))
grid3 <- plot_grid(plot_CC_Rich, blankplot, legends2, blankplot,
                   rel_heights = c(1, .2, .4, .02),
                   scale = .95,
                   ncol = 1, 
                   labels = c("D", '', '', ''), vjust = .3, hjust = -1)

finalgrid2 <- plot_grid(grid1, grid3, 
                        nrow = 1, rel_widths = c(2.9,1),
                        scale = 1) +
  theme(plot.margin = unit(c(.5, .1, .5, .1), "cm")) #T, R, B, L
#finalgrid2

ggsave(finalgrid2, device = 'pdf', width = 173, height = 110, dpi = 600, units = 'mm', filename = 'figures/nestedness_competence.pdf')

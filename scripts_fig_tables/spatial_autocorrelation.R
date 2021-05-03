#1. morans I correlogram on the mean residuals of the best performing plot-level model (Fig S2)
#2. visualizing estimaetd covariance bewteen plots using the GP bernoulli models (Fig A1, A2)

library(raster)
library(gstat)
library(ncf)
library(tidyverse)
theme_set(theme_classic())

#functions----------------------------------------

#pearsons resid = Pobs - Ppred / sqrt(Ppred(1-Ppred))
fresid <- function(I, n, post){
  obs = I / n
  pred = colMeans(post$y_rep) / n
  (obs - pred) / sqrt(pred*(1-pred))
} 



#get the data-------------------------------------
plot_level <- read_csv('data/plot_level_data.csv')

#posteriors from prevalence models 
fit <- read_rds('../../../Box/Stan_model_outputs/Big_Sur/post_plotlevel_all.RDS')
post <- rethinking::extract.samples(fit[[2]])


#calculate residuals------------------------------
resid <- fresid(plot_level$I_all, plot_level$n_all, post)
#resid <- fdev(plot_level$I_all, plot_level$n_all, post)
plot_level$resid <- resid


#check out spatial correlogram.-------------------
fit1 <- correlog(x = plot_level$Easting/1000, y = plot_level$Northing/1000, z = resid, increment = .2, 
                 resamp = 500, quiet = F)

tmp <- data.frame(x = fit1$mean.of.class, y = fit1$correlation, p = fit1$p)
p1 <- ggplot(tmp %>% filter(x<5), aes(x, y)) +
  geom_point(aes(pch = p<.05), size = 3) +
  geom_line() +
  scale_shape_manual(values = c(1, 16)) +
  labs(x = 'Distance (km)', y = 'Correlation') +
  theme(legend.position = 'none')
p1 #no positive spatial autocorrelaiton, i.e. no clustering.


pdf('figures/correlogram.pdf', width = 4.5, height = 3.5)
p1 #units are km
dev.off()


#look at distribution of pairwise distances---------
km <- as.vector(dist(cbind(plot_level$Easting, plot_level$Northing)))/1000

data.frame(d = km) %>% 
  ggplot(., aes(d)) +
  geom_histogram(aes(fill = d > 1.25),  
                 binwidth = .5, color = NA) +
  scale_fill_manual(values = c('firebrick4', 'grey80'),
                    name = 'Distance',
                    labels = c(expression(""<=" 1.25 km"), '> 1.25 km')) +
  labs(x = "Distance (km)", y = 'Frequency', title = "Pairwise distances between plots") 


ggsave('figures/histogram_of_distances.pdf', width = 5, height = 3)

sum(km < 1.25) / length(km) # 2% of plots are within 1.25 km



#plot covariance vs distance-------
post_indlevel <- read_rds('../../../Box/Stan_model_outputs/Big_Sur/post_indlevel_GP.RDS')
dmat <- read_rds('data/dmat.RDS')

#plot posterior of covariance
postGP <- rethinking::extract.samples(post_indlevel)

# compute posterior median covariance among plots
Kmedian <- matrix(0,nrow=151,ncol=151)
for ( i in 1:151){
  for ( j in 1:151){
    Kmedian[i,j] <- median(postGP$eta_sq) *
      exp( -.5/median(postGP$rho_sq) * dmat[i,j]^2 )
  }
}
diag(Kmedian) <- median(postGP$eta_sq) + median(postGP$sig_sq)
#convert to correlation matrix
Kmedian <- cov2cor(Kmedian)
cormedian <- data.frame(cor = as.vector(Kmedian), d = as.vector(dmat))

#repeat, but for individual draws
ndraws <- 100
K <- array(dim = c(ndraws, 151, 151)) 
for(z in 1:ndraws){
  for ( i in 1:151){
    for ( j in 1:151){
      K[z,i,j] <- postGP$eta_sq[z] *
        exp( -.5/postGP$rho_sq[z] * dmat[i,j]^2 )
    }
  }
}
for(z in 1:ndraws){
  diag(K[z,,]) <- postGP$eta_sq[z] + postGP$sig_sq[z]
}

#convert to correlation matrix
Kcor <- array(dim = c(ndraws, 151, 151)) 
for(z in 1:ndraws){
  Kcor[z,,] <- cov2cor(K[z,,])
}

tmp <- lapply(1:ndraws, function(z) 
  cbind(draw = z, cor = as.vector(Kcor[z,,]), d = as.vector(dmat)))
tmp2 <- as.data.frame(do.call(rbind, tmp))

tmp2 %>% filter(d < 5, d> 0) %>% 
  ggplot(., aes(d, cor)) +
  geom_line(lwd = .2, aes(group = draw), alpha = .2) +
  geom_line(data = cormedian %>% filter(d<5, d>0), 
            aes(d, cor), color = 'firebrick4') +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = 'Distance (km)', y = 'Correlation', title = "Gaussian process posterior")

ggsave('figures/gaussian_process_correlation.pdf', width = 5, height = 3)


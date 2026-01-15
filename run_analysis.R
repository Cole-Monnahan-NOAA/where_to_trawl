
## Main script to run the analysis

library(tidyr)
library(dplyr)
library(ggplot2)
theme_set(theme_bw())
library(fmesher)
library(fields)




# loads the spatial simulator function
source("code/simulate_correlated_fields.R")


## ---------------- POLLOCK SCENARIO ---------------------------
# category = unique combination of class, spp, size, etc. so a
# single row
categories <- expand.grid(scenario='EBS-pollock', 
                          class='PLK1', # acoustic class
                          size=c('15cm', '25cm', '35cm', '>35cm'),
                          species=c('pollock', 'other-fish'))
categories <- mutate(categories, category=paste(species,class,size, sep='_'))
# assumed population level mean of log abundance (I guess billions for pollock?)
categories$logmean <- log(c(10,15,5,3, 3, 1, .5, .1))
Nc <- nrow(categories)

# spatial properties
range <- 5 # decorrelation range = how many km away to get 10% correlation in space
spatial.var <- .1 # spatial variation controls range of simulated log-abundance
mu <- categories$logmean
var <- rep(1,Nc) # variance of category, assumed 1 for now
# build the covariance of our categories
Sigma <- diag(Nc)
dimnames(Sigma) <- list(categories$category, categories$category)
# correlations for pollock
for(i in 2:4) for(j in 1:(i-1)) Sigma[i,j] <- .7
# correlations for other-fish
for(i in 6:8) for(j in 5:(i-1)) Sigma[i,j] <- .3
# make symmetric
Sigma[upper.tri(Sigma)] <- t(Sigma)[upper.tri(Sigma)]
solve(Sigma) # ensure it's positive definite (a valid correlation)
# check it
#corrplot::corrplot(Sigma, is.corr=TRUE)
# make it a covariance by scaling by the variances
Sigma <- Sigma * tcrossprod(sqrt(var))
#corrplot::corrplot(Sigma, is.corr=FALSE)

# interval locations, for now one transect
xy <- expand.grid(x=1:300, y=1) # 300 1km intervals
# simulate spatial fields of log abundance
sim <- simulate_fields_gp(xy=xy, mu=mu, Sigma=Sigma,
                          range=range, spatial.var=spatial.var,
                          nreps=3)
# add on the category information
sim <- merge(categories, sim, by='category')

# quick visual checks, !! only works for a single transect!!
ggplot(sim, aes(x, y=logA, color=size)) + geom_line() + facet_grid(species~replicate)

write.csv(sim, file='results/sim_test.csv')

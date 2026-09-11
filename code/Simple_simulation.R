
## Main script to run the analysis
# based on run_analysis - modfied by ADR

library(tidyr)
library(dplyr)
library(ggplot2)
theme_set(theme_bw())
library(fmesher)
library(fields)

# loads the spatial simulator function
source("code/simulate_correlated_fields.R")
source("code/simulate_AT_survey.R")

# 1. Setup paramters for pollock survey ----
# Cole's parameters
#
# category = unique combination of class, spp, size, etc. so a
# single row
categories <- expand.grid(scenario='EBS-pollock', 
                          class='PK1', # acoustic class
                          size=c('15_cm', '25_cm', '35_cm', '50_cm'),
                          species=c('pollock', 'other-fish'))
categories <- mutate(categories, category=paste(species,class,size, sep='_'))
# assumed population level mean of log abundance (I guess billions for pollock?)
categories$logmean <- log(c(10,15,5,3, 3, 1, .5, .1)) # log mean abundance
categories$FL=c(15,25,35,50, 15,25,35,50) # Fork length for each category
categories$sigma_bs=10^((20*log10(categories$FL)-66)/10) # add sigma BS for a single fish - for now assume same TS relationship

#
# spatial properties
range <- 10 # decorrelation range = how many km away to get 10% correlation in space # initial was 50
spatial.var <- .2 # spatial variation controls range of simulated log-abundance

#
# Alex's paramters
# Let's just look at the first replicate 
param=list()
param$n_edsu = 500 # number of EDSU's to simulate
param$n_hauls=2  # number of hauls to sample
param$replicate =1 # which replicate to process
param$nreps=1 # number of replicates to simulate

 



# 2. setup for  simulation ----
Nc <- nrow(categories)
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
xy <- expand.grid(x=1:param$n_edsu, y=1) # n_edsu 1km intervals
xy$interval=seq_len(nrow(xy)) # add a interval number for convenience


#3. simulate spatial fields of log abundance ----
sim <- simulate_fields_gp(xy=xy, mu=mu, Sigma=Sigma,
                          range=range, spatial.var=spatial.var,
                          nreps=param$nreps)
sim <- merge(categories, sim, by='category')# # add on the category information
sim <- merge(xy, sim, by=c('x','y')) # add the interval data

# now add quantities needed for simulation to/from acoustic units
sim<-mutate(sim,A=exp(logA)) # linear abundance
sim<-mutate(sim, sa=sigma_bs*A) # sa or total backscatter in unit of m^2

# compute fraction of sa for each combination of interval and replicate
sim <- sim %>%
  group_by(interval, replicate) %>%
  mutate(prop_sa = sa / sum(sa, na.rm = TRUE)) %>%
  ungroup()


#4. now compute survey abundance ----

result=tibble() # start empty tibble for results

replicate_on=1 # sim can produce many replicates at one time.  not suing for now

n_hauls=5 # number of hauls
thresh_dist=round(length(sim$sa)/(n_hauls*4)) # min # EDSU's apart

result_rand<-simulate_survey(sim, replicate_on, sampling_method='random',n_hauls, thresh_dist)
result_max<-simulate_survey(sim, replicate_on, sampling_method='max',n_hauls, thresh_dist)
result_systematic<-simulate_survey(sim, replicate_on, sampling_method='systematic',n_hauls, thresh_dist)
result_cumsum<-simulate_survey(sim, replicate_on, sampling_method='cumsum',n_hauls, thresh_dist)

# combine results for all methods.
result<-bind_rows(result_rand,result_max,result_systematic,result_cumsum)


# 5. Report results --- 


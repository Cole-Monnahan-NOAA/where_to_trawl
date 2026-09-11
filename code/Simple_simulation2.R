
## Main script to run the analysis
# based on run_analysis - modfied by ADR

rm(list = ls()) # clear all workspace variables

library(tidyverse)
library(here)
theme_set(theme_bw())
library(fmesher)
library(fields)

# loads the spatial simulator function
source("code/simulate_correlated_fields.R") # Cole's base function
source("code/simulate_abundance.R") # set up correlated fields
source("code/simulate_AT_survey.R") # survey calculations

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



# Simulation paramters
param=list()
param$n_edsu = 500 # number of EDSU's to simulate
param$n_hauls=100  # number of hauls to sample
param$thresh_dist=round(param$n_edsu/(param$n_hauls*4)) # min # EDSU's apart for hauls
param$replicate_sims =1 # which replicate to process within sim [typically 1]
param$n_iter=50 # number of replicates to simulate
param$save_name = "results/delme" # path & name to save output
# spatial properties
param$range <- 25 # decorrelation range = how many km away to get 10% correlation in space # initial was 50
param$spatial_var <- .4 # spatial variation controls range of simulated log-abundance
param$categories<-categories



# Pre-allocate a list to store results across iterations
results_list <- vector("list", param$n_iter)

for (i in seq_len(param$n_iter)) {
  message("Starting iteration ", i, " of ", param$n_iter)
  # 1. Setup simulation
  sim <- simulate_abundance(categories=param$categories, n_edsu = param$n_edsu, range=param$range, 
    spatial.var=param$spatial_var, nreps = param$replicate_sims)
  
  # make a plot for first time
  if (i==1) {
    ggplot(sim, aes(x = interval, y = sa,color=category)) +
      geom_line() +
      labs(x = "interval", y = "sa") +
      theme_minimal()
    ggsave(here(paste0(param$save_name,'-fig.png')), width = 10, height = 6, dpi = 150)
  }
  
  
  # 2. Compute survey abundance across trawling strategies
  result_rand <- simulate_survey(sim, replicate_on, sampling_method = 'random', param$n_hauls, param$thresh_dist)
  result_max <- simulate_survey(sim, replicate_on, sampling_method = 'max', param$n_hauls, param$thresh_dist)
  result_systematic <- simulate_survey(sim, replicate_on, sampling_method = 'systematic', param$n_hauls, param$thresh_dist)
  result_cumsum <- simulate_survey(sim, replicate_on, sampling_method = 'cumsum', param$n_hauls, param$thresh_dist)
  
  # 3. Combine methods for iteration i and tag with iteration index
  results_list[[i]] <- bind_rows(
    result_rand,
    result_max,
    result_systematic,
    result_cumsum
  ) %>%
    mutate(iter = i)
}

# 4. Stack all iteration outputs into a single dataframe
result <- bind_rows(results_list)

# clean up a bit
rm(result_cumsum, result_max,result_rand,result_systematic,results_list)
# 5. Report results --- 
write_csv(result, here(paste0(param$save_name,'-results.csv'))) # results
capture.output(print(param), file = here(paste0(param$save_name,'-params.txt'))) # parameters

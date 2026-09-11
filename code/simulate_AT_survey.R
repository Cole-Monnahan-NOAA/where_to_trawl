# Functions to simulate survey

# survey calculations
# simulate_AT_survey 
#
#Trawling strategy 
# get_trawl_locs_max - highest values, but not too close 
# get_trawl_locs_random - random, but not too close
# get_trawl_locs_cum - cumulative sum
# get_trawl_locs_systematic - systematic sampling


# simulate_AT_survey.R
#
#' simulate survey - compute survey abundance given simulated abundance field
#'
#' @param sim simulated data field
#' @param replicate replicate to process (should be 1 unless expanded)
#' @param sampling_method  method used for trawl sampling 
#'       options: "max", "random", "cumsum", "systematic"
#' @param n_hauls  - number of hauls
#' @param thresh_dist - min # EDSUs between trawl samples
#' 
#'  @return result tibble comparing realized and actual results over th survey area
#' class, category, replicate A_survey, A_true, percent error
#'
#'sample call 
#' tmp<-simulate_survey(sim, replicate_on, sampling_method='max',n_hauls=5, thresh_dist=50)

simulate_survey <- function(sim, replicate_on, sampling_method, n_hauls, thresh_dist) {
  
  replicate_on=1 # only process the first replicate
  
  # now, colllapse to the EDSU-level acoustic observations
  survey<- sim %>%
    group_by(x, y, class,replicate) %>%
    summarise(sa = sum(sa, na.rm = TRUE), .groups = c("drop"))
  # quick diagnostic plot assuming that x defines unique points
 # ggplot(EDSUs, aes(x=x, y=sa, color=class)) + geom_line() + facet_grid(~replicate)
  
  
  # now for a single replicate
  survey <- survey%>% filter(replicate == replicate_on )  # Echosounder observations
  
  # haul locs - get indices into trawl locations
  # "max", "random", "cumsum", "systematic"
if (sampling_method == 'max') {
  # trawl at high densities
  ind <- get_trawl_locs_max(survey$sa, n_hauls, thresh_dist)
} else if (sampling_method == 'random') {
  ind <- get_trawl_locs_random(survey$sa, n_hauls, thresh_dist)
} else if (sampling_method == 'cumsum') {
  ind <- get_trawl_locs_cum(survey$sa, n_hauls)
} else if (sampling_method == 'systematic') {
  ind <- get_trawl_locs_systematic(survey$sa, n_hauls)
} else {
  print("Undefined sampling method")
  return()
}
  
  # define the nearest haul
  hauls <- sim %>%
    slice(ind$index) %>%
    select(x, y, interval)
  
  hauls$haul_id=seq(1,nrow(hauls)) # assign haul_ID
  
  # add nearest haul indices to EDSU_obs
  survey <- survey %>%
    mutate(haul_id = apply(
      outer(x, hauls$x, `-`)^2 + outer(y, hauls$y, `-`)^2, 
      1, 
      which.min
    ),
    haul_x = hauls$x[haul_id],
    haul_y = hauls$y[haul_id]
    )
  
  survey <- survey %>%
    rename(sa_total=sa) 
  
  # get information from trawl locations from sim
  catch_obs <- sim %>% 
    filter(replicate == param$replicate) %>%
    inner_join(hauls, by = c("x", "y")) %>%
    rename(interval=interval.x, prop_sa_haul=prop_sa)%>%
    select(haul_id, interval, category, class, FL, sigma_bs, prop_sa_haul)
  
  # now, compute biomass from acoustic observations
  survey<- survey %>%
    left_join(catch_obs, by = "haul_id", relationship = "many-to-many")%>%
    rename(class=class.x)%>%
    mutate(sa=prop_sa_haul*sa_total) %>%
    mutate(A=sa/sigma_bs)%>%  
    select(-class.y)
  
  survey_result<-survey %>% group_by(class,category,replicate) %>%
    summarise(A=sum(A, na.rm = TRUE), .groups="drop")
  
  true_result<-sim %>% group_by(class,category,replicate) %>%
    summarise(A=sum(A, na.rm = TRUE), .groups="drop")
  
  result<- survey_result %>%
    left_join(true_result, by = c("replicate","category"), relationship = "many-to-many")%>%
    rename(A_survey=A.x, A_true=A.y, class=class.x)%>%
    select(-class.y)%>%
    mutate(percent_error=((A_survey-A_true)/A_true)*100, method=sampling_method)
  
  rm (true_result, survey_result) # clean up  
  return(result)
}

############################## Functions to get Trawl locations ########################

# get_trawl_locs_max
# 
# identify trawl locations based on maximum observed values further than n_thresh indicies apart
#
# inputs
# Get trawl locations from sa vector
# n_hauls - number of hauls
# thresh dist - maxima must be at least thresh_dist indicies  apart (i.e. if 50 then no maxima within 50)
# 
# outputs
# data frame with index [index into sa vector]
# 
# usage ind=get_trawl_locs_max(sa, n_hauls = 5, thresh_dist =  round(length(sa)/(n_hauls*4)))
get_trawl_locs_max <- function(sa, n_hauls, thresh_dist = 0) {
  len <- length(sa)
  if (len == 0) {
    return(data.frame(index = integer(0), value = numeric(0)))
  }
  
  # 1. Order all indices by value in descending order (no peak requirement)
  sorted_indices <- order(sa, decreasing = TRUE)
  
  # 2. Greedily select top values separated by > thresh_dist samples
  selected <- integer(0)
  for (idx in sorted_indices) {
    if (length(selected) == 0 || all(abs(idx - selected) > thresh_dist)) {
      selected <- c(selected, idx)
    }
    if (length(selected) == n_hauls) break
  }
  
  #3. Return data frame sorted chronologically by index
  selected <- sort(selected)
  return(data.frame(index = selected, value = sa[selected]))
}

# get_trawl_locs_random
# 
# identify trawl locations based on random locations further than n_thresh indicies apart
#
# inputs
# Get trawl locations from sa vector
# n_hauls - number of hauls
# thresh dist - maxima must be at least thresh_dist indicies  apart (i.e. if 50 then no maxima within 50)
# 
# outputs
# data frame with index [index into sa vector]
# useage ind=get_trawl_locs_random(sa, n_hauls = 5, thresh_dist =  round(length(sa)/(n_hauls*4)))
get_trawl_locs_random <- function(sa, n_hauls, thresh_dist = 0) {
  len <- length(sa)
  
  if (len == 0) {
    return(data.frame(index = integer(0), value = numeric(0)))
  }
  
  # 1. Randomize all indices
  random_indices <- sample(seq_along(sa))
  
  # 2. Greedily select valid indices separated by > thresh_dist samples
  selected <- integer(0)
  for (idx in random_indices) {
    if (length(selected) == 0 || all(abs(idx - selected) > thresh_dist)) {
      selected <- c(selected, idx)
    }
    if (length(selected) == n_hauls) break
  }
  
  # Return data frame sorted chronologically by index
  selected <- sort(selected)
  return(data.frame(index = selected, value = sa[selected]))
}


# get_trawl_locs_cum
# 
# identify trawl locations based on cumulative distance apart
# e.g. for 2 samples would be .25, .75 of cumulative sum
#
# inputs
# Get trawl locations from sa vector
# n_hauls - number of hauls
# 
# outputs
# data frame with index [index into sa vector]
# useage ind=get_trawl_locs_cum(sa, n_hauls)
## pick the cumulative values (i.e. if n_hauls 5, pick 10, 30, 50, 70, 90)
get_trawl_locs_cum <- function(sa, n_hauls = 5) {
  len <- length(sa)
  
  if (len == 0) {
    return(data.frame(index = integer(0), value = numeric(0)))
  }
  
  #  cumulative sum  
  cum_sa <- cumsum(sa)/sum(sa)
  # Centered target fractions:  
  targets <- (seq_len(n_hauls) - 0.5) / n_hauls
  # Pick index closest to each target cumulative value
  selected <- sapply(targets, function(t) {
    which.min(abs(cum_sa - t))
  })
  selected <- sort(unique(selected))
  return(data.frame(index = selected, value = sa[selected]))
}



# get_trawl_locs_systematic
# 
# identify trawl locations based on random locations further than n_thresh indicies apart
#
# inputs
# Get trawl locations from sa vector
# n_hauls - number of hauls
# thresh dist - maxima must be at least thresh_dist indicies  apart (i.e. if 50 then no maxima within 50)
# 
# outputs
# data frame with index [index into sa vector]
# useage ind=get_trawl_locs_systematic(sa, n_hauls = 5)
get_trawl_locs_systematic <- function(sa, n_hauls) {
  len <- length(sa)
  
  if (len == 0) {
    return(data.frame(index = integer(0), value = numeric(0)))
  }
  
  # Calculate evenly centered indices across vector length
  target_indices <- round(((seq_len(n_hauls) - 0.5) / n_hauls) * len)
  
  # Clamp indices within range [1, len] and keep unique values
  selected <- pmax(1, pmin(len, target_indices))
  selected <- sort(unique(selected))
  
  return(data.frame(index = selected, value = sa[selected]))
}

#'simulate_abundance
#' 
#' #' Simulates spatial pollock abundance fields across EDSU intervals  
#' uses Cole's simulate_fields_gp 
#' requires simulate_correlated_fields.R function
#' 
#' @param categories Data frame containing category metadata  
#' @param n_edsu Integer specifying the number of spatial intervals.
#' @param range Spatial correlation range parameter.
#' @param spatial.var Spatial variance parameter.
#' @param nreps Integer number of simulation repetitions (generally 1)
#'
#' @return sim A merged frame containing simulated spatial abundances for each edsu/category/size_class
#' 
#' @examples
#'sim<-simulate_abundance(categories, n_edsu, range, spatial.var, nreps)
simulate_abundance <- function(categories, n_edsu, range, spatial.var, nreps) {
  # build up the parameters
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

  
  #3. simulate spatial fields of log abundance using simpulate_fields_gp ----
  sim <- simulate_fields_gp(xy=xy, mu=mu, Sigma=Sigma,
                            range=range, spatial.var=spatial.var,
                            nreps=nreps)
  
  #4. information needed
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
  return(sim)
}



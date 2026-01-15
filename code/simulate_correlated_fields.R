
# Very early version of this that does not leverage sparsity


#' Gaussian process simulation of two correlated 2d spatial
#' fields using a Matern covariance function
#'
#' @param xy A data.frame of spatial locations
#' @param mu A vector of mean log abundance for each category
#' @param Sigma A covariance matrix of the categories
#' @param range The spatial decorrelation range
#' @param spatial.var The spatial variation
#' @param nreps Number of replicate simulations to do
#' @param seed An optional random number seed for reproducibility
#' @return A data.frame with locations and simulated
#'   log-abundance for each category for each replicate
simulate_fields_gp <- function(xy, mu,  Sigma, range, 
                               spatial.var, nreps, seed=NULL) {
  stopifnot(all.equal(length(mu), ncol(Sigma), nrow(Sigma)))
  # build spatial covariance matrix from Matern function 
  dists <- as.matrix(dist(xy))
  cov_s <- spatial.var*Matern(dists, range=range, smoothness=1)
  # combine space and categories
  cov_total <- kronecker(Sigma, cov_s)
  
  # Simulate the data using cholesky of the total covariance
  L <- t(chol(cov_total)) # This returns the lower triangular matrix
  # uncorrelated white noise: N(0,1)
  set.seed(seed)
  out <- list()
  for(rr in 1:nreps){
    epsilon <- rnorm(nrow(cov_total))
    # Correlate values (assuming mean = 0) based on cov_total
    z <- L %*% epsilon
    z <- matrix(z, ncol=ncol(Sigma))
    # add means by category
    z <- as.data.frame(z + mu[col(z)])
    names(z) <- dimnames(Sigma)[[1]]#paste0('z', 1:ncol(Sigma))
    out[[rr]] <- cbind(replicate=rr, x=xy$x, y=xy$y, z)
  }
  out <- bind_rows(out) |>
    pivot_longer(cols=-c(replicate,x,y), names_to='category', values_to = 'logA')
  return(out)
}



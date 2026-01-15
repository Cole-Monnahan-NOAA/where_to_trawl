
# Very early version of this that does not leverage sparsity
library(tidyr)
library(ggplot2)
library(fmesher)
library(fields)

# Gaussian process simulation of two correlated 2d spatial fields
# using a Matern covariance function
simulate_fields_gp <- function(xy, mean=c(0,0), alpha = 1, 
                                sigma2 = 1, rho=0) {
   dist_matrix <- as.matrix(dist(xy))
   # alpha: range parameter (inverse of rho)
   # sigma2: partial sill (variance)
   # nu: smoothness parameter
   nu <- 1
   cov_s <- sigma2 * fields::Matern(dist_matrix, range = alpha, smoothness = nu)
   cov_v <- matrix(c(1,rho,rho,1), nrow=2)
   cov_matrix <- kronecker(cov_v, cov_s)
   #corrplot::corrplot(cov2cor(cov_matrix))
  
   # Simulate the data using Cholesky cholesky of the total covariance
   # We use the formula: y = Mu + L %*% epsilon
   # where L is the lower triangular Cholesky factor (L L^T = Cov)
  n <- nrow(cov_matrix)
  L <- chol(cov_matrix) # This returns the upper triangular matrix
  epsilon <- rnorm(n)
  
  # Generate the spatially correlated values (assuming mean = 0)
  z <- t(L) %*% epsilon
  z <- matrix(z, ncol=2) 
  out <- data.frame(x=xy$x, y=xy$y, z1=z[,1]+mean[1], z2=z[,2]+mean[2])
  return(out)
}

xy <- expand.grid(x=1:100, y=seq(1:5)*2)
sim <- simulate_fields_gp(xy=xy, mean=c(8,0), sigma2=5, alpha = 1, rho=.9)
sim.long <- tidyr::pivot_longer(sim, cols=c('z1', 'z2'))
ggplot(sim.long, aes(x, y=value, color=name)) + geom_point() + geom_line() +
  facet_wrap('y') # y is a proxy for transect number for now
ggplot(sim.long, aes(x, y=y, size=value, color=value)) + geom_point() +
  facet_wrap('name')

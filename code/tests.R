





simulate_spde_matern <- function(xy, alpha = 0.2, sigma2 = 1, nu = 1) {
  # nu = 1 in SPDE corresponds to a Matern smoothness of alpha_spde - d/2
  # For 2D (d=2), alpha_spde = nu + 1. We'll fix alpha_spde = 2 for simplicity.
  
  if (!requireNamespace("fmesher", quietly = TRUE)) {
    stop("Package 'fmesher' is required. Please install it.")
  }
  
  # 1. Create a Delaunay triangulation (mesh) over the spatial domain
  # The mesh allows us to define the sparsity pattern based on vertex neighbors
  mesh <- fmesher::fm_mesh_2d_inla(loc = xy, max.n.strict=nrow(xy), max.edge=.01)
  
  # 2. Define SPDE parameters
  # kappa relates to the range (alpha): kappa = sqrt(8*nu)/range
  kappa <- sqrt(8 * 1) / alpha
  
  # 3. Construct the sparse Finite Element matrices
  # C is the mass matrix, G is the stiffness matrix
  fem <- fmesher::fm_fem(mesh)
  C <- fem$c0
  G <- fem$g1
  
  # 4. Build the Precision Matrix (Q)
  # For alpha_spde = 2: Q = (kappa^2 * C + G)
  # This matrix is sparse!
  Q <- (kappa^4 * C + 2 * kappa^2 * G + fem$g2) * (1 / (4 * pi * kappa^2 * sigma2))
  Q1 <- fem$c0*(1+kappa^2)^2 +
       fem$g1*(1+kappa^2)*kappa +
       fem$g2 *kappa^2
  # 5. Simulate the field on the mesh vertices
  # Since Q*z = e (where e ~ N(0, Q)), we solve L'z = epsilon
  # But a simpler way for sparse Q is using Matrix::chol
  L <- Matrix::Cholesky(Q, Imult = 0.0001) # Precision factor
  z_mesh <- Matrix::solve(L, rnorm(nrow(Q)), system = "Lt")
  ggplot(data.frame(x=mesh$loc[,1], y=mesh$loc[,2],z=z_mesh), 
         aes(x,y, size=z_mesh, color=z_mesh)) + geom_point()
  # 6. Project the values from mesh vertices back to your original 'xy' locations
  A <- fmesher::fm_basis(mesh, loc = as.matrix(xy))
  z_out <- as.vector(A %*% z_mesh)
  
  return(z_out)
}

xy <- expand.grid(x=1:100, y=seq(1:1)*10)
xy <- xy[sample(1:nrow(xy), size=10, replace=FALSE),]
xy <- data.frame(x=runif(1000), y=runif(1000))
set.seed(1)
sim <- cbind(xy, d=simulate_spde_matern(xy=xy,  sigma2=.1, alpha = 100))
ggplot(sim, aes(x,y=d, color=y)) + geom_point()
ggplot(sim, aes(x,y=y, size=d, color=d)) + geom_point()



simulate_multivariate_spde <- function(xy, alpha = 0.2, sigma2 = 1, nu = 1, 
                                       rho = 0.9, v_scales = c(1, 1), mean=c(5,5)) {
  if (!requireNamespace("fmesher", quietly = TRUE)) stop("Install 'fmesher'")
  if (!requireNamespace("Matrix", quietly = TRUE)) stop("Install 'Matrix'")
  # 1. Spatial Component (SPDE)
  mesh <- fmesher::fm_mesh_2d_inla(loc = xy, max.n.strict=nrow(xy), max.edge=.01)
  #plot(mesh)
  mesh$n
  kappa <- sqrt(8 * 1) / alpha
  fem <- fmesher::fm_fem(mesh)
  
  # Sparse Spatial Precision Matrix (Qs)
  Qs <- (kappa^4 * fem$c0 + 2 * kappa^2 * fem$g1 + fem$g2) * (1 / (4 * pi * kappa^2 * sigma2))
  #Matrix::image(Qs)
  # 2. Multivariate Component (2x2)
  # Create a correlation matrix for the 2 variables
  Corr_v <- matrix(c(1, rho, rho, 1), nrow = 2)
  # Scale by marginal variances (v_scales) to get Covariance
  Sigma_v <- diag(v_scales) %*% Corr_v %*% diag(v_scales)
  # Invert to get the Variable Precision Matrix (Qv)
  Qv <- solve(Sigma_v)
  
  # 3. Kronecker Expansion (The Precision Matrix Trick)
  # Q_total = Qs (spatial) %x% Qv (variables)
  # Note: %x% works with sparse Matrix objects to keep the result sparse
 # browser()
  Q_total <- Matrix::kronecker(Qs, Qv)
 # Matrix::image(Q_total)
# corrplot::corrplot(cov2cor(solve(Q_total)))
  # 4. Simulate using Sparse Cholesky
  # We solve L'z = e where e ~ N(0, I) and L L' = Q_total
  L <- Matrix::Cholesky(Q_total)
 # Matrix::image(L)
  z_mesh_all <- Matrix::solve(L, rnorm(nrow(Q_total)), system = "Lt")
  
  # 5. Reorganize and Project back to 'xy'
  # The simulation gives a vector where [Var1_mesh, Var2_mesh, Var1_mesh...] 
  # We split it into two vectors (one for each variable)
  n_mesh <- nrow(Qs)
  z_mesh_v1 <- z_mesh_all[seq(1, by = 2, length.out = n_mesh)]
  z_mesh_v2 <- z_mesh_all[seq(2, by = 2, length.out = n_mesh)]
 plot(cbind(z_mesh_v1, z_mesh_v2))
  # Project to original xy locations
  A <- fmesher::fm_basis(mesh, loc = as.matrix(xy))
  z_out_v1 <- as.vector(A %*% z_mesh_v1) + mean[1]
  z_out_v2 <- as.vector(A %*% z_mesh_v2) + mean[2]
  return(data.frame(x = xy[,1], y = xy[,2], Var1 = z_out_v1, Var2 = z_out_v2))
}


xy <- expand.grid(x=1:100, y=seq(1:1)/2)
#xy <- data.frame(x=rnorm(5), y=rnorm(5))
set.seed(1)
sim <- simulate_multivariate_spde(xy=xy, sigma2=.1, mean=c(10,-10), alpha = 1000, rho=-.99)
sim.long <- tidyr::pivot_longer(sim, cols = c('Var1', 'Var2'))
ggplot(sim.long, aes(x, y=value, color=name)) + geom_point()
ggplot(sim, aes(Var1, Var2)) + geom_point()
ggplot(sim.long, aes(x,y=y, color=value, size=value)) + geom_point()  +
  facet_wrap('name')



# 
# simulate_multivariate_gmrf <- function(n_x, n_y, rho_spatial = 0.95, rho_var = 0.7, sigma2 = 1) {
#   require(Matrix)
#   
#   # 1. Create Spatial Structure (Q_s)
#   # For a grid, we use a 1st order neighborhood (rook's case)
#   N <- n_x * n_y
#   adj <- Matrix::sparseMatrix(
#     i = c(rep(1:n_x, each = n_y - 1), rep(1:(n_x - 1), each = n_y)),
#     j = c(rep(1:n_x, each = n_y - 1) + rep(1:(n_y - 1), n_x), 
#           rep(1:(n_x - 1), each = n_y) + n_y),
#     x = 1, dims = c(N, N), symmetric = TRUE
#   )
#   
#   # Degree matrix (number of neighbors)
#   D <- Diagonal(x = rowSums(adj))
#   
#   # Precision for a Proper CAR model: Qs = D - rho * W
#   # (Ensures positive definiteness by scaling rho)
#   Qs <- D - rho_spatial * adj
#   
#   # 2. Create Variable Structure (Q_v)
#   # Precision matrix for 2 variables with correlation rho_var
#   Sigma_v <- matrix(c(1, rho_var, rho_var, 1), 2, 2) * sigma2
#   Qv <- solve(Sigma_v)
#   
#   # 3. Combine using Kronecker
#   # This results in a sparse (2N) x (2N) matrix
#   Q_total <- kronecker(Qs, Qv)
#   
#   # 4. Simulate
#   # Standard GMRF simulation via Cholesky: Q = LL' -> Solve L'z = epsilon
#   # We use the Matrix-specific Cholesky for speed
#   chol_Q <- Matrix::Cholesky(Q_total, LDL = FALSE)
#   z_raw <- Matrix::solve(chol_Q, rnorm(2 * N), system = "Lt")
#   
#   # 5. Format Output
#   # Extract interleaved variables
#   var1 <- z_raw[seq(1, 2 * N, by = 2)]
#   var2 <- z_raw[seq(2, 2 * N, by = 2)]
#   
#   grid <- expand.grid(y = 1:n_y, x = 1:n_x)
#   return(cbind(grid, Var1 = as.vector(var1), Var2 = as.vector(var2)))
# }
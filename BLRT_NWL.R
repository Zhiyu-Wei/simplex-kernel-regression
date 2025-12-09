# -----------------------------------------------------------
# Lognormal kernel utilities for BLRT simulations
#
# Contains:
#   - lognormal_kernel_matrix()
#   - lognormal_kernel_reg_cv()
#
# All comments and error messages are in English.
# -----------------------------------------------------------

# -----------------------------------------------------------
# lognormal_kernel_matrix
#
# Input:
#   X    : n × d matrix of compositional predictors (strictly > 0)
#   h    : bandwidth
#   S    : covariance matrix of log(X) (optional; estimated if NULL)
#   ridge: small ridge added to S for numerical stability
#
# Output:
#   A list with components:
#     K    : n × n lognormal kernel matrix
#     Z    : n × d log-transformed data
#     S    : d × d covariance of Z (after ridge)
#     H    : d × d bandwidth matrix = h^2 * S
#     Hinv : inverse of H
# -----------------------------------------------------------
lognormal_kernel_matrix <- function(X, h, S = NULL, ridge = 1e-6) {
  X <- as.matrix(X)
  n <- nrow(X)
  d <- ncol(X)
  
  if (any(X <= 0)) {
    stop("All entries of X must be strictly positive for the log transform.")
  }
  
  Z <- log(X)  # n × d
  
  if (is.null(S)) {
    S <- stats::cov(Z)
  }
  # Ridge regularization for numerical stability
  S <- S + ridge * diag(d)
  
  H    <- (h^2) * S
  Hinv <- solve(H)
  
  # Quadratic form:
  # Q_ij = (z_i - z_j)^T Hinv (z_i - z_j)
  #      = d_i + d_j - 2 * z_i^T Hinv z_j
  A    <- Z %*% Hinv          # n × d
  dvec <- rowSums(A * Z)      # length n, d_i = z_i^T Hinv z_i
  
  Q <- outer(dvec, dvec, "+") - 2 * (A %*% t(Z))  # n × n
  
  # Gaussian-type kernel (normalizing constant omitted)
  K <- exp(-0.5 * Q)
  
  list(K = K, Z = Z, S = S, H = H, Hinv = Hinv)
}


# -----------------------------------------------------------
# lognormal_kernel_reg_cv
#
# Nadaraya–Watson regression with a lognormal kernel and
# least-squares cross-validation to choose h.
#
# Input:
#   X     : n × d compositional predictors (strictly > 0)
#   y     : response vector of length n
#   h_seq : optional grid of bandwidths; if NULL, an automatic
#           search is performed (coarse grid + optimize)
#   ridge : ridge term used in covariance matrix of log(X)
#
# Output:
#   An object of class "lognorm_kreg_cv" with components:
#     h_opt  : selected bandwidth
#     cv_min : minimum CV score
#     grid   : data.frame(h, cv) for the explored grid
#     cv_fun : underlying CV function (h → CV(h))
# -----------------------------------------------------------
lognormal_kernel_reg_cv <- function(X, y,
                                    h_seq = NULL,
                                    ridge = 1e-6) {
  X <- as.matrix(X)
  n <- nrow(X)
  d <- ncol(X)
  
  if (length(y) != n) {
    stop("length(y) must be equal to nrow(X).")
  }
  if (any(X <= 0)) {
    stop("All entries of X must be strictly positive for the log transform.")
  }
  
  Z <- log(X)
  S <- stats::cov(Z) + ridge * diag(d)
  
  # Given h, build K(h) and compute LOOCV MSE
  cv_fun <- function(h) {
    kmat <- lognormal_kernel_matrix(X, h, S = S, ridge = ridge)
    K    <- kmat$K        # n × n
    
    num <- K %*% y        # n × 1, includes self
    den <- rowSums(K)     # length n
    
    diagK <- diag(K)
    
    num_loo <- as.numeric(num) - diagK * y
    den_loo <- den - diagK
    
    # Guard against zero or negative denominators
    bad   <- den_loo <= .Machine$double.eps
    y_hat <- num_loo / den_loo
    if (any(bad)) {
      y_hat[bad] <- mean(y)   # fallback: global mean
    }
    
    mean((y - y_hat)^2)
  }
  
  if (is.null(h_seq)) {
    # Rule-of-thumb initial scale
    h0    <- n^(-1 / (d + 4))
    h_min <- h0 / 5
    h_max <- h0 * 5
    
    # Coarse grid search
    grid_h  <- exp(seq(log(h_min), log(h_max), length.out = 15))
    grid_cv <- sapply(grid_h, cv_fun)
    
    h_init <- grid_h[which.min(grid_cv)]
    
    # Refine around the best grid value using optimize()
    lower <- h_init / 2
    upper <- h_init * 2
    
    opt <- stats::optimize(cv_fun, interval = c(lower, upper))
    
    h_opt  <- opt$minimum
    cv_min <- opt$objective
    
    res <- list(
      h_opt  = h_opt,
      cv_min = cv_min,
      grid   = data.frame(h = grid_h, cv = grid_cv),
      cv_fun = cv_fun
    )
  } else {
    grid_cv <- sapply(h_seq, cv_fun)
    k       <- which.min(grid_cv)
    h_opt   <- h_seq[k]
    
    res <- list(
      h_opt  = h_opt,
      cv_min = grid_cv[k],
      grid   = data.frame(h = h_seq, cv = grid_cv),
      cv_fun = cv_fun
    )
  }
  
  class(res) <- "lognorm_kreg_cv"
  return(res)
}

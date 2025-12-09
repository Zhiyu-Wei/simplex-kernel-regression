# -------------------------------------------------------
# Kernel Regression Utilities for Simplex Predictors
#
# Contains:
#   - d.kernel_fast
#   - make_simplex_grid
#   - build_K_dirichlet_fast
#   - build_K_dirichlet_fast_UX
#   - build_S_locallinear_dirichlet
#   - build_S_locallinear_dirichlet_UX
#   - lognormal_kernel_matrix
#   - lognormal_kernel_matrix_new
#
# Only functions used in the simulation are included.
# All comments and error messages are in English.
# -------------------------------------------------------

# -------------------------------------------------------
# Fast Dirichlet kernel matrix K(U, X)
# U: nU × p
# X: nX × p
# -------------------------------------------------------
d.kernel_fast <- function(U, X, alpha, h) {
  U <- as.matrix(U)
  X <- as.matrix(X)
  nU <- nrow(U)
  nX <- nrow(X)
  p  <- ncol(U)
  
  UA     <- U / h + alpha
  Utilde <- UA - 1
  
  logX <- log(X)
  
  const <- lgamma(1 / h + p * alpha)
  C_vec <- const - rowSums(lgamma(UA))
  
  B <- Utilde %*% t(logX)
  
  K_log <- sweep(B, 1, C_vec, "+")
  exp(K_log)
}

# -------------------------------------------------------
# Generate a grid on the 3D simplex:
#   x1 + x2 + x3 = 1, xj ≥ 0
# k controls the resolution (step size = 1/k)
# -------------------------------------------------------
make_simplex_grid <- function(k) {
  pts <- list()
  idx <- 1
  for (i in 0:k) {
    for (j in 0:(k - i)) {
      x <- i / k
      y <- j / k
      z <- 1 - x - y
      pts[[idx]] <- c(x, y, z)
      idx <- idx + 1
    }
  }
  do.call(rbind, pts)
}

# -------------------------------------------------------
# Dirichlet kernel matrix K(X, X) for local linear smoother
# X: n × p
# -------------------------------------------------------
build_K_dirichlet_fast <- function(X, alpha, h) {
  X <- as.matrix(X)
  n <- nrow(X)
  p <- ncol(X)
  
  XA    <- X / h + alpha
  const <- lgamma(1 / h + p * alpha)
  C_vec <- const - rowSums(lgamma(XA))
  
  Xtilde <- XA - 1
  logX   <- log(X)
  
  B <- Xtilde %*% t(logX)
  K_log <- sweep(B, 1, C_vec, "+")
  exp(K_log)
}

# -------------------------------------------------------
# Dirichlet kernel matrix K(U, X) for local linear smoother
# U: nU × p
# X: nX × p
# -------------------------------------------------------
build_K_dirichlet_fast_UX <- function(U, X, alpha, h) {
  U <- as.matrix(U)
  X <- as.matrix(X)
  nU <- nrow(U)
  p  <- ncol(U)
  
  UA     <- U / h + alpha
  Utilde <- UA - 1
  logX   <- log(X)
  
  const <- lgamma(1 / h + p * alpha)
  C_vec <- const - rowSums(lgamma(UA))
  
  B <- Utilde %*% t(logX)
  
  K_log <- sweep(B, 1, C_vec, "+")
  exp(K_log)
}

# -------------------------------------------------------
# Local linear Dirichlet smoother S (train | train)
#
# X_kern : n × p_kern (kernel coordinates, e.g. x1,x2,x3)
# X_loc  : n × q      (local-linear coordinates, e.g. x1,x2)
#
# Returns:
#   S : n × n smoother matrix
# -------------------------------------------------------
build_S_locallinear_dirichlet <- function(X_kern, X_loc, alpha, h) {
  X_kern <- as.matrix(X_kern)
  X_loc  <- as.matrix(X_loc)
  n  <- nrow(X_kern)
  q  <- ncol(X_loc)
  
  K <- build_K_dirichlet_fast(X_kern, alpha = alpha, h = h)
  S <- matrix(NA_real_, n, n)
  
  for (i in 1:n) {
    w <- K[i, ]
    W <- diag(w)
    
    Xs <- cbind(1, sweep(X_loc, 2, X_loc[i, ], "-"))
    e1 <- c(1, rep(0, q))
    
    A <- solve(t(Xs) %*% W %*% Xs)
    S[i, ] <- as.numeric(t(e1) %*% A %*% t(Xs) %*% W)
  }
  
  S
}

# -------------------------------------------------------
# Local linear Dirichlet smoother S_test (test | train)
#
# U_kern : nU × p_kern (test kernel coordinates)
# U_loc  : nU × q      (test local-linear coordinates)
# X_kern : nX × p_kern (train kernel coordinates)
# X_loc  : nX × q      (train local-linear coordinates)
#
# Returns:
#   S_test : nU × nX smoother matrix
# -------------------------------------------------------
build_S_locallinear_dirichlet_UX <- function(U_kern, U_loc,
                                             X_kern, X_loc,
                                             alpha, h) {
  U_kern <- as.matrix(U_kern)
  U_loc  <- as.matrix(U_loc)
  X_kern <- as.matrix(X_kern)
  X_loc  <- as.matrix(X_loc)
  
  nU <- nrow(U_kern)
  q  <- ncol(X_loc)
  
  K_UX   <- build_K_dirichlet_fast_UX(U_kern, X_kern, alpha = alpha, h = h)
  nX     <- nrow(X_kern)
  S_test <- matrix(NA_real_, nU, nX)
  
  for (i in 1:nU) {
    w <- K_UX[i, ]
    W <- diag(w)
    
    Xs <- cbind(1, sweep(X_loc, 2, U_loc[i, ], "-"))
    e1 <- c(1, rep(0, q))
    
    A <- solve(t(Xs) %*% W %*% Xs)
    S_test[i, ] <- as.numeric(t(e1) %*% A %*% t(Xs) %*% W)
  }
  
  S_test
}

# -------------------------------------------------------
# Lognormal kernel matrix K(X, X) in log-space
#
# X    : n × d compositional predictors (strictly positive)
# h    : bandwidth
# S    : covariance matrix of log(X) (optional, estimated if NULL)
# ridge: small ridge added to S for numerical stability
#
# Returns:
#   list(K, Z, S, H, Hinv)
# -------------------------------------------------------
lognormal_kernel_matrix <- function(X, h, S = NULL, ridge = 1e-6) {
  X <- as.matrix(X)
  n <- nrow(X)
  d <- ncol(X)
  
  if (any(X <= 0)) {
    stop("All entries of X must be strictly positive for log transformation.")
  }
  
  Z <- log(X)
  
  if (is.null(S)) {
    S <- stats::cov(Z)
  }
  S <- S + ridge * diag(d)
  
  H    <- (h^2) * S
  Hinv <- solve(H)
  
  A    <- Z %*% Hinv
  dvec <- rowSums(A * Z)
  
  Q <- outer(dvec, dvec, "+") - 2 * (A %*% t(Z))
  
  K <- exp(-0.5 * Q)
  
  list(K = K, Z = Z, S = S, H = H, Hinv = Hinv)
}

# -------------------------------------------------------
# Lognormal kernel matrix K_new(X_new, X_train) in log-space
#
# X_train : n × d, training compositional predictors
# X_new   : m × d, new compositional predictors
# h       : bandwidth
# S       : covariance of log(X_train) (optional, estimated if NULL)
# ridge   : ridge term added to S
#
# Returns:
#   K_new : m × n kernel matrix
# -------------------------------------------------------
lognormal_kernel_matrix_new <- function(X_train, X_new, h,
                                        S = NULL, ridge = 1e-6) {
  X_train <- as.matrix(X_train)
  X_new   <- as.matrix(X_new)
  
  if (any(X_train <= 0) || any(X_new <= 0)) {
    stop("All entries of X_train and X_new must be strictly positive for log transformation.")
  }
  
  Zt <- log(X_train)
  Zn <- log(X_new)
  d  <- ncol(Zt)
  
  if (is.null(S)) {
    S <- stats::cov(Zt)
  }
  S <- S + ridge * diag(d)
  
  H    <- (h^2) * S
  Hinv <- solve(H)
  
  At <- Zt %*% Hinv
  An <- Zn %*% Hinv
  
  dt <- rowSums(At * Zt)
  dn <- rowSums(An * Zn)
  
  cross <- An %*% t(Zt)
  Q <- outer(dn, dt, "+") - 2 * cross
  
  K_new <- exp(-0.5 * Q)
  
  K_new
}

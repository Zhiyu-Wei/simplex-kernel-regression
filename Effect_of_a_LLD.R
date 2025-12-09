# -------------------------------------------------------
# Dirichlet Kernel + Local Linear Smoother Utilities
#
# This file contains:
#   - d.kernel, d.kernel_fast
#   - build_K_dirichlet_fast, build_K_dirichlet_fast_UX
#   - build_S_locallinear_dirichlet
#   - build_S_locallinear_dirichlet_UX
#   - fit_dirichlet_plm_ll_dirichlet (partial linear model)
#
# All functions are independent utilities and can be sourced
# into any simulation or analysis script.
# -------------------------------------------------------

# -------------------------------------------------------
# Basic Dirichlet kernel for a single pair (x, u)
# -------------------------------------------------------
d.kernel <- function(x, u, alpha, h, p = 3) {
  if (all(x > 0 & x < 1) & all(u > 0 & u < 1)) {
    result <- lgamma(1 / h + p * alpha) -
      sum(lgamma(x / h + alpha)) +
      sum((x / h + alpha - 1) * log(u))
    return(exp(result))
  } else {
    stop("All entries of x and u must lie in (0, 1).")
  }
}

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
  logX   <- log(X)
  
  const <- lgamma(1 / h + p * alpha)
  C_vec <- const - rowSums(lgamma(UA))
  
  B <- Utilde %*% t(logX)
  
  K_log <- sweep(B, 1, C_vec, "+")
  exp(K_log)
}

# -------------------------------------------------------
# Generate a simplex grid of size (k+1)(k+2)/2
# Each point satisfies x1 + x2 + x3 = 1
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
# Dirichlet kernel matrix K(U, X) for test | train
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
# Local linear smoother S (train|train)
#
# X_kern : n × p_kern (kernel coordinates)
# X_loc  : n × q      (local linear coordinates)
# -------------------------------------------------------
build_S_locallinear_dirichlet <- function(X_kern, X_loc, alpha, h) {
  X_kern <- as.matrix(X_kern)
  X_loc  <- as.matrix(X_loc)
  n  <- nrow(X_kern)
  q  <- ncol(X_loc)
  
  K <- build_K_dirichlet_fast(X_kern, alpha, h)
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
# Local linear smoother S_test (test|train)
#
# U_kern : nU × p_kern
# U_loc  : nU × q
# X_kern : nX × p_kern
# X_loc  : nX × q
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
  
  K_UX <- build_K_dirichlet_fast_UX(U_kern, X_kern, alpha, h)
  S_test <- matrix(NA_real_, nU, nrow(X_kern))
  
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
# Partial linear model with local linear Dirichlet kernel
#
# Model:
#   y = β1 z + m(x) + ε
#
# Returns:
#   list(beta_hat, se_beta1, mse)
# -------------------------------------------------------
fit_dirichlet_plm_ll_dirichlet <- function(
    X_train, z_train, y_train,
    X_test,  z_test,  y_test,
    alpha, h) {
  
  X_train <- as.matrix(X_train)
  X_test  <- as.matrix(X_test)
  
  n <- nrow(X_train)
  I <- diag(n)
  
  # Kernel coordinates = (x1, x2, x3)
  X_kern_train <- X_train
  X_kern_test  <- X_test
  
  # Local coordinates = (x1, x2)
  X_loc_train <- X_train[, 1:2, drop = FALSE]
  X_loc_test  <- X_test[, 1:2, drop = FALSE]
  
  # Smoother matrices
  S_train <- build_S_locallinear_dirichlet(
    X_kern = X_kern_train,
    X_loc  = X_loc_train,
    alpha  = alpha,
    h      = h
  )
  
  S_test <- build_S_locallinear_dirichlet_UX(
    U_kern = X_kern_test,
    U_loc  = X_loc_test,
    X_kern = X_kern_train,
    X_loc  = X_loc_train,
    alpha  = alpha,
    h      = h
  )
  
  # Explicit PLM solution
  W   <- matrix(z_train, ncol = 1)
  A   <- I - S_train
  Xt  <- t(W)
  
  XtAW <- Xt %*% A %*% W
  XtAy <- Xt %*% A %*% y_train
  beta_hat <- solve(XtAW, XtAy)
  
  # Nonparametric component
  mu_hat_train <- S_train %*% (y_train - W %*% beta_hat)
  
  # Predictions
  y_lin_test    <- z_test * as.numeric(beta_hat)
  y_nonlin_test <- S_test %*% mu_hat_train
  y_hat_test    <- as.numeric(y_lin_test + y_nonlin_test)
  
  mse_test <- mean((y_test - y_hat_test)^2)
  
  # Variance estimate for beta
  residual   <- (I - S_train) %*% (y_train - W %*% beta_hat)
  edf        <- sum(diag(S_train))
  sigma2_hat <- sum(residual^2) / (n - edf)
  cov_beta   <- sigma2_hat * solve(XtAW)
  se_beta1   <- sqrt(cov_beta[1, 1])
  
  list(
    beta_hat = as.numeric(beta_hat),
    se_beta1 = se_beta1,
    mse      = mse_test
  )
}

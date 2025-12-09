# --------------------------------------------------------
# Bandwidth selection for Dirichlet kernel regression
# - Nadaraya–Watson Dirichlet (NWD)
# - Local Linear Dirichlet (LLD)
# --------------------------------------------------------

# This file only defines generic utilities and CV-based
# bandwidth selectors. Parallel simulation or data
# generation can be written separately.
#
# All functions assume compositional predictors on the
# simplex (each row lies in (0,1)^p and sums to 1).

# --------------------------------------------------------
# Dirichlet kernel utilities
# --------------------------------------------------------

# Single-point Dirichlet kernel:
#   K_h(u | x) = C(x; h, alpha) * prod_j u_j^{x_j/h + alpha - 1}
# x, u: numeric vectors of length p in (0,1)
dirichlet_kernel <- function(x, u, alpha, h) {
  x <- as.numeric(x)
  u <- as.numeric(u)
  p <- length(x)
  
  if (!all(is.finite(x)) || !all(is.finite(u))) {
    stop("Non-finite values in x or u.")
  }
  if (any(x <= 0) || any(x >= 1) || any(u <= 0) || any(u >= 1)) {
    stop("All entries of x and u must lie in (0, 1).")
  }
  if (!is.finite(h) || h <= 0) {
    stop("Bandwidth h must be positive and finite.")
  }
  
  xa    <- x / h + alpha                  # length p
  const <- lgamma(1 / h + p * alpha) - sum(lgamma(xa))
  val   <- const + sum((xa - 1) * log(u))
  
  exp(val)
}

# Kernel matrix K(X, X): n x n, K_ij = K_h(x_j | x_i)
# X: n x p matrix of simplex predictors
dirichlet_kernel_matrix <- function(X, alpha, h) {
  X <- as.matrix(X)     # n x p
  n <- nrow(X)
  p <- ncol(X)
  
  if (!is.finite(h) || h <= 0) {
    stop("Bandwidth h must be positive and finite.")
  }
  
  XA    <- X / h + alpha                  # n x p
  const <- lgamma(1 / h + p * alpha)
  C_vec <- const - rowSums(lgamma(XA))    # length n
  
  Xtilde <- XA - 1                        # n x p
  logX   <- log(X)                        # n x p
  
  # B_ij = sum_k (x_i,k/h + alpha - 1) log x_jk
  B <- Xtilde %*% t(logX)                 # n x n
  
  K_log <- sweep(B, 1, C_vec, "+")       # add C(u_i) rowwise
  K     <- exp(K_log)
  K[!is.finite(K)] <- 0
  
  K
}

# Kernel matrix K(U, X): n_U x n_X, K_ij = K_h(x_j | u_i)
# Useful for prediction at new points U
dirichlet_kernel_UX <- function(U, X, alpha, h) {
  U <- as.matrix(U)   # n_U x p
  X <- as.matrix(X)   # n_X x p
  nU <- nrow(U)
  nX <- nrow(X)
  p  <- ncol(U)
  
  if (!is.finite(h) || h <= 0) {
    stop("Bandwidth h must be positive and finite.")
  }
  
  UA     <- U / h + alpha
  Utilde <- UA - 1
  logX   <- log(X)
  
  const <- lgamma(1 / h + p * alpha)
  C_vec <- const - rowSums(lgamma(UA))    # length n_U
  
  B <- Utilde %*% t(logX)                # n_U x n_X
  
  K_log <- sweep(B, 1, C_vec, "+")
  K     <- exp(K_log)
  K[!is.finite(K)] <- 0
  
  K
}

# --------------------------------------------------------
# 1. Nadaraya–Watson Dirichlet (NWD):
#    CV-based bandwidth selection
# --------------------------------------------------------

# Leave-one-out CV for NWD:
# x   : n x p matrix of simplex predictors
# y   : length-n response vector
# alpha: centering parameter for Dirichlet kernel
# interval: search interval for bandwidth h (c(lower, upper))
cv_bw_dirichlet_nw <- function(x, y, alpha, interval = c(1e-4, 1)) {
  x <- as.matrix(x)
  y <- as.numeric(y)
  n <- nrow(x)
  
  if (length(y) != n) {
    stop("Length of y must match nrow(x).")
  }
  
  I_n <- diag(n)
  
  cv_fun <- function(h) {
    # penalize non-positive or non-finite bandwidths
    if (!is.finite(h) || h <= 0) return(Inf)
    
    S <- dirichlet_kernel_matrix(x, alpha, h)
    if (any(!is.finite(S))) return(Inf)
    
    # Leave-one-out: zero out diagonal
    S0 <- S
    diag(S0) <- 0
    
    # W^{-1} = diag(1 / w_i), where w_i = sum_j S_ij - S_ii
    row_sum <- rowSums(S)
    w_diag  <- row_sum - diag(S)
    
    eps <- 1e-8
    w_diag[!is.finite(w_diag) | w_diag <= eps] <- eps
    
    W_inv <- diag(1 / w_diag, nrow = n, ncol = n)
    
    # A = I - W^{-1} S0
    A  <- I_n - W_inv %*% S0
    Ay <- A %*% y
    
    mean(Ay^2)
  }
  
  interval <- sort(interval)
  opt <- optimize(cv_fun, interval = interval)
  opt$minimum
}

# Safe wrapper: return NA_real_ if CV fails
safe_cv_bw_dirichlet_nw <- function(x, y, alpha, interval = c(1e-4, 1)) {
  out <- try(
    cv_bw_dirichlet_nw(x, y, alpha, interval = interval),
    silent = TRUE
  )
  if (inherits(out, "try-error") || !is.finite(out)) NA_real_ else out
}

# --------------------------------------------------------
# 2. Local Linear Dirichlet (LLD):
#    CV-based bandwidth selection
# --------------------------------------------------------

# Local linear Dirichlet estimator with LOOCV for bandwidth.
#
# x       : n x p matrix of simplex predictors
# y       : length-n response vector
# alpha   : centering parameter for Dirichlet kernel
# interval: search interval for h
# ll_index: which columns of x enter the local-linear expansion
#           (by default, first p-1 coordinates)
# ridge   : small ridge added to X'WX for numerical stability
cv_bw_dirichlet_lld <- function(x, y, alpha,
                                interval = c(1e-3, 1),
                                ll_index = NULL,
                                ridge = 1e-8) {
  x <- as.matrix(x)
  y <- as.numeric(y)
  n <- nrow(x)
  p <- ncol(x)
  
  if (length(y) != n) {
    stop("Length of y must match nrow(x).")
  }
  
  if (is.null(ll_index)) {
    if (p < 2) stop("x must have at least 2 columns for local linear LLD.")
    ll_index <- 1:(p - 1)
  }
  
  X_ll    <- x[, ll_index, drop = FALSE]
  interval <- sort(interval)
  
  cv_fun <- function(h) {
    if (!is.finite(h) || h <= 0) return(1e12)
    
    # Kernel weights S(i, j) = K_h(x_j | x_i)
    S <- dirichlet_kernel_matrix(x, alpha, h)
    if (all(S == 0) || any(!is.finite(S))) return(1e12)
    
    pred <- numeric(n)
    
    for (i in seq_len(n)) {
      w <- S[i, ]
      w[i] <- 0   # leave-one-out
      
      # If all weights degenerate, fall back to y_i
      if (all(w <= 0) || all(!is.finite(w))) {
        pred[i] <- y[i]
        next
      }
      
      # Weighted local linear fit at x_i
      W  <- diag(w, nrow = n, ncol = n)
      Xs <- cbind(1, sweep(X_ll, 2, X_ll[i, ], "-"))  # (n x q), q = 1 + dim(ll_index)
      e1 <- c(1, rep(0, ncol(Xs) - 1))                # picks intercept
      
      XtW  <- t(Xs) %*% W
      XtWX <- XtW %*% Xs
      XtWX <- XtWX + diag(ridge, ncol(XtWX))          # ridge stabilization
      
      beta_hat <- tryCatch(
        solve(XtWX, XtW %*% y),
        error = function(e) rep(NA_real_, ncol(Xs))
      )
      
      if (any(!is.finite(beta_hat))) {
        pred[i] <- y[i]
      } else {
        pred[i] <- sum(e1 * beta_hat)
      }
    }
    
    mean((y - pred)^2)
  }
  
  opt <- optimize(cv_fun, interval = interval)
  opt$minimum
}

# Safe wrapper: return NA_real_ if CV fails
safe_cv_bw_dirichlet_lld <- function(x, y, alpha,
                                     interval = c(1e-3, 1),
                                     ll_index = NULL,
                                     ridge = 1e-8) {
  out <- tryCatch(
    cv_bw_dirichlet_lld(
      x        = x,
      y        = y,
      alpha    = alpha,
      interval = interval,
      ll_index = ll_index,
      ridge    = ridge
    ),
    error   = function(e) NA_real_,
    warning = function(w) {
      # Sometimes warnings still lead to a usable solution; try once more
      suppressWarnings(
        tryCatch(
          cv_bw_dirichlet_lld(
            x        = x,
            y        = y,
            alpha    = alpha,
            interval = interval,
            ll_index = ll_index,
            ridge    = ridge
          ),
          error = function(e2) NA_real_
        )
      )
    }
  )
  
  if (!is.finite(out)) NA_real_ else out
}

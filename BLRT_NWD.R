# -------------------------------------------------------
# Dirichlet Kernel Utilities for Nadaraya–Watson Regression
#
# Contains:
#   - d.kernel_fast()           : K(U, X) for general U, X
#   - build_K_dirichlet_fast()  : K(X, X) for training data
#   - build_K_dirichlet_fast_UX(): K(U, X) for new vs train
#   - CV.Best.h()               : bandwidth selection by CV
#
# Notes:
#   * All functions assume compositional predictors in (0,1)
#     with rows summing to 1 (simplex points).
#   * CV.Best.h() requires build_K_dirichlet_fast().
#   * Comments and messages are in English.
# -------------------------------------------------------

# -------------------------------------------------------
# Fast Dirichlet kernel matrix K(U, X)
#
# U : nU × p matrix of simplex points
# X : nX × p matrix of simplex points
# alpha : Dirichlet parameter (> 0)
# h     : bandwidth (> 0)
#
# Returns:
#   nU × nX kernel matrix K(U_i, X_j)
# -------------------------------------------------------
d.kernel_fast <- function(U, X, alpha, h) {
  U <- as.matrix(U)   # nU × p
  X <- as.matrix(X)   # nX × p
  nU <- nrow(U)
  nX <- nrow(X)
  p  <- ncol(U)
  
  # u/h + alpha (nU × p)
  UA     <- U / h + alpha
  Utilde <- UA - 1
  
  # log X  (nX × p)
  logX <- log(X)
  
  # C(u) = lgamma(1/h + p*alpha) - sum_k lgamma(UA_ik)
  const <- lgamma(1 / h + p * alpha)
  C_vec <- const - rowSums(lgamma(UA))   # length nU
  
  # B_ij = sum_k (u_k/h + alpha - 1) * log x_{jk}
  #      = Utilde[i, ] · logX[j, ]'
  # → (nU × p) %*% (p × nX) = (nU × nX)
  B <- Utilde %*% t(logX)
  
  # log K = C(u_i) + B_ij
  K_log <- sweep(B, 1, C_vec, "+")
  
  exp(K_log)
}

# -------------------------------------------------------
# Dirichlet kernel matrix K(X, X)
#
# X     : n × p matrix of simplex points
# alpha : Dirichlet parameter (> 0)
# h     : bandwidth (> 0)
#
# Returns:
#   n × n kernel matrix K(X_i, X_j)
# -------------------------------------------------------
build_K_dirichlet_fast <- function(X, alpha, h) {
  X <- as.matrix(X)     # n × p
  n <- nrow(X)
  p <- ncol(X)
  
  XA    <- X / h + alpha           # n × p
  const <- lgamma(1 / h + p * alpha)
  C_vec <- const - rowSums(lgamma(XA))  # length n
  
  Xtilde <- XA - 1                 # n × p
  logX   <- log(X)                 # n × p
  
  # B_ij = sum_k (x_ik/h + alpha - 1) * log x_jk
  B <- Xtilde %*% t(logX)          # n × n
  
  K_log <- sweep(B, 1, C_vec, "+")
  exp(K_log)                       # n × n
}

# -------------------------------------------------------
# Dirichlet kernel matrix K(U, X) for new vs train
#
# U     : nU × p matrix (e.g., new points)
# X     : nX × p matrix (training points)
# alpha : Dirichlet parameter (> 0)
# h     : bandwidth (> 0)
#
# Returns:
#   nU × nX kernel matrix K(U_i, X_j)
# -------------------------------------------------------
build_K_dirichlet_fast_UX <- function(U, X, alpha, h) {
  U <- as.matrix(U)   # nU × p
  X <- as.matrix(X)   # nX × p
  nU <- nrow(U)
  p  <- ncol(U)
  
  UA     <- U / h + alpha          # nU × p
  Utilde <- UA - 1                 # nU × p
  logX   <- log(X)                 # nX × p
  
  const <- lgamma(1 / h + p * alpha)
  C_vec <- const - rowSums(lgamma(UA))   # length nU
  
  B <- Utilde %*% t(logX)          # nU × nX
  
  K_log <- sweep(B, 1, C_vec, "+")
  exp(K_log)                       # nU × nX
}

# -------------------------------------------------------
# Bandwidth selection for Nadaraya–Watson Dirichlet kernel
# via cross-validation.
#
# data         : n × p matrix of simplex predictors
# pred         : length-n response vector
# alp          : Dirichlet parameter (> 0)
# opt.interval : search interval for h (lower, upper)
#
# Returns:
#   h_hat : CV-selected bandwidth
# -------------------------------------------------------
CV.Best.h <- function(data, pred, alp,
                      opt.interval = c(0.00001, 0.1)) {
  data <- as.matrix(data)
  pred <- as.vector(pred)
  
  CV1 <- function(h) {
    n <- nrow(data)
    
    # Dirichlet kernel matrix K(X, X)
    S <- build_K_dirichlet_fast(data, alp, h)
    
    # Identity matrix and auxiliary objects
    I   <- diag(n)
    one <- matrix(1, n, 1)
    
    # S_0: kernel matrix with zero diagonal
    S_0 <- S
    diag(S_0) <- 0
    
    # Diagonal weight matrix W
    # W_ii = sum_j K_ij - K_ii
    W <- diag(n)
    diag(W) <- apply(I, 2, function(i) (t(i) %*% S %*% one - S[i, i]))
    Winv <- solve(W)
    
    y_i <- matrix(pred, n, 1)
    
    # CV criterion: (1/n) * || (I - W^{-1} S_0) y ||^2
    res <- (I - Winv %*% S_0) %*% y_i
    as.numeric((1 / n) * t(res) %*% res)
  }
  
  result <- optimize(CV1, interval = opt.interval)
  result$minimum
}

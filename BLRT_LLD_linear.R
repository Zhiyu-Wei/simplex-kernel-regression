## ============================================================
##  Dirichlet Kernel + Local Linear Smoother Utility Functions
##  (English-annotated version)
## ============================================================

library(MCMCpack)   # used only if rdirichlet is needed in simulation


## ============================================================
## 1. Fast Dirichlet kernel K(U, X)
##    Returns an nU × nX kernel matrix
## ============================================================
d.kernel_fast <- function(U, X, alpha, h) {
  U <- as.matrix(U)   # nU × p
  X <- as.matrix(X)   # nX × p
  p <- ncol(U)
  
  # Compute U/h + α
  UA <- U / h + alpha
  Utilde <- UA - 1
  
  # log(X)
  logX <- log(X)
  
  # Constant term: lgamma(1/h + pα) − Σ lgamma(UA)
  const <- lgamma(1/h + p * alpha)
  C_vec <- const - rowSums(lgamma(UA))   # length nU
  
  # B_ij = Σ_k (u_k/h + α − 1) log(x_jk)
  # Result is nU × nX
  B <- Utilde %*% t(logX)
  
  # log K = C(u_i) + B_ij
  K_log <- sweep(B, 1, C_vec, "+")
  
  exp(K_log)
}


## ============================================================
## 2. Dirichlet kernel K(X, X)
##    Square n × n kernel matrix used for smoothing
## ============================================================
build_K_dirichlet_fast <- function(X, alpha, h) {
  X <- as.matrix(X)
  p <- ncol(X)
  
  XA <- X / h + alpha
  const <- lgamma(1/h + p * alpha)
  C_vec <- const - rowSums(lgamma(XA))     # length n
  
  Xtilde <- XA - 1
  logX <- log(X)
  
  # B_ij = Σ_k (x_ik/h + α − 1) log(x_jk)
  B <- Xtilde %*% t(logX)
  
  K_log <- sweep(B, 1, C_vec, "+")
  exp(K_log)
}


## ============================================================
## 3. Dirichlet kernel K(U, X)
##    Similar to above but allows U ≠ X (e.g., prediction points)
## ============================================================
build_K_dirichlet_fast_UX <- function(U, X, alpha, h) {
  U <- as.matrix(U)
  X <- as.matrix(X)
  p <- ncol(U)
  
  UA <- U / h + alpha
  Utilde <- UA - 1
  logX <- log(X)
  
  const <- lgamma(1/h + p * alpha)
  C_vec <- const - rowSums(lgamma(UA))
  
  B <- Utilde %*% t(logX)
  K_log <- sweep(B, 1, C_vec, "+")
  exp(K_log)
}


## ============================================================
## 4. Bandwidth selection via LOOCV for Local Linear Smoother
##    Returns optimal h ∈ (interval)
## ============================================================
WLS_best_h <- function(data, obj, alpha, interval = c(0.001, 1)) {
  data <- as.matrix(data)
  n <- nrow(data)
  
  # Local-linear uses first two coordinates of X (your convention)
  X <- data[, 1:2, drop = FALSE]
  
  # ----------------------------------------
  # LOOCV objective
  # ----------------------------------------
  WLS_cv <- function(h) {
    pred <- numeric(n)
    
    for (i in 1:n) {
      # Compute all kernel weights K(u_i, X)
      K_mat <- d.kernel_fast(U = data,
                             X = data[i, , drop = FALSE],
                             alpha = alpha,
                             h = h)
      
      w <- as.numeric(K_mat[, 1])
      w[i] <- 0  # Remove self-weight (LOOCV)
      
      W <- diag(w)
      
      # Local linear design: [1, X − X_i]
      Xs <- cbind(1, sweep(X, 2, X[i, ], "-"))
      e1 <- c(1, rep(0, ncol(Xs) - 1))
      
      # Weighted least squares
      beta_hat <- solve(t(Xs) %*% W %*% Xs,
                        t(Xs) %*% W %*% obj)
      
      pred[i] <- drop(e1 %*% beta_hat)
    }
    
    mean((obj - pred)^2)
  }
  
  # Optimize h in the given interval
  opt <- optimize(WLS_cv, interval = interval)
  opt$minimum
}


## ============================================================
## 5. Local Linear Dirichlet Smoother Matrix S
##    Returns S (n × n)
## ============================================================
build_S_locallinear_dirichlet <- function(X_kern, X_loc, alpha, h) {
  X_kern <- as.matrix(X_kern)   # for kernel weights
  X_loc  <- as.matrix(X_loc)    # for local linear component
  n <- nrow(X_kern)
  
  # Precompute kernel matrix K(X, X)
  K <- build_K_dirichlet_fast(X_kern, alpha, h)
  
  S <- matrix(NA, n, n)
  q <- ncol(X_loc)
  
  for (i in 1:n) {
    w <- K[i, ]                 # weights for point i
    W <- diag(w)
    
    # Local-linear design matrix
    Xs <- cbind(1, sweep(X_loc, 2, X_loc[i, ], "-"))
    e1 <- c(1, rep(0, q))
    
    # A = (X'WX)^(-1)
    A <- solve(t(Xs) %*% W %*% Xs)
    
    # s_i' = e1' A X' W
    S[i, ] <- as.numeric(t(e1) %*% A %*% t(Xs) %*% W)
  }
  
  S
}

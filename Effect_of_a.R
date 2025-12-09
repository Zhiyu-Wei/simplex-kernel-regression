# -------------------------------------------------------
# Kernel Regression with Simplex Predictors
# Utility functions for Dirichlet kernel on the simplex
# -------------------------------------------------------

#' Compute the Dirichlet kernel K_h(u | x) at a single pair (x, u)
#'
#' @param x Numeric vector of length p, compositional predictor (each in (0,1))
#' @param u Numeric vector of length p, evaluation point (each in (0,1))
#' @param alpha Centring parameter (>0)
#' @param h Bandwidth (h > 0)
#'
#' @return A scalar, Dirichlet kernel value K_h(u | x).
#'
dirichlet_kernel_scalar <- function(x, u, alpha, h) {
  x <- as.numeric(x)
  u <- as.numeric(u)
  p <- length(x)
  
  if (!all(is.finite(x)) || !all(is.finite(u))) {
    stop("x and u must be finite.")
  }
  if (!all(x > 0 & x < 1) || !all(u > 0 & u < 1)) {
    stop("All entries of x and u must lie in (0, 1).")
  }
  if (h <= 0) {
    stop("Bandwidth h must be positive.")
  }
  
  xa     <- x / h + alpha
  const  <- lgamma(1 / h + p * alpha) - sum(lgamma(xa))
  val_log <- const + sum((xa - 1) * log(u))
  
  exp(val_log)
}

#' Fast matrix version of the Dirichlet kernel
#'
#' Computes an n_U x n_X matrix with entries K_h(U_i | X_j),
#' where rows of U, X are points on the simplex.
#'
#' @param U Matrix (n_U x p) of evaluation points
#' @param X Matrix (n_X x p) of design points
#' @param alpha Centring parameter
#' @param h Bandwidth (h > 0)
#'
#' @return An n_U x n_X matrix of kernel values.
#'
dirichlet_kernel_matrix <- function(U, X, alpha, h) {
  U <- as.matrix(U)   # n_U x p
  X <- as.matrix(X)   # n_X x p
  
  nU <- nrow(U)
  nX <- nrow(X)
  p  <- ncol(U)
  
  if (h <= 0) {
    stop("Bandwidth h must be positive.")
  }
  
  if (any(U <= 0) || any(U >= 1) || any(X <= 0) || any(X >= 1)) {
    stop("All entries of U and X must lie in (0, 1).")
  }
  
  # U / h + alpha (n_U x p)
  UA     <- U / h + alpha
  Utilde <- UA - 1
  
  # log X (n_X x p)
  logX <- log(X)
  
  # C(u) = lgamma(1/h + p*alpha) - sum_k lgamma(UA_ik)
  const <- lgamma(1 / h + p * alpha)
  C_vec <- const - rowSums(lgamma(UA))   # length n_U
  
  # B_ij = sum_k (u_k/h + alpha - 1) * log x_{jk}
  #      = Utilde[i, ] %*% logX[j, ]'
  # (n_U x p) %*% (p x n_X) = (n_U x n_X)
  B <- Utilde %*% t(logX)
  
  # log K = C(u_i) + B_ij
  K_log <- sweep(B, 1, C_vec, "+")
  
  K <- exp(K_log)
  K[!is.finite(K)] <- 0
  K
}

#' Generate a regular grid on the 2-simplex with mesh size 1/k
#'
#' @param k Positive integer, number of subdivisions.
#'
#' @return A matrix with rows (x1, x2, x3) such that x1 + x2 + x3 = 1.
#'
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

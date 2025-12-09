# -------------------------------------------------------
# Kernel Regression with Simplex Predictors
# Simulation: partial linear model with Dirichlet kernel
# alpha = 0, 0.5, 1 and fixed bandwidths
# Produces summary tables for beta1 and prediction MSE.
# -------------------------------------------------------

library(foreach)
library(doParallel)
library(MCMCpack)   # for rdirichlet

# Source Dirichlet kernel utilities (in R branch / R folder)
source("https://raw.githubusercontent.com/Zhiyu-Wei/simplex-kernel-regression/R/Effect_of_a.R")


# -------------------------------------------------------
# True regression components on the simplex
# m_1, ..., m_4
# -------------------------------------------------------
r1 <- 1
r2 <- 1.5
r3 <- -0.5

m_list <- list(
  function(s1, s2) r1 * s1 + r2 * s2 + r3 * (1 - s1 - s2),
  function(s1, s2) log(1 + s1 + s2),
  function(s1, s2) s1 * s2,
  function(s1, s2) cos(6 * pmin(s1, s2, 1 - s1 - s2))
)

# -------------------------------------------------------
# Helper: fit partial linear model with Dirichlet kernel
# y = beta1 * z1 + m(x) + error
#
# Inputs:
#   X_train  : n x 3 matrix of simplex covariates (training)
#   z_train  : length-n vector, linear covariate (training)
#   y_train  : length-n vector, response (training)
#   X_test   : n_test x 3 matrix, test covariates
#   z_test   : length-n_test vector, test linear covariate
#   y_test   : length-n_test vector, test response
#   alpha, h : Dirichlet kernel parameters
#
# Output:
#   list(beta_hat, se_beta1, mse)
#   mse is test prediction MSE (used as accuracy criterion).
# -------------------------------------------------------
fit_dirichlet_plm <- function(X_train, z_train, y_train,
                              X_test,  z_test,  y_test,
                              alpha, h) {
  X_train <- as.matrix(X_train)
  X_test  <- as.matrix(X_test)
  
  n <- nrow(X_train)
  I <- diag(n)
  
  # S_ij = K_h(x_j | x_i)
  S <- dirichlet_kernel_matrix(X_train, X_train, alpha = alpha, h = h)
  row_sum <- rowSums(S)
  # Avoid division by zero (should not happen, but keep safe)
  row_sum[row_sum == 0] <- 1e-12
  Sij <- S / row_sum
  
  # S_test: kernel between test x and train x
  S_test <- dirichlet_kernel_matrix(X_test, X_train, alpha = alpha, h = h)
  row_sum_test <- rowSums(S_test)
  row_sum_test[row_sum_test == 0] <- 1e-12
  Sij_test <- S_test / row_sum_test
  
  W   <- matrix(z_train, ncol = 1)  # n x 1
  obj <- y_train
  
  A <- I - Sij
  Xt <- t(W)
  
  # Explicit partial linear solution:
  # beta_hat = (W^T (I - S) W)^{-1} W^T (I - S) y
  XtAW <- Xt %*% A %*% W
  XtAy <- Xt %*% A %*% obj
  beta_hat <- solve(XtAW, XtAy)
  
  # Nonparametric component at training points
  mu_hat_train <- Sij %*% (y_train - W %*% beta_hat)
  
  # Predictions at test points
  y_lin_test    <- z_test * as.numeric(beta_hat)
  y_nonlin_test <- Sij_test %*% mu_hat_train
  y_hat_test    <- as.numeric(y_lin_test + y_nonlin_test)
  
  mse_test <- mean((y_test - y_hat_test)^2)
  
  # Estimate variance of beta1 using EDF-based sigma^2
  residual   <- (I - Sij) %*% (y_train - W %*% beta_hat)
  edf        <- sum(diag(Sij))             # effective degrees of freedom
  sigma2_hat <- sum(residual^2) / (n - edf)
  
  cov_beta   <- sigma2_hat * solve(XtAW)
  se_beta1   <- sqrt(cov_beta[1, 1])
  
  list(
    beta_hat = as.numeric(beta_hat),
    se_beta1 = se_beta1,
    mse      = mse_test
  )
}

# -------------------------------------------------------
# Simulation settings
# -------------------------------------------------------

B      <- 100   # number of Monte Carlo replications
beta1  <- 2
N      <- 100    # sample size
shape  <- c(5, 5, 5)  # Dirichlet concentration for covariates

# Bandwidths for each m_k and each alpha
# h_alpha0[k], h_alpha05[k], h_alpha1[k]
h_alpha0  <- c(0.01639899, 0.02058882, 0.02634324, 0.01274597)
h_alpha05 <- c(0.02876762, 0.03786284, 0.04176682, 0.01426423)
h_alpha1  <- c(0.02018845, 0.02500206, 0.03002680, 0.01280003)

# Fixed 1000 test grid on the simplex (interior points)
grid_all <- as.data.frame(make_simplex_grid(47))
colnames(grid_all) <- c("x1", "x2", "x3")
grid_all <- subset(
  grid_all,
  x1 > 0.01 & x1 < 1 &
    x2 > 0.01 & x2 < 1 &
    x3 > 0.01 & x3 < 1
)
# Optionally thin out (as in original code)
grid_all <- grid_all[-round(seq(1, nrow(grid_all), length.out = 35), 0), ]
test_n   <- nrow(grid_all)

# -------------------------------------------------------
# Parallel setup
# -------------------------------------------------------

cores <- max(1, parallel::detectCores() - 2)
cl    <- parallel::makeCluster(cores)
doParallel::registerDoParallel(cl)

dir.create("results", showWarnings = FALSE)

start.time <- Sys.time()

# -------------------------------------------------------
# Main simulation over m_1, ..., m_4
# For each true function, compute summaries for alpha=0, 0.5, 1
# -------------------------------------------------------

for (k in seq_along(m_list)) {
  f_k     <- m_list[[k]]
  h0      <- h_alpha0[k]
  h05     <- h_alpha05[k]
  h1      <- h_alpha1[k]
  
  result_mat <- foreach(
    rep = 1:B,
    .combine  = rbind,
    .packages = c("MCMCpack")
  ) %dopar% {
    # --- Generate training data on simplex ---
    X_train <- MCMCpack::rdirichlet(N, shape)
    x1 <- X_train[, 1]
    x2 <- X_train[, 2]
    
    z_train <- rnorm(N, mean = 0, sd = 0.01)
    eps     <- rnorm(N, mean = 0, sd = 0.01)
    y_train <- beta1 * z_train + f_k(x1, x2) + eps
    
    # --- Generate test data on fixed grid ---
    X_test <- as.matrix(grid_all)
    z_test <- rnorm(test_n, mean = 0, sd = 0.01)
    eps_t  <- rnorm(test_n, mean = 0, sd = 0.01)
    y_test <- beta1 * z_test + f_k(X_test[, 1], X_test[, 2]) + eps_t
    
    # alpha = 0
    fit0 <- fit_dirichlet_plm(
      X_train, z_train, y_train,
      X_test,  z_test,  y_test,
      alpha = 0, h = h0
    )
    
    # alpha = 0.5
    fit05 <- fit_dirichlet_plm(
      X_train, z_train, y_train,
      X_test,  z_test,  y_test,
      alpha = 0.5, h = h05
    )
    
    # alpha = 1
    fit1 <- fit_dirichlet_plm(
      X_train, z_train, y_train,
      X_test,  z_test,  y_test,
      alpha = 1, h = h1
    )
    
    c(
      fit0$beta_hat,  fit0$se_beta1,  fit0$mse,
      fit05$beta_hat, fit05$se_beta1, fit05$mse,
      fit1$beta_hat,  fit1$se_beta1,  fit1$mse
    )
  }
  
  # Extract blocks for each alpha
  res0  <- result_mat[, 1:3]
  res05 <- result_mat[, 4:6]
  res1  <- result_mat[, 7:9]
  
  # Column summaries: mean(beta_hat), mean(se_beta1), mean(mse)
  mean0  <- colMeans(res0)
  mean05 <- colMeans(res05)
  mean1  <- colMeans(res1)
  
  # SD of beta_hat and mse (for table reporting)
  sd_beta0  <- sd(res0[, 1])
  sd_mse0   <- sd(res0[, 3])
  sd_beta05 <- sd(res05[, 1])
  sd_mse05  <- sd(res05[, 3])
  sd_beta1  <- sd(res1[, 1])
  sd_mse1   <- sd(res1[, 3])
  
  table_res <- rbind(
    c(mean0[1],  sd_beta0,  mean0[2], mean0[3], sd_mse0),
    c(mean05[1], sd_beta05, mean05[2], mean05[3], sd_mse05),
    c(mean1[1],  sd_beta1,  mean1[2], mean1[3], sd_mse1),
    c(beta1,     NA,        NA,       NA,       NA)
  )
  
  rownames(table_res) <- c("alpha=0", "alpha=0.5", "alpha=1", "True beta1")
  colnames(table_res) <- c(
    "Mean(beta1_hat)",
    "SD(beta1_hat)",
    "Mean(SE(beta1_hat))",
    "MSE",            # prediction MSE on test grid
    "SD(MSE)"
  )
  
  # Write result table for m_k
  outfile <- file.path("results",
                       sprintf("dirichlet_plm_m%d_alpha_0_0.5_1.csv", k))
  write.csv(round(table_res, 6), outfile, row.names = TRUE)
}

parallel::stopCluster(cl)

end.time <- Sys.time()
print(end.time - start.time)

# -------------------------------------------------------
# Kernel Regression with Simplex Predictors
# Simulation: coefficient accuracy (boundary design)
#
# Model:
#   y = beta1 * z1 + m_k(x) + error
#
# x   ~ Dirichlet(0.6, 0.6, 0.6)    (boundary-concentrated design)
# z1  ~ N(0, 0.01)
# eps ~ N(0, 0.01)
#
# Methods compared:
#   1. NW(Dirichlet)    (alpha = 0, bandwidth h)
#   2. NW(Lognormal)    (Gaussian kernel in log-space, bandwidth h_logn)
#   3. Local linear Dirichlet (alpha = 1, bandwidth h_wls)
#   4. Linear model (ordinary least squares)
#
# For each method and m_k, we summarize:
#   - Mean(beta1_hat)
#   - SD(beta1_hat)
#   - Mean(SE(beta1_hat))
#   - ISE (approximated by MSE on a fixed simplex test grid)
#   - SD(ISE)
#
# Output:
#   CSV files in ./results/coef_accuracy_boundary/
# -------------------------------------------------------

library(foreach)
library(doParallel)
library(MCMCpack)  # for rdirichlet

# Source Dirichlet local linear utilities
source("https://raw.githubusercontent.com/Zhiyu-Wei/simplex-kernel-regression/R/Kernel_method_comparison.R")

# -------------------------------------------------------
# True regression components m_1, ..., m_4
# -------------------------------------------------------
r1 <-  1
r2 <-  1.5
r3 <- -0.5

m_list <- list(
  function(s1, s2) r1 * s1 + r2 * s2 + r3 * (1 - s1 - s2),
  function(s1, s2) log(1 + s1 + s2),
  function(s1, s2) s1 * s2,
  function(s1, s2) cos(6 * pmin(s1, s2, 1 - s1 - s2))
)

# -------------------------------------------------------
# Parallel setup
# -------------------------------------------------------
cores <- 12
cl    <- makeCluster(cores)
registerDoParallel(cl)

start.time <- Sys.time()

# -------------------------------------------------------
# Simulation settings
# -------------------------------------------------------
B       <- 1000       # number of Monte Carlo replications
h       <- 0.05       # bandwidth for NW(Dirichlet)
h_wls   <- 0.05       # bandwidth for local linear Dirichlet
h_logn  <- 0.5        # bandwidth for lognormal kernel
beta1   <- 2
N       <- 100        # sample size

# -------------------------------------------------------
# Fixed test grid on the simplex (interior points)
# -------------------------------------------------------
data1_test <- as.data.frame(make_simplex_grid(50))
colnames(data1_test) <- c("x1", "x2", "x3")

# keep interior points only
data1_test <- subset(
  data1_test,
  x1 > 0.03 & x1 < 1 &
    x2 > 0.03 & x2 < 1 &
    x3 > 0.03 & x3 < 1
)

# thin out grid (mimic original code)
data1_test <- data1_test[-round(seq(1, 1035, length.out = 35), 0), ]
test_n     <- nrow(data1_test)

# -------------------------------------------------------
# Output directory
# -------------------------------------------------------
out_dir <- file.path("results", "coef_accuracy_boundary")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# -------------------------------------------------------
# Main simulation loop over m_k
# -------------------------------------------------------
for (k in 1:4) {
  f <- m_list[[k]]
  cat(sprintf(">>> Simulating for m_%d ...\n", k))
  
  result_list <- foreach(
    b = 1:B,
    .combine  = rbind,
    .packages = c("MCMCpack")
  ) %dopar% {
    
    # ---------------------------------------------------
    # Generate training data
    # ---------------------------------------------------
    alpha_dir <- 0
    X <- rdirichlet(N, c(0.6, 0.6, 0.6))
    x1 <- X[, 1]
    x2 <- X[, 2]
    x3 <- X[, 3]
    
    z1 <- rnorm(N, 0, 0.01)
    eps <- rnorm(N, 0, 0.01)
    y   <- beta1 * z1 + f(x1, x2) + eps
    
    data1_train <- data.frame(
      x1 = x1,
      x2 = x2,
      x3 = x3,
      z1 = z1,
      y  = y
    )
    
    # ---------------------------------------------------
    # Generate test data on fixed grid
    # ---------------------------------------------------
    z1_test <- rnorm(test_n, 0, 0.01)
    eps_t   <- rnorm(test_n, 0, 0.01)
    
    y_test <- beta1 * z1_test +
      f(data1_test$x1, data1_test$x2) +
      eps_t
    
    data1_test_full <- data.frame(
      x1 = data1_test$x1,
      x2 = data1_test$x2,
      x3 = data1_test$x3,
      z1 = z1_test,
      y  = y_test
    )
    
    data2_train <- as.matrix(data1_train[, c("x1", "x2", "x3")])
    data2_test  <- as.matrix(data1_test_full[, c("x1", "x2", "x3")])
    
    obj <- data1_train$y
    n   <- nrow(data1_train)
    
    Wmatrix      <- matrix(data1_train$z1, ncol = 1)
    Wmatrix_test <- matrix(data1_test_full$z1, ncol = 1)
    
    I <- diag(n)
    
    # ---------------------------------------------------
    # 1. Linear model (OLS)
    # ---------------------------------------------------
    lm.result <- lm(y ~ . - 1, data = data1_train)
    lm_coef   <- coef(lm.result)
    
    beta1_hat_lm <- lm_coef["z1"]
    se_beta1_lm  <- summary(lm.result)$coefficients["z1", "Std. Error"]
    
    pred_lm <- predict(lm.result, newdata = data1_test_full)
    MSE_lm  <- mean((data1_test_full$y - pred_lm)^2)
    
    # ---------------------------------------------------
    # 2. NW(Dirichlet) with alpha = 0
    # ---------------------------------------------------
    S <- d.kernel_fast(
      U     = data2_train,
      X     = data2_train,
      alpha = alpha_dir,
      h     = h
    )
    Sij <- S / rowSums(S)
    
    S_test <- d.kernel_fast(
      U     = data2_test,
      X     = data2_train,
      alpha = alpha_dir,
      h     = h
    )
    Sij_test <- S_test / rowSums(S_test)
    
    A_dir    <- I - Sij
    XtAW_dir <- t(Wmatrix) %*% A_dir %*% Wmatrix
    XtAy_dir <- t(Wmatrix) %*% A_dir %*% obj
    
    ipBetahat <- solve(XtAW_dir, XtAy_dir)
    
    ipmuhat_full <- Sij %*% (data1_train$y - Wmatrix %*% ipBetahat)
    
    ip_Y_linear_hat_full <- Wmatrix_test %*% ipBetahat
    ip_Y_nonlin_hat_full <- Sij_test %*% ipmuhat_full
    ipYhat_full          <- ip_Y_linear_hat_full + ip_Y_nonlin_hat_full
    
    MSE_Imp <- mean((data1_test_full$y - ipYhat_full)^2)
    
    residual <- (I - Sij) %*% (data1_train$y - Wmatrix %*% ipBetahat)
    edf      <- sum(diag(Sij))
    sigma2_hat <- sum(residual^2) / (n - edf)
    
    cov_beta    <- sigma2_hat * solve(XtAW_dir)
    ip.se_beta1 <- sqrt(cov_beta[1, 1])
    
    # ---------------------------------------------------
    # 3. Local linear Dirichlet (alpha = 1)
    # ---------------------------------------------------
    alpha_ll <- 1
    
    X_kern_train <- data2_train
    X_loc_train  <- data2_train[, 1:2, drop = FALSE]
    
    X_kern_test <- data2_test
    X_loc_test  <- data2_test[, 1:2, drop = FALSE]
    
    Sij_LL <- build_S_locallinear_dirichlet(
      X_kern = X_kern_train,
      X_loc  = X_loc_train,
      alpha  = alpha_ll,
      h      = h_wls
    )
    
    Sij_test_LL <- build_S_locallinear_dirichlet_UX(
      U_kern = X_kern_test,
      U_loc  = X_loc_test,
      X_kern = X_kern_train,
      X_loc  = X_loc_train,
      alpha  = alpha_ll,
      h      = h_wls
    )
    
    A_ll    <- I - Sij_LL
    XtAW_ll <- t(Wmatrix) %*% A_ll %*% Wmatrix
    XtAy_ll <- t(Wmatrix) %*% A_ll %*% obj
    
    ipBetahat_LL <- solve(XtAW_ll, XtAy_ll)
    
    ipmuhat_full_LL <- Sij_LL %*% (data1_train$y - Wmatrix %*% ipBetahat_LL)
    
    ip_Y_linear_hat_full_LL <- Wmatrix_test %*% ipBetahat_LL
    ip_Y_nonlin_hat_full_LL <- Sij_test_LL %*% ipmuhat_full_LL
    ipYhat_full_LL          <- ip_Y_linear_hat_full_LL + ip_Y_nonlin_hat_full_LL
    
    MSE_Imp_LL <- mean((data1_test_full$y - ipYhat_full_LL)^2)
    
    residual_LL   <- (I - Sij_LL) %*% (data1_train$y - Wmatrix %*% ipBetahat_LL)
    edf_LL        <- sum(diag(Sij_LL))
    sigma2_hat_LL <- sum(residual_LL^2) / (n - edf_LL)
    
    cov_beta_LL    <- sigma2_hat_LL * solve(XtAW_ll)
    ip.se_beta1_LL <- sqrt(cov_beta_LL[1, 1])
    
    # ---------------------------------------------------
    # 4. NW(Lognormal)
    # ---------------------------------------------------
    S_logn <- lognormal_kernel_matrix(data2_train, h_logn)$K
    Sij_logn <- S_logn / rowSums(S_logn)
    
    S_test_logn <- lognormal_kernel_matrix_new(
      data2_train, data2_test, h_logn, ridge = 1e-4
    )
    Sij_test_logn <- S_test_logn / rowSums(S_test_logn)
    
    A_logn    <- I - Sij_logn
    XtAW_logn <- t(Wmatrix) %*% A_logn %*% Wmatrix
    XtAy_logn <- t(Wmatrix) %*% A_logn %*% obj
    
    ipBetahat_logn <- solve(XtAW_logn, XtAy_logn)
    
    ipmuhat_full_logn <- Sij_logn %*% (data1_train$y - Wmatrix %*% ipBetahat_logn)
    
    ip_Y_linear_hat_full_logn <- Wmatrix_test %*% ipBetahat_logn
    ip_Y_nonlin_hat_full_logn <- Sij_test_logn %*% ipmuhat_full_logn
    ipYhat_full_logn          <- ip_Y_linear_hat_full_logn + ip_Y_nonlin_hat_full_logn
    
    MSE_Imp_logn <- mean((data1_test_full$y - ipYhat_full_logn)^2)
    
    residual_logn   <- (I - Sij_logn) %*% (data1_train$y - Wmatrix %*% ipBetahat_logn)
    edf_logn        <- sum(diag(Sij_logn))
    sigma2_hat_logn <- sum(residual_logn^2) / (n - edf_logn)
    
    cov_beta_logn    <- sigma2_hat_logn * solve(XtAW_logn)
    ip.se_beta1_logn <- sqrt(cov_beta_logn[1, 1])
    
    # ---------------------------------------------------
    # Return all results for this replication
    # ---------------------------------------------------
    c(
      ipBetahat,         ip.se_beta1,         MSE_Imp,          # NW(Dirichlet)
      lm_coef,           se_beta1_lm,        MSE_lm,           # Linear model
      ipBetahat_LL,      ip.se_beta1_LL,     MSE_Imp_LL,       # Local linear
      ipBetahat_logn,    ip.se_beta1_logn,   MSE_Imp_logn      # NW(Lognormal)
    )
  }
  
  # -------------------------------------------------------
  # Extract blocks
  # -------------------------------------------------------
  ip.beta.mse      <- result_list[, 1:3]    # NW(Dirichlet)
  lm.beta.mse      <- result_list[, 4:9]    # linear model (x1,x2,x3,z1, se_z1, MSE)
  ip.beta.mse_LL   <- result_list[, 10:12]  # local linear
  ip.beta.mse_logn <- result_list[, 13:15]  # NW(Lognormal)
  
  ip.result        <- round(colMeans(ip.beta.mse), 5)
  lmresult         <- round(colMeans(lm.beta.mse), 5)
  ip.result_LL     <- round(colMeans(ip.beta.mse_LL), 5)
  ip.result_logn   <- round(colMeans(ip.beta.mse_logn, na.rm = TRUE), 5)
  
  result <- rbind(
    c(ip.result[1],       round(sd(ip.beta.mse[, 1]), 5),
      ip.result[2],       ip.result[3],
      round(sd(ip.beta.mse[, 3]), 5)),
    c(ip.result_logn[1],  round(sd(ip.beta.mse_logn[, 1]), 5),
      ip.result_logn[2],  ip.result_logn[3],
      round(sd(ip.beta.mse_logn[, 3], na.rm = TRUE), 5)),
    c(ip.result_LL[1],    round(sd(ip.beta.mse_LL[, 1]), 5),
      ip.result_LL[2],    ip.result_LL[3],
      round(sd(ip.beta.mse_LL[, 3]), 5)),
    c(lmresult[4],        round(sd(lm.beta.mse[, 4]), 5),
      lmresult[5],        lmresult[6],
      round(sd(lm.beta.mse[, 3]), 5)),
    c(beta1,              ".", ".", ".", ".")
  )
  
  rownames(result) <- c(
    "NW(Dirichlet)",
    "NW(Lognormal)",
    "Local linear",
    "Linear model",
    "TRUE"
  )
  colnames(result) <- c(
    "Mean(beta1_hat)",
    "SD(beta1_hat)",
    "Mean(SE(beta1_hat))",
    "ISE",
    "SD(ISE)"
  )
  
  outfile <- file.path(
    out_dir,
    sprintf("m%d_linear_accuracy_boundary.csv", k)
  )
  write.csv(result, outfile)
}

stopCluster(cl)

Sys.time() - start.time

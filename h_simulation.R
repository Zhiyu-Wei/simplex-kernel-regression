library(MCMCpack)
library(doParallel)
library(foreach)

source("https://raw.githubusercontent.com/Zhiyu-Wei/simplex-kernel-regression/main/bandwidth_selection.r")

r1 <- 1
r2 <- 1.5
r3 <- -0.5

m_list <- list(
  function(s1, s2) r1 * s1 + r2 * s2 + r3 * (1 - s1 - s2),
  function(s1, s2) log(1 + s1 + s2),
  function(s1, s2) s1 * s2,
  function(s1, s2) cos(10 * pmin(s1, s2, 1 - s1 - s2))
)

N         <- 100
beta1     <- 2
loc       <- 0.6
pos       <- c(loc, loc, loc)
alpha_vec <- c(0, 0.5, 1)
K         <- length(m_list)
M         <- 1000

set.seed(123)

nCores <- parallel::detectCores() - 1
cl     <- parallel::makeCluster(nCores)
doParallel::registerDoParallel(cl)

# Example: summarize NW-bandwidths only
h_summary_nw <- array(NA_real_,
                      dim = c(length(alpha_vec), K, 4),
                      dimnames = list(
                        paste0("alpha=", alpha_vec),
                        paste0("m", 1:K),
                        c("mean", "median", "IQR", "na_rate")
                      )
)

for (ai in seq_along(alpha_vec)) {
  alp <- alpha_vec[ai]
  
  for (k in seq_len(K)) {
    f_k <- m_list[[k]]
    
    h_vec <- foreach(
      rep = 1:M,
      .combine  = "c",
      .packages = "MCMCpack"
    ) %dopar% {
      X  <- MCMCpack::rdirichlet(N, pos)
      x1 <- X[, 1]
      x2 <- X[, 2]
      z1 <- rnorm(N, 0, 0.01)
      y  <- beta1 * z1 + f_k(x1, x2) + rnorm(N, 0, 0.01)
      
      safe_cv_bw_dirichlet_nw(X, y, alpha = alp, interval = c(1e-4, 1))
    }
    
    is_na   <- is.na(h_vec)
    h_clean <- h_vec[!is_na]
    
    h_summary_nw[ai, k, "mean"]    <- mean(h_clean)
    h_summary_nw[ai, k, "median"]  <- median(h_clean)
    h_summary_nw[ai, k, "IQR"]     <- IQR(h_clean)
    h_summary_nw[ai, k, "na_rate"] <- mean(is_na)
  }
}

parallel::stopCluster(cl)

h_summary_nw

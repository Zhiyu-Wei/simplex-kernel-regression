d.kernel_fast<- function(U, X, alpha, h) {
  U <- as.matrix(U)   # nU × p
  X <- as.matrix(X)   # nX × p
  nU <- nrow(U)
  nX <- nrow(X)
  p  <- ncol(U)
  
  # u/h + alpha (nU × p)
  UA <- U / h + alpha
  Utilde <- UA - 1
  
  # log X  (nX × p)
  logX <- log(X)
  
  # C(u) = lgamma(1/h + p*alpha) - sum_k lgamma(UA_ik)
  const <- lgamma(1/h + p * alpha)
  C_vec <- const - rowSums(lgamma(UA))   # 長度 nU
  
  # B_ij = sum_k (u_k/h + alpha - 1) * log x_{jk}
  #  = Utilde[i, ] · logX[j, ]'
  # → (nU × p) × (p × nX) = (nU × nX)
  B <- Utilde %*% t(logX)
  
  # log K = outer(C_vec, rep(1,nX)) + B
  K_log <- sweep(B, 1, C_vec, "+")
  
  exp(K_log)
}
build_K_dirichlet_fast <- function(X, alpha, h) {
  X <- as.matrix(X)     # n × p
  n <- nrow(X)
  p <- ncol(X)
  
  XA <- X / h + alpha           # n×p
  const <- lgamma(1/h + p*alpha)
  C_vec <- const - rowSums(lgamma(XA))  # n
  
  Xtilde <- XA - 1              # n×p
  logX   <- log(X)              # n×p
  
  B <- Xtilde %*% t(logX)       # n×n, B_ij = sum_k (x_i/h+α-1) log x_jk
  
  K_log <- sweep(B, 1, C_vec, "+")
  exp(K_log)                    # K (n×n)
}
build_K_dirichlet_fast_UX <- function(U, X, alpha, h) {
  U <- as.matrix(U)   # nU × p
  X <- as.matrix(X)   # nX × p
  nU <- nrow(U)
  nX <- nrow(X)
  p  <- ncol(U)
  
  UA <- U / h + alpha           # nU×p
  Utilde <- UA - 1              # nU×p
  logX <- log(X)                # nX×p
  
  const <- lgamma(1/h + p * alpha)
  C_vec <- const - rowSums(lgamma(UA))   # nU
  
  B <- Utilde %*% t(logX)       # nU×nX
  
  K_log <- sweep(B, 1, C_vec, "+")
  exp(K_log)                    # nU×nX
}
build_S_locallinear_dirichlet <- function(X_kern, X_loc, alpha, h) {
  X_kern <- as.matrix(X_kern)   # n × p_kern (e.g., 3)
  X_loc  <- as.matrix(X_loc)    # n × q     (e.g., 2)
  n  <- nrow(X_kern)
  
  # 先一次算好 Dirichlet kernel 矩陣 K(X, X)
  K <- build_K_dirichlet_fast(X_kern, alpha = alpha, h = h)  # n×n
  
  # 準備 output
  S <- matrix(NA_real_, n, n)
  
  # local linear 的設計矩陣部分
  q <- ncol(X_loc)
  ones <- rep(1, n)
  
  for (i in 1:n) {
    # 權重：用第 i 列 kernel 對所有 training 資料
    w <- K[i, ]                     # 長度 n 向量
    W_s_i <- diag(w)                # n×n (可以改成更省記憶體版本，但先寫清楚)
    
    # centered design：對 X_loc 以 X_loc[i,] 為中心
    Xs <- cbind(1, sweep(X_loc, 2, X_loc[i, ], "-"))  # n×(q+1)
    e1 <- c(1, rep(0, q))                             # (q+1)×1
    
    # A = (Xs' W Xs)^{-1}
    A <- solve(t(Xs) %*% W_s_i %*% Xs)
    
    # smoother row: s_i' = e1' A Xs' W
    S[i, ] <- as.numeric(t(e1) %*% A %*% t(Xs) %*% W_s_i)
  }
  
  S
}
build_S_locallinear_dirichlet_UX <- function(U_kern, U_loc,
                                             X_kern, X_loc,
                                             alpha, h) {
  U_kern <- as.matrix(U_kern)   # nU × p_kern
  U_loc  <- as.matrix(U_loc)    # nU × q
  X_kern <- as.matrix(X_kern)   # nX × p_kern
  X_loc  <- as.matrix(X_loc)    # nX × q
  
  nU <- nrow(U_kern)
  nX <- nrow(X_kern)
  q  <- ncol(X_loc)
  
  # 先算 Dirichlet kernel：K(U, X)
  K_UX <- build_K_dirichlet_fast_UX(U_kern, X_kern, alpha = alpha, h = h)  # nU × nX
  
  S_test <- matrix(NA_real_, nU, nX)
  
  for (i in 1:nU) {
    w <- K_UX[i, ]                 # 長度 nX
    W_i <- diag(w)                 # nX×nX
    
    # local linear 的 centered design (以 U_loc[i,] 為中心)
    Xs <- cbind(1, sweep(X_loc, 2, U_loc[i, ], "-"))  # nX×(q+1)
    e1 <- c(1, rep(0, q))                              # (q+1)×1
    
    A <- solve(t(Xs) %*% W_i %*% Xs)
    
    # row i: s_i'(·) = e1' A Xs' W_i
    S_test[i, ] <- as.numeric(t(e1) %*% A %*% t(Xs) %*% W_i)
  }
  
  S_test
}



draw_composition_plot <- function(data, pred, lab.title = "", original.data=NULL,show_boundary = FALSE,
                                  boundary.color="red",
                                  scale_color ="D") {
  # Convert input to matrix and set column names
  datapoints <- as.matrix(data)
  if (!is.null(original.data)) {
    observed.data<- as.matrix(original.data)
  }
  pred <- as.numeric(pred)
  colnames(datapoints) <- colnames(data)
  
  # Create barycentric transformation matrix
  Bmatrix <- matrix(c(0,1, 0.5,0, 0,sqrt(3)/2), 3,2)
  
  # Minimum values for each dimension (used for axis setup)
  Lmin <- min(datapoints[,1])
  Rmin <- min(datapoints[,2])
  Tmin <- min(datapoints[,3])
  sep <- (1 - Lmin - Rmin - Tmin) / 10
  
  # Triangle border in composition space
  bord.matrix <- matrix(c(1-Rmin-Tmin, Lmin, Lmin,
                          Rmin, 1-Lmin-Tmin, Rmin,
                          Tmin, Tmin, 1-Lmin-Rmin), 3, 3)
  
  # Limits based on observed data range (used for boundary lines)
  if (!is.null(original.data)) {
    limit.matrix <- matrix(c(
      max(original.data[,1]), Rmin, 1-max(original.data[,1])-Rmin,
      max(original.data[,1]), 1-max(original.data[,1])-Tmin, Tmin,
      min(original.data[,1]), Rmin, 1-min(original.data[,1])-Rmin,
      min(original.data[,1]), 1-min(original.data[,1])-Tmin, Tmin,
      Lmin, max(original.data[,2]), 1-Lmin-max(original.data[,2]),
      1-max(original.data[,2])-Tmin, max(original.data[,2]), Tmin,
      Lmin, min(original.data[,2]), 1-Lmin-min(original.data[,2]),
      1-min(original.data[,2])-Tmin, min(original.data[,2]), Tmin,
      Lmin, 1-Lmin-max(original.data[,3]), max(original.data[,3]),
      1-Rmin-max(original.data[,3]), Rmin, max(original.data[,3]),
      Lmin, 1-Lmin-min(original.data[,3]), min(original.data[,3]),
      1-Rmin-min(original.data[,3]), Rmin, min(original.data[,3])
    ), 12, 3, byrow = TRUE)
  }
  
  # Project compositions to 2D coordinates
  Bbord <- bord.matrix %*% Bmatrix
  if (!is.null(original.data)) {
    Blimit <- limit.matrix %*% Bmatrix
  }
  Bdata <- datapoints %*% Bmatrix
  
  # Main triangle and expanded outer triangle for plotting
  tar.triangle <- data.frame(
    x = c(Bbord[1,1], Bbord[2,1], Bbord[3,1], Bbord[1,1]),
    y = c(Bbord[1,2], Bbord[2,2], Bbord[3,2], Bbord[1,2])
  )
  exp.triangle <- data.frame(
    x = c(Bbord[1,1]-sep/2*sqrt(3)*2, Bbord[2,1]+sep/2*sqrt(3)*2, Bbord[3,1], Bbord[1,1]-sep/2*sqrt(3)*2),
    y = c(Bbord[1,2]-sep, Bbord[2,2]-sep, Bbord[3,2]+sep*2, Bbord[1,2]-sep)
  )
  
  # Data points projected to 2D
  point_df <- data.frame(x = Bdata[,1], y = Bdata[,2], pred = pred)
  
  # Helper function for axis tick labels
  make_axis_labels <- function(axis, coord_mat, axis_label) {
    data.frame(
      x = coord_mat[,1],
      y = coord_mat[,2],
      label = round(axis_label, 2)
    )
  }
  
  # Generate axis grid positions (L, R, T axes)
  Lpos <- matrix(c(seq(1-Rmin-Tmin, Lmin, by=-sep),
                   rep(Rmin, 11),
                   seq(Tmin, 1-Lmin-Rmin, by=sep)), 11, 3)
  Rpos <- matrix(c(seq(1-Rmin-Tmin, Lmin, by=-sep),
                   seq(Rmin, 1-Lmin-Tmin, by=sep),
                   rep(Tmin, 11)), 11, 3)
  Tpos <- matrix(c(rep(Lmin, 11),
                   seq(1-Lmin-Tmin, Rmin, by=-sep),
                   seq(Tmin, 1-Lmin-Rmin, by=sep)), 11, 3)
  
  B.Lpos <- Lpos %*% Bmatrix
  B.Rpos <- Rpos %*% Bmatrix
  B.Tpos <- Tpos %*% Bmatrix
  
  Lmarkpos <- cbind(B.Lpos[,1] - sep/10, B.Lpos[,2] + sep/10 * sqrt(3))
  Rmarkpos <- cbind(B.Rpos[,1] - sep/10, B.Rpos[,2] - sep/10 * sqrt(3))
  Tmarkpos <- cbind(B.Tpos[,1] + sep/10 * 2, B.Tpos[,2])
  
  # Tick labels for each axis
  Lname <- make_axis_labels("L", Lmarkpos, Lpos[,1])
  Rname <- make_axis_labels("R", Rmarkpos, Rpos[,2])
  Tname <- make_axis_labels("T", Tmarkpos, Tpos[,3])
  
  # Axis names (e.g., sleep / exercise / other)
  axis_labels <- data.frame(
    x = c(Bbord[1,1], Bbord[2,1], Bbord[3,1]),
    y = c(Bbord[1,2], Bbord[2,2], Bbord[3,2]),
    label = colnames(data),
    hjust = c(-1, 1, -1),
    vjust = c(-4, 4, -4),
    angle = c(54, 0, -60)
  )
  
  # Begin plot construction
  p <- ggplot() +
    geom_polygon(data = tar.triangle, aes(x, y), fill = NA, color = "steelblue1", alpha = 0.3) +
    geom_polygon(data = exp.triangle, aes(x, y), fill = NA, color = "white") +
    geom_point(data = point_df, aes(x = x, y = y, color = pred), size = 2, alpha = 0.8) +
    scale_color_viridis_c(option = scale_color) +
    guides(color = guide_colorbar(title = "rate(%)")) +
    theme_minimal() +
    labs(title = lab.title, x = "", y = "") +
    theme(axis.text = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank()) +
    geom_text(data = axis_labels, aes(x = x, y = y, label = label),
              size = 5, hjust = axis_labels$hjust,
              vjust = axis_labels$vjust, angle = axis_labels$angle) +
    geom_segment(aes(x = B.Lpos[,1], y = B.Lpos[,2], xend = Lmarkpos[,1], yend = Lmarkpos[,2]), color = "lightblue") +
    geom_segment(aes(x = B.Rpos[,1], y = B.Rpos[,2], xend = Rmarkpos[,1], yend = Rmarkpos[,2]), color = "lightblue") +
    geom_segment(aes(x = B.Tpos[,1], y = B.Tpos[,2], xend = Tmarkpos[,1], yend = Tmarkpos[,2]), color = "lightblue") +
    geom_text(data = Lname, aes(x, y, label = label), size = 3, hjust = 1, vjust = 0, angle = -60) +
    geom_text(data = Rname, aes(x, y, label = label), size = 3, hjust = 1, vjust = 0, angle = 60) +
    geom_text(data = Tname, aes(x, y, label = label), size = 3, hjust = 0, vjust = 0, angle = 0) +
    geom_segment(aes(x = B.Lpos[2:10,1], y = B.Lpos[2:10,2], xend = B.Rpos[2:10,1], yend = B.Rpos[2:10,2]),
                 color = "lightblue", linetype = 2, alpha = 0.3) +
    geom_segment(aes(x = B.Rpos[2:10,1], y = B.Rpos[2:10,2], xend = B.Tpos[10:2,1], yend = B.Tpos[10:2,2]),
                 color = "lightblue", linetype = 2, alpha = 0.3) +
    geom_segment(aes(x = B.Tpos[2:10,1], y = B.Tpos[2:10,2], xend = B.Lpos[2:10,1], yend = B.Lpos[2:10,2]),
                 color = "lightblue", linetype = 2, alpha = 0.3)
  
  # Optionally show data bounds
  if (!is.null(original.data) && show_boundary) {
    boundary_df <- data.frame(
      x = Blimit[seq(1, 11, 2), 1],
      y = Blimit[seq(1, 11, 2), 2],
      xend = Blimit[seq(2, 12, 2), 1],
      yend = Blimit[seq(2, 12, 2), 2]
    )
    p <- p + geom_segment(data = boundary_df,
                          aes(x = x, y = y, xend = xend, yend = yend),
                          color = boundary.color)
  }
  
  return(p)
}

composition_LL_kernel_smoother=function(data2, obj, resolution, 
                                        h, alp = 0, show_boundary = T, 
                                        original.data = NULL, boundary.color="purple",
                                        scale_color ="D",
                                        lab.title = ""){
  # Minimum values for each dimension (used for axis setup)
  Lmin <- min(data2[,1])
  Rmin <- min(data2[,2])
  Tmin <- min(data2[,3])
  
  
  
  grid_points <- expand.grid(
    x = seq(0, 1, by = resolution),
    y = seq(0, 1, by = resolution)
  )
  grid_points$z <- 1 - grid_points$x - grid_points$y
  
  grid_points <- subset(grid_points,x > Lmin & x < 1-Rmin-Tmin&y > Rmin & y < 1-Lmin-Tmin&z > Tmin & z < 1-Lmin-Rmin)
  
  
  X_kern_train <- as.matrix(data2[, 1:3])   # x1,x2,x3
  X_loc_train  <- as.matrix(data2[, 1:2])   # x1,x2
  # test 部分
  X_kern_test <- as.matrix(grid_points[, 1:3])
  X_loc_test  <- as.matrix(grid_points[, 1:2])
  
  
  Sij <- build_S_locallinear_dirichlet_UX(   # test–train smoother
    U_kern = X_kern_test,
    U_loc  = X_loc_test,
    X_kern = X_kern_train,
    X_loc  = X_loc_train,
    alpha  = alp,
    h      = h
  )
  
  
  output <- Sij%*%obj
  colnames(grid_points) <- colnames(data2)
  draw_composition_plot(grid_points, output, lab.title = lab.title,original.data = original.data ,show_boundary = show_boundary,
                        boundary.color=boundary.color,scale_color = scale_color)
  
}


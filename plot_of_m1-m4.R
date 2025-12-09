# ---------------------------------------------
# Kernel Regression with Simplex Predictors
# Utility functions: Dirichlet kernel
# ternary plot for compositional predictors
# ---------------------------------------------

library(ggplot2)
library(Cairo)

# ---------------------------------------------
# Dirichlet kernel functions
# ---------------------------------------------

# Dirichlet kernel at a single point K_h(u | x)
# x, u: numeric vectors of length p in (0, 1)
# alpha: centering parameter (scalar)
# h: bandwidth
# p: dimension of the simplex (default 3)
dirichlet_kernel <- function(x, u, alpha, h, p = 3) {
  x <- as.numeric(x)
  u <- as.numeric(u)
  
  if (!all(x > 0 & x < 1) || !all(u > 0 & u < 1)) {
    stop("All entries of x and u must lie in (0, 1).")
  }
  
  result <- lgamma(1 / h + p * alpha) -
    sum(lgamma(x / h + alpha)) +
    sum((x / h + alpha - 1) * log(u))
  
  exp(result)
}

# Fast matrix version:
#   U: n_U x p matrix of evaluation points
#   X: n_X x p matrix of design points
# Returns: n_U x n_X matrix with K_h(U_i | X_j)
dirichlet_kernel_matrix <- function(U, X, alpha, h) {
  U <- as.matrix(U)   # n_U x p
  X <- as.matrix(X)   # n_X x p
  nU <- nrow(U)
  nX <- nrow(X)
  p  <- ncol(U)
  
  # U / h + alpha  (n_U x p)
  UA <- U / h + alpha
  Utilde <- UA - 1
  
  # log X (n_X x p)
  logX <- log(X)
  
  # C(u) = lgamma(1/h + p*alpha) - sum_k lgamma(UA_ik)
  const <- lgamma(1 / h + p * alpha)
  C_vec <- const - rowSums(lgamma(UA))   # length n_U
  
  # B_ij = sum_k (u_k/h + alpha - 1) * log x_{jk}
  #     = Utilde[i, ] %*% logX[j, ]'
  # (n_U x p) %*% (p x n_X) = (n_U x n_X)
  B <- Utilde %*% t(logX)
  
  # log K = C(u_i) + B_ij
  K_log <- sweep(B, 1, C_vec, "+")
  
  exp(K_log)
}
# ---------------------------------------------
# Ternary plot for compositional predictors
# ---------------------------------------------

# data: n x 3 matrix/data frame of simplex covariates
# pred: length-n vector (e.g., fitted values)
# lab.title: plot title
# original.data: optional original design to draw data bounds
# show_boundary: whether to draw boundary segments
# boundary.color: color of the boundary lines
# bar.title: title of colorbar (can be expression)
# scale_color: two-color gradient (low, high)
draw_composition_plot <- function(
    data,
    pred,
    lab.title = "Title",
    original.data = NULL,
    show_boundary = FALSE,
    boundary.color = "red",
    bar.title = "Y",
    scale_color = c("yellow", "blue")
) {
  datapoints <- as.matrix(data)
  if (!is.null(original.data)) {
    observed.data <- as.matrix(original.data)
  }
  pred <- as.numeric(pred)
  colnames(datapoints) <- colnames(data)
  
  # Barycentric transformation matrix for ternary plot
  Bmatrix <- matrix(
    c(0, 1,
      0.5, 0,
      0, sqrt(3) / 2),
    nrow = 3, byrow = TRUE
  )
  
  # Minimum values for each component
  Lmin <- min(datapoints[, 1])
  Rmin <- min(datapoints[, 2])
  Tmin <- min(datapoints[, 3])
  sep  <- (1 - Lmin - Rmin - Tmin) / 10
  
  # Triangle border in composition space
  bord.matrix <- matrix(
    c(
      1 - Rmin - Tmin, Lmin,               Lmin,
      Rmin,            1 - Lmin - Tmin,    Rmin,
      Tmin,            Tmin,               1 - Lmin - Rmin
    ),
    nrow = 3, byrow = TRUE
  )
  
  # Data range boundaries (optional)
  if (!is.null(original.data)) {
    limit.matrix <- matrix(
      c(
        max(original.data[, 1]), Rmin,                    1 - max(original.data[, 1]) - Rmin,
        max(original.data[, 1]), 1 - max(original.data[, 1]) - Tmin, Tmin,
        min(original.data[, 1]), Rmin,                    1 - min(original.data[, 1]) - Rmin,
        min(original.data[, 1]), 1 - min(original.data[, 1]) - Tmin, Tmin,
        Lmin,                    max(original.data[, 2]), 1 - Lmin - max(original.data[, 2]),
        1 - max(original.data[, 2]) - Tmin, max(original.data[, 2]), Tmin,
        Lmin,                    min(original.data[, 2]), 1 - Lmin - min(original.data[, 2]),
        1 - min(original.data[, 2]) - Tmin, min(original.data[, 2]), Tmin,
        Lmin,                    1 - Lmin - max(original.data[, 3]), max(original.data[, 3]),
        1 - Rmin - max(original.data[, 3]), Rmin,         max(original.data[, 3]),
        Lmin,                    1 - Lmin - min(original.data[, 3]), min(original.data[, 3]),
        1 - Rmin - min(original.data[, 3]), Rmin,         min(original.data[, 3])
      ),
      nrow = 12, byrow = TRUE
    )
  }
  
  # Map compositions to 2D coordinates
  Bbord <- bord.matrix %*% Bmatrix
  if (!is.null(original.data)) {
    Blimit <- limit.matrix %*% Bmatrix
  }
  Bdata <- datapoints %*% Bmatrix
  
  # Main and expanded triangles for plotting
  tar.triangle <- data.frame(
    x = c(Bbord[1, 1], Bbord[2, 1], Bbord[3, 1], Bbord[1, 1]),
    y = c(Bbord[1, 2], Bbord[2, 2], Bbord[3, 2], Bbord[1, 2])
  )
  exp.triangle <- data.frame(
    x = c(
      Bbord[1, 1] - sep / 2 * sqrt(3) * 2,
      Bbord[2, 1] + sep / 2 * sqrt(3) * 2,
      Bbord[3, 1],
      Bbord[1, 1] - sep / 2 * sqrt(3) * 2
    ),
    y = c(
      Bbord[1, 2] - sep,
      Bbord[2, 2] - sep,
      Bbord[3, 2] + sep * 2,
      Bbord[1, 2] - sep
    )
  )
  
  # Data points (in 2D)
  point_df <- data.frame(
    x    = Bdata[, 1],
    y    = Bdata[, 2],
    pred = pred
  )
  
  # Helper to create axis labels
  make_axis_labels <- function(coord_mat, axis_label) {
    data.frame(
      x     = coord_mat[, 1],
      y     = coord_mat[, 2],
      label = round(axis_label, 2)
    )
  }
  
  # Axis grid positions for the three components
  Lpos <- matrix(
    c(
      seq(1 - Rmin - Tmin, Lmin, by = -sep),
      rep(Rmin, 11),
      seq(Tmin, 1 - Lmin - Rmin, by = sep)
    ),
    nrow = 11
  )
  Rpos <- matrix(
    c(
      seq(1 - Rmin - Tmin, Lmin, by = -sep),
      seq(Rmin, 1 - Lmin - Tmin, by = sep),
      rep(Tmin, 11)
    ),
    nrow = 11
  )
  Tpos <- matrix(
    c(
      rep(Lmin, 11),
      seq(1 - Lmin - Tmin, Rmin, by = -sep),
      seq(Tmin, 1 - Lmin - Rmin, by = sep)
    ),
    nrow = 11
  )
  
  B.Lpos <- Lpos %*% Bmatrix
  B.Rpos <- Rpos %*% Bmatrix
  B.Tpos <- Tpos %*% Bmatrix
  
  Lmarkpos <- cbind(B.Lpos[, 1] - sep / 10, B.Lpos[, 2] + sep / 10 * sqrt(3))
  Rmarkpos <- cbind(B.Rpos[, 1] - sep / 10, B.Rpos[, 2] - sep / 10 * sqrt(3))
  Tmarkpos <- cbind(B.Tpos[, 1] + sep / 10 * 2, B.Tpos[, 2])
  
  Lname <- make_axis_labels(Lmarkpos, Lpos[, 1])
  Rname <- make_axis_labels(Rmarkpos, Rpos[, 2])
  Tname <- make_axis_labels(Tmarkpos, Tpos[, 3])
  
  # Names of the three components (e.g., x1, x2, x3)
  axis_labels <- data.frame(
    x     = c(Bbord[1, 1], Bbord[2, 1], Bbord[3, 1]),
    y     = c(Bbord[1, 2], Bbord[2, 2], Bbord[3, 2]),
    label = colnames(data),
    hjust = c(-1, 1, -1),
    vjust = c(-4, 4, -4),
    angle = c(54, 0, -60)
  )
  
  p <- ggplot() +
    geom_polygon(
      data = tar.triangle,
      aes(x, y),
      fill = NA,
      color = "steelblue1",
      alpha = 0.3
    ) +
    geom_polygon(
      data = exp.triangle,
      aes(x, y),
      fill = NA,
      color = "white"
    ) +
    geom_point(
      data = point_df,
      aes(x = x, y = y, color = pred),
      size = 2,
      alpha = 0.8
    ) +
    scale_color_gradient(low = scale_color[1], high = scale_color[2]) +
    guides(color = guide_colorbar(title = bar.title)) +
    theme_minimal() +
    labs(title = lab.title, x = "", y = "") +
    theme(
      axis.text  = element_blank(),
      axis.ticks = element_blank(),
      panel.grid = element_blank()
    ) +
    geom_text(
      data = axis_labels,
      aes(x = x, y = y, label = label),
      size  = 5,
      hjust = axis_labels$hjust,
      vjust = axis_labels$vjust,
      angle = axis_labels$angle
    ) +
    geom_segment(
      aes(
        x = B.Lpos[, 1], y = B.Lpos[, 2],
        xend = Lmarkpos[, 1], yend = Lmarkpos[, 2]
      ),
      color = "lightblue"
    ) +
    geom_segment(
      aes(
        x = B.Rpos[, 1], y = B.Rpos[, 2],
        xend = Rmarkpos[, 1], yend = Rmarkpos[, 2]
      ),
      color = "lightblue"
    ) +
    geom_segment(
      aes(
        x = B.Tpos[, 1], y = B.Tpos[, 2],
        xend = Tmarkpos[, 1], yend = Tmarkpos[, 2]
      ),
      color = "lightblue"
    ) +
    geom_text(
      data = Lname,
      aes(x, y, label = label),
      size = 3, hjust = 1, vjust = 0, angle = -60
    ) +
    geom_text(
      data = Rname,
      aes(x, y, label = label),
      size = 3, hjust = 1, vjust = 0, angle = 60
    ) +
    geom_text(
      data = Tname,
      aes(x, y, label = label),
      size = 3, hjust = 0, vjust = 0, angle = 0
    ) +
    geom_segment(
      aes(
        x = B.Lpos[2:10, 1], y = B.Lpos[2:10, 2],
        xend = B.Rpos[2:10, 1], yend = B.Rpos[2:10, 2]
      ),
      color = "lightblue",
      linetype = 2,
      alpha = 0.3
    ) +
    geom_segment(
      aes(
        x = B.Rpos[2:10, 1], y = B.Rpos[2:10, 2],
        xend = B.Tpos[10:2, 1], yend = B.Tpos[10:2, 2]
      ),
      color = "lightblue",
      linetype = 2,
      alpha = 0.3
    ) +
    geom_segment(
      aes(
        x = B.Tpos[2:10, 1], y = B.Tpos[2:10, 2],
        xend = B.Lpos[2:10, 1], yend = B.Lpos[2:10, 2]
      ),
      color = "lightblue",
      linetype = 2,
      alpha = 0.3
    )
  
  # Optional boundary segments based on original data
  if (!is.null(original.data) && show_boundary) {
    boundary_df <- data.frame(
      x    = Blimit[seq(1, 11, 2), 1],
      y    = Blimit[seq(1, 11, 2), 2],
      xend = Blimit[seq(2, 12, 2), 1],
      yend = Blimit[seq(2, 12, 2), 2]
    )
    p <- p +
      geom_segment(
        data = boundary_df,
        aes(x = x, y = y, xend = xend, yend = yend),
        color = boundary.color
      )
  }
  
  p
}

# ---------------------------------------------
# Simplex grid generator (not used below, but kept as utility)
# ---------------------------------------------

# Generate a regular grid on the 2-simplex with step 1/k
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

# ---------------------------------------------
# Example: data generation on a simplex grid
# ---------------------------------------------

resolution <- 0.01

grid_points <- expand.grid(
  x = seq(0, 1, by = resolution),
  y = seq(0, 1, by = resolution)
)
grid_points$z <- 1 - grid_points$x - grid_points$y

# Keep interior points
grid_points_1 <- subset(
  grid_points,
  x > 0 & x < 1 & y > 0 & y < 1 & z > 0 & z < 1
)

# Optionally add boundary points on one edge
add_x <- seq(0, 1, by = resolution)
add_y <- 1 - add_x
add_z <- rep(0, length(add_x))
add_matrix <- cbind(x = add_x, y = add_y, z = add_z)

grid_points_all <- rbind(grid_points_1, as.data.frame(add_matrix))
grid_points_all <- subset(
  grid_points_all,
  x > 0 & x < 1 & y > 0 & y < 1 & z > 0 & z < 1
)

# Linear covariate and noise
set.seed(123)
X1    <- rnorm(nrow(grid_points_all), mean = 1, sd = 0.01)
beta1 <- 2
error <- rnorm(nrow(grid_points_all), mean = 0, sd = 0.01)

# Example regression functions m_1, ..., m_4
Y1 <- 1  * grid_points_all$x +
  1.5 * grid_points_all$y -
  0.5 * grid_points_all$z +
  beta1 * X1 + error

Y2 <- log(1 + grid_points_all$x + grid_points_all$y) +
  beta1 * X1 + error

Y3 <- grid_points_all$x * grid_points_all$y +
  beta1 * X1 + error

Y4 <- cos(6 * pmin(
  grid_points_all$x,
  grid_points_all$y,
  1 - grid_points_all$x - grid_points_all$y
)) +
  beta1 * X1 + error

colnames(grid_points_all) <- c("x1", "x2", "x3")

# NOTE: If you also have Y5–Y8, define them here before binding.
simu_data_all <- cbind(
  grid_points_all,
  Y1 = Y1,
  Y2 = Y2,
  Y3 = Y3,
  Y4 = Y4
)

# ---------------------------------------------
# Save ternary plots of the regression functions
# ---------------------------------------------

# Use a relative path for GitHub (e.g., "figures/m1.png")
dir.create("figures", showWarnings = FALSE)

for (i in 1:4) {
  CairoPNG(
    file   = sprintf("figures/m%d.png", i),
    width  = 2140,
    height = 1605,
    dpi    = 300
  )
  
  print(
    draw_composition_plot(
      simu_data_all[, 1:3],
      simu_data_all[, i + 3],
      lab.title = "",
      bar.title = expression(m(s))
    )
  )
  
  dev.off()
}

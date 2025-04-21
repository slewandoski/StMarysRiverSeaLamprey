library(sf)
library(sfnetworks)
library(Matrix)
library(RTMB)
library(ggplot2)

#-------------------------------------------------------------------------------
# 1. read and transform sampling data
#-------------------------------------------------------------------------------
data <- read.csv("all_data_through_2024/2024-10-24-CatchALL.csv")
data <- st_as_sf(data, coords = c("longitude", "latitude"), crs = 4326)
data <- st_transform(data, crs = 26916)
data <- subset(data, year >= 2000)
data$n <- data$sl.larv.n

# read and transform shapefiles
regions <- st_read("region_polys/Regions_2013.shp")
regions <- st_transform(regions, crs = 26916)

shore <- st_read("AOC_MI_StMarys_2022/AOC_MI_StMarys_2022.shp")
shore <- st_transform(shore, crs = 26916)
shore_region <- shore[regions, ]

#-------------------------------------------------------------------------------
# 2. build grid over river polygon
#-------------------------------------------------------------------------------
grid <- st_as_sf(st_make_grid(shore_region, cellsize = 2000,
                              what = "polygons", square = FALSE))
grid_in_river <- grid[shore_region, ]
grid_clipped <- st_intersection(grid_in_river, shore_region)
centroids <- st_centroid(grid_clipped)

# plot grid and centroids for checking
ggplot() +
  geom_sf(data = grid_clipped, fill = "white", color = "gray") +
  geom_sf(data = centroids, color = "darkorchid2", size = 1.5) +
  geom_sf(data = data, aes(size = n), alpha = 0.6) +
  labs(title = "grid over st. marys river with centroids and observations") +
  theme_minimal()

#-------------------------------------------------------------------------------
# 3. build adjacency-based precision matrix Q
#-------------------------------------------------------------------------------
adj <- st_intersects(grid_clipped, remove_self = TRUE)
n <- length(adj)
Q <- Matrix(0, n, n, sparse = TRUE)
for (i in 1:n) {
  Q[i, i] <- length(adj[[i]])
  Q[i, adj[[i]]] <- -1
}

# quick diagnostic: show sparsity pattern
image(Q, main = "sparsity of Q (neighborhood structure)")

#-------------------------------------------------------------------------------
# 4. manually build projection matrix A_is using nearest centroids
#-------------------------------------------------------------------------------
nearest <- st_nearest_feature(data, centroids)
A_is <- Matrix(0, nrow(data), length(adj), sparse = TRUE)
for (i in seq_along(nearest)) A_is[i, nearest[i]] <- 1

stopifnot(ncol(A_is) == nrow(Q))
stopifnot(nrow(A_is) == nrow(data))

#-------------------------------------------------------------------------------
# 5. prepare RTMB inputs
#-------------------------------------------------------------------------------
year_set <- min(data$year):max(data$year)
data_list <- list(
  "c_i" = data$n,
  "n_t" = length(year_set),
  "n_i" = nrow(data),
  "t_i" = data$year - min(data$year) + 1,
  "A_is" = A_is,
  "Q" = Q
)

par <- list(
  "beta0" = 0,
  "ln_tauO" = 0,
  "ln_tauE" = 0,
  "ln_kappa" = 0,
  "trans_rho" = 0,
  "omega_s" = numeric(ncol(A_is)),
  "epsilon_st" = matrix(0.0, ncol(A_is), length(year_set))
)

to_cor <- function(x) {
  2 / (1 + exp(-2 * x)) - 1
}

#-------------------------------------------------------------------------------
# 6. objective function
#-------------------------------------------------------------------------------
f <- function(par) {
  getAll(data_list, par, warn = FALSE)
  c_i <- OBS(c_i)
  rho <- to_cor(trans_rho)
  Q_spde <- exp(4 * ln_kappa) * Diagonal(ncol(A_is)) + Q
  jnll <- 0
  jnll <- jnll - dgmrf(omega_s, 0.0, Q_spde, TRUE, scale = 1 / exp(ln_tauO))
  jnll <- jnll - dgmrf(epsilon_st[, 1],
                       mu = 0, Q = Q_spde,
                       log = TRUE,
                       scale = 1 / exp(ln_tauE) / sqrt(1 - rho^2)
  )
  for (t in 2:n_t) {
    jnll <- jnll - dgmrf(epsilon_st[, t],
                         mu = rho * epsilon_st[, t - 1],
                         Q = Q_spde, log = TRUE, scale = 1 / exp(ln_tauE)
    )
  }
  omega_i <- A_is %*% omega_s
  epsilon_it <- A_is %*% epsilon_st
  ln_pred_i <- numeric(length(c_i))
  for (i in seq_along(c_i)) {
    ln_pred_i[i] <- beta0 + omega_i[i] + epsilon_it[i, t_i[i]]
  }
  jnll <- jnll - sum(dpois(c_i, exp(ln_pred_i), TRUE))
  range <- sqrt(8) / exp(ln_kappa)
  ADREPORT(range)
  jnll
}

#-------------------------------------------------------------------------------
# 7. fit model
#-------------------------------------------------------------------------------
obj <- MakeADFun(f, par, random = c("omega_s", "epsilon_st"))
obj$fn()
obj$gr()
opt <- nlminb(obj$par, obj$fn, obj$gr)
sdr <- sdreport(obj, bias.correct = TRUE)
sdr

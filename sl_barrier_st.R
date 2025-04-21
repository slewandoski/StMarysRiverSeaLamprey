#-------------------------------------------------------------------------------
# st. mary's larval lamprey in space-time
# cahill + lewandoski
#-------------------------------------------------------------------------------

library(sf)
library(fmesher)
library(INLAspacetime)
library(RTMB)

# read csv and convert to sf with lat/lon in degrees
data <- read.csv("all_data_through_2024/2024-10-24-CatchALL.csv")
data <- st_as_sf(data, coords = c("longitude", "latitude"), crs = 4326)

# transform to utm zone 16n (nad83)
data <- st_transform(data, crs = 26916)

# keep data from 2000 onward
data <- subset(data, year >= 2000)

# rename n for simplicity
data$n <- data$sl.larv.n

#-------------------------------------------------------------------------------
# set up spde approximation for barrier model
#-------------------------------------------------------------------------------

# load land shapefiles
regions <- st_read("region_polys/Regions_2013.shp")
shoreline <- st_read("AOC_MI_StMarys_2022/AOC_MI_StMarys_2022.shp")
shoreline_f <- st_transform(shoreline, st_crs(regions))

# create coarse water polygon
bbox_sf <- st_as_sfc(st_bbox(shoreline_f))
st_crs(bbox_sf) <- st_crs(shoreline_f)
bbox_sf_clipped <- st_difference(bbox_sf, st_union(shoreline_f))

# limit domain to detection buffer
obs_buffer <- st_buffer(st_union(data), dist = 50)
obs_buffer <- st_transform(obs_buffer, st_crs(bbox_sf_clipped))
water_clip <- st_intersection(bbox_sf_clipped, obs_buffer)
land_clip <- st_intersection(shoreline_f, obs_buffer)

# combine water and land for barrier domain
water <- st_sf(geometry = water_clip)
land <- st_sf(geometry = st_geometry(land_clip))
water$barrier <- FALSE
land$barrier <- TRUE
domain <- rbind(water["barrier"], land["barrier"])

# build mesh
mesh <- fm_mesh_2d(
  boundary = domain,
  max.edge = c(8000, 70000),
  cutoff = 250
)

# identify land triangles
tris <- mesh$graph$tv
tl <- nrow(tris)
tri_centroids <- matrix(NA, nrow = tl, ncol = 2)
for (i in seq_len(tl)) {
  tri_centroids[i, ] <- colMeans(mesh$loc[tris[i, ], 1:2, drop = FALSE])
}
tri_sf <- st_as_sf(as.data.frame(tri_centroids),
  coords = c(1, 2), crs = st_crs(land_clip)
)
land_simple <- st_simplify(land, dTolerance = 100)
tri_on_land <- lengths(st_intersects(tri_sf, land_simple)) > 0

# construct spde precision components
spde <- INLAspacetime::mesh2fem.barrier(
  mesh,
  barrier.triangles = which(tri_on_land)
)
M0 <- as(Diagonal(x = as.numeric(spde$C[[1]])), "CsparseMatrix")
M1 <- spde$D[[1]]
M2 <- spde$D[[2]]

# build A matrix for observed points
coords_obs <- st_coordinates(data)
A_is <- fm_evaluator(mesh, loc = coords_obs)$proj$A

# create prediction grid over water
tiles <- st_make_grid(water, cellsize = 50, what = "polygons")
tiles <- st_as_sf(tiles)
tiles <- st_intersection(tiles, water)
centroids <- st_centroid(tiles)
coords_pred <- st_coordinates(centroids)
A_pred <- fm_evaluator(mesh, loc = coords_pred)$proj$A
areas <- as.numeric(st_area(tiles))

# plot the water background
plot(st_geometry(water), col = "blue", border = NA, main = "Prediction Grid Locations")

# overlay prediction points (centroids)
plot(st_geometry(centroids), add = TRUE, pch = 16, col = "darkorchid2", cex = 0.6)

# plot mesh points as reference
points(mesh$loc[,1], mesh$loc[,2], pch = ".", col = "gray")

#-------------------------------------------------------------------------------
# visualize mesh with land barrier
#-------------------------------------------------------------------------------

# plot mesh with domain
plot(mesh, asp = 1, main = "mesh with domain boundary")
plot(st_geometry(domain), add = TRUE, border = "black", lwd = 2)

# plot triangle centroids with land overlay
plot(st_geometry(land_simple),
  col = "grey80",
  border = NA, main = "land overlap"
)
plot(tri_sf,
  col = ifelse(tri_on_land, "red", "blue"),
  pch = 16, cex = 0.5, add = TRUE
)

# plot observation locations with mesh
plot(mesh, asp = 1, main = "mesh and observation locations")
points(coords_obs, pch = 16, cex = 0.4, col = "darkorchid4")

#-------------------------------------------------------------------------------
# fit spatial-temporal Poisson GLMM using SPDE approximation via RTMB
#-------------------------------------------------------------------------------

year_set <- min(data$year):max(data$year)
data <- list(
  "c_i" = data$n,
  "n_t" = length(year_set),
  "n_i" = nrow(A_is),
  "t_i" = data$year - min(data$year) + 1,
  "A_is" = A_is,
  "A_pred" = A_pred,
  "area_pred" = areas,
  "M0" = M0, "M1" = M1, "M2" = M2
)

par <- list(
  "beta0" = 0,
  "ln_tauO" = 0,
  "ln_tauE" = 0,
  "ln_kappa" = 0,
  "trans_rho" = 0,
  "omega_s" = numeric(nrow(M0)),
  "epsilon_st" = matrix(0.0, nrow(M0), data$n_t)
)

to_cor <- function(x) {
  2 / (1 + exp(-2 * x)) - 1
}

f <- function(par) {
  getAll(data, par, warn = FALSE)
  c_i <- OBS(c_i)
  rho <- to_cor(trans_rho)
  Q <- exp(4 * ln_kappa) * M0 + 2 * exp(2 * ln_kappa) * M1 + M2
  jnll <- 0
  # spatial random effect:
  jnll <- jnll - dgmrf(omega_s, 0.0, Q, TRUE, scale = 1 / exp(ln_tauO))
  # temporally evolving AR-1 st field
  # initialize first time slice
  jnll <- jnll - dgmrf(epsilon_st[, 1],
    mu = 0, Q = Q,
    log = TRUE,
    scale = 1 / exp(ln_tauE) / sqrt(1 - rho^2)
  )
  # auto-regress remaining slices
  for (t in 2:n_t) {
    jnll <- jnll - dgmrf(epsilon_st[, t],
      mu = rho * epsilon_st[, t - 1],
      Q = Q, log = TRUE, scale = 1 / exp(ln_tauE)
    )
  }
  omega_i <- A_is %*% omega_s
  epsilon_it <- A_is %*% epsilon_st
  ln_pred_i <- numeric(length(c_i))
  for (i in 1:length(c_i)) {
    ln_pred_i[i] <- beta0 + omega_i[i] + epsilon_it[i, t_i[i]]
  }
  jnll <- jnll - sum(dpois(c_i, exp(ln_pred_i), TRUE))
  range <- sqrt(8) / exp(ln_kappa)
  omega_pred <- A_pred %*% omega_s
  epsilon_pred <- A_pred %*% epsilon_st
  N <- numeric(n_t)
  for (t in 1:n_t) {
    log_lambda <- beta0 + omega_pred + epsilon_pred[, t]
    lambda_pred <- exp(log_lambda)
    N[t] <- sum(lambda_pred * area_pred)
  }
  ADREPORT(range)
  ADREPORT(N)
  jnll
}

obj <- MakeADFun(f, par, random = c("omega_s", "epsilon_st"))
opt <- nlminb(obj$par, obj$fn, obj$gr)
sdr <- sdreport(obj, bias.correct = TRUE)

# Extract values
i <- which(names(sdr$value) == "N")
N <- sdr$value[i]
SE <- sdr$sd[i]

par(mfrow=c(2,1))
# Plot points
plot(N ~ year_set,
  type = "b", pch = 16, cex = 0.5,
  ylim = range(N - 1.96 * SE, N + 1.96 * SE),
  col = "darkorchid4", xlab = "Year",
  ylab = "Larval abundance", 
  main = "Poisson"
)

# CI lines aligned with year_set
segments(
  x0 = year_set, y0 = N - 1.96 * SE,
  x1 = year_set, y1 = N + 1.96 * SE,
  col = "darkorchid4"
)

#-------------------------------------------------------------------------------
# fit spatial-temporal Negative Binomial GLMM using SPDE approximation via RTMB
#-------------------------------------------------------------------------------

par <- list(
  "beta0" = 0,
  "ln_tauO" = 0,
  "ln_tauE" = 0,
  "ln_kappa" = 0,
  "trans_rho" = 0,
  "ln_theta" = 0,
  "omega_s" = numeric(nrow(M0)),
  "epsilon_st" = matrix(0.0, nrow(M0), data$n_t)
)

f <- function(par) {
  getAll(data, par, warn = FALSE)
  c_i <- OBS(c_i)
  rho <- to_cor(trans_rho)
  Q <- exp(4 * ln_kappa) * M0 + 2 * exp(2 * ln_kappa) * M1 + M2
  jnll <- 0
  # spatial random effect:
  jnll <- jnll - dgmrf(omega_s, 0.0, Q, TRUE, scale = 1 / exp(ln_tauO))
  # temporally evolving AR-1 st field
  # initialize first time slice
  jnll <- jnll - dgmrf(epsilon_st[, 1],
    mu = 0, Q = Q,
    log = TRUE,
    scale = 1 / exp(ln_tauE) / sqrt(1 - rho^2)
  )
  # auto-regress remaining slices
  for (t in 2:n_t) {
    jnll <- jnll - dgmrf(epsilon_st[, t],
      mu = rho * epsilon_st[, t - 1],
      Q = Q, log = TRUE, scale = 1 / exp(ln_tauE)
    )
  }
  omega_i <- A_is %*% omega_s
  epsilon_it <- A_is %*% epsilon_st
  ln_pred_i <- numeric(length(c_i))
  for (i in 1:length(c_i)) {
    ln_pred_i[i] <- beta0 + omega_i[i] + epsilon_it[i, t_i[i]]
  }
  mu_i <- exp(ln_pred_i)
  theta <- exp(ln_theta)
  prob_i <- theta / (theta + mu_i)
  jnll <- jnll - sum(dnbinom(c_i, size = theta, prob = prob_i, log = TRUE))

  range <- sqrt(8) / exp(ln_kappa)
  omega_pred <- A_pred %*% omega_s
  epsilon_pred <- A_pred %*% epsilon_st
  N <- numeric(n_t)
  for (t in 1:n_t) {
    log_lambda <- beta0 + omega_pred + epsilon_pred[, t]
    lambda_pred <- exp(log_lambda)
    N[t] <- sum(lambda_pred * area_pred)
  }
  ADREPORT(range)
  ADREPORT(N)
  jnll
}

obj <- MakeADFun(f, par, random = c("omega_s", "epsilon_st"))
opt <- nlminb(obj$par, obj$fn, obj$gr)
sdr <- sdreport(obj, bias.correct = TRUE)

# Extract values
i <- which(names(sdr$value) == "N")
N <- sdr$value[i]
SE <- sdr$sd[i]

# Plot points
plot(N ~ year_set,
  type = "b", pch = 16, cex = 0.5,
  ylim = range(N - 1.96 * SE, N + 1.96 * SE),
  col = "darkorchid4", xlab = "Year",
  ylab = "Larval abundance", 
  main = "Negative binomial"
)

# CI lines aligned with year_set
segments(
  x0 = year_set, y0 = N - 1.96 * SE,
  x1 = year_set, y1 = N + 1.96 * SE,
  col = "darkorchid4"
)

#NOT FINISHED#

library(sf)
library(fmesher)
library(INLAspacetime)
library(RTMB)

# Step 1: Load and project sample data
data <- read.csv("all_data_through_2024/2024-10-24-CatchALL.csv")
data <- st_as_sf(data, coords = c("longitude", "latitude"), crs = 4326)
data <- st_transform(data, crs = 26916)
data$n <- data$sl.larv.n
data <- subset(data, year == 2024)

# Step 2: Load land shapefiles
regions <- st_read('region_polys/Regions_2013.shp')
shoreline <- st_read('AOC_MI_StMarys_2022/AOC_MI_StMarys_2022.shp')
shoreline_f <- st_transform(shoreline, st_crs(regions))

# Step 3: Create coarse water polygon
bbox_sf <- st_as_sfc(st_bbox(shoreline_f))
st_crs(bbox_sf) <- st_crs(shoreline_f)
bbox_sf_clipped <- st_difference(bbox_sf, st_union(shoreline_f))

# Step 4: Limit domain to detection buffer
obs_buffer <- st_buffer(st_union(data), dist = 50)
obs_buffer <- st_transform(obs_buffer, st_crs(bbox_sf_clipped))
water_clip <- st_intersection(bbox_sf_clipped, obs_buffer)
land_clip <- st_intersection(shoreline_f, obs_buffer)

# Step 5: Combine water and land for barrier domain
water <- st_sf(geometry = water_clip)
land <- st_sf(geometry = st_geometry(land_clip))
water$barrier <- FALSE
land$barrier <- TRUE
domain <- rbind(water["barrier"], land["barrier"])

# Step 6: Build mesh
mesh <- fm_mesh_2d(
  boundary = domain,
  max.edge = c(8000, 70000),
  cutoff = 250
)

# Step 7: Identify land triangles
tris <- mesh$graph$tv
tl <- nrow(tris)
tri_centroids <- matrix(NA, nrow = tl, ncol = 2)
for (i in seq_len(tl)) {
  tri_centroids[i, ] <- colMeans(mesh$loc[tris[i, ], 1:2, drop = FALSE])
}
tri_sf <- st_as_sf(as.data.frame(tri_centroids), coords = c(1, 2), crs = st_crs(land_clip))
land_simple <- st_simplify(land, dTolerance = 100)
tri_on_land <- lengths(st_intersects(tri_sf, land_simple)) > 0

# Step 8: Construct SPDE precision components
spde <- INLAspacetime::mesh2fem.barrier(mesh, barrier.triangles = which(tri_on_land))
M0 <- as(Diagonal(x = as.numeric(spde$C[[1]])), "CsparseMatrix")
M1 <- spde$D[[1]]
M2 <- spde$D[[2]]

# Step 9: Build A matrix for observed points
coords_obs <- st_coordinates(data)
A_is <- fm_evaluator(mesh, loc = coords_obs)$proj$A

# Step 10: Create prediction grid over water
tiles <- st_make_grid(water, cellsize = 50, what = "polygons")
tiles <- st_as_sf(tiles)
tiles <- st_intersection(tiles, water)
centroids <- st_centroid(tiles)
coords_pred <- st_coordinates(centroids)
A_pred <- fm_evaluator(mesh, loc = coords_pred)$proj$A
areas <- as.numeric(st_area(tiles))

# Step 11: Fit barrier SPDE model via RTMB
data <- list(
  c_i = data$n,
  A_is = A_is,
  A_pred = A_pred,
  M0 = M0,
  M1 = M1,
  M2 = M2,
  area_pred = areas
)

par <- list(
  beta0 = 0,
  ln_tau = 0,
  ln_kappa = 0,
  omega_s = numeric(nrow(M0))
)

f <- function(par) {
  getAll(data, par, warn = FALSE)
  c_i <- OBS(c_i)
  sigE <- 1 / sqrt(4 * pi * exp(2 * ln_tau) * exp(2 * ln_kappa))
  Q <- exp(4 * ln_kappa) * M0 + 2 * exp(2 * ln_kappa) * M1 + M2
  jnll <- 0
  jnll <- jnll - dgmrf(omega_s, 0.0, Q, TRUE, scale = 1 / exp(ln_tau))
  omega_i <- A_is %*% omega_s
  jnll <- jnll - sum(dpois(c_i, exp(beta0 + omega_i), TRUE))
  
  # Derived quantities
  range <- sqrt(8) / exp(ln_kappa)
  omega_pred <- A_pred %*% omega_s
  log_lambda <- beta0 + omega_pred
  lambda_pred <- exp(log_lambda)
  N <- sum(lambda_pred * area_pred)
  
  ADREPORT(N)
  ADREPORT(range)
  REPORT(sigE)
  REPORT(lambda_pred)
  jnll
}

obj <- MakeADFun(f, par, random = "omega_s")
opt <- nlminb(obj$par, obj$fn, obj$gr)
sdr <- sdreport(obj, bias.correct = TRUE)

sdr$value
sdr$sd

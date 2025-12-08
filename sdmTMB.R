# fit some index standardization models using sdmTMB
#
# see https://sdmtmb.github.io/sdmTMB/
#
# note, this takes about 66 minutes to run on my desktop
#
# todo:
# - compare abundance predictions among models

library(sdmTMB)
library(sf)
library(ggplot2)
library(raster)
library(patchwork)

# read CSV and convert to sf with lat/lon in degrees
data <- read.csv("all_data_through_2024/2024-10-24-CatchALL.csv")
data <- data[which(data$hab.type != 3), ]
data <- st_as_sf(data, coords = c("longitude", "latitude"), crs = 4326) # WGS84

# transform to UTM Zone 16N (NAD83)
data <- st_transform(data, crs = 26916) # EPSG:26916

# extract UTM coordinates in km
coords <- st_coordinates(data) / 1000
data$X <- coords[, 1]
data$Y <- coords[, 2]
data <- sf::st_drop_geometry(data)

# rename for simplicity
data$c <- data$sl.larv.n
data$hab_type <- data$hab.type   

# plot non-zero counts in color, zero in grey
ggplot() +
  # zeros first, as faint background
  geom_point(
    data = subset(data, c == 0),
    aes(X, Y),
    colour = "grey90",
    alpha = 0.95,
    size = 0.6
  ) +
  # non-zeros on top, with colour scale
  geom_point(
    data = subset(data, c > 0),
    aes(X, Y, colour = c),
    size = 1
  ) +
  coord_fixed() +
  scale_colour_viridis_c(trans = "asinh") +
  labs(colour = "count", title = "Surveys with counts > 0")

mesh <- make_mesh(data, xy_cols = c("X", "Y"), cutoff = 0.25)
# mesh$mesh$n
# plot(mesh)

# no space or year effect
m1 <- sdmTMB(
  c ~ as.factor(hab_type) + s(depth, k = 5), # smoother on depth
  data = data,
  mesh = mesh,
  family = nbinom2(link = "log"),
  spatial = "off",
  silent = FALSE
)

data$resids <- residuals(m1, type = "mle-mvn") # randomized quantile residuals
qqnorm(data$resids)
abline(0, 1)

# space effect shared among years
m2 <- sdmTMB(
  c ~ 1 + as.factor(hab_type) + s(depth, k = 5),
  data = data,
  mesh = mesh,
  family = nbinom2(link = "log"),
  spatial = "on",
  silent = FALSE
)

data$resids <- residuals(m2, type = "mle-mvn") # randomized quantile residuals
qqnorm(data$resids)
abline(0, 1)

# intercept factor by year, no space or time effects
m3 <- sdmTMB(
  c ~ 0 + as.factor(year) + as.factor(hab_type) + s(depth, k = 5),
  data = data,
  mesh = mesh,
  family = nbinom2(link = "log"),
  spatial = "off",
  silent = FALSE
)

# intercept factor by year, spatial effect but not time
m4 <- sdmTMB(
  c ~ 0 + as.factor(year) + as.factor(hab_type) + s(depth, k = 5),
  time = "year",
  data = data,
  mesh = mesh,
  family = nbinom2(link = "log"),
  spatial = "on",
  silent = FALSE
)

# one intercept all years, but IID spatial-temporal fields
m5 <- sdmTMB(
  c ~ 1 + as.factor(hab_type) + s(depth, k = 5),
  time = "year",
  data = data,
  mesh = mesh,
  family = nbinom2(link = "log"),
  spatial = "on",
  spatiotemporal = "iid",
  silent = FALSE
)

# one intercept all years, but spatiotemporal random walk among years
m6 <- sdmTMB(
  c ~ 1 + as.factor(hab_type) + s(depth, k = 5),
  time = "year",
  data = data,
  mesh = mesh,
  family = nbinom2(link = "log"),
  spatial = "on",
  spatiotemporal = "rw",
  extra_time = c(1997, 1998, 2020), # ensure regular spacing
  silent = FALSE
)

# on intercept all years, spatiotemporal AR1 among years
m7 <- sdmTMB(
  c ~ 1 + as.factor(hab_type) + s(depth, k = 5),
  time = "year",
  data = data,
  mesh = mesh,
  family = nbinom2(link = "log"),
  spatial = "on",
  spatiotemporal = "ar1",
  extra_time = c(1997, 1998, 2020), # ensure regular spacing
  silent = FALSE
)

# intercept factor by year, spatio-temporal AR1 among years
data$year_fac <- as.factor(data$year)
m8 <- sdmTMB(
  c ~ 0 + year_fac + as.factor(hab_type) + s(depth, k = 5),
  time = "year",
  data = data,
  mesh = mesh,
  family = nbinom2(link = "log"),
  spatial = "on",
  spatiotemporal = "ar1",
  extra_time = c(1997, 1998, 2020), # ensure regular spacing
  silent = FALSE
)

print(AIC(m1, m2, m3, m4, m5, m6, m7, m8))

# results
#    df      AIC
# m1  5 30375.93
# m2  7 26822.42
# m3 33 29484.83
# m4 36 26232.10
# m5  8 26292.35
# m6  8 26399.01
# m7  9 26279.90
# m8 37 26229.16

# general plots and diagnostics, can (should) be run for each model
# m4
# tidy(m4, conf.int = TRUE)
# tidy(m4, effects = "ran_pars", conf.int = TRUE)
# sanity(m4)
#  plot depth effect
# ggeffects::ggpredict(m4, terms="depth[0:40, by = 2]") |> plot()
# data$resids <- residuals(m4) # randomized quantile residuals
# hist(data$resids)
# ggplot(data, aes(X, Y, col = resids)) + scale_colour_gradient2() +
#     geom_point() + facet_wrap(~year, nrow = 3) + coord_fixed()
# set.seed(19283)
# s <- simulate(m4, nsim = 1000, type = "mle-mvn")
# dharma_residuals(s, m4)
# abline(0,1)
#
# ggplot(data, aes(X, Y, col = resids)) +
#   scale_colour_gradient2() +
#   geom_point() +
#   facet_wrap(~year, nrow = 4) +
#   coord_fixed()

#########################################################################
#########################################################################
#########################################################################
# All plotting/prediction grid stuff below here
#------------------------------------------------------------------------
# set up a prediction grid -- note this part is interactive to remove
# a portion of the prediction grid
# 
# target_crs <- 26916 # utm zone 16n, nad83
# grid_spacing <- 250 # meters between grid points
# bbox_buffer <- 100 # how far beyond survey extent to keep grid (m)
# crop_pad <- 100 # extra padding around grid+points in plots (m)
# 
# # ---------------------------
# # read shoreline -> utm 26916
# # ---------------------------
# 
# shore_raw <- st_read("AOC_MI_StMarys_2022/AOC_MI_StMarys_2022.shp")
# shore <- st_transform(shore_raw, target_crs)
# shore_union <- st_union(shore)
# 
# # ---------------------------
# # read survey data -> utm 26916 with X,Y
# # ---------------------------
# 
# dat <- read.csv("all_data_through_2024/2024-10-24-CatchALL.csv")
# dat <- dat[dat$hab.type != 3, ]
# 
# dat <- subset(dat, !is.na(longitude) & !is.na(latitude))
# 
# dat_sf <- st_as_sf(dat, coords = c("longitude", "latitude"), crs = 4326)
# dat_sf <- st_transform(dat_sf, target_crs)
# 
# coords <- st_coordinates(dat_sf)
# dat_sf$X <- coords[, 1]
# dat_sf$Y <- coords[, 2]
# dat_df <- st_drop_geometry(dat_sf)
# 
# # ---------------------------
# # bounding box of survey points (utm) with buffer
# # ---------------------------
# 
# bb <- st_bbox(dat_sf)
# 
# xmin_t <- bb["xmin"] - bbox_buffer
# xmax_t <- bb["xmax"] + bbox_buffer
# ymin_t <- bb["ymin"] - bbox_buffer
# ymax_t <- bb["ymax"] + bbox_buffer
# 
# # ---------------------------
# # build prediction grid over shoreline, then trim to survey extent
# # ---------------------------
# 
# grid_all <- st_make_grid(
#   shore_union,
#   cellsize = grid_spacing,
#   what     = "centers"
# )
# 
# grid_sf <- st_as_sf(grid_all)
# 
# inside_shore <- st_within(grid_sf, shore_union, sparse = FALSE)[, 1]
# grid_sf <- grid_sf[inside_shore, ]
# 
# grid_coords <- st_coordinates(grid_sf)
# keep_bbox <- grid_coords[, 1] >= xmin_t &
#   grid_coords[, 1] <= xmax_t &
#   grid_coords[, 2] >= ymin_t &
#   grid_coords[, 2] <= ymax_t
# 
# grid_sf_sub <- grid_sf[keep_bbox, ]
# grid_coords_sub <- grid_coords[keep_bbox, , drop = FALSE]
# 
# pred_grid <- data.frame(
#   x       = grid_coords_sub[, 1],
#   y       = grid_coords_sub[, 2],
#   area_m2 = grid_spacing^2,
#   area_ha = (grid_spacing^2) / 10000
# )
# 
# # ---------------------------
# # common crop limits from grid + points (for original grid)
# # ---------------------------
# 
# x_min <- min(c(pred_grid$x, dat_df$X)) - crop_pad
# x_max <- max(c(pred_grid$x, dat_df$X)) + crop_pad
# y_min <- min(c(pred_grid$y, dat_df$Y)) - crop_pad
# y_max <- max(c(pred_grid$y, dat_df$Y)) + crop_pad
# 
# # ---------------------------
# # p1: original grid; p2: sampling points
# # ---------------------------
# 
# p1 <- ggplot() +
#   geom_sf(data = shore_union, fill = "grey90", colour = NA) +
#   geom_point(
#     data  = pred_grid,
#     aes(x = x, y = y),
#     size  = 0.1,
#     alpha = 0.5
#   ) +
#   coord_sf(
#     xlim = c(x_min, x_max),
#     ylim = c(y_min, y_max),
#     expand = FALSE
#   ) +
#   labs(
#     title = "prediction grid (original)",
#     subtitle = paste0(
#       "grid spacing = ", grid_spacing,
#       " m, bbox buffer = ", bbox_buffer, " m"
#     )
#   )
# 
# p2 <- ggplot() +
#   geom_sf(data = shore_union, fill = "grey90", colour = NA) +
#   geom_point(
#     data = dat_df,
#     aes(x = X, y = Y),
#     size = 0.1
#   ) +
#   coord_sf(
#     xlim = c(x_min, x_max),
#   ) +
#   labs(title = "survey locations")
# 
# # ---------------------------------------------------------------------
# # hacky way to remove that unsampled lake with poly_click
# # ---------------------------------------------------------------------
# 
# # work in lon/lat so axes are intuitive
# shore_ll <- st_read("AOC_MI_StMarys_2022/AOC_MI_StMarys_2022.shp") # original NAD83
# dat_ll <- read.csv("all_data_through_2024/2024-10-24-CatchALL.csv")
# 
# plot(st_geometry(shore_ll), col = "lightblue")
# points(dat_ll$longitude, dat_ll$latitude, pch = 16, cex = 0.5)
# 
# # now click around the region you want to remove, then right-click / esc to stop
# poly_click <- locator(type = "l") # draws as you click
# 
# # combine locator points into matrix (lon, lat)
# bad_coords_ll <- cbind(poly_click$x, poly_click$y)
# 
# # close the polygon by repeating the first point
# bad_coords_ll <- rbind(bad_coords_ll, bad_coords_ll[1, , drop = FALSE])
# 
# bad_poly_ll <- st_sfc(
#   st_polygon(list(bad_coords_ll)),
#   crs = st_crs(shore_ll) # same crs as when you clicked
# )
# 
# # transform bad polygon to utm 26916 (same as pred_grid)
# bad_poly <- st_transform(bad_poly_ll, target_crs)
# 
# # mask pred_grid: drop cells inside bad polygon
# grid_sf_pred <- st_as_sf(pred_grid, coords = c("x", "y"), crs = target_crs)
# 
# inside_bad <- st_within(grid_sf_pred, bad_poly, sparse = FALSE)[, 1]
# 
# grid_sf_keep <- grid_sf_pred[!inside_bad, ]
# grid_coords_keep <- st_coordinates(grid_sf_keep)
# 
# pred_grid_fix <- data.frame(
#   x       = grid_coords_keep[, 1],
#   y       = grid_coords_keep[, 2],
#   area_m2 = grid_spacing^2,
#   area_ha = (grid_spacing^2) / 10000
# )
# 
# # recompute crop limits from fixed grid + points (optional but cleaner)
# x_min_fix <- min(c(pred_grid_fix$x, dat_df$X)) - crop_pad
# x_max_fix <- max(c(pred_grid_fix$x, dat_df$X)) + crop_pad
# y_min_fix <- min(c(pred_grid_fix$y, dat_df$Y)) - crop_pad
# y_max_fix <- max(c(pred_grid_fix$y, dat_df$Y)) + crop_pad
# 
# # ---------------------------
# # p3: fixed prediction grid
# # ---------------------------
# 
# p3 <- ggplot() +
#   geom_sf(data = shore_union, fill = "grey90", colour = NA) +
#   geom_point(
#     data  = pred_grid_fix,
#     aes(x = x, y = y),
#     size  = 0.1,
#     alpha = 0.5
#   ) +
#   coord_sf(
#     xlim = c(x_min_fix, x_max_fix),
#     ylim = c(y_min_fix, y_max_fix),
#     expand = FALSE
#   ) +
#   labs(
#     title = "prediction grid (masked)"
#   )
# 
# # ---------------------------
# # three-panel layout: original grid, surveys, fixed grid
# # ---------------------------
# 
# fig <- p1 + p2 + p3
# ggsave(
#   filename = "prediction_grid_diagnostics.pdf",
#   plot     = fig,
#   width    = 14,
#   height   = 11,
# )
# 
# grid_sf_out <- st_as_sf(
#   pred_grid_fix,
#   coords = c("x", "y"),
#   crs    = target_crs
# )
# 
# st_write(
#   grid_sf_out,
#   "prediction_grid_utm26916_masked.gpkg",
#   delete_dsn = TRUE
# )

library(sf)
library(terra)
library(gstat)
library(ggplot2)
library(patchwork)

#------------------------------------------------------------------------
# set up a prediction grid -- note this part is interactive to remove
# a portion of the prediction grid
#
target_crs <- 26916 # utm zone 16n, nad83
grid_spacing <- 100 # meters between grid points
bbox_buffer <- 100 # how far beyond survey extent to keep grid (m)
crop_pad <- 100 # extra padding around grid+points in plots (m)

# ---------------------------
# read shoreline -> utm 26916
# ---------------------------

shore_raw <- st_read("AOC_MI_StMarys_2022/AOC_MI_StMarys_2022.shp")
shore <- st_transform(shore_raw, target_crs)
shore_union <- st_union(shore)

# ---------------------------
# read survey data -> utm 26916 with X,Y
# ---------------------------

dat <- read.csv("all_data_through_2024/2024-10-24-CatchALL.csv")
dat <- dat[dat$hab.type != 3, ]
dat <- subset(dat, !is.na(longitude) & !is.na(latitude))
dat_sf <- st_as_sf(dat, coords = c("longitude", "latitude"), crs = 4326)
dat_sf <- st_transform(dat_sf, target_crs)

coords <- st_coordinates(dat_sf)
dat_sf$X <- coords[, 1]
dat_sf$Y <- coords[, 2]
dat_df <- st_drop_geometry(dat_sf)

# ---------------------------
# bounding box of survey points (utm) with buffer
# ---------------------------

bb <- st_bbox(dat_sf)

xmin_t <- bb["xmin"] - bbox_buffer
xmax_t <- bb["xmax"] + bbox_buffer
ymin_t <- bb["ymin"] - bbox_buffer
ymax_t <- bb["ymax"] + bbox_buffer

# ---------------------------
# build prediction grid over shoreline, then trim to survey extent
# ---------------------------

grid_all <- st_make_grid(
  shore_union,
  cellsize = grid_spacing,
  what     = "centers"
)

grid_sf <- st_as_sf(grid_all)

inside_shore <- st_within(grid_sf, shore_union, sparse = FALSE)[, 1]
grid_sf <- grid_sf[inside_shore, ]

grid_coords <- st_coordinates(grid_sf)
keep_bbox <- grid_coords[, 1] >= xmin_t &
  grid_coords[, 1] <= xmax_t &
  grid_coords[, 2] >= ymin_t &
  grid_coords[, 2] <= ymax_t

grid_sf_sub <- grid_sf[keep_bbox, ]
grid_coords_sub <- grid_coords[keep_bbox, , drop = FALSE]

pred_grid <- data.frame(
  x       = grid_coords_sub[, 1],
  y       = grid_coords_sub[, 2],
  area_m2 = grid_spacing^2,
  area_ha = (grid_spacing^2) / 10000
)

# ---------------------------
# common crop limits from grid + points (for original grid)
# ---------------------------

x_min <- min(c(pred_grid$x, dat_df$X)) - crop_pad
x_max <- max(c(pred_grid$x, dat_df$X)) + crop_pad
y_min <- min(c(pred_grid$y, dat_df$Y)) - crop_pad
y_max <- max(c(pred_grid$y, dat_df$Y)) + crop_pad

# ---------------------------
# p1: original grid; p2: sampling points
# ---------------------------

p1 <- ggplot() +
  geom_sf(data = shore_union, fill = "grey90", colour = NA) +
  geom_point(
    data  = pred_grid,
    aes(x = x, y = y),
    size  = 0.1,
    alpha = 0.5
  ) +
  coord_sf(
    xlim = c(x_min, x_max),
    ylim = c(y_min, y_max),
    expand = FALSE
  ) +
  labs(
    title = "prediction grid (original)",
    subtitle = paste0(
      "grid spacing = ", grid_spacing,
      " m, bbox buffer = ", bbox_buffer, " m"
    )
  )

p2 <- ggplot() +
  geom_sf(data = shore_union, fill = "grey90", colour = NA) +
  geom_point(
    data = dat_df,
    aes(x = X, y = Y),
    size = 0.1
  ) +
  coord_sf(
    xlim = c(x_min, x_max),
  ) +
  labs(title = "survey locations")

# ---------------------------------------------------------------------
# hacky way to remove that unsampled lake with poly_click
# ---------------------------------------------------------------------

# work in lon/lat so axes are intuitive
shore_ll <- st_read("AOC_MI_StMarys_2022/AOC_MI_StMarys_2022.shp") # original NAD83
dat_ll <- read.csv("all_data_through_2024/2024-10-24-CatchALL.csv")

plot(st_geometry(shore_ll), col = "lightblue")
points(dat_ll$longitude, dat_ll$latitude, pch = 16, cex = 0.5)

# now click around the region you want to remove, then right-click / esc to stop
poly_click <- locator(type = "l") # draws as you click

# combine locator points into matrix (lon, lat)
bad_coords_ll <- cbind(poly_click$x, poly_click$y)

# close the polygon by repeating the first point
bad_coords_ll <- rbind(bad_coords_ll, bad_coords_ll[1, , drop = FALSE])

bad_poly_ll <- st_sfc(
  st_polygon(list(bad_coords_ll)),
  crs = st_crs(shore_ll) # same crs as when you clicked
)

# transform bad polygon to utm 26916 (same as pred_grid)
bad_poly <- st_transform(bad_poly_ll, target_crs)

# mask pred_grid: drop cells inside bad polygon
grid_sf_pred <- st_as_sf(pred_grid, coords = c("x", "y"), crs = target_crs)

inside_bad <- st_within(grid_sf_pred, bad_poly, sparse = FALSE)[, 1]

grid_sf_keep <- grid_sf_pred[!inside_bad, ]
grid_coords_keep <- st_coordinates(grid_sf_keep)

pred_grid_fix <- data.frame(
  x       = grid_coords_keep[, 1],
  y       = grid_coords_keep[, 2],
  area_m2 = grid_spacing^2,
  area_ha = (grid_spacing^2) / 10000
)

# recompute crop limits from fixed grid + points (optional but cleaner)
x_min_fix <- min(c(pred_grid_fix$x, dat_df$X)) - crop_pad
x_max_fix <- max(c(pred_grid_fix$x, dat_df$X)) + crop_pad
y_min_fix <- min(c(pred_grid_fix$y, dat_df$Y)) - crop_pad
y_max_fix <- max(c(pred_grid_fix$y, dat_df$Y)) + crop_pad

# ---------------------------
# p3: fixed prediction grid
# ---------------------------

p3 <- ggplot() +
  geom_sf(data = shore_union, fill = "grey90", colour = NA) +
  geom_point(
    data  = pred_grid_fix,
    aes(x = x, y = y),
    size  = 0.1,
    alpha = 0.5
  ) +
  coord_sf(
    xlim = c(x_min_fix, x_max_fix),
    ylim = c(y_min_fix, y_max_fix),
    expand = FALSE
  ) +
  labs(
    title = "prediction grid (masked)"
  )

# ---------------------------
# three-panel layout: original grid, surveys, fixed grid
# ---------------------------

fig <- p1 + p2 + p3
ggsave(
  filename = "prediction_grid.pdf",
  plot     = fig,
  width    = 14,
  height   = 11,
)

grid_sf_out <- st_as_sf(
  pred_grid_fix,
  coords = c("x", "y"),
  crs    = target_crs
)

st_write(
  grid_sf_out,
  "prediction_grid_utm26916_masked.gpkg",
  delete_dsn = TRUE
)

# -------------------------------------------------------------------
# idw depth + habitat probabilities, extract to pred grid, write pred_grid.gpkg
# -------------------------------------------------------------------

# read and prep spatial data for rasters
regions <- st_read("region_polys/Regions_2013.shp") |>
  st_transform(26916)

shoreline <- st_read("AOC_MI_StMarys_2022/AOC_MI_StMarys_2022.shp") |>
  st_transform(26916)

shoreline_f <- st_intersection(shoreline, regions)

# catch and habitat data
catch <- read.csv("all_data_through_2024/2024-10-24-CatchALL.csv")

catch_sf <- st_as_sf(
  catch,
  coords = c("longitude", "latitude"),
  crs = 4326
) |>
  st_transform(26916)

coords <- st_coordinates(catch_sf)
catch_sf$x <- coords[, 1]
catch_sf$y <- coords[, 2]

catch_sf$Type1 <- as.integer(catch_sf$hab.type == 1)
catch_sf$Type2 <- as.integer(catch_sf$hab.type == 2)
catch_sf$Type3 <- as.integer(catch_sf$hab.type == 3)

# gstat idw models
depth_gs <- gstat(
  id = "depth",
  formula = depth ~ 1,
  locations = ~ x + y,
  data = as.data.frame(catch_sf),
  nmax = 7,
  set = list(idp = 2)
)

habTypeI_gs <- gstat(
  id = "prob_type_I",
  formula = Type1 ~ 1,
  locations = ~ x + y,
  data = as.data.frame(catch_sf),
  nmax = 7,
  set = list(idp = 2)
)

habTypeII_gs <- gstat(
  id = "prob_type_II",
  formula = Type2 ~ 1,
  locations = ~ x + y,
  data = as.data.frame(catch_sf),
  nmax = 7,
  set = list(idp = 2)
)

habTypeIII_gs <- gstat(
  id = "prob_type_III",
  formula = Type3 ~ 1,
  locations = ~ x + y,
  data = as.data.frame(catch_sf),
  nmax = 7,
  set = list(idp = 2)
)

# prediction rasters
pred_raster_temp <- rast(
  nrows = 1000,
  ncols = 500,
  extent = ext(st_bbox(shoreline_f)),
  crs = "epsg:26916"
)

depth_raster <- interpolate(pred_raster_temp, depth_gs, debug.level = 0, index = 1)
TypeIhab_rast <- interpolate(pred_raster_temp, habTypeI_gs, debug.level = 0, index = 1)
TypeIIhab_rast <- interpolate(pred_raster_temp, habTypeII_gs, debug.level = 0, index = 1)
TypeIIIhab_rast <- interpolate(pred_raster_temp, habTypeIII_gs, debug.level = 0, index = 1)

cellArea <- cellSize(pred_raster_temp)

habitat_stack <- c(
  depth_raster,
  TypeIhab_rast,
  TypeIIhab_rast,
  TypeIIIhab_rast,
  cellArea
)

habitat_stack2river <- mask(habitat_stack, shoreline_f)

# extract depth and probabilities to prediction grid
grid_sf_out <- st_read("prediction_grid_utm26916_masked.gpkg")
grid_v <- vect(grid_sf_out)

hab_prob_rast <- habitat_stack2river[[c(
  "depth.pred",
  "prob_type_I.pred",
  "prob_type_II.pred",
  "prob_type_III.pred"
)]]

vals <- terra::extract(hab_prob_rast, grid_v)

grid_sf_out <- cbind(grid_sf_out, vals[, -1])

names(grid_sf_out)[names(grid_sf_out) == "depth.pred"] <- "depth"
names(grid_sf_out)[names(grid_sf_out) == "prob_type_I.pred"] <- "prob_I"
names(grid_sf_out)[names(grid_sf_out) == "prob_type_II.pred"] <- "prob_II"
names(grid_sf_out)[names(grid_sf_out) == "prob_type_III.pred"] <- "prob_III"

grid_sf_out <- subset(grid_sf_out, !is.na(depth))

summary(grid_sf_out$prob_I + grid_sf_out$prob_II + grid_sf_out$prob_III)

# plotting to compare rasters and point extractions
par(mfrow = c(2, 2))

plot(habitat_stack2river$depth.pred, main = "depth raster")
plot(st_geometry(grid_sf_out), add = TRUE, pch = 16, cex = 0.4)

plot(habitat_stack2river$prob_type_I.pred, main = "type I raster")
points(
  st_coordinates(grid_sf_out),
  pch = 16,
  cex = 0.4,
  col = gray(1 - grid_sf_out$prob_I)
)

plot(habitat_stack2river$prob_type_II.pred, main = "type II raster")
points(
  st_coordinates(grid_sf_out),
  pch = 16,
  cex = 0.4,
  col = gray(1 - grid_sf_out$prob_II)
)

plot(habitat_stack2river$prob_type_III.pred, main = "type III raster")
points(
  st_coordinates(grid_sf_out),
  pch = 16,
  cex = 0.4,
  col = gray(1 - grid_sf_out$prob_III)
)

par(mfrow = c(1, 1))

# write prediction grid for sdmTMB

# deal with habitat factor -- based on max probability
df_probs <- st_drop_geometry(grid_sf_out)

idx_max <- max.col(
  df_probs[, c("prob_I", "prob_II", "prob_III")],
  ties.method = "first"
)

grid_sf_out$hab_type <- idx_max
# could do this other ways...

nrow(grid_sf_out)
grid_sf_out <- grid_sf_out[which(grid_sf_out$hab_type != 3), ]
nrow(grid_sf_out)

# save the final grid
# st_write(
#  grid_sf_out,
#  "pred_grid.gpkg",
#  layer = "pred_grid",
#  delete_layer = TRUE
# )

# quick comparison plot: sampled data vs prediction grid

coords_grid <- st_coordinates(grid_sf_out)

x_min_all <- min(c(coords_grid[,1], dat_df$X))
x_max_all <- max(c(coords_grid[,1], dat_df$X))
y_min_all <- min(c(coords_grid[,2], dat_df$Y))
y_max_all <- max(c(coords_grid[,2], dat_df$Y))

p_sample <- ggplot() +
  geom_sf(data = shore_union, fill = "grey90", colour = NA) +
  geom_point(
    data = dat_df,
    aes(x = X, y = Y),
    size = 0.1
  ) +
  coord_sf(
    xlim = c(x_min_all, x_max_all),
    ylim = c(y_min_all, y_max_all),
    expand = FALSE
  ) +
  labs(title = "sampled data")

p_predgrid <- ggplot() +
  geom_sf(data = shore_union, fill = "grey90", colour = NA) +
  geom_sf(
    data = grid_sf_out,
    size = 0.1
  ) +
  coord_sf(
    xlim = c(x_min_all, x_max_all),
    ylim = c(y_min_all, y_max_all),
    expand = FALSE
  ) +
  labs(title = "prediction grid")

fig_pred <- p_sample + p_predgrid + plot_layout(ncol = 2)

# ggsave(
#   filename = "pred_grid.pdf",
#   plot = fig_pred,
#   width = 11,
#   height = 8
# )
# 

library(sf)
library(terra)
library(gstat)
library(ggplot2)

# -------------------------------------------------------------------
# idw depth + habitat probabilities, extract values and write pred_grid.gpkg
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

# set raster resolution based on intended cell area
shoreline_f_bbox <- st_bbox(shoreline_f)
width_m <- shoreline_f_bbox["xmax"] - shoreline_f_bbox["xmin"]
height_m <- shoreline_f_bbox["ymax"] - shoreline_f_bbox["ymin"]

cell_height <- 100
cell_width <- 100

nrows_pred_raster <- round(width_m / cell_width)
ncols_pred_raster <- round(height_m / cell_height)
# prediction rasters
pred_raster_temp <- rast(
  nrows = nrows_pred_raster,
  ncols = ncols_pred_raster,
  extent = ext(shoreline_f_bbox),
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

grid_out <- as.data.frame(habitat_stack2river, xy = TRUE)
grid_sf_out <- st_as_sf(grid_out, coords = c("x", "y"), crs = 26916)

names(grid_sf_out)[names(grid_sf_out) == "depth.pred"] <- "depth"
names(grid_sf_out)[names(grid_sf_out) == "prob_type_I.pred"] <- "prob_I"
names(grid_sf_out)[names(grid_sf_out) == "prob_type_II.pred"] <- "prob_II"
names(grid_sf_out)[names(grid_sf_out) == "prob_type_III.pred"] <- "prob_III"

summary(grid_sf_out$prob_I + grid_sf_out$prob_II + grid_sf_out$prob_III)

# plotting to compare rasters and point extractions
par(mfrow = c(2, 2))

plot(habitat_stack2river$depth.pred, main = "depth raster")

plot(habitat_stack2river$prob_type_I.pred, main = "type I raster")

plot(habitat_stack2river$prob_type_II.pred, main = "type II raster")

plot(habitat_stack2river$prob_type_III.pred, main = "type III raster")

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

ggplot() +
  geom_sf(
    data = grid_sf_out,
    aes(color = as.factor(hab_type)),
    size = 0.3
  ) + 
  guides(colour = guide_legend(override.aes = list(size = 5)))

# remove type three habitat (no samples)
nrow(grid_sf_out)
grid_sf_out <- grid_sf_out[which(grid_sf_out$hab_type != 3), ]
nrow(grid_sf_out)

xy <- st_coordinates(grid_sf_out)

# add X and Y columns in Km like analysis
grid_sf_out$X <- xy[, 1]/1000
grid_sf_out$Y <- xy[, 2]/1000

grid_sf_out <- as.data.frame(grid_sf_out)
grid_sf_out$geom <- NULL
grid_sf_out <- grid_sf_out[,c("X", "Y", "depth", "hab_type", "area")]

# save the final grid
save(grid_sf_out, file = "pred_grid.RData")

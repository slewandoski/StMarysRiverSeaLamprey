# fit some index standardization models using sdmTMB
# each sampled unit is 2.44 m^2
# see https://sdmtmb.github.io/sdmTMB/
#
# note, this takes about 66 minutes to run on my desktop
#
# todo:
# - compare abundance predictions among models

library(sdmTMB)
library(sf)
library(ggplot2)
library(dplyr)
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
plot(mesh)

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
# note these will change if changing mesh...

# general plots and diagnostics, can (should) be run for each model
sanity(m8)
tidy(m8, conf.int = TRUE)
tidy(m8, effects = "ran_pars", conf.int = TRUE)
# plot depth effect
ggeffects::ggpredict(m8, terms="depth[0:40, by = 2]") |> plot()
data$resids <- residuals(m8) # randomized quantile residuals
hist(data$resids)
ggplot(data, aes(X, Y, col = resids)) + scale_colour_gradient2() +
    geom_point() + facet_wrap(~year, nrow = 3) + coord_fixed()
set.seed(19283)
s <- simulate(m8, nsim = 1000, type = "mle-mvn")
dharma_residuals(s, m8)
abline(0,1)

ggplot(data, aes(X, Y, col = resids)) +
  scale_colour_gradient2() +
  geom_point() +
  facet_wrap(~year, nrow = 4) +
  coord_fixed()

#----------------------------------------------------------------------
# predictions 
#----------------------------------------------------------------------

data$year_fac <- as.factor(data$year)
pred_grid <- st_read("pred_grid.gpkg")
# extract coordinates
xy <- st_coordinates(pred_grid)

# add X and Y columns
pred_grid$X <- xy[, 1]/1000
pred_grid$Y <- xy[, 2]/1000

pred_grid <- as.data.frame(pred_grid)
pred_grid$geom <- NULL
pred_grid <- pred_grid[,c("X", "Y", "depth", "hab_type")]
grid_yrs <- replicate_df(pred_grid, "year", unique(data$year))
grid_yrs$year_fac <- as.factor(grid_yrs$year)
predictions <- predict(m8, newdata = grid_yrs, return_tmb_object = TRUE)

# function to make maps
plot_map <- function(dat, column) {
  ggplot(dat, aes(X, Y, fill = {{ column }})) +
    geom_raster() +
    facet_wrap(~year, nrow = 3) +
    coord_fixed()
}

p1 <- plot_map(predictions$data, exp(est)) +
  scale_fill_viridis_c(trans = "sqrt") +
  ggtitle("Prediction (fixed effects + all random effects)")

p2 <- plot_map(predictions$data, exp(est_non_rf)) +
  ggtitle("Prediction (fixed effects only)") +
  scale_fill_viridis_c(trans = "sqrt")

p3 <- plot_map(predictions$data, omega_s) +
  ggtitle("Spatial random effects only") +
  scale_fill_gradient2()

p4 <- plot_map(predictions$data, epsilon_st) +
  ggtitle("Spatiotemporal random effects only") +
  scale_fill_gradient2()

# area of each grid is 100 m by 100 m = 0.01 km^2
index <- get_index(predictions, area = 0.01, bias_correct = TRUE)
p5 <- ggplot(index, aes(year, est)) + geom_line() +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.4) +
  xlab('Year') + ylab('Number of Larvae')

# multi-page pdf"
pdf("ar1st_idx_std.pdf", width = 15, height = 10)
print(p1)
print(p2)
print(p3)
print(p4)
print(p5)
dev.off()

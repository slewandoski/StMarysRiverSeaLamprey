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


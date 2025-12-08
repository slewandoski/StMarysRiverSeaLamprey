# fit some index standardization models using sdmTMB
#
# see https://sdmtmb.github.io/sdmTMB/ 
#
# note, this takes about 66 minutes to run on my desktop
# 
# todo: 
# - make a prediction grid during duration of study
# - think about incorporating the habitat factor sean talked about
# - create abundance indeces based on count modeling from negative binomial 
# - compare abundance predictions among models

library(sdmTMB)
library(sf)
library(ggplot2)

# read CSV and convert to sf with lat/lon in degrees
data <- read.csv("all_data_through_2024/2024-10-24-CatchALL.csv")
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

mesh <- make_mesh(data, xy_cols = c("X", "Y"), cutoff = 0.25)
# mesh$mesh$n
# plot(mesh$mesh,
#     main = "study area with mesh",
#     xlab = "Easting", ylab = "Northing"
# )
# points(coords, cex = 0.1, pch = 1, col = "steelblue4")

# no space or year effect
m1 <- sdmTMB(
             c ~ s(depth, k = 5), # smoother on depth
             data = data,
             mesh = mesh,
             family = nbinom2(link = "log"),
             spatial = "off",
             silent = FALSE
)

# space effect shared among years
m2 <- sdmTMB(
             c ~ 1 + s(depth, k = 5),
             data = data,
             mesh = mesh,
             family = nbinom2(link = "log"),
             spatial = "on",
             silent = FALSE
)

# intercept factor by year, no space or time effects
m3 <- sdmTMB(
             c ~ 0 + as.factor(year) + s(depth, k = 5),
             time = "year",
             data = data,
             mesh = mesh,
             family = nbinom2(link = "log"),
             spatial = "off",
             silent = FALSE
)

# intercept factor by year, spatial effect but not time
m4 <- sdmTMB(
             c ~ 0 + as.factor(year) + s(depth, k = 5),
             time = "year",
             data = data,
             mesh = mesh,
             family = nbinom2(link = "log"),
             spatial = "on",
             silent = FALSE
)

# one intercept all years, but IID spatial-temporal fields
m5 <- sdmTMB(
             c ~ 1 + s(depth, k = 5),
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
             c ~ 1 + s(depth, k = 5),
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
             c ~ 1 + s(depth, k = 5),
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
             c ~ 0 + year_fac + s(depth, k = 5),
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
#   df      AIC
#m1  4 31714.61
#m2  6 28159.21
#m3 34 28177.78 
#m4 35 27447.36 fairly "basic" idx standardization model
#m5  7 27524.55
#m6  7 27667.27
#m7  8 27516.50
#m8 36 27447.02 fairly complicated model, interpolating missing yrs

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

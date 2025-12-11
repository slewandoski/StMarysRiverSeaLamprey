# fit some index standardization models using sdmTMB
# see https://sdmtmb.github.io/sdmTMB/
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
# test mesh
mesh <- make_mesh(data, xy_cols = c("X", "Y"), cutoff = 1)
# mesh$mesh$n
plot(mesh)

# intercept factor by year, spatio-temporal AR1 among years
data$year_fac <- as.factor(data$year)
m1 <- sdmTMB(
             c ~ 1 + year_fac + as.factor(hab_type) + s(depth, k = 5),
             time = "year",
             data = data,
             mesh = mesh,
             family = nbinom2(link = "log"),
             spatial = "on",
             spatiotemporal = "off",
             extra_time = c(1997, 1998, 2020), # ensure regular spacing
             silent = FALSE
)

sanity(m1)

m2 <- sdmTMB(
             c ~ 1 + year_fac + as.factor(hab_type) + s(depth, k = 5),
             time = "year",
             data = data,
             mesh = mesh,
             family = nbinom2(link = "log"),
             spatial = "on",
             spatiotemporal = "iid",
             extra_time = c(1997, 1998, 2020), # ensure regular spacing
             silent = FALSE
)
sanity(m2)

m3 <- sdmTMB(
             c ~ 1 + year_fac + as.factor(hab_type) + s(depth, k = 5),
             time = "year",
             data = data,
             mesh = mesh,
             family = nbinom2(link = "log"),
             spatial = "on",
             spatiotemporal = "ar1",
             extra_time = c(1997, 1998, 2020), # ensure regular spacing
             silent = FALSE
)
sanity(m3)

m4 <- sdmTMB(
             c ~ 1 + year_fac + as.factor(hab_type) + s(depth, k = 5),
             time = "year",
             data = data,
             mesh = mesh,
             family = poisson(link = "log"),
             spatial = "on",
             spatiotemporal = "off",
             extra_time = c(1997, 1998, 2020), # ensure regular spacing
             silent = FALSE
)
sanity(m4)

m5 <- sdmTMB(
             c ~ 1 + year_fac + as.factor(hab_type) + s(depth, k = 5),
             time = "year",
             data = data,
             mesh = mesh,
             family = poisson(link = "log"),
             spatial = "on",
             spatiotemporal = "iid",
             extra_time = c(1997, 1998, 2020), # ensure regular spacing
             silent = FALSE
)
sanity(m5)

m6 <- sdmTMB(
             c ~ 1 + year_fac + as.factor(hab_type) + s(depth, k = 5),
             time = "year",
             data = data,
             mesh = mesh,
             family = poisson(link = "log"),
             spatial = "on",
             spatiotemporal = "iid",
             extra_time = c(1997, 1998, 2020), # ensure regular spacing
             silent = FALSE
)
sanity(m6)

print(AIC(m1, m2, m3, m4, m5, m6))

dir.create("output", showWarnings = FALSE)

save(
  m1, m2, m3, m4, m5, m6,
  mesh,
  data,
  file = "output/sdmTMB_models.RData"
)

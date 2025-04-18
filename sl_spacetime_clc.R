#-------------------------------------------------------------------------------
# St. Mary's larval lamprey in space-time
# Cahill 18 April 2025
#
# TODO:
# figure out a way to get rid of that overdispersion
#
#-------------------------------------------------------------------------------
library(fmesher)
library(sf)
library(RTMB)

# read CSV and convert to sf with lat/lon in degrees
data <- read.csv("all_data_through_2024/2024-10-24-CatchALL.csv")
data <- st_as_sf(data, coords = c("longitude", "latitude"), crs = 4326) # WGS84

# Transform to UTM Zone 16N (NAD83)
data <- st_transform(data, crs = 26916) # EPSG:26916

# start with 2024 data
data <- subset(data, year >= 2000)

# extract UTM coordinates in km
coords <- st_coordinates(data) / 1000
data$easting_km <- coords[, 1]
data$northing_km <- coords[, 2]

# rename n for simplicity
data$n <- data$sl.larv.n

# filter for 2024 and plot
data_2024 <- subset(data, year == 2024)
plot(data_2024$northing_km ~ data_2024$easting_km,
  las = 1,
  xlab = "Easting (km)",
  ylab = "Northing (km)",
  cex = 0.5
)

# proportion of zeros
sum(data$n == 0) / nrow(data)

# simple glm
m <- glm(data$n ~ 1, family = poisson)
-logLik(m)

#-------------------------------------------------------------------------------
# set up spde approximation
#-------------------------------------------------------------------------------
mesh <- fm_mesh_2d(coords, refine = TRUE, cutoff = 0.9)
mesh$n # mesh nodes --> random effects we must estimate, so coarse to start
plot(mesh,
  main = "study area with mesh",
  xlab = "Easting", ylab = "Northing"
)
points(coords, cex = 0.1, pch = 1, col = "steelblue4")

# create matrices in fmesher / INLA
spde <- fm_fem(mesh, order = 2)

# create projection matrix from vertices to sample locations
A_is <- fm_evaluator(mesh, loc = st_coordinates(data) / 1000)$proj$A

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
  "M0" = spde$c0, "M1" = spde$g1, "M2" = spde$g2
)

par <- list(
  "beta0" = 0,
  "ln_tauO" = 0,
  "ln_tauE" = 0,
  "ln_kappa" = 0,
  "trans_rho" = 0,
  "omega_s" = numeric(nrow(spde$c0)),
  "epsilon_st" = matrix(0.0, nrow(spde$c0), data$n_t)
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
  ADREPORT(range)
  jnll
}

obj <- MakeADFun(f, par, random = c("omega_s", "epsilon_st"))
obj$fn()
obj$gr()

opt <- nlminb(obj$par, obj$fn, obj$gr)
opt
sdr <- sdreport(obj, bias.correct = TRUE)
sdr

# extract spatial range
sdr$value
sdr$sd

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
  "omega_s" = numeric(nrow(spde$c0)),
  "epsilon_st" = matrix(0.0, nrow(spde$c0), data$n_t)
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
  ln_var_minus_mu_c <- 2 * ln_pred_i - ln_theta
  jnll <- jnll - sum(dnbinom_robust(c_i, ln_pred_i,
    log_var_minus_mu = ln_var_minus_mu_c,
    TRUE
  ))
  range <- sqrt(8) / exp(ln_kappa)
  # compute avg overdispersion across observations:
  # Overdispersion = 1 + mu / theta, where mu = exp(ln_pred_i)
  od <- mean(1 + exp(ln_pred_i - ln_theta))
  # Predicted mean and variance for each observation
  mu_i <- exp(ln_pred_i)
  theta <- exp(ln_theta)
  var_i <- mu_i + mu_i^2 / theta
  
  # Compute Pearson residuals
  resid_pearson <- (c_i - mu_i) / sqrt(var_i)
  
  # Report individual residuals (e.g., for spatial plotting)
  REPORT(mu_i)
  REPORT(resid_pearson)
  ADREPORT(od)
  ADREPORT(range)
  jnll
}

obj <- MakeADFun(f, par, random = c("omega_s", "epsilon_st"))
obj$fn()
obj$gr()

opt <- nlminb(obj$par, obj$fn, obj$gr)
opt

sdr <- sdreport(obj, bias.correct = TRUE)
sdr

sdr$value
sdr$sd

plot(obj$report()$'mu_i', obj$report()$'resid_pearson',
     pch = 20, main = "Pearson residuals vs Fitted", 
     xlab = "Fitted values (mu_i)", ylab = "Pearson residuals")
abline(h = 0, col = "gray")

# start of processing script

# general plots and diagnostics, can (should) be run for each model
tidy(m1, conf.int = TRUE)
tidy(m1, effects = "ran_pars", conf.int = TRUE)
# plot depth effect
ggeffects::ggpredict(m1, terms="depth[0:40, by = 2]") |> plot()
data$resids <- residuals(m1) # randomized quantile residuals
hist(data$resids)
ggplot(data, aes(X, Y, col = resids)) + scale_colour_gradient2() +
    geom_point() + facet_wrap(~year, nrow = 3) + coord_fixed()
set.seed(19283)
s <- simulate(m1, nsim = 1000, type = "mle-mvn")
dharma_residuals(s, m1)
abline(0,1)

ggplot(data, aes(X, Y, col = resids)) +
    scale_colour_gradient2() +
    geom_point() +
    facet_wrap(~year, nrow = 4) +
    coord_fixed()

# install.packages("remotes")
# remotes::install_github("pbs-assess/sdmTMBextra", dependencies = TRUE)
#set.seed(123)
#samps <- sdmTMBextra::predict_mle_mcmc(m1, mcmc_iter = 800, mcmc_warmup = 400)
#mcmc_res <- residuals(m1, type = "mle-mcmc", mcmc_samples = samps)
#qqnorm(mcmc_res)
#abline(0, 1)

# simulation based stuff - number of zeros
s_nb2 <- simulate(m1, nsim = 500, type = "mle-mvn")
sum(data$c == 0) / length(data$c) # observed
sum(s_nb2 == 0)/length(s_nb2) # predicted

r_nbinom <- dharma_residuals(s_nb2, m1, return_DHARMa = TRUE)
plot(r_nbinom)
DHARMa::testResiduals(r_nbinom)

#----------------------------------------------------------------------
# predictions 
#----------------------------------------------------------------------
load("pred_grid.RData")
data$year_fac <- as.factor(data$year)
grid_yrs <- replicate_df(pred_grid, "year", unique(data$year))
grid_yrs$year_fac <- as.factor(grid_yrs$year)

predictions <- predict(m2, newdata = grid_yrs, return_tmb_object = TRUE)

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

# area of each grid is 100 m by 100 m, and quadrat is 2.44m^2
index <- get_index(predictions, area = grid_yrs$area/2.44, bias_correct = TRUE)
p5 <- ggplot(index, aes(year, est)) + geom_line() +
    geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.4) +
    xlab('Year') + ylab('Number of juvenile lamprey')

# multi-page pdf"
pdf("ar1st_idx_std.pdf", width = 15, height = 11)
print(p1)
print(p2)
print(p3)
print(p4)
print(p5)
dev.off()

# play
# center of gravity:
#cog <- get_cog(predictions, format = "wide")
#cog
#ggplot(cog, aes(est_x, est_y, colour = year)) +
#  geom_point() +
#  geom_linerange(aes(xmin = lwr_x, xmax = upr_x)) +
#  geom_linerange(aes(ymin = lwr_y, ymax = upr_y)) +
#  scale_colour_viridis_c()

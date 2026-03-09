setwd("~/School/STAT-S799/IncomesDataAnalysis/")
source("classes.R")

n <- 100000
mu <- 10
log_data <- data.frame(
  x = rlnorm(n, mu, 2),
  weights = rlnorm(n, 2, 0.5) * rbinom(n, 1, 0.5),
  time = 1
)

for (i in 1:39) {
  mu <- c(mu, mu[i] * 0.98 + rnorm(1, 0, 0.01))
  new_log_data <- data.frame(
    x = rlnorm(n, mu[i + 1], 1),
    weights = rlnorm(n, 2, 0.5) * rbinom(n, 1, 0.5),
    time = i + 1
  )
  log_data <- rbind(log_data, new_log_data)
}

analysis_obj <- DensityTimeSeries$new(
  log_data$x,
  log_data$time,
  log_data$weights
)

# If n=4096, we get errors from dens2quantile.
# Presumably this is due to numerical precison issues
analysis_obj$create_dens_grid(h = 10, n = 2048)

# Find a better way to pick h.
# Using a CV technique would be best, but the custom function
# is quite slow, so this would take too long.
analysis_obj$create_dens(h = 10)

# getBinnedData in fdapace starts binning data if we
# make this increment too small.
analysis_obj$create_quants(seq(0, 1, 0.01))

# FDA AR1 Model

FDA_AR_fits <- lapply(20:40, analysis_obj$fda_ar)
FDA_AR_distances <- sapply(FDA_AR_fits, function(x) x[[1]])
mean(FDA_AR_distances)

# Wasserstein AR model

WAR_fits <- lapply(20:40, analysis_obj$wasserstein_ar, order = 1)
WAR_distances <- sapply(WAR_fits, function(x) x[[1]])
mean(WAR_distances)

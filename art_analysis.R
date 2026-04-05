setwd("~/School/STAT-S799/IncomesDataAnalysis/")
source("classes.R")

create_data <- function(n = 1000, mu = 10, times = 40) {
  log_data <- data.frame(
    x = rlnorm(n, mu, 2),
    weights = rlnorm(n, 2, 0.5) * rbinom(n, 1, 0.5),
    time = 1
  )

  for (i in 1:(times - 1)) {
    mu <- c(mu, mu[i] * 0.98 + rnorm(1, 0, 0.01))
    new_log_data <- data.frame(
      x = rlnorm(n, mu[i + 1], 1),
      weights = rlnorm(n, 2, 0.5) * rbinom(n, 1, 0.5),
      time = i + 1
    )
    log_data <- rbind(log_data, new_log_data)
  }

  log_data
}

create_analaysis_obj <- function(data) {
  analysis_obj <- DensityTimeSeries$new(
    data$x,
    data$time,
    data$weights
  )

  # If n=4096, we get errors from dens2quantile.
  # Presumably this is due to numerical precison issues
  analysis_obj$create_dens_grid(n = 1024)

  # getBinnedData in fdapace starts binning data if we
  # make this increment too small.
  analysis_obj$create_quants(seq(0, 1, 0.01))

  k <- analysis_obj$select_knn_bandwidth(16, verbose = TRUE)

  analysis_obj$create_dens_knn(k)

  analysis_obj
}

log_data <- create_data()
analysis_obj <- create_analaysis_obj(log_data)

# Function to test models
test_model <- function(times, func, ...) {
  ar_fits <- lapply(times, func, ...)
  ar_distances <- sapply(ar_fits, function(x) analysis_obj$wass_dist_ar(x))
  mean(ar_distances)
}

test_all_models <- function(times, obj) {
  c(
    "FDA" = test_model(times, analysis_obj$fda_ar),
    "Bayes" = test_model(times, analysis_obj$bayes_ar),
    "LQD" = test_model(times, analysis_obj$lqd_ar),
    "Wasserstein" = mean(sapply(times, function(x) {
      ar_obj <- analysis_obj$wasserstein_ar(x, order = 1)
      ar_obj[[1]]
    }))
  )
}

test_all_models(20:40, analysis_obj)

setwd("~/School/STAT-S799/IncomesDataAnalysis/")
source("classes.R")
source("plot.R")

set.seed(1234)

create_analaysis_obj <- function(data, ...) {
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

  k <- analysis_obj$select_knn_bandwidth(16, verbose = TRUE, ...)

  analysis_obj$create_dens_knn(k)

  analysis_obj
}

# Creates a list of models to plot with
models_func <- function(analysis_obj) {
  list(
    "FDA" = function(t) analysis_obj$fda_ar(t, ensure_positive = TRUE),
    "Bayes" = analysis_obj$bayes_ar,
    "LQD" = analysis_obj$lqd_ar,
    "Wasserstein" = analysis_obj$wasserstein_ar
  )
}

# Function to test models
test_model <- function(times, analysis_obj, func) {
  ar_fits <- lapply(times, func)
  ar_distances <- sapply(ar_fits, function(x) analysis_obj$wass_dist_ar(x))
  mean(ar_distances)
}

test_all_models <- function(times, analysis_obj, models) {
  wass_dists <- numeric(length(models))
  names(wass_dists) <- names(models)
  for (i in seq_along(models)) {
    wass_dists[i] <- test_model(times, analysis_obj, models[[i]])
  }

  wass_dists
}

# Create animation
create_anim <- function(title, analysis_obj, times, models, asinh_scale = FALSE) {
  anim <- animate_all_models(title, analysis_obj, times, models, asinh_scale = asinh_scale)
  animate(anim, nframes = 60, fps = 3, width = 2000, height = 2000, res = 200)
}

# Log data
create_data_log <- function(n = 1000, mu = 10, times = 40) {
  data <- data.frame(
    x = rlnorm(n, mu, 1),
    weights = 1,
    time = 1
  )

  for (i in 1:(times - 1)) {
    mu <- c(mu, mu[i] * 0.98 + rnorm(1, 0, 0.01))
    new_data <- data.frame(
      x = rlnorm(n, mu[i + 1], 1),
      weights = 1,
      time = i + 1
    )
    data <- rbind(data, new_data)
  }

  data
}
data <- create_data_log()
analysis_obj <- create_analaysis_obj(data)
models <- models_func(analysis_obj)
test_all_models(20:40, analysis_obj, models)

create_anim("Log Normal", analysis_obj, 20:40, models, asinh_scale = TRUE)
anim_save("log_norm.gif")

create_data_norm <- function(n = 1000, mu = 0, times = 40) {
  data <- data.frame(
    x = rnorm(n, mu, 1),
    weights = 1,
    time = 1
  )

  for (i in 1:(times - 1)) {
    mu <- c(mu, mu[i] * 1.1 + rnorm(1, 0, 0.10))
    new_data <- data.frame(
      x = rnorm(n, mu[i + 1], 1),
      weights = 1,
      time = i + 1
    )
    data <- rbind(data, new_data)
  }

  data
}
data <- create_data_norm()
analysis_obj <- create_analaysis_obj(data)
models <- models_func(analysis_obj)
test_all_models(20:40, analysis_obj, models)

create_anim("Normal", analysis_obj, 20:40, models)
anim_save("norm.gif")

create_data_unif <- function(n = 1000, a = 0, b = 1, times = 40) {
  data <- data.frame(
    x = runif(n, a, b),
    weights = 1,
    time = 1
  )

  for (i in 1:(times - 1)) {
    a <- a * 1.05 + rnorm(1, 0, 0.10)
    b <- b * 1.05 + rnorm(1, 0, 0.10)
    new_data <- data.frame(
      x = runif(n, a, b),
      weights = 1,
      time = i + 1
    )
    data <- rbind(data, new_data)
  }

  data
}
data <- create_data_unif()
analysis_obj <- create_analaysis_obj(data)
models <- models_func(analysis_obj)
test_all_models(20:40, analysis_obj, models)

create_anim("Uniform", analysis_obj, 20:40, models)
anim_save("unif.gif")

setwd("~/School/STAT-S799/IncomesDataAnalysis/")
source("classes.R")

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
    "Wasserstein" = test_model(times, analysis_obj$wasserstein_ar)
  )
}

source("plot.R")

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
test_all_models(20:40, analysis_obj)

models <- list(
  "FDA" = analysis_obj$fda_ar,
  "Bayes" = analysis_obj$bayes_ar,
  "LQD" = analysis_obj$lqd_ar,
  "Wasserstein" = function(t) analysis_obj$wasserstein_ar(t, order = 1)
)

plot_all_models_vs_actual(analysis_obj, 40, models, log_scale = TRUE)
anim <- animate_all_models(analysis_obj, 20:40, models, log_scale = TRUE)
animate(anim, nframes = 50, fps = 5, width = 1000, height = 1000)

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
analysis_obj <- create_analaysis_obj(data, lower_mul = 1 / 50, upper_mul = 1 / 2)
test_all_models(20:40, analysis_obj)

models <- list(
  "FDA" = analysis_obj$fda_ar,
  "Bayes" = analysis_obj$bayes_ar,
  "LQD" = analysis_obj$lqd_ar,
  "Wasserstein" = function(t) analysis_obj$wasserstein_ar(t, order = 1)
)

# Generate static plot for timestep 40
plot_all_models_vs_actual(analysis_obj, 30, models)
# Generate animation for timesteps 20 through 40
anim <- animate_all_models(analysis_obj, 20:40, models)
# Render the animation (adjust frames/fps as needed for execution speed)
animate(anim, nframes = 50, fps = 5, width = 1000, height = 1000)

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
analysis_obj <- create_analaysis_obj(data, lower_mul = 1 / 20, upper_mul = 1 / 2)
test_all_models(20:40, analysis_obj)

models <- list(
  "FDA" = analysis_obj$fda_ar,
  "Bayes" = analysis_obj$bayes_ar,
  "LQD" = analysis_obj$lqd_ar,
  "Wasserstein" = function(t) analysis_obj$wasserstein_ar(t, order = 1)
)

# Generate static plot for timestep 40
plot_all_models_vs_actual(analysis_obj, 40, models)
# Generate animation for timesteps 20 through 40
anim <- animate_all_models(analysis_obj, 20:40, models)
# Render the animation (adjust frames/fps as needed for execution speed)
animate(anim, nframes = 50, fps = 5, width = 1000, height = 1000)

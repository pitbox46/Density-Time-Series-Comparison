library(patchwork)
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

# Create plot showing raw FPCA
fda_obj <- analysis_obj$fda_ar(40)
fda_obj <- data.frame(
  x = asinh(fda_obj$forecast_dens$mean$x),
  y = fda_obj$forecast_pdf
)
plott <- ggplot(fda_obj, aes(x = x, y = y)) +
  geom_line() +
  geom_hline(yintercept = 0, color = "grey") +
  theme_minimal() +
  labs(
    title = "Log Normal - FPCA Forecast (t = 40)",
    x = "Inverse Hyperbolic Sine Value",
    y = "Probability"
  )
ggsave(
  "media/log_norm40fpca.png",
  plot = plott,
  width = 1600,
  height = 1600,
  units = "px",
  dpi = 200
)

# Create plot showing raw Bayes
bayes_obj <- analysis_obj$bayes_ar(40)
bayes_obj <- data.frame(
  x = asinh(bayes_obj$forecast_dens$mean$x),
  y = bayes_obj$forecast_raw,
  z = bayes_obj$forecast_pdf
)
plott1 <- ggplot(bayes_obj, aes(x = x, y = y)) +
  geom_line() +
  geom_hline(yintercept = 0, color = "grey") +
  theme_minimal() +
  labs(
    title = "Bayes Space - Untransformed Forecast (t = 40)",
    x = "Inverse Hyperbolic Sine Value",
    y = "Bayes Space Value"
  )
plott2 <- ggplot(bayes_obj, aes(x = x, y = z)) +
  geom_line() +
  geom_hline(yintercept = 0, color = "grey") +
  theme_minimal() +
  labs(
    title = "Bayes Space - PDF Forecast (t = 40)",
    x = "Inverse Hyperbolic Sine Value",
    y = "Probability"
  )
plott <- plott1 + plott2
ggsave(
  "media/log_norm40bayes.png",
  plot = plott,
  width = 3200,
  height = 1600,
  units = "px",
  dpi = 200
)

plot_title <- "Log Normal"

save_plot(
  plot_title,
  "media/log_norm40.png",
  analysis_obj,
  40,
  models,
  asinh_scale = TRUE
)

create_anim(plot_title, analysis_obj, 20:40, models, asinh_scale = TRUE)
anim_save("media/log_norm.mp4")

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

plot_title <- "Normal"

save_plot(
  plot_title,
  "media/norm40.png",
  analysis_obj,
  40,
  models,
  asinh_scale = FALSE
)

create_anim(plot_title, analysis_obj, 20:40, models)
anim_save("media/norm.mp4")

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

plot_title <- "Uniform"

save_plot(
  plot_title,
  "media/unif40.png",
  analysis_obj,
  40,
  models,
  asinh_scale = FALSE
)

create_anim(plot_title, analysis_obj, 20:40, models)
anim_save("media/unif.mp4")

# Running uniform again with higher order for ftsm

analysis_obj <- create_analaysis_obj(data)
analysis_obj$ftsm_order <- 6
models <- models_func(analysis_obj)
test_all_models(20:40, analysis_obj, models)

plot_title <- "Uniform"

save_plot(
  plot_title,
  "media/unif40_order6.png",
  analysis_obj,
  40,
  models,
  asinh_scale = FALSE
)

create_anim(plot_title, analysis_obj, 20:40, models)
anim_save("media/unif_order6.mp4")

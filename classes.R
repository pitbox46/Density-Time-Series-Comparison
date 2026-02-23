library(plyr)
library(Hmisc)
library(fdadensity)
library(ftsa)
library(WRI)

# Static functions

# A mask for the approx() function
cdf2quant <- function(cdf, cdf_grid, quant_grid) {
  approx(cdf, cdf_grid, quant_grid)
}

wass_dist <- function(quant1, quant2, quant_grid) {
  ret <- 0
  for (i in 1:(length(quant_grid) - 1)) {
    step <- quant_grid[i + 1] - quant_grid[i]

    ret <- ret + (quant2[i] - quant1[i])^2 * step
  }
  sqrt(ret)
}

# Computes a KDE using the normal kernel from a supplied grid of points
density_from_grid <- function(data, h, grid, weights = NULL) {
  if (is.null(weights)) {
    n <- length(data)
    weights <- rep(1, n)
  } else {
    n <- sum(weights)
    data <- data[weights != 0]
    weights <- weights[weights != 0]
  }

  grid |>
    sapply(function(x) weights * dnorm((data - x) / h, 0, 1)) |>
    apply(2, sum) / n / h
}

library(R6)
DensityTimeSeries <- R6Class(
  "DensityTimeSeries",
  public = list(
    data = NULL,
    dens_grid = NULL,
    dens_mat = NULL,
    quant_grid = NULL,
    quant_mat = NULL,
    initialize = function(data_x = NA, data_time = NA, data_weights = NA) {
      self$data <- data.frame(
        time = data_time,
        x = data_x,
        weights = ifelse(is.na(data_weights), 1, data_weights)
      )
    },
    create_dens = function(h, dens_grid) {
      self$dens_grid <- dens_grid
      densities <- daply(
        self$data,
        .(time),
        function(xx) {
          dens <- density_from_grid(
            xx$x,
            h = h,
            grid = dens_grid,
            weights = xx$weights
          )
          dens
        }
      )
      self$dens_mat <- t(densities)
    },
    create_quants = function(quant_grid) {
      self$quant_grid <- quant_grid
      self$quant_mat <- daply(
        self$data,
        .(time),
        \(xx) wtd.quantile(xx$x, weights = xx$weights, probs = quant_grid)
      ) |> t()
    },
    get_data = function(target_time = NA) {
      if (is.na(target_time)) {
        self$data
      } else {
        subset(self$data, time == target_time)
      }
    },
    get_quant = function(target_time = NA) {
      if (is.na(target_time)) {
        self$quant_mat
      } else {
        self$quant_mat[, which(colnames(self$quant_mat) == target_time)]
      }
    },
    # Predicts the target year using a model built from all years prior
    # Uses a naive FPCA approach
    fda_ar = function(target_year) {
      dens_fts <- fts(
        self$dens_grid,
        self$dens_mat[, which(colnames(self$dens_mat) < target_year)]
      )
      dens_ftsm <- ftsm(dens_fts, order = 3)
      # Not sure if this is AR(1) or some other model
      forecast_dens <- forecast(
        dens_ftsm,
        method = "arima",
        level = 95,
        h = 1
      )
      forecast_pdf <- forecast_dens$mean$y[, 1]

      # Sometimes pdf_forecast is less than zero, which is problematic
      wasserstein_dist <- wass_dist(
        self$get_quant(target_year),
        dens2quantile(
          matrix(ifelse(forecast_pdf < 0, 0, forecast_pdf), nrow = 1),
          forecast_dens$mean$x,
          self$quant_grid
        ),
        self$quant_grid
      )

      list(wasserstein_dist, forecast_dens)
    },
    # Predicts the target year using a model built from all years prior
    # Uses WARp model
    wasserstein_ar = function(target_year, order) {
      data_WAR1 <- WARp(
        self$quant_mat[, which(colnames(self$quant_mat) < target_year)],
        self$quant_grid,
        order
      )
      forecast_war <- predict(data_WAR1, self$dens_grid, self$dens_grid)

      forecast_quant <- cdf2quant(
        forecast_war$pred.cdf,
        self$dens_grid[-length(self$dens_grid)],
        self$quant_grid
      )

      wasserstein_dist <- wass_dist(
        self$get_quant(target_year),
        forecast_quant$y,
        self$quant_grid
      )

      list(wasserstein_dist, forecast_quant)
    }
  ),
  private = list()
)

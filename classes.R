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
Rcpp::sourceCpp("density_from_grid.cpp")
density_from_grid <- function(data, h, grid, weights = NULL, cutoff = 6) {
  if (is.null(weights)) {
    weights <- rep_len(1, length(data))
  }
  keep <- weights != 0
  density_from_grid_cpp(data[keep], h, grid, weights[keep], cutoff)
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
      # Sort data for easier computations later
      self$data <- self$data[order(self$data$time, self$data$x), ]
    },
    # Creates a grid for KDE evaultion based on the quanitiles of all data
    create_dens_grid = function(h, n) {
      self$dens_grid <- Hmisc::wtd.quantile(
        self$data$x,
        weights = self$data$weights,
        probs = seq(0, 1, length.out = n - 2)
      )
      self$dens_grid <- c(
        min(self$dens_grid) - 3 * h,
        self$dens_grid,
        max(self$dens_grid) + 3 * h
      )
    },
    # A normal density is non-zero everywhere, but that means we must compute
    # several values which are essentially zero.
    # Adding a cutoff limits the number of data points we must iterate over.
    # The actual cutoff "radius" is h * cutoff
    create_dens = function(h, cutoff = 6) {
      densities <- daply(
        self$data,
        .(time),
        function(xx) {
          dens <- density_from_grid(
            xx$x,
            h = h,
            grid = self$dens_grid,
            weights = xx$weights,
            cutoff = cutoff
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
    fda_ar = function(target_time) {
      dens_fts <- fts(
        self$dens_grid,
        self$dens_mat[, which(as.numeric(colnames(self$dens_mat)) < target_time)]
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
        self$get_quant(target_time),
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
    wasserstein_ar = function(target_time, order) {
      data_WAR1 <- WARp(
        self$quant_mat[, which(colnames(self$quant_mat) < target_time)],
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
        self$get_quant(target_time),
        forecast_quant$y,
        self$quant_grid
      )

      list(wasserstein_dist, forecast_quant)
    }
  ),
  private = list()
)

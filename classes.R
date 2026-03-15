library(plyr)
library(Hmisc)
library(fdadensity)
library(ftsa)
library(WRI)
library(parallel)

# Configuration

THREADS <- 16

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

# Takes a sorted x and k to get adaptive bandwidths
# Change minimum_dist to prevent h = 0
knn_bandwidth <- function(x, k, minimum_dist = 0.1) {
  n <- length(x)

  left <- pmax(seq_len(n) - k, 1)
  right <- pmin(seq_len(n) + k, n)

  h <- pmax(
    abs(x - x[left]),
    abs(x[right] - x)
  )

  pmax(h, minimum_dist)
}

# Computes a KDE using the normal kernel from a supplied grid of points
Rcpp::sourceCpp("density_from_grid.cpp")
density_from_grid <- function(data, h, grid, weights = NULL, cutoff = 6) {
  if (is.null(weights)) {
    weights <- rep_len(1, length(data))
  }
  if (length(h) < length(data)) {
    h <- rep_len(h, length(data))
  }
  keep <- weights != 0
  density_from_grid_cpp(data[keep], h, grid, weights[keep], cutoff)
}

library(R6)
DensityTimeSeries <- R6Class(
  "DensityTimeSeries",
  public = list(
    data = NULL,
    data_split = NULL,
    times = NULL,
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

      # Remove 0 or less weighted entries
      self$data <- self$data[self$data$weights > 0, ]

      # Create splits
      self$data_split <- split(self$data, self$data$time)
      self$times <- as.numeric(names(self$data_split))
    },
    # Creates a grid for KDE evaultion based on the quanitiles of all data
    create_dens_grid = function(n) {
      self$dens_grid <- Hmisc::wtd.quantile(
        self$data$x,
        weights = self$data$weights,
        probs = seq(0, 1, length.out = n - 2)
      )
      self$dens_grid <- c(
        min(self$dens_grid) - 3,
        self$dens_grid,
        max(self$dens_grid) + 3
      )
    },
    # A normal density is non-zero everywhere, but that means we must compute
    # several values which are essentially zero.
    # Adding a cutoff limits the number of data points we must iterate over.
    # The actual cutoff "radius" is h * cutoff
    create_dens = function(h, cutoff = 6) {
      densities <- mclapply(
        self$data_split,
        function(xx) {
          density_from_grid(
            xx$x,
            h = h,
            grid = self$dens_grid,
            weights = xx$weights,
            cutoff = cutoff
          )
        },
        mc.cores = THREADS
      )
      self$dens_mat <- do.call(cbind, densities)
    },
    create_dens_knn = function(k, cutoff = 6) {
      densities <- mclapply(
        self$data_split,
        function(xx) {
          density_from_grid(
            xx$x,
            h = knn_bandwidth(xx$x, k),
            grid = self$dens_grid,
            weights = xx$weights,
            cutoff = cutoff
          )
        },
        mc.cores = THREADS
      )
      self$dens_mat <- do.call(cbind, densities)
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
    # KNN selection stuff
    compute_k_grid = function(length.out = 16) {
      sqrt_n_avg <- nrow(self$data) / length(self$times)
      round(seq(sqrt_n_avg / 4, sqrt_n_avg, length.out = length.out))
    },
    compute_knn = function(data_split, k_grid) {
      mclapply(data_split, function(xx) {
        lapply(k_grid, function(k) {
          knn_bandwidth(xx$x, k)
        })
      },
      mc.cores = THREADS
      )
    },
    # Calculates the Wasserstein score for a particular density matrix
    score_density = function(dens_mat, times, start_times) {
      scores <- numeric(length(times) - start_times)

      idx <- 1

      for (t in times[(start_times + 1):length(times)]) {
        result <- self$fda_ar(t, dens_mat = dens_mat)

        scores[idx] <- result[[1]]

        idx <- idx + 1
      }

      mean(scores)
    },
    select_knn_bandwidth = function(k_grid_length,
                                    start_times = 5,
                                    cutoff = 6,
                                    verbose = FALSE) {
      k_grid <- self$compute_k_grid(k_grid_length)
      data_split <- self$data_split
      times <- self$times

      h_all <- self$compute_knn(data_split, k_grid)

      scores <- numeric(length(k_grid))

      scores <- unlist(mclapply(
        seq_along(k_grid),
        function(ii) {
          h_list <- lapply(h_all, `[[`, ii)

          dens <- mclapply(
            seq_along(data_split),
            function(i) {
              xx <- data_split[[i]]
              h <- h_list[[i]]

              density_from_grid(
                xx$x,
                h,
                grid = self$dens_grid,
                weights = xx$weights,
                cutoff = cutoff
              )
            },
            mc.cores = THREADS
          )

          dens_mat <- do.call(cbind, dens)
          colnames(dens_mat) <- self$times

          self$score_density(dens_mat, times, start_times)
        },
        mc.cores = THREADS
      ))

      if (verbose) {
        print(cbind(k = k_grid, wasserstein = scores))
      }

      k_grid[which.min(scores)]
    },
    # Predicts the target year using a model built from all years prior
    # Uses a naive FPCA approach
    fda_ar = function(target_time,
                      transformation = \(x, y) x,
                      inv_transformation = \(x, y) x,
                      dens_grid = self$dens_grid,
                      dens_mat = self$dens_mat) {
      dens_mat_bayes <- transformation(dens_mat, dens_grid)
      dens_fts <- fts(
        dens_grid,
        dens_mat_bayes[, which(as.numeric(colnames(dens_mat)) < target_time)]
      )
      dens_ftsm <- ftsm(dens_fts, order = 3)
      # Not sure if this is AR(1) or some other model
      forecast_dens <- forecast(
        dens_ftsm,
        method = "arima",
        level = 95,
        h = 1
      )
      forecast_bayes <- forecast_dens$mean$y[, 1]

      forecast_pdf <- inv_transformation(forecast_bayes, dens_grid)

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
    },
    # Bayes Space stuff
    bayes_ar = function(target_time, dens_grid = self$dens_grid, dens_mat = self$dens_mat) {
      self$fda_ar(
        target_time,
        transformation = function(dens, grid) {
          log_dens <- log(pmax(dens, 1e-12))
          mean_log <- mean(log_dens)
          log_dens - mean_log
        },
        inv_transformation = function(h, grid) {
          f <- exp(h)
          dx <- diff(grid)
          dx <- c(dx, dx[length(dx)])
          f / sum(f * dx)
        },
        dens_grid = dens_grid,
        dens_mat = dens_mat
      )
    }
  ),
  private = list()
)

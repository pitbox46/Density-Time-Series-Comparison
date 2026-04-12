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
  suppressWarnings(approx(cdf, cdf_grid, quant_grid))
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
        ar_obj <- self$fda_ar(t, dens_mat = dens_mat)

        scores[idx] <- self$wass_dist_ar(ar_obj)

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

          dens <- lapply(
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
            }
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
    # Takes an object from fda_ar
    wass_dist_ar = function(obj) {
      # Sometimes pdf_forecast is less than zero, which is problematic
      wass_dist(
        self$get_quant(obj$target_time),
        suppressWarnings(dens2quantile(
          matrix(ifelse(obj$forecast_pdf < 0, 0, obj$forecast_pdf), nrow = 1),
          obj$forecast_dens$mean$x,
          self$quant_grid
        )),
        self$quant_grid
      )
    },
    # Predicts the target year using a model built from all years prior
    # Uses a naive FPCA approach
    fda_ar = function(target_time,
                      dens_grid = self$dens_grid,
                      dens_mat = self$dens_mat) {
      dens_fts <- fts(
        dens_grid,
        dens_mat[, which(as.numeric(colnames(dens_mat)) < target_time)]
      )
      dens_ftsm <- ftsm(dens_fts, order = 3)
      # Not sure if this is AR(1) or some other model
      forecast_dens <- forecast(
        dens_ftsm,
        method = "arima",
        level = 95,
        h = 1,
        max.p = 1, max.d = 0, max.q = 0
      )
      forecast_pdf <- forecast_dens$mean$y[, 1]

      list(
        target_time = target_time,
        dens_grid = dens_grid,
        dens_mat = dens_mat,
        forecast_pdf = forecast_pdf,
        forecast_dens = forecast_dens
      )
    },
    # Predicts the target year using a model built from all years prior
    # Uses WARp model
    wasserstein_ar = function(target_time, order = 1, dens_grid = self$dens_grid, dens_mat = self$dens_mat) {
      # Fit the WARp model
      data_WAR1 <- WARp(
        self$quant_mat[, which(colnames(self$quant_mat) < target_time)],
        self$quant_grid,
        order
      )

      # Predict on the density grid
      forecast_war <- predict(data_WAR1, dens_grid, dens_grid)

      # Extract the predicted density from the WARp predict output
      forecast_pdf <- as.numeric(forecast_war$pred.pdf)
      forecast_pdf <- c(0, forecast_pdf)

      # Ensure non-negativity and normalize it to integrate to 1
      forecast_pdf <- pmax(forecast_pdf, 1e-12)
      dx <- diff(dens_grid)
      dx <- c(dx, dx[length(dx)])
      forecast_pdf <- forecast_pdf / sum(forecast_pdf * dx)

      # Mock the forecast_dens structure so wass_dist_ar() can successfully
      # access obj$forecast_dens$mean$x
      forecast_dens <- list(
        mean = list(
          x = dens_grid
        )
      )

      list(
        target_time = target_time,
        dens_grid = dens_grid,
        dens_mat = dens_mat,
        forecast_pdf = forecast_pdf,
        forecast_dens = forecast_dens
      )
    },
    # Bayes Space stuff
    bayes_ar = function(target_time, dens_grid = self$dens_grid, dens_mat = self$dens_mat) {
      # Transformation
      log_dens <- log(pmax(dens_mat, 1e-12))
      mean_log <- mean(log_dens)
      dens_mat <- log_dens - mean_log

      ar_obj <- self$fda_ar(
        target_time,
        dens_grid = dens_grid,
        dens_mat = dens_mat
      )

      # Inverse transformation
      f <- exp(ar_obj$forecast_pdf)
      dx <- diff(dens_grid)
      dx <- c(dx, dx[length(dx)])
      ar_obj$forecast_pdf <- f / sum(f * dx)

      ar_obj
    },
    lqd_ar = function(target_time,
                      dens_grid = self$dens_grid,
                      dens_mat = self$dens_mat) {
      idx <- which(as.numeric(colnames(dens_mat)) < target_time)

      # Transform
      lqd_mat <- sapply(idx, function(j) {
        dens_at_q <- approx(
          dens_grid,
          pmax(dens_mat[, j], 1e-12),
          xout = self$quant_mat[, j],
          rule = 2
        )$y

        -log(pmax(dens_at_q, 1e-12))
      })

      colnames(lqd_mat) <- colnames(dens_mat)[idx]

      # ---- FPCA ----
      ar_obj <- self$fda_ar(
        target_time,
        dens_grid = self$quant_grid,
        dens_mat = lqd_mat
      )

      g <- ar_obj$forecast_pdf

      # Inverse transform
      q <- exp(g)
      dt <- diff(self$quant_grid)

      Q <- cumsum(c(0, dt * (head(q, -1) + tail(q, -1)) / 2))

      # shift using median (anchor)
      anchor_idx <- which.min(abs(self$quant_grid - 0.5))
      Q <- Q - Q[anchor_idx] + self$quant_mat[anchor_idx, idx[length(idx)]]

      # convert back to density
      F <- approx(Q, self$quant_grid, xout = dens_grid, rule = 2)$y
      qF <- approx(self$quant_grid, q, xout = F, rule = 2)$y

      f <- 1 / pmax(qF, 1e-12)

      # normalize
      dx <- diff(dens_grid)
      dx <- c(dx, dx[length(dx)])
      ar_obj$forecast_pdf <- f / sum(f * dx)

      # fix grid for wass_dist_ar
      ar_obj$forecast_dens$mean$x <- dens_grid

      ar_obj
    }
  ),
  private = list()
)

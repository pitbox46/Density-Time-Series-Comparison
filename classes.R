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
  # 1. Compute the squared differences at every grid point
  sq_diff <- (quant1 - quant2)^2

  # 2. Calculate the step sizes (dp)
  dp <- diff(quant_grid)

  # 3. Integrate using the Trapezoidal Rule
  # Area of a trapezoid = height * (left_edge + right_edge) / 2
  integral <- sum(dp * (head(sq_diff, -1) + tail(sq_diff, -1)) / 2)

  sqrt(integral)
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
    ftsm_order = 3,
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
    # Creates a grid for KDE evaluation based on the quantiles of all data,
    # enforcing a maximum distance between adjacent points.
    create_dens_grid = function(n, sd_multiplier = 4, max_gap = NULL) {
      base_grid <- Hmisc::wtd.quantile(
        self$data$x,
        weights = self$data$weights,
        probs = seq(0, 1, length.out = n)
      )

      # Default max_gap to a fraction of standard deviation if not provided
      if (is.null(max_gap)) {
        max_gap <- sd(self$data$x) / 4
      }

      # 1. Fill sparse gaps in the quantile grid
      diffs <- diff(base_grid)
      inserts <- pmax(1, ceiling(diffs / max_gap))

      base_grid <- unlist(lapply(1:(length(base_grid) - 1), function(i) {
        # Generate sequence and drop the last point to avoid duplicates
        seq(base_grid[i], base_grid[i + 1], length.out = inserts[i] + 1)[-(inserts[i] + 1)]
      }))

      # Append the final boundary point
      base_grid <- c(base_grid, max(self$data$x))

      # 2. Tail padding
      pad_amt <- sd_multiplier * sd(self$data$x)

      # Use the smaller of the average step or max_gap for the tails
      tail_step <- max_gap
      pad_amt <- max(pad_amt, tail_step)

      left_tail <- seq(min(base_grid) - pad_amt, min(base_grid) - tail_step, by = tail_step)
      right_tail <- seq(max(base_grid) + tail_step, max(base_grid) + pad_amt, by = tail_step)

      self$dens_grid <- unname(c(left_tail, base_grid, right_tail))
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
    compute_k_grid = function(length.out = 16, lower_mul = 1 / 4, upper_mul = 2) {
      sqrt_n_avg <- sqrt(nrow(self$data)) / length(self$times)
      round(seq(sqrt_n_avg * lower_mul, sqrt_n_avg * upper_mul, length.out = length.out))
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
                                    verbose = FALSE,
                                    ...) {
      k_grid <- self$compute_k_grid(k_grid_length, ...)
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
          matrix(pmax(0, obj$forecast_pdf), nrow = 1),
          obj$forecast_dens$mean$x,
          self$quant_grid,
          useSplines = FALSE
        )),
        self$quant_grid
      )
    },
    # Predicts the target year using a model built from all years prior
    # Uses a naive FPCA approach
    fda_ar = function(target_time,
                      dens_grid = self$dens_grid,
                      dens_mat = self$dens_mat,
                      ensure_positive = FALSE) {
      dens_fts <- fts(
        dens_grid,
        dens_mat[, which(as.numeric(colnames(dens_mat)) < target_time)]
      )
      dens_ftsm <- ftsm(dens_fts, order = self$ftsm_order)
      # Uses auto.arima() under the hood
      forecast_dens <- forecast(
        dens_ftsm,
        method = "arima",
        level = 95,
        h = 1,
        max.p = 1, max.d = 2, max.q = 0
      )
      forecast_pdf <- forecast_dens$mean$y[, 1]

      if (ensure_positive) {
        # Ensures that the output PDF is positive
        forecast_pdf <- pmax(forecast_pdf, 1e-12)
      }

      list(
        target_time = target_time,
        dens_grid = dens_grid,
        dens_mat = dens_mat,
        forecast_pdf = forecast_pdf,
        forecast_dens = forecast_dens,
        forecast_raw = forecast_pdf
      )
    },
    # Predicts the target year using a model built from all years prior
    # Uses WARp model
    wasserstein_ar = function(target_time,
                              order = 1,
                              dens_grid = self$dens_grid,
                              dens_mat = self$dens_mat) {
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
    bayes_ar = function(target_time,
                        dens_grid = self$dens_grid,
                        dens_mat = self$dens_mat) {
      # Transformation: Proper Centered Log-Ratio (CLR) per column
      log_dens <- log(pmax(dens_mat, 1e-12))

      dx <- diff(dens_grid)
      dx <- c(dx, dx[length(dx)])
      domain_length <- sum(dx)

      # Calculate the integral of the log-density for EACH time step
      col_integrals <- apply(log_dens, 2, function(y) sum(y * dx))
      col_means <- col_integrals / domain_length

      # Sweep subtracts the corresponding mean from each column individually
      dens_mat <- sweep(log_dens, 2, col_means, "-")

      ar_obj <- self$fda_ar(
        target_time,
        dens_grid = dens_grid,
        dens_mat = dens_mat
      )

      g <- ar_obj$forecast_pdf

      # Inverse transformation: Log-Sum-Exp trick for Overflow Protection
      # Shift the log-density by its maximum before exponentiating.
      g_shifted <- g - max(g)
      f <- exp(g_shifted)

      # normalize
      ar_obj$forecast_pdf <- f / sum(f * dx)

      ar_obj
    },
    lqd_ar = function(target_time,
                      dens_grid = self$dens_grid,
                      dens_mat = self$dens_mat) {
      idx <- which(as.numeric(colnames(dens_mat)) < target_time)

      # ---------------------------------------------------------
      # DYNAMIC KDE QUANTILES: Calculate locally for LQD only
      # ---------------------------------------------------------
      dx <- diff(dens_grid)
      dx <- c(dx, dx[length(dx)])

      kde_quant_mat <- sapply(idx, function(j) {
        pdf <- pmax(dens_mat[, j], 1e-12)
        cdf <- cumsum(pdf * dx)
        cdf <- cdf / tail(cdf, 1) # Normalize to 1
        approx(x = cdf, y = dens_grid, xout = self$quant_grid, rule = 2)$y
      })
      # Note: kde_quant_mat columns now map 1-to-1 with the `idx` vector

      # Transform
      lqd_mat <- sapply(seq_along(idx), function(k) {
        j <- idx[k] # The actual column index in dens_mat

        dens_at_q <- approx(
          dens_grid,
          pmax(dens_mat[, j], 1e-12),
          xout = kde_quant_mat[, k], # Use the dynamically generated KDE quantile
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
      # Use the most recent KDE quantile (the last column) to anchor
      Q <- Q - Q[anchor_idx] + kde_quant_mat[anchor_idx, ncol(kde_quant_mat)]

      # convert back to density
      F <- approx(Q, self$quant_grid, xout = dens_grid, rule = 2)$y
      qF <- approx(self$quant_grid, q, xout = F, rule = 2)$y

      f <- 1 / pmax(qF, 1e-12)

      # Enforce zero density strictly outside of Q's bounds to prevent extrapolation flatlines
      f[dens_grid < min(Q) | dens_grid > max(Q)] <- 1e-12

      # normalize
      ar_obj$forecast_pdf <- f / sum(f * dx)

      # fix grid for wass_dist_ar
      ar_obj$forecast_dens$mean$x <- dens_grid
      ar_obj$dens_grid <- dens_grid

      ar_obj
    }
  ),
  private = list()
)

options(error = function() traceback(3))

# Load data file or get the data from the CPS API
import_cps <- function(load_file, years, data_file = "./data.RData") {
  # If we have one saved, we use the given file instead
  if (load_file) {
    load(data_file)
    return(cps)
  }

  # Get CPI data
  cpi_data <- data.frame(
    Year = years,
    CPI = c(163.0, 166.6, 172.2, 177.1, 179.9, 184.0, 188.9, 195.3, 201.6, 207.3, 215.3, 214.5, 218.1, 224.9, 229.6, 233.0, 236.7, 237.0, 240.0, 245.1, 251.1, 255.7, 258.8, 271.0, 292.7, 304.7, 314.4)
  )
  cpi_data$Inflation <- cpi_data$CPI[nrow(cpi_data)] / cpi_data$CPI

  # Get CPS data and combine with CPI data
  library(cpsR)
  cps <- data.frame(
    ptotval = numeric(),
    a_ernlwt = numeric(),
    year = numeric()
  )
  cols_select <- c(colnames(cps)[-ncol(cps)])

  for (year in years) {
    cps_year <- get_asec(year, cols_select) # Requires API key. See cpsR documentation
    cps_year$ptotval <- cps_year$ptotval * cpi_data[cpi_data$Year == year, "Inflation"]
    cps_year$Year <- year
    cps <- rbind(cps, cps_year)
  }

  save(cps, file = data_file)
  return(cps)
}

setwd("~/School/STAT-S799/IncomesDataAnalysis/")
cps <- import_cps(load_file = TRUE, years = NA)

cps$ptotval <- cps$ptotval + runif(length(cps$ptotval), -10, 10)

colnames(cps) <- c("x", "weights", "time")

# Initialize
create_analaysis_obj <- function(data) {
  source("classes.R")
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

  # KNN bandwidths
  k <- analysis_obj$select_knn_bandwidth(16, verbose = TRUE)

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

analysis_obj <- create_analaysis_obj(cps)
test_all_models(2010:2024, analysis_obj)

models <- list(
  "FDA" = analysis_obj$fda_ar,
  "Bayes" = analysis_obj$bayes_ar,
  "LQD" = analysis_obj$lqd_ar,
  "Wasserstein" = function(t) analysis_obj$wasserstein_ar(t, order = 1)
)

# Generate static plot for timestep 40
plot_all_models_vs_actual(analysis_obj, 2020, models)

# Generate animation for timesteps 20 through 40
anim <- animate_all_models(analysis_obj, 2010:2024, models, log_scale = TRUE)
# Render the animation (adjust frames/fps as needed for execution speed)
animate(anim, nframes = 50, fps = 5, width = 1000, height = 1000)

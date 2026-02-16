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

setwd("~/School/Now/STAT-S799/IncomesDataAnalysis/")
cps <- import_cps(load_file = TRUE, years = NA)

# Create Densities

library(plyr)

data_densities <- cps |>
  daply(
    .(Year),
    function(x) {
      dens <- density(
        x$ptotval,
        weights = x$a_ernlwt,
        bw = 50000,
        from = -200000,
        to = 3000000,
        n = 1024
      )
      densities_grid <<- dens$x
      dens$y
    }
  ) |>
  t()

# Create Quantiles

library(Hmisc)
library(fdadensity)

quantile_grid <- seq(0, 1, 0.005)
data_quantiles <- cps |>
  daply(
    .(Year),
    \(x) wtd.quantile(x$ptotval, weights = x$a_ernlwt, probs = quantile_grid)
  ) |>
  t()

get_data_by_year <- function(target_year) {
  subset(cps, Year == target_year)
}

get_quantiles_by_year <- function(target_year_data, quantile_grid) {
  wtd.quantile(
    target_year_data$ptotval,
    weights = target_year_data$a_ernlwt,
    probs = quantile_grid
  )
}

# Wasserstein Distance for model comparison

wass_dist <- function(quantile1, quantile2, quantile_grid) {
  ret <- 0
  for (i in 1:(length(quantile_grid) - 1)) {
    step <- quantile_grid[i + 1] - quantile_grid[i]

    ret <- ret + (quantile2[i] - quantile1[i])^2 * step
  }
  sqrt(ret)
}

# FDA AR1 Model

library(ftsa)

fda_ar <- function(target_year, densities_grid, quantile_grid) {
  dens_fts <- fts(
    densities_grid,
    data_densities[, which(colnames(data_densities) < target_year)]
  )
  dens_ftsm <- ftsm(dens_fts, order = 3)
  # Not sure if this is AR(1) or some other model
  dens_forecast <- forecast(
    dens_ftsm,
    method = "arima",
    level = 95,
    h = 1
  )
  pdf_forecast <- dens_forecast$mean$y[, 1]

  target_year_data <- get_data_by_year(target_year)
  target_year_quantiles <- get_quantiles_by_year(target_year_data, quantile_grid)

  wasserstein_dist <- wass_dist(
    target_year_quantiles,
    dens2quantile(
      matrix(ifelse(pdf_forecast < 0, 0, pdf_forecast), nrow = 1),
      dens_forecast$mean$x,
      quantile_grid
    ),
    quantile_grid
  )

  list(
    target_year_data = target_year_data,
    target_year = target_year,
    target_year_quantiles = target_year_quantiles,
    wasserstein_dist = wasserstein_dist
  )
}

FDA_AR_fits <- lapply(2010:2024, fda_ar, densities_grid = densities_grid, quantile_grid = quantile_grid)
FDA_AR_distances <- sapply(FDA_AR_fits, function(x) x$wasserstein_dist)

# Wasserstein AR model

library(WRI)

# Computes a prediction for the target year using WAR1
wasserstein_ar <- function(target_year, densities_grid, quantile_grid, order) {
  data_WAR1 <- WARp(data_quantiles, quantile_grid, order)
  war_pred <- predict(data_WAR1, densities_grid, densities_grid)

  target_year_data <- get_data_by_year(target_year)
  target_year_quantiles <- get_quantiles_by_year(target_year_data, quantile_grid)
  wasserstein_dist <- wass_dist(war_pred$dSup, target_year_quantiles, quantile_grid)

  list(
    WARp_obj = data_WAR1,
    WARp_pred = war_pred,
    target_year_data = target_year_data,
    target_year = target_year,
    target_year_quantiles = target_year_quantiles,
    wasserstein_dist = wasserstein_dist
  )
}

WAR_fits <- lapply(2010:2024, wasserstein_ar, quantile_grid = quantile_grid, order = 1)
WAR_distances <- sapply(WAR_fits, function(x) x$wasserstein_dist)

WAR_fits <- lapply(2010:2024, wasserstein_ar, quantile_grid = quantile_grid, order = 2)
WAR_distances <- sapply(WAR_fits, function(x) x$wasserstein_dist)

# Plots

plot(log(target_year_quantiles), quantile_grid, type = "l", lwd = 2)
lines(log(war_pred$dSup), war_pred$pred.cdf, col = "red")

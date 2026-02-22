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

setwd("~/School/Now/STAT-S799/IncomesDataAnalysis/")
cps <- import_cps(load_file = TRUE, years = NA)

# Initialize

analysis_obj <- DensityTimeSeries$new(cps$ptotval, cps$Year, cps$a_ernlwt)

analysis_obj$create_dens(
  bw = 50000,
  from = -200000,
  to = 3000000,
  n = 4096
)
analysis_obj$create_quants(seq(0, 1, 0.001))

# FDA AR1 Model

analysis_obj$fda_ar(2010)
FDA_AR_fits <- lapply(2010:2024, analysis_obj$fda_ar)
FDA_AR_distances <- sapply(FDA_AR_fits, function(x) x$wasserstein_dist)
mean(FDA_AR_distances)

# Wasserstein AR model

WAR_fits <- lapply(2010:2024, wasserstein_ar, densities_grid = densities_grid, quantile_grid = quantile_grid, order = 1)
WAR_distances <- sapply(WAR_fits, function(x) x$wasserstein_dist)
mean(WAR_distances)

WAR_fits <- lapply(2010:2024, wasserstein_ar, densities_grid = densities_grid, quantile_grid = quantile_grid, order = 2)
WAR_distances <- sapply(WAR_fits, function(x) x$wasserstein_dist)
mean(WAR_distances)

WAR_fits <- lapply(2010:2024, wasserstein_ar, densities_grid = densities_grid, quantile_grid = quantile_grid, order = 3)
WAR_distances <- sapply(WAR_fits, function(x) x$wasserstein_dist)
mean(WAR_distances)

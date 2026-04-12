library(ggplot2)
library(reshape2)
library(viridis)

# 1. Visualize the evolution of the density distributions over time
plot_density_evolution <- function(analysis_obj, title_suffix = "") {
  dens_mat <- analysis_obj$dens_mat
  grid <- analysis_obj$dens_grid

  # Ensure columns are properly identified as time steps
  df <- melt(dens_mat)
  colnames(df) <- c("GridIndex", "Time", "Density")
  df$x <- grid[df$GridIndex]

  ggplot(df, aes(x = x, y = Density, group = Time, color = Time)) +
    geom_line(alpha = 0.7) +
    scale_color_viridis(option = "plasma") +
    theme_minimal() +
    labs(
      title = paste("Evolution of Density Over Time", title_suffix),
      x = "Value", y = "Density"
    )
}

# 2. Compare Actual vs Predicted Density for a specific time step
plot_actual_vs_predicted <- function(analysis_obj, target_time, model_func, model_name = "Model", log_scale = FALSE) {
  # Get the forecast object
  ar_obj <- model_func(target_time)

  if (log_scale) {
    grid <- log(ar_obj$dens_grid)
  } else {
    grid <- ar_obj$dens_grid
  }
  predicted_pdf <- ar_obj$forecast_pdf

  # Get actual density from the matrix (requires column names to be time values)
  actual_pdf <- analysis_obj$dens_mat[, as.character(target_time)]

  df <- data.frame(
    x = grid,
    Actual = actual_pdf,
    Predicted = predicted_pdf
  )

  df_melt <- melt(df, id.vars = "x", variable.name = "Type", value.name = "Density")

  ggplot(df_melt, aes(x = x, y = Density, color = Type, linetype = Type)) +
    geom_line(linewidth = 1) +
    scale_color_manual(values = c("Actual" = "black", "Predicted" = "blue")) +
    theme_minimal() +
    labs(
      title = sprintf("%s: Actual vs Predicted Density (t = %s)", model_name, target_time),
      x = ifelse(log_scale, "Log Value", "Value"), y = "Density"
    )
}

# 3. Plot the residual difference (Actual - Predicted)
plot_prediction_error <- function(analysis_obj, target_time, model_func, model_name = "Model") {
  ar_obj <- model_func(target_time)
  grid <- ar_obj$dens_grid
  actual_pdf <- analysis_obj$dens_mat[, as.character(target_time)]

  # Calculate residuals
  residuals <- actual_pdf - ar_obj$forecast_pdf

  df <- data.frame(x = grid, Error = residuals)

  ggplot(df, aes(x = x, y = Error)) +
    geom_line(color = "red", linewidth = 0.8) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
    theme_minimal() +
    labs(
      title = sprintf("%s: Density Prediction Error (t = %s)", model_name, target_time),
      x = "Value", y = "Actual - Predicted Density"
    )
}

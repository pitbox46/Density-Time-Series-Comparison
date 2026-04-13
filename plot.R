library(ggplot2)
library(reshape2)
library(viridis)

# 1. Visualize the evolution of the quantile functions over time
plot_quantile_evolution <- function(analysis_obj, title_suffix = "") {
  quant_mat <- analysis_obj$quant_mat
  prob_grid <- analysis_obj$quant_grid

  # melt() on a matrix returns Var1 (row index), Var2 (col name), and value
  df <- melt(quant_mat)
  colnames(df) <- c("GridIndex", "Time", "Value")

  # Replace row indices with actual probabilities (0 to 1)
  df$Probability <- prob_grid[df$GridIndex]

  ggplot(df, aes(x = Probability, y = Value, group = Time, color = Time)) +
    geom_line(alpha = 0.7) +
    scale_color_viridis(option = "plasma") +
    theme_minimal() +
    labs(
      title = paste("Evolution of Quantile Functions", title_suffix),
      x = "Probability", y = "Domain Value"
    )
}

# 2. Compare Actual vs Predicted Quantiles for a specific time step
plot_actual_vs_predicted <- function(analysis_obj, target_time, model_func, model_name = "Model", log_scale = FALSE) {
  # Get actual quantiles from the matrix
  actual_quant <- analysis_obj$quant_mat[, as.character(target_time)]
  prob_grid <- analysis_obj$quant_grid

  # Get the forecast object (PDF)
  ar_obj <- model_func(target_time)
  grid <- ar_obj$dens_grid
  predicted_pdf <- ar_obj$forecast_pdf

  # Numerically integrate PDF to get CDF
  dx <- diff(grid)
  dx <- c(dx, dx[length(dx)])
  predicted_cdf <- cumsum(predicted_pdf * dx)
  predicted_cdf <- predicted_cdf / tail(predicted_cdf, 1) # Ensure max is exactly 1

  # Interpolate to find predicted quantiles at the exact probability grid points
  predicted_quant <- approx(x = predicted_cdf, y = grid, xout = prob_grid, rule = 2)$y

  if (log_scale) {
    actual_quant <- log(actual_quant)
    predicted_quant <- log(predicted_quant)
  }

  df <- data.frame(
    Probability = prob_grid,
    Actual = actual_quant,
    Predicted = predicted_quant
  )

  df_melt <- melt(df, id.vars = "Probability", variable.name = "Type", value.name = "Value")

  ggplot(df_melt, aes(x = Probability, y = Value, color = Type, linetype = Type)) +
    geom_line(linewidth = 1) +
    scale_color_manual(values = c("Actual" = "black", "Predicted" = "blue")) +
    theme_minimal() +
    labs(
      title = sprintf("%s: Actual vs Predicted Quantiles (t = %s)", model_name, target_time),
      x = "Probability", y = ifelse(log_scale, "Log Value", "Value")
    )
}

# 3. Plot the residual difference (Actual - Predicted)
plot_prediction_error <- function(analysis_obj, target_time, model_func, model_name = "Model") {
  actual_quant <- analysis_obj$quant_mat[, as.character(target_time)]
  prob_grid <- analysis_obj$quant_grid

  ar_obj <- model_func(target_time)
  grid <- ar_obj$dens_grid
  predicted_pdf <- ar_obj$forecast_pdf

  # Convert predicted PDF to CDF
  dx <- diff(grid)
  dx <- c(dx, dx[length(dx)])
  predicted_cdf <- cumsum(predicted_pdf * dx)
  predicted_cdf <- predicted_cdf / tail(predicted_cdf, 1)

  # Extract predicted quantiles
  predicted_quant <- approx(x = predicted_cdf, y = grid, xout = prob_grid, rule = 2)$y

  # Calculate residuals in Wasserstein space (difference in quantiles)
  residuals <- actual_quant - predicted_quant

  df <- data.frame(Probability = prob_grid, Error = residuals)

  ggplot(df, aes(x = Probability, y = Error)) +
    geom_line(color = "red", linewidth = 0.8) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
    theme_minimal() +
    labs(
      title = sprintf("%s: Quantile Prediction Error (t = %s)", model_name, target_time),
      x = "Probability", y = "Actual - Predicted Quantile"
    )
}

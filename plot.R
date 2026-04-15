library(ggplot2)
library(reshape2)
library(viridis)
library(gganimate)

# --- Helper Function ---
# Extracts and interpolates the predicted quantiles for a given model
get_predicted_quantiles <- function(analysis_obj, target_time, model_func) {
  ar_obj <- model_func(target_time)

  grid <- ar_obj$forecast_dens$mean$x
  predicted_pdf <- ar_obj$forecast_pdf

  # Convert predicted PDF to CDF
  dx <- diff(grid)
  dx <- c(dx, dx[length(dx)])
  predicted_cdf <- cumsum(predicted_pdf * dx)
  predicted_cdf <- predicted_cdf / tail(predicted_cdf, 1)

  # Extract predicted quantiles
  approx(x = predicted_cdf, y = grid, xout = analysis_obj$quant_grid, rule = 2)$y
}

# --- 1. Static Plot: All Models vs Actual ---
plot_all_models_vs_actual <- function(analysis_obj, target_time, models_list, log_scale = FALSE) {
  prob_grid <- analysis_obj$quant_grid
  actual_quant <- analysis_obj$quant_mat[, as.character(target_time)]

  # Initialize data frame with actual data
  df <- data.frame(Probability = prob_grid, Actual = actual_quant)

  # Append predictions from each model
  for (model_name in names(models_list)) {
    df[[model_name]] <- get_predicted_quantiles(analysis_obj, target_time, models_list[[model_name]])
  }

  if (log_scale) {
    df[, -1] <- log(pmax(df[, -1], 1e-12))
  }

  df_melt <- melt(df, id.vars = "Probability", variable.name = "Model", value.name = "Value")

  ggplot(df_melt, aes(x = Probability, y = Value, color = Model, group = Model)) +
    geom_line(aes(linewidth = Model == "Actual", linetype = Model == "Actual")) +
    scale_linewidth_manual(values = c("TRUE" = 1.2, "FALSE" = 0.8), guide = "none") +
    scale_linetype_manual(values = c("TRUE" = "solid", "FALSE" = "dashed"), guide = "none") +
    scale_color_manual(values = c(
      "Actual" = "black",
      "FDA" = "#E69F00",
      "Bayes" = "#56B4E9",
      "LQD" = "#009E73",
      "Wasserstein" = "#CC79A7"
    )) +
    theme_minimal() +
    labs(
      title = sprintf("Model Comparison (t = %s)", target_time),
      x = "Probability",
      y = ifelse(log_scale, "Log Value", "Value")
    )
}

# --- 2. Animated Plot: Evolution Across Time ---
animate_all_models <- function(analysis_obj, times, models_list, log_scale = FALSE) {
  all_dfs <- list()

  # Build a combined dataframe for all specified times
  for (t in times) {
    actual_quant <- analysis_obj$quant_mat[, as.character(t)]
    df_t <- data.frame(Time = t, Probability = analysis_obj$quant_grid, Actual = actual_quant)

    for (model_name in names(models_list)) {
      df_t[[model_name]] <- get_predicted_quantiles(analysis_obj, t, models_list[[model_name]])
    }

    if (log_scale) {
      cols_to_log <- !(names(df_t) %in% c("Time", "Probability"))
      df_t[, cols_to_log] <- log(df_t[, cols_to_log])
    }

    all_dfs[[as.character(t)]] <- df_t
  }

  final_df <- do.call(rbind, all_dfs)
  df_melt <- melt(final_df, id.vars = c("Time", "Probability"), variable.name = "Model", value.name = "Value")

  # Create the base plot
  p <- ggplot(df_melt, aes(x = Probability, y = Value, color = Model, group = Model)) +
    geom_line(aes(linewidth = Model == "Actual", linetype = Model == "Actual")) +
    scale_linewidth_manual(values = c("TRUE" = 1.2, "FALSE" = 0.8), guide = "none") +
    scale_linetype_manual(values = c("TRUE" = "solid", "FALSE" = "dashed"), guide = "none") +
    scale_color_manual(values = c(
      "Actual" = "black",
      "FDA" = "#E69F00",
      "Bayes" = "#56B4E9",
      "LQD" = "#009E73",
      "Wasserstein" = "#CC79A7"
    )) +
    theme_minimal() +
    labs(
      title = "Model Comparison (t = {frame_time})",
      x = "Probability",
      y = ifelse(log_scale, "Log Value", "Value")
    ) +
    transition_time(Time) +
    ease_aes("linear")

  return(p)
}

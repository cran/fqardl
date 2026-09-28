# =============================================================================
# Visualization Functions for FQARDL
# Publication-ready plots
# =============================================================================

#' Plot FQARDL Results
#'
#' @description
#' Creates various diagnostic and result plots for FQARDL models.
#'
#' @param x An object of class "fqardl"
#' @param type Type of plot: "coefficients", "multipliers", "3d", "heatmap", "residuals"
#' @param variable Variable name for coefficient plots
#' @param ... Additional arguments passed to plotting functions
#'
#' @return A ggplot object or plotly object for 3D plots
#'
#' @export
plot.fqardl <- function(x, type = c("coefficients", "multipliers", "3d", 
                                     "heatmap", "residuals"), 
                        variable = NULL, ...) {
  
  type <- match.arg(type)
  
  switch(type,
    "coefficients" = plot_coefficients(x, variable),
    "multipliers" = plot_multipliers(x),
    "3d" = plot_3d_surface(x, variable),
    "heatmap" = plot_heatmap(x),
    "residuals" = plot_residuals(x)
  )
}


#' Plot Coefficients Across Quantiles
#'
#' @param obj FQARDL object
#' @param variable Variable name (NULL for all)
#'
#' @return ggplot object
#'
#' @keywords internal
plot_coefficients <- function(obj, variable = NULL) {tau <- obj$tau
  
  # Extract coefficients
  coef_list <- lapply(seq_along(obj$qardl_results), function(i) {
    res <- obj$qardl_results[[i]]
    data.frame(
      tau = res$tau,
      variable = names(res$coefficients),
      estimate = res$coefficients,
      se = res$std_errors,
      stringsAsFactors = FALSE
    )
  })
  
  coef_df <- do.call(rbind, coef_list)
  
  # Filter variable if specified
  if (!is.null(variable)) {
    coef_df <- coef_df[grep(variable, coef_df$variable), ]
  }
  
  # Remove intercept and Fourier terms for cleaner plot
  coef_df <- coef_df[!grepl("^intercept$|^sin_|^cos_", coef_df$variable), ]
  
  # Create plot
  p <- ggplot(coef_df, aes(x = tau, y = estimate, group = variable)) +
    geom_ribbon(aes(ymin = estimate - 1.96 * se, 
                    ymax = estimate + 1.96 * se,
                    fill = variable), alpha = 0.2) +
    geom_line(aes(color = variable), linewidth = 1) +
    geom_point(aes(color = variable), size = 2) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    facet_wrap(~ variable, scales = "free_y") +
    labs(
      title = "Coefficient Estimates Across Quantiles",
      subtitle = "With 95% Confidence Intervals",
      x = expression(tau~"(Quantile)"),
      y = "Coefficient Estimate"
    ) +
    theme_minimal() +
    theme(
      legend.position = "bottom",
      strip.background = element_rect(fill = "gray90"),
      strip.text = element_text(face = "bold")
    )
  
  return(p)
}


#' Plot Long-run and Short-run Multipliers
#'
#' @param obj FQARDL object
#'
#' @return ggplot object
#'
#' @keywords internal
plot_multipliers <- function(obj) {tau <- obj$tau
  
  # Long-run multipliers
  lr_df <- as.data.frame(obj$long_run)
  lr_df$tau <- tau
  lr_df$type <- "Long-run"
  
  lr_long <- pivot_longer(lr_df, 
                          cols = -c(tau, type),
                          names_to = "variable",
                          values_to = "multiplier")
  
  # Short-run multipliers
  sr_df <- as.data.frame(obj$short_run)
  sr_df$tau <- tau
  sr_df$type <- "Short-run"
  
  sr_long <- pivot_longer(sr_df,
                          cols = -c(tau, type),
                          names_to = "variable",
                          values_to = "multiplier")
  
  # Combine
  mult_df <- rbind(lr_long, sr_long)
  
  # Plot
  p <- ggplot(mult_df, aes(x = tau, y = multiplier, color = variable)) +
    geom_line(linewidth = 1) +
    geom_point(size = 2) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    facet_wrap(~ type, scales = "free_y") +
    labs(
      title = "Dynamic Multipliers Across Quantiles",
      subtitle = "Long-run and Short-run Effects",
      x = expression(tau~"(Quantile)"),
      y = "Multiplier"
    ) +
    scale_color_brewer(palette = "Set1") +
    theme_minimal() +
    theme(
      legend.position = "bottom",
      strip.background = element_rect(fill = "gray90"),
      strip.text = element_text(face = "bold", size = 12)
    )
  
  return(p)
}


#' Plot 3D Surface of Coefficients
#'
#' @param obj FQARDL object
#' @param variable Variable to plot
#'
#' @return plotly object
#'
#' @keywords internal
plot_3d_surface <- function(obj, variable = NULL) {
  
  if (!requireNamespace("plotly", quietly = TRUE)) {
    stop("Package 'plotly' is required for 3D plots")
  }
  tau <- obj$tau
  
  # If variable not specified, use first independent variable
  if (is.null(variable)) {
    variable <- obj$x_names[1]
  }
  
  # Extract coefficients for the variable across quantiles
  coef_name <- paste0("x1_lag1")  # Simplified - adjust based on variable index
  
  coefs <- sapply(obj$qardl_results, function(x) {
    idx <- grep(coef_name, names(x$coefficients))[1]
    if (length(idx) > 0) x$coefficients[idx] else NA
  })
  
  # For a true 3D surface, we need another dimension
  # Use time or simulate different scenarios
  time_points <- 1:obj$n
  
  # Create matrix for surface
  # This is a simplified version - full implementation would use bootstrap
  z_matrix <- outer(coefs, seq(-2, 2, length.out = 20), function(c, m) c * m)
  
  fig <- plotly::plot_ly(
    x = tau,
    y = seq(-2, 2, length.out = 20),
    z = t(z_matrix)
  )
  fig <- plotly::add_surface(fig,
      colorscale = "Viridis",
      contours = list(
        z = list(show = TRUE, usecolormap = TRUE, highlightcolor = "#ff0000", project = list(z = TRUE))
      )
    )
  fig <- plotly::layout(fig,
      title = paste("3D Surface:", variable),
      scene = list(
        xaxis = list(title = "Quantile (tau)"),
        yaxis = list(title = "Effect Size"),
        zaxis = list(title = "Coefficient")
      )
    )
  
  return(fig)
}


#' Plot Heatmap of Coefficients
#'
#' @param obj FQARDL object
#'
#' @return ggplot object
#'
#' @keywords internal
plot_heatmap <- function(obj) {# Create coefficient matrix
  coef_matrix <- t(sapply(obj$qardl_results, function(x) {
    x$coefficients[grep("_lag1$", names(x$coefficients))]
  }))
  
  rownames(coef_matrix) <- paste0("tau_", obj$tau)
  
  # Convert to long format
  coef_df <- as.data.frame(coef_matrix)
  coef_df$tau <- obj$tau
  
  coef_long <- pivot_longer(coef_df,
                            cols = -tau,
                            names_to = "variable",
                            values_to = "coefficient")
  
  # Create heatmap
  p <- ggplot(coef_long, aes(x = variable, y = factor(tau), fill = coefficient)) +
    geom_tile(color = "white") +
    geom_text(aes(label = round(coefficient, 3)), size = 3) +
    scale_fill_gradient2(
      low = "#2166AC", mid = "white", high = "#B2182B",
      midpoint = 0,
      name = "Coefficient"
    ) +
    labs(
      title = "Coefficient Heatmap Across Quantiles",
      x = "Variable",
      y = expression(tau~"(Quantile)")
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid = element_blank()
    )
  
  return(p)
}


#' Plot Residual Diagnostics
#'
#' @param obj FQARDL object
#'
#' @return ggplot object
#'
#' @keywords internal
plot_residuals <- function(obj) {# Use median quantile
  median_idx <- which(obj$tau == 0.5)
  if (length(median_idx) == 0) median_idx <- 1
  
  res <- obj$qardl_results[[median_idx]]
  residuals <- res$residuals
  fitted <- res$fitted
  
  # Create data frame
  diag_df <- data.frame(
    fitted = fitted,
    residuals = residuals,
    index = 1:length(residuals)
  )
  
  # Residuals vs Fitted
  p1 <- ggplot(diag_df, aes(x = fitted, y = residuals)) +
    geom_point(alpha = 0.5) +
    geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
    geom_smooth(method = "loess", se = FALSE, color = "blue") +
    labs(title = "Residuals vs Fitted", x = "Fitted Values", y = "Residuals") +
    theme_minimal()
  
  # Residuals over time
  p2 <- ggplot(diag_df, aes(x = index, y = residuals)) +
    geom_line(alpha = 0.7) +
    geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
    labs(title = "Residuals Over Time", x = "Observation", y = "Residuals") +
    theme_minimal()
  
  # Histogram
  p3 <- ggplot(diag_df, aes(x = residuals)) +
    geom_histogram(aes(y = after_stat(density)), bins = 30, fill = "steelblue", alpha = 0.7) +
    geom_density(color = "red", linewidth = 1) +
    labs(title = "Residual Distribution", x = "Residuals", y = "Density") +
    theme_minimal()
  
  # Q-Q plot
  p4 <- ggplot(diag_df, aes(sample = residuals)) +
    stat_qq() +
    stat_qq_line(color = "red") +
    labs(title = "Q-Q Plot", x = "Theoretical Quantiles", y = "Sample Quantiles") +
    theme_minimal()
  
  # Combine
  gridExtra::grid.arrange(p1, p2, p3, p4, ncol = 2)
}


#' Plot Persistence Profile
#'
#' @description
#' Plots the persistence profile showing the adjustment path
#' to long-run equilibrium after a shock.
#'
#' @param obj FQARDL object
#' @param horizons Number of periods for persistence profile
#'
#' @return ggplot object
#'
#' @export
plot_persistence <- function(obj, horizons = 20) {
  tau <- obj$tau
  
  # Extract ECT coefficients (speed of adjustment)
  ect <- sapply(obj$qardl_results, function(x) x$coefficients["y_lag1"])
  
  # Generate persistence profiles
  profiles <- lapply(seq_along(tau), function(i) {
    phi <- ect[i]
    profile <- (1 + phi)^(0:horizons)
    data.frame(
      horizon = 0:horizons,
      persistence = profile,
      tau = tau[i]
    )
  })
  
  profile_df <- do.call(rbind, profiles)
  profile_df$tau <- factor(profile_df$tau)
  
  # Plot
  p <- ggplot(profile_df, aes(x = horizon, y = persistence, color = tau)) +
    geom_line(linewidth = 1) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "gray50") +
    geom_hline(yintercept = 0, linetype = "solid", color = "black") +
    labs(
      title = "Persistence Profile",
      subtitle = "Adjustment to Long-run Equilibrium After Unit Shock",
      x = "Horizon (Periods)",
      y = "Cumulative Effect",
      color = expression(tau)
    ) +
    scale_color_viridis_d() +
    theme_minimal() +
    theme(legend.position = "right")
  
  return(p)
}

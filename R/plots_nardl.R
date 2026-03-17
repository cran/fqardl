#' =============================================================================
#' Visualization Functions for FNARDL
#' Dynamic Multiplier Plots and Asymmetry Analysis
#' =============================================================================

#' Plot FNARDL Results
#'
#' @param x An object of class "fnardl"
#' @param type Type of plot
#' @param variable Variable to plot
#' @param horizon Horizon for dynamic multipliers
#' @param ... Additional arguments
#'
#' @export
plot.fnardl <- function(x, type = c("asymmetry", "dynamic", "cumulative", "comparison"),
                        variable = NULL, horizon = 20, ...) {
  
  type <- match.arg(type)
  
  if (is.null(variable)) {
    variable <- x$decompose[1]
  }
  
  switch(type,
    "asymmetry" = plot_asymmetry(x, variable),
    "dynamic" = plot_dynamic_multipliers(x, variable, horizon),
    "cumulative" = plot_cumulative_multipliers(x, variable, horizon),
    "comparison" = plot_pos_neg_comparison(x)
  )
}


#' Plot Asymmetry Bar Chart
#'
#' @param obj FNARDL object
#' @param variable Variable to plot
#'
#' @keywords internal
plot_asymmetry <- function(obj, variable = NULL) {if (is.null(variable)) {
    variable <- obj$decompose
  }
  
  # Prepare data
  lr_pos <- obj$multipliers$long_run_positive
  lr_neg <- obj$multipliers$long_run_negative
  sr_pos <- obj$multipliers$short_run_positive
  sr_neg <- obj$multipliers$short_run_negative
  
  plot_data <- data.frame(
    variable = rep(names(lr_pos), 4),
    type = rep(c("Long-run (+)", "Long-run (-)", "Short-run (+)", "Short-run (-)"), 
               each = length(lr_pos)),
    value = c(lr_pos, lr_neg, sr_pos, sr_neg),
    direction = rep(c("Positive", "Negative", "Positive", "Negative"), 
                    each = length(lr_pos)),
    horizon = rep(c("Long-run", "Long-run", "Short-run", "Short-run"), 
                  each = length(lr_pos))
  )
  
  # Filter if variable specified
  if (!is.null(variable) && length(variable) < length(obj$decompose)) {
    plot_data <- plot_data[plot_data$variable %in% variable, ]
  }
  
  # Create plot
  p <- ggplot(plot_data, aes(x = variable, y = value, fill = type)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    facet_wrap(~ horizon, scales = "free_y") +
    scale_fill_manual(values = c(
      "Long-run (+)" = "#2166AC",
      "Long-run (-)" = "#B2182B",
      "Short-run (+)" = "#4393C3",
      "Short-run (-)" = "#D6604D"
    )) +
    labs(
      title = "Asymmetric Multipliers",
      subtitle = "Positive vs Negative Effects",
      x = "Variable",
      y = "Multiplier",
      fill = "Effect"
    ) +
    theme_minimal() +
    theme(
      legend.position = "bottom",
      axis.text.x = element_text(angle = 45, hjust = 1),
      strip.background = element_rect(fill = "gray90"),
      strip.text = element_text(face = "bold")
    )
  
  return(p)
}


#' Plot Dynamic Multipliers
#'
#' @param obj FNARDL object
#' @param variable Variable to plot
#' @param horizon Number of periods
#'
#' @export
plot_dynamic_multipliers <- function(obj, variable, horizon = 20) {
  # Get coefficients
  phi <- obj$multipliers$ect
  
  lr_pos <- obj$multipliers$long_run_positive[variable]
  lr_neg <- obj$multipliers$long_run_negative[variable]
  sr_pos <- obj$multipliers$short_run_positive[variable]
  sr_neg <- obj$multipliers$short_run_negative[variable]
  
  # Compute dynamic multipliers
  # Using simplified dynamic path: m_h = SR + (LR - SR) * (1 - (1+phi)^h)
  
  h <- 0:horizon
  
  dm_pos <- sr_pos + (lr_pos - sr_pos) * (1 - (1 + phi)^h)
  dm_neg <- sr_neg + (lr_neg - sr_neg) * (1 - (1 + phi)^h)
  
  plot_data <- data.frame(
    horizon = rep(h, 2),
    multiplier = c(dm_pos, dm_neg),
    type = rep(c("Positive shock", "Negative shock"), each = length(h))
  )
  
  p <- ggplot(plot_data, aes(x = horizon, y = multiplier, color = type)) +
    geom_line(linewidth = 1.2) +
    geom_point(size = 2) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    geom_hline(yintercept = lr_pos, linetype = "dotted", color = "#2166AC", alpha = 0.7) +
    geom_hline(yintercept = lr_neg, linetype = "dotted", color = "#B2182B", alpha = 0.7) +
    scale_color_manual(values = c("Positive shock" = "#2166AC", "Negative shock" = "#B2182B")) +
    labs(
      title = paste("Dynamic Multipliers:", variable),
      subtitle = "Response to Unit Positive and Negative Shocks",
      x = "Horizon (Periods)",
      y = "Multiplier",
      color = "Shock Type"
    ) +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  return(p)
}


#' Plot Cumulative Multipliers
#'
#' @param obj FNARDL object
#' @param variable Variable to plot
#' @param horizon Number of periods
#'
#' @export
plot_cumulative_multipliers <- function(obj, variable, horizon = 20) {
  phi <- obj$multipliers$ect
  
  lr_pos <- obj$multipliers$long_run_positive[variable]
  lr_neg <- obj$multipliers$long_run_negative[variable]
  sr_pos <- obj$multipliers$short_run_positive[variable]
  sr_neg <- obj$multipliers$short_run_negative[variable]
  
  h <- 0:horizon
  
  # Dynamic multipliers
  dm_pos <- sr_pos + (lr_pos - sr_pos) * (1 - (1 + phi)^h)
  dm_neg <- sr_neg + (lr_neg - sr_neg) * (1 - (1 + phi)^h)
  
  # Cumulative (sum of dynamic)
  cm_pos <- cumsum(dm_pos)
  cm_neg <- cumsum(dm_neg)
  
  # Asymmetry (difference)
  asym <- cm_pos - cm_neg
  
  plot_data <- data.frame(
    horizon = rep(h, 3),
    value = c(cm_pos, cm_neg, asym),
    type = rep(c("Cumulative (+)", "Cumulative (-)", "Asymmetry"), each = length(h))
  )
  
  p <- ggplot(plot_data, aes(x = horizon, y = value, color = type)) +
    geom_line(linewidth = 1.2) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    scale_color_manual(values = c(
      "Cumulative (+)" = "#2166AC",
      "Cumulative (-)" = "#B2182B",
      "Asymmetry" = "#1B7837"
    )) +
    labs(
      title = paste("Cumulative Dynamic Multipliers:", variable),
      subtitle = "Cumulative Response and Asymmetry Over Time",
      x = "Horizon",
      y = "Cumulative Multiplier",
      color = ""
    ) +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  return(p)
}


#' Plot Positive vs Negative Comparison
#'
#' @param obj FNARDL object
#'
#' @keywords internal
plot_pos_neg_comparison <- function(obj) {
  lr_pos <- obj$multipliers$long_run_positive
  lr_neg <- obj$multipliers$long_run_negative
  
  plot_data <- data.frame(
    variable = names(lr_pos),
    positive = lr_pos,
    negative = lr_neg
  )
  
  # Add asymmetry test results
  plot_data$asymmetric <- sapply(plot_data$variable, function(v) {
    obj$asymmetry_tests[[v]]$decision == "Asymmetric"
  })
  
  p <- ggplot(plot_data, aes(x = positive, y = negative)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
    geom_point(aes(color = asymmetric, shape = asymmetric), size = 4) +
    geom_text(aes(label = variable), vjust = -1, size = 3) +
    scale_color_manual(values = c("TRUE" = "#B2182B", "FALSE" = "#2166AC"),
                       labels = c("TRUE" = "Asymmetric", "FALSE" = "Symmetric")) +
    scale_shape_manual(values = c("TRUE" = 17, "FALSE" = 16),
                       labels = c("TRUE" = "Asymmetric", "FALSE" = "Symmetric")) +
    labs(
      title = "Long-run Multiplier Comparison",
      subtitle = "45? line indicates symmetry",
      x = "Positive Shock Multiplier",
      y = "Negative Shock Multiplier",
      color = "Asymmetry",
      shape = "Asymmetry"
    ) +
    coord_fixed() +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  return(p)
}


#' Generate FNARDL Report
#'
#' @description
#' Generates a comprehensive report with all plots and tables.
#'
#' @param obj FNARDL object
#' @param file Output file path (HTML or PDF)
#' @param horizon Horizon for dynamic multipliers
#' @param verbose Logical. Print completion message (default: TRUE)
#'
#' @export
generate_fnardl_report <- function(obj, file = "fnardl_report.html", horizon = 20,
                                   verbose = TRUE) {
  # Create temporary Rmd file
  rmd_content <- sprintf('
---
title: "FNARDL Analysis Report"
date: "`r Sys.Date()`"
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = FALSE, warning = FALSE, message = FALSE)
```

## Model Summary

```{r}
summary(obj)
```
## Asymmetric Multipliers

```{r fig.width=10, fig.height=6}
plot(obj, type = "asymmetry")
```

## Dynamic Multipliers

```{r fig.width=10, fig.height=6}
for (var in obj$decompose) {
  print(plot_dynamic_multipliers(obj, var, %d))
}
```

## Cumulative Multipliers

```{r fig.width=10, fig.height=6}
for (var in obj$decompose) {
  print(plot_cumulative_multipliers(obj, var, %d))
}
```
', horizon, horizon)
  
  # Write and render
  tmp_rmd <- tempfile(fileext = ".Rmd")
  writeLines(rmd_content, tmp_rmd)
  
  rmarkdown::render(tmp_rmd, output_file = file, envir = new.env())
  
  if (verbose) message(sprintf("Report saved to: %s", file))
}

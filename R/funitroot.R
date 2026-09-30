# =============================================================================
# Fourier Unit Root Tests
# Based on Enders & Lee (2012) and Becker, Enders & Lee (2006)
# R implementation: Muhammad Alkhalaf (Rufyq Elngeh)
# =============================================================================

#' Fourier ADF Test
#'
#' @description
#' Tests for unit roots allowing for smooth structural breaks using
#' Fourier approximation. Implements Enders & Lee (2012) methodology.
#'
#' @param y Numeric vector of time series data
#' @param model Model specification: "c" (constant), "ct" (constant + trend)
#' @param max_freq Maximum Fourier frequency to test (default: 3)
#' @param max_lag Maximum lag for ADF (default: NULL, auto-select)
#' @param criterion Lag selection criterion ("AIC", "BIC", "t-sig")
#'
#' @return Object of class "fadf" with test results
#'
#' @examples
#' \donttest{
#' set.seed(123)
#' y <- cumsum(rnorm(200))  # Random walk
#' result <- fourier_adf_test(y, model = "c", max_freq = 3)
#' print(result)
#' }
#'
#' @references
#' Enders, W., & Lee, J. (2012). The flexible Fourier form and Dickey-Fuller
#' type unit root tests. Economics Letters, 117(1), 196-199.
#'
#' @export
fourier_adf_test <- function(y, model = c("c", "ct"), 
                              max_freq = 3, 
                              max_lag = NULL,
                              criterion = c("AIC", "BIC", "t-sig"),
                              verbose = TRUE) {
  
  model <- match.arg(model)
  criterion <- match.arg(criterion)
  if (max_freq < 1 || max_freq > 5)
    stop("'max_freq' must be between 1 and 5 (critical values are tabulated for k <= 5).")
  
  n <- length(y)
  
  # Default max lag (Schwert rule)
  if (is.null(max_lag)) {
    max_lag <- floor(12 * (n / 100)^0.25)
  }
  
  if (verbose) {
    message("=================================================================")
    message("   Fourier ADF Unit Root Test")
    message("   Enders & Lee (2012) Economics Letters")
    message("=================================================================\n")
  }
  
  # Test for each frequency and find optimal k
  results_by_k <- list()
  ssr_values <- numeric(max_freq)
  
  for (k in 1:max_freq) {
    results_by_k[[k]] <- estimate_fadf(y, k, model, max_lag, criterion)
    ssr_values[k] <- results_by_k[[k]]$ssr
  }
  
  # Optimal k minimizes SSR
  optimal_k <- which.min(ssr_values)
  best_result <- results_by_k[[optimal_k]]
  
  if (verbose) {
    message(sprintf("Model: %s", if (model == "c") "Constant" else "Constant + Trend"))
    message(sprintf("Sample size: %d", n))
    message(sprintf("Maximum lag tested: %d", max_lag))
    message(sprintf("Maximum frequency: %d", max_freq))
    message(sprintf("Lag selection: %s\n", criterion))
    message("-----------------------------------------------------------------")
    message(sprintf("Optimal frequency (k): %d", optimal_k))
    message(sprintf("Optimal lag (p): %d", best_result$optimal_lag))
    message(sprintf("ADF statistic: %.4f", best_result$adf_stat))
    message("-----------------------------------------------------------------")
  }
  
  # Critical values
  cv <- get_fadf_critical_values(n, model, optimal_k)
  if (verbose) {
    message("\nCritical values:")
    message(sprintf("   1%%  : %.4f %s", cv[1], 
                if (best_result$adf_stat < cv[1]) "*" else ""))
    message(sprintf("   5%%  : %.4f %s", cv[2], 
                if (best_result$adf_stat < cv[2]) "*" else ""))
    message(sprintf("   10%% : %.4f %s", cv[3], 
                if (best_result$adf_stat < cv[3]) "*" else ""))
  }
  
  # F-test for linearity
  f_test <- fadf_f_test(y, model, optimal_k, best_result$optimal_lag)
  if (verbose) {
    message("-----------------------------------------------------------------")
    message("F-test for linearity (H0: no Fourier terms needed):")
    message(sprintf("   F-statistic: %.4f", f_test$f_stat))
    message(sprintf("   Critical values 1%%/5%%/10%%: %.2f / %.2f / %.2f",
                    f_test$critical_values[1], f_test$critical_values[2],
                    f_test$critical_values[3]))
    message(sprintf("   -> %s", 
                if (f_test$reject) "Reject linearity: Fourier terms ARE significant" 
                else "Cannot reject linearity: Consider standard ADF"))
  }
  
  # Conclusion
  reject <- best_result$adf_stat < cv[2]  # 5% level
  if (verbose) {
    message("-----------------------------------------------------------------")
    message(sprintf("Conclusion: %s null hypothesis of unit root",
                if (reject) "Reject" else "Cannot reject"))
    message("            at 5% significance level")
    message("=================================================================")
  }
  
  result <- list(
    statistic = best_result$adf_stat,
    p_value = best_result$p_value,
    optimal_frequency = optimal_k,
    optimal_lag = best_result$optimal_lag,
    critical_values = cv,
    f_test = f_test,
    reject_null = reject,
    model = model,
    n = n,
    ssr_by_k = ssr_values,
    all_results = results_by_k
  )
  
  class(result) <- "fadf"
  return(result)
}


#' Estimate Fourier ADF for given k
#'
#' @keywords internal
estimate_fadf <- function(y, k, model, max_lag, criterion) {
  
  n <- length(y)
  t_index <- 1:n
  
  # Fourier terms
  sin_term <- sin(2 * pi * k * t_index / n)
  cos_term <- cos(2 * pi * k * t_index / n)
  
  # First difference
  dy <- diff(y)
  y_lag1 <- y[1:(n-1)]
  sin_adj <- sin_term[2:n]
  cos_adj <- cos_term[2:n]
  
  # Find optimal lag
  best_aic <- Inf
  best_lag <- 0
  
  for (p in 0:max_lag) {
    if (p >= (n - 10)) next
    
    tryCatch({
      result <- fit_fadf_model(dy, y_lag1, sin_adj, cos_adj, p, model)
      
      if (criterion == "AIC" && result$aic < best_aic) {
        best_aic <- result$aic
        best_lag <- p
      } else if (criterion == "BIC" && result$bic < best_aic) {
        best_aic <- result$bic
        best_lag <- p
      }
    }, error = function(e) NULL)
  }
  
  # Final estimation with optimal lag
  final <- fit_fadf_model(dy, y_lag1, sin_adj, cos_adj, best_lag, model)
  
  return(list(
    adf_stat = final$adf_stat,
    p_value = final$p_value,
    optimal_lag = best_lag,
    ssr = final$ssr,
    coefficients = final$coefficients,
    std_errors = final$std_errors
  ))
}


#' Fit Fourier ADF Model
#'
#' @keywords internal
fit_fadf_model <- function(dy, y_lag1, sin_term, cos_term, p, model) {
  
  n <- length(dy)
  max_lag <- p
  n_eff <- n - max_lag
  
  # Adjust series for lags
  dy_adj <- dy[(max_lag + 1):n]
  y_lag1_adj <- y_lag1[(max_lag + 1):n]
  sin_adj <- sin_term[(max_lag + 1):n]
  cos_adj <- cos_term[(max_lag + 1):n]
  
  # Build design matrix
  if (model == "c") {
    X <- cbind(1, y_lag1_adj, sin_adj, cos_adj)
    col_names <- c("const", "y_lag1", "sin", "cos")
  } else {
    trend <- 1:n_eff
    X <- cbind(1, trend, y_lag1_adj, sin_adj, cos_adj)
    col_names <- c("const", "trend", "y_lag1", "sin", "cos")
  }
  
  # Add lagged differences
  if (p > 0) {
    for (j in 1:p) {
      lag_idx <- (max_lag + 1 - j):(n - j)
      dy_lag <- dy[lag_idx]
      X <- cbind(X, dy_lag)
      col_names <- c(col_names, paste0("dy_lag", j))
    }
  }
  
  colnames(X) <- col_names
  
  # OLS estimation
  model_fit <- lm(dy_adj ~ X - 1)
  summary_fit <- summary(model_fit)
  
  # Get delta (coefficient on y_lag1)
  delta_idx <- which(col_names == "y_lag1")
  delta <- coef(model_fit)[delta_idx]
  se_delta <- summary_fit$coefficients[delta_idx, 2]
  t_stat <- delta / se_delta
  
  # The distribution is non-standard and only critical values are
  # tabulated (Enders and Lee, 2012), so no p-value is reported
  p_value <- NA_real_
  
  # Information criteria
  k <- length(coef(model_fit))
  ssr <- sum(residuals(model_fit)^2)
  aic <- n_eff * log(ssr / n_eff) + 2 * k
  bic <- n_eff * log(ssr / n_eff) + k * log(n_eff)
  
  return(list(
    adf_stat = t_stat,
    p_value = p_value,
    coefficients = coef(model_fit),
    std_errors = summary_fit$coefficients[, 2],
    ssr = ssr,
    aic = aic,
    bic = bic,
    residuals = residuals(model_fit)
  ))
}


#' F-test for Linearity in Fourier ADF
#'
#' @description
#' Tests H0: gamma1 = gamma2 = 0 (no Fourier terms needed)
#'
#' @param y Time series
#' @param model Model specification
#' @param k Fourier frequency
#' @param p Number of lags
#'
#' @return List with the F statistic, its critical values (Enders and Lee,
#'   2012, Table 1b) and p_value = NA (the distribution is non-standard)
#'
#' @export
fadf_f_test <- function(y, model, k, p) {
  
  n <- length(y)
  t_index <- 1:n
  
  # Fourier terms
  sin_term <- sin(2 * pi * k * t_index / n)
  cos_term <- cos(2 * pi * k * t_index / n)
  
  dy <- diff(y)
  y_lag1 <- y[1:(n-1)]
  
  n_adj <- length(dy) - p
  dy_adj <- dy[(p + 1):length(dy)]
  y_lag1_adj <- y_lag1[(p + 1):length(y_lag1)]
  sin_adj <- sin_term[(p + 2):n]
  cos_adj <- cos_term[(p + 2):n]
  
  # Unrestricted model (with Fourier)
  if (model == "c") {
    X_unrestricted <- cbind(1, y_lag1_adj, sin_adj, cos_adj)
  } else {
    trend <- 1:n_adj
    X_unrestricted <- cbind(1, trend, y_lag1_adj, sin_adj, cos_adj)
  }
  
  # Add lagged diffs
  if (p > 0) {
    for (j in 1:p) {
      dy_lag <- dy[(p + 1 - j):(length(dy) - j)]
      X_unrestricted <- cbind(X_unrestricted, dy_lag)
    }
  }
  
  model_unrestricted <- lm(dy_adj ~ X_unrestricted - 1)
  ssr_unrestricted <- sum(residuals(model_unrestricted)^2)
  
  # Restricted model (without Fourier)
  if (model == "c") {
    X_restricted <- cbind(1, y_lag1_adj)
  } else {
    trend <- 1:n_adj
    X_restricted <- cbind(1, trend, y_lag1_adj)
  }
  
  if (p > 0) {
    for (j in 1:p) {
      dy_lag <- dy[(p + 1 - j):(length(dy) - j)]
      X_restricted <- cbind(X_restricted, dy_lag)
    }
  }
  
  model_restricted <- lm(dy_adj ~ X_restricted - 1)
  ssr_restricted <- sum(residuals(model_restricted)^2)
  
  # F-statistic
  q <- 2  # Number of restrictions (sin and cos)
  df2 <- n_adj - ncol(X_unrestricted)
  f_stat <- ((ssr_restricted - ssr_unrestricted) / q) / (ssr_unrestricted / df2)
  # Non-standard distribution under the unit root null: critical values
  # from Enders and Lee (2012, Table 1b), no p-value
  p_value <- NA_real_
  cv <- get_fadf_f_critical_values(n, model)
  
  return(list(
    f_stat = f_stat,
    p_value = p_value,
    critical_values = cv,
    reject = f_stat > cv[2]  # Reject at 5%
  ))
}


#' Get Fourier ADF Critical Values
#'
#' @keywords internal
get_fadf_critical_values <- function(n, model, k) {
  # Enders and Lee (2012), Table 1a, as tabulated in TSPDLIB (Nazlioglu,
  # _getFourierADFCrit). Rows k = 1..5; columns 1%, 5%, 10%; four
  # sample-size blocks (T <= 150, <= 350, <= 500, > 500).
  tabs <- if (model == "c") list(
    c(-4.42, -3.81, -3.49, -3.97, -3.27, -2.91, -3.77, -3.07, -2.71,
      -3.64, -2.97, -2.64, -3.58, -2.93, -2.60),
    c(-4.37, -3.78, -3.47, -3.93, -3.26, -2.92, -3.74, -3.06, -2.72,
      -3.62, -2.98, -2.65, -3.55, -2.94, -2.62),
    c(-4.35, -3.76, -3.46, -3.91, -3.26, -2.91, -3.70, -3.06, -2.72,
      -3.62, -2.97, -2.66, -3.56, -2.94, -2.62),
    c(-4.31, -3.75, -3.45, -3.89, -3.25, -2.90, -3.69, -3.05, -2.71,
      -3.61, -2.96, -2.64, -3.53, -2.93, -2.61)
  ) else list(
    c(-4.95, -4.35, -4.05, -4.69, -4.05, -3.71, -4.45, -3.78, -3.44,
      -4.29, -3.65, -3.29, -4.20, -3.56, -3.22),
    c(-4.87, -4.31, -4.02, -4.62, -4.01, -3.69, -4.38, -3.77, -3.43,
      -4.27, -3.63, -3.31, -4.18, -3.56, -3.24),
    c(-4.81, -4.29, -4.01, -4.57, -3.99, -3.67, -4.38, -3.76, -3.43,
      -4.25, -3.64, -3.31, -4.18, -3.56, -3.25),
    c(-4.80, -4.27, -4.00, -4.58, -3.98, -3.67, -4.38, -3.75, -3.43,
      -4.24, -3.63, -3.30, -4.16, -3.55, -3.24))
  if (k > 5) return(rep(NA_real_, 3))
  blk <- if (n <= 150) 1 else if (n <= 350) 2 else if (n <= 500) 3 else 4
  matrix(tabs[[blk]], 5, 3, byrow = TRUE)[k, ]
}


# F test for the Fourier terms: Enders and Lee (2012), Table 1b, as
# tabulated in TSPDLIB (fstatADFcv). Columns 1%, 5%, 10%.
get_fadf_f_critical_values <- function(n, model) {
  tab <- if (model == "c") rbind(c(10.35, 7.58, 6.35), c(10.02, 7.41, 6.25),
                                 c(9.78, 7.29, 6.16), c(9.72, 7.25, 6.11))
         else rbind(c(12.21, 9.14, 7.78), c(11.70, 8.88, 7.62),
                    c(11.52, 8.76, 7.53), c(11.35, 8.71, 7.50))
  blk <- if (n <= 150) 1 else if (n <= 350) 2 else if (n <= 500) 3 else 4
  tab[blk, ]
}


#' @export
print.fadf <- function(x, ...) {
  cat("\nFourier ADF Test (Enders & Lee, 2012)\n")
  cat("=====================================\n")
  cat(sprintf("ADF Statistic: %.4f\n", x$statistic))
  cat(sprintf("Critical values 1%%/5%%/10%%: %.4f / %.4f / %.4f\n",
              x$critical_values[1], x$critical_values[2], x$critical_values[3]))
  cat(sprintf("Optimal frequency: k = %d\n", x$optimal_frequency))
  cat(sprintf("Optimal lag: p = %d\n", x$optimal_lag))
  cat(sprintf("Conclusion: %s\n", 
              if (x$reject_null) "Stationary" else "Unit root"))
  invisible(x)
}


#' Fourier KPSS Test
#'
#' @description
#' Tests for stationarity allowing for smooth structural breaks.
#' Implements Becker, Enders & Lee (2006) methodology.
#'
#' @param y Numeric vector of time series data
#' @param model Model specification: "c" (constant), "ct" (constant + trend)
#' @param max_freq Maximum Fourier frequency (default: 3)
#'
#' @return Object of class "fkpss" with test results
#'
#' @references
#' Becker, R., Enders, W., & Lee, J. (2006). A stationarity test in the presence
#' of an unknown number of smooth breaks. Journal of Time Series Analysis, 27(3), 381-409.
#'
#' @param verbose Logical. Print progress messages (default: TRUE)
#' @export
fourier_kpss_test <- function(y, model = c("c", "ct"), max_freq = 3, verbose = TRUE) {
  
  model <- match.arg(model)
  if (max_freq < 1 || max_freq > 5)
    stop("'max_freq' must be between 1 and 5 (critical values are tabulated for k <= 5).")
  n <- length(y)
  t_index <- 1:n
  
  if (verbose) {
    message("=================================================================")
    message("   Fourier KPSS Stationarity Test")
    message("   Becker, Enders & Lee (2006)")
    message("=================================================================\n")
  }
  
  # Test for each frequency
  results_by_k <- list()
  ssr_values <- numeric(max_freq)
  
  for (k in 1:max_freq) {
    results_by_k[[k]] <- estimate_fkpss(y, k, model)
    ssr_values[k] <- results_by_k[[k]]$ssr
  }
  
  # Optimal k minimizes SSR
  optimal_k <- which.min(ssr_values)
  best_result <- results_by_k[[optimal_k]]
  
  if (verbose) {
    message(sprintf("Model: %s", if (model == "c") "Level" else "Trend"))
    message(sprintf("Sample size: %d", n))
    message(sprintf("Optimal frequency: k = %d\n", optimal_k))
    message("-----------------------------------------------------------------")
    message(sprintf("KPSS statistic: %.6f", best_result$kpss_stat))
    message("-----------------------------------------------------------------")
  }
  
  # Critical values
  cv <- get_fkpss_critical_values(model, optimal_k, n)
  if (verbose) {
    message("\nCritical values:")
    message(sprintf("   1%%  : %.4f %s", cv[1], 
                if (best_result$kpss_stat > cv[1]) "*" else ""))
    message(sprintf("   5%%  : %.4f %s", cv[2], 
                if (best_result$kpss_stat > cv[2]) "*" else ""))
    message(sprintf("   10%% : %.4f %s", cv[3], 
                if (best_result$kpss_stat > cv[3]) "*" else ""))
  }
  
  # Conclusion
  reject <- best_result$kpss_stat > cv[2]
  if (verbose) {
    message("-----------------------------------------------------------------")
    message(sprintf("Conclusion: %s null hypothesis of stationarity",
                if (reject) "Reject" else "Cannot reject"))
    message(sprintf("            -> Series is %s",
                if (reject) "NON-STATIONARY" else "STATIONARY"))
    message("=================================================================")
  }
  
  result <- list(
    statistic = best_result$kpss_stat,
    p_value = best_result$p_value,
    optimal_frequency = optimal_k,
    critical_values = cv,
    reject_null = reject,
    model = model,
    n = n
  )
  
  class(result) <- "fkpss"
  return(result)
}


#' Estimate Fourier KPSS
#'
#' @keywords internal
estimate_fkpss <- function(y, k, model) {
  
  n <- length(y)
  t_index <- 1:n
  
  # Fourier terms
  sin_term <- sin(2 * pi * k * t_index / n)
  cos_term <- cos(2 * pi * k * t_index / n)
  
  # Regression
  if (model == "c") {
    X <- cbind(1, sin_term, cos_term)
  } else {
    X <- cbind(1, t_index, sin_term, cos_term)
  }
  
  model_fit <- lm(y ~ X - 1)
  resid <- residuals(model_fit)
  ssr <- sum(resid^2)
  
  # Partial sums
  S <- cumsum(resid)
  
  # Long-run variance (Newey-West)
  sigma2_lr <- newey_west_variance(resid, floor(4 * (n/100)^(2/9)))
  
  # KPSS statistic
  kpss_stat <- sum(S^2) / (n^2 * sigma2_lr)
  
  # Only critical values are tabulated (Becker, Enders and Lee, 2006)
  p_value <- NA_real_

  return(list(
    kpss_stat = kpss_stat,
    p_value = p_value,
    ssr = ssr
  ))
}


#' Newey-West Variance Estimator
#'
#' @keywords internal
newey_west_variance <- function(resid, bandwidth) {
  n <- length(resid)
  gamma0 <- sum(resid^2) / n
  
  gamma_sum <- 0
  for (j in 1:bandwidth) {
    weight <- 1 - j / (bandwidth + 1)
    gamma_j <- sum(resid[1:(n-j)] * resid[(j+1):n]) / n
    gamma_sum <- gamma_sum + 2 * weight * gamma_j
  }
  
  return(gamma0 + gamma_sum)
}


#' Get Fourier KPSS Critical Values
#'
#' @keywords internal
get_fkpss_critical_values <- function(model, k, n = 100) {
  # Becker, Enders and Lee (2006), as tabulated in TSPDLIB
  # (_getFourierKPSSCrit). Rows k = 1..5; columns 1%, 5%, 10%; sample-size
  # blocks T <= 250, <= 500, > 500.
  tabs <- if (model == "c") list(
    c(0.2699, 0.1720, 0.1318, 0.6671, 0.4152, 0.3150, 0.7182, 0.4480, 0.3393,
      0.7222, 0.4592, 0.3476, 0.7386, 0.4626, 0.3518),
    c(0.2709, 0.1696, 0.1294, 0.6615, 0.4075, 0.3053, 0.7046, 0.4424, 0.3309,
      0.7152, 0.4491, 0.3369, 0.7344, 0.4571, 0.3415),
    c(0.2706, 0.1704, 0.1295, 0.6526, 0.4047, 0.3050, 0.7086, 0.4388, 0.3304,
      0.7163, 0.4470, 0.3355, 0.7297, 0.4525, 0.3422)
  ) else list(
    c(0.0716, 0.0546, 0.0471, 0.2022, 0.1321, 0.1034, 0.2103, 0.1423, 0.1141,
      0.2170, 0.1478, 0.1189, 0.2177, 0.1484, 0.1201),
    c(0.0720, 0.0539, 0.0463, 0.1968, 0.1278, 0.0995, 0.2091, 0.1404, 0.1123,
      0.2111, 0.1441, 0.1155, 0.2178, 0.1465, 0.1178),
    c(0.0718, 0.0538, 0.0461, 0.1959, 0.1275, 0.0994, 0.2081, 0.1398, 0.1117,
      0.2139, 0.1436, 0.1149, 0.2153, 0.1451, 0.1163))
  if (k > 5) return(rep(NA_real_, 3))
  blk <- if (n <= 250) 1 else if (n <= 500) 2 else 3
  matrix(tabs[[blk]], 5, 3, byrow = TRUE)[k, ]
}


print.fkpss <- function(x, ...) {
  cat("\nFourier KPSS Test (Becker, Enders & Lee, 2006)\n")
  cat("==============================================\n")
  cat(sprintf("KPSS Statistic: %.6f\n", x$statistic))
  cat(sprintf("Critical values 1%%/5%%/10%%: %.4f / %.4f / %.4f\n",
              x$critical_values[1], x$critical_values[2], x$critical_values[3]))
  cat(sprintf("Optimal frequency: k = %d\n", x$optimal_frequency))
  cat(sprintf("Conclusion: %s\n", 
              if (x$reject_null) "Non-stationary" else "Stationary"))
  invisible(x)
}


#' Complete Unit Root Analysis
#'
#' @description
#' Performs both Fourier ADF and Fourier KPSS tests for comprehensive
#' unit root analysis.
#'
#' @param y Time series
#' @param name Optional name for the series
#' @param max_freq Maximum Fourier frequency
#' @param verbose Logical. Print progress messages (default: TRUE)
#'
#' @return List with results from both tests and joint conclusion
#'
#' @export
fourier_unit_root_analysis <- function(y, name = "Series", max_freq = 3, verbose = TRUE) {
  
  if (verbose) {
    message("")
    message("###############################################################")
    message(sprintf("   Unit Root Analysis: %s", name))
    message("###############################################################\n")
  }
  
  # Fourier ADF
  if (verbose) message(">>> Fourier ADF Test <<<\n")
  adf_result <- fourier_adf_test(y, model = "c", max_freq = max_freq, verbose = verbose)
  
  if (verbose) message("\n>>> Fourier KPSS Test <<<\n")
  kpss_result <- fourier_kpss_test(y, model = "c", max_freq = max_freq, verbose = verbose)
  
  # Joint interpretation
  if (adf_result$reject_null && !kpss_result$reject_null) {
    conclusion <- "STATIONARY (I(0)) - Both tests agree"
    integration_order <- 0
  } else if (!adf_result$reject_null && kpss_result$reject_null) {
    conclusion <- "UNIT ROOT (I(1)) - Both tests agree"
    integration_order <- 1
  } else if (!adf_result$reject_null && !kpss_result$reject_null) {
    conclusion <- "INCONCLUSIVE (ADF: unit root, KPSS: stationary)"
    integration_order <- NA
  } else {
    conclusion <- "CONTRADICTORY (needs further analysis)"
    integration_order <- NA
  }
  
  if (verbose) {
    message("")
    message("===============================================================")
    message("   JOINT CONCLUSION")
    message("===============================================================")
    message(sprintf("\n%s\n", conclusion))
    
    if (adf_result$f_test$reject || kpss_result$optimal_frequency > 1) {
      message("Note: Fourier terms ARE significant -> Structural breaks present")
    } else {
      message("Note: Fourier terms not significant -> Standard tests may suffice")
    }
    message("===============================================================")
  }
  
  return(list(
    adf = adf_result,
    kpss = kpss_result,
    conclusion = conclusion,
    integration_order = integration_order,
    structural_breaks = adf_result$f_test$reject
  ))
}

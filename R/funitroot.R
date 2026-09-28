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
    message(sprintf("P-value: %.4f", best_result$p_value))
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
    message(sprintf("   P-value: %.4f", f_test$p_value))
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
  
  # P-value from Fourier ADF distribution (approximation)
  p_value <- fadf_pvalue(t_stat, length(dy_adj), model)
  
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
#' @return List with F-statistic and p-value
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
  p_value <- 1 - pf(f_stat, q, df2)
  
  # Critical values from Enders & Lee (2012)
  cv <- c(10.02, 7.41, 6.25)  # 1%, 5%, 10%
  
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
  # Critical values from Enders & Lee (2012) Table 1
  # Rows: k=1,2,3; Columns: 1%, 5%, 10%
  
  if (model == "c") {
    cv_table <- matrix(c(
      -4.37, -3.78, -3.47,  # k=1
      -4.20, -3.62, -3.32,  # k=2
      -4.11, -3.54, -3.25   # k=3
    ), nrow = 3, byrow = TRUE)
  } else {  # ct
    cv_table <- matrix(c(
      -4.81, -4.25, -3.96,  # k=1
      -4.67, -4.11, -3.82,  # k=2
      -4.58, -4.02, -3.73   # k=3
    ), nrow = 3, byrow = TRUE)
  }
  
  k_adj <- min(k, 3)
  return(cv_table[k_adj, ])
}


#' Approximate P-value for Fourier ADF
#'
#' @keywords internal
fadf_pvalue <- function(t_stat, n, model) {
  # Approximation based on response surface
  # Uses interpolation between critical values
  
  cv <- get_fadf_critical_values(n, model, 1)
  
  if (t_stat < cv[1]) {
    return(0.001)
  } else if (t_stat < cv[2]) {
    return(0.01 + (t_stat - cv[1]) / (cv[2] - cv[1]) * 0.04)
  } else if (t_stat < cv[3]) {
    return(0.05 + (t_stat - cv[2]) / (cv[3] - cv[2]) * 0.05)
  } else {
    return(0.10 + (t_stat - cv[3]) * 0.1)
  }
}


#' @export
print.fadf <- function(x, ...) {
  cat("\nFourier ADF Test (Enders & Lee, 2012)\n")
  cat("=====================================\n")
  cat(sprintf("ADF Statistic: %.4f\n", x$statistic))
  cat(sprintf("P-value: %.4f\n", x$p_value))
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
    message(sprintf("P-value: %.4f", best_result$p_value))
    message("-----------------------------------------------------------------")
  }
  
  # Critical values
  cv <- get_fkpss_critical_values(model, optimal_k)
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
  
  # P-value approximation
  cv <- get_fkpss_critical_values(model, k)
  if (kpss_stat > cv[1]) {
    p_value <- 0.001
  } else if (kpss_stat > cv[2]) {
    p_value <- 0.05 - (cv[1] - kpss_stat) / (cv[1] - cv[2]) * 0.04
  } else if (kpss_stat > cv[3]) {
    p_value <- 0.10 - (cv[2] - kpss_stat) / (cv[2] - cv[3]) * 0.05
  } else {
    p_value <- 0.10 + (cv[3] - kpss_stat) * 0.5
  }
  
  return(list(
    kpss_stat = kpss_stat,
    p_value = min(1, max(0, p_value)),
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
get_fkpss_critical_values <- function(model, k) {
  # From Becker, Enders & Lee (2006) Table 1
  
  if (model == "c") {
    cv_table <- matrix(c(
      0.2699, 0.1720, 0.1318,  # k=1
      0.2022, 0.1321, 0.1034,  # k=2
      0.1669, 0.1117, 0.0886   # k=3
    ), nrow = 3, byrow = TRUE)
  } else {
    cv_table <- matrix(c(
      0.1295, 0.0912, 0.0756,  # k=1
      0.1084, 0.0787, 0.0661,  # k=2
      0.0958, 0.0709, 0.0602   # k=3
    ), nrow = 3, byrow = TRUE)
  }
  
  k_adj <- min(k, 3)
  return(cv_table[k_adj, ])
}


#' @export
print.fkpss <- function(x, ...) {
  cat("\nFourier KPSS Test (Becker, Enders & Lee, 2006)\n")
  cat("==============================================\n")
  cat(sprintf("KPSS Statistic: %.6f\n", x$statistic))
  cat(sprintf("P-value: %.4f\n", x$p_value))
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

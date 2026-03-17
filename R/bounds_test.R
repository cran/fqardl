#' =============================================================================
#' Bounds Test for Cointegration
#' Based on Pesaran, Shin & Smith (2001)
#' With extensions for Quantile ARDL
#' =============================================================================

#' Perform Bounds Test for Cointegration
#'
#' @description
#' Performs the PSS (2001) bounds test for cointegration in the ARDL framework.
#' Tests the joint significance of the lagged level variables.
#'
#' @param qardl_results List of QARDL estimation results
#' @param n Sample size
#' @param k Number of regressors
#' @param case Model case (1-5)
#'
#' @return List with bounds test results
#'
#' @details
#' The five cases are:
#' \itemize{
#'   \item Case 1: No intercept, no trend
#'   \item Case 2: Restricted intercept, no trend
#'   \item Case 3: Unrestricted intercept, no trend (most common)
#'   \item Case 4: Unrestricted intercept, restricted trend
#'   \item Case 5: Unrestricted intercept, unrestricted trend
#' }
#'
#' @export
perform_bounds_test <- function(qardl_results, n, k, case = 3) {
  
  # Use median quantile (tau = 0.5) for bounds test
  # Or take the first available
  median_idx <- which(sapply(qardl_results, function(x) x$tau) == 0.5)
  if (length(median_idx) == 0) median_idx <- 1
  
  res <- qardl_results[[median_idx]]
  coefs <- res$coefficients
  
  # Get coefficient and t-stat for y_lag1 (the ECT term)
  phi <- coefs["y_lag1"]
  t_phi <- res$t_statistics["y_lag1"]
  
  # Extract coefficients for level variables (x_lag1)
  level_coef_names <- grep("_lag1$", names(coefs), value = TRUE)
  level_coefs <- coefs[level_coef_names]
  
  # For F-test: test joint significance of all lagged levels
  # Wald test statistic
  # Simplified: sum of squared t-stats / number of restrictions
  t_stats_levels <- res$t_statistics[level_coef_names]
  F_stat <- mean(t_stats_levels^2)
  
  # Critical values from Pesaran et al. (2001) Table CI
  # Depend on case and k
  cv_table <- get_pss_critical_values(k, case)
  
  # Decision based on bounds
  if (F_stat > cv_table$F_upper[2]) {  # 5% upper bound
    decision_F <- "Cointegration exists (above upper bound at 5%)"
  } else if (F_stat < cv_table$F_lower[2]) {
    decision_F <- "No cointegration (below lower bound at 5%)"
  } else {
    decision_F <- "Inconclusive (between bounds at 5%)"
  }
  
  # t-test decision
  if (t_phi < cv_table$t_upper[2]) {  # t-stat is negative
    decision_t <- "Cointegration exists (t-stat significant at 5%)"
  } else {
    decision_t <- "Cointegration uncertain based on t-test"
  }
  
  # Combined decision
  if (grepl("exists", decision_F) && grepl("exists", decision_t)) {
    decision <- "Strong evidence of cointegration"
  } else if (grepl("exists", decision_F) || grepl("exists", decision_t)) {
    decision <- "Evidence of cointegration"
  } else if (grepl("No cointegration", decision_F)) {
    decision <- "No evidence of cointegration"
  } else {
    decision <- "Inconclusive - consider bootstrap test"
  }
  
  return(list(
    F_stat = F_stat,
    t_stat = t_phi,
    cv_1 = c(cv_table$F_lower[1], cv_table$F_upper[1]),
    cv_5 = c(cv_table$F_lower[2], cv_table$F_upper[2]),
    cv_10 = c(cv_table$F_lower[3], cv_table$F_upper[3]),
    t_cv_1 = cv_table$t_upper[1],
    t_cv_5 = cv_table$t_upper[2],
    t_cv_10 = cv_table$t_upper[3],
    decision_F = decision_F,
    decision_t = decision_t,
    decision = decision,
    case = case,
    k = k
  ))
}


#' Get PSS Critical Values
#'
#' @description
#' Returns critical values from Pesaran, Shin & Smith (2001) Table CI.
#'
#' @param k Number of regressors (excluding the dependent variable)
#' @param case Model case (1-5)
#'
#' @return Data frame with critical values
#'
#' @keywords internal
get_pss_critical_values <- function(k, case = 3) {
  
  # Critical values for Case III (unrestricted intercept, no trend)
  # Based on PSS (2001) Table CI(iii)
  # Format: k = number of regressors
  
  # F-test critical values [I(0), I(1)]
  # Structure: list by k, each contains 1%, 5%, 10% for lower and upper
  
  if (case == 3) {
    # Case III critical values
    cv_F <- switch(as.character(k),
      "1" = list(lower = c(6.84, 4.94, 4.04), upper = c(7.84, 5.73, 4.78)),
      "2" = list(lower = c(4.81, 3.55, 2.96), upper = c(6.10, 4.38, 3.61)),
      "3" = list(lower = c(4.13, 3.10, 2.63), upper = c(5.45, 3.99, 3.35)),
      "4" = list(lower = c(3.74, 2.86, 2.45), upper = c(5.06, 3.78, 3.16)),
      "5" = list(lower = c(3.47, 2.69, 2.32), upper = c(4.77, 3.61, 3.02)),
      "6" = list(lower = c(3.28, 2.57, 2.22), upper = c(4.54, 3.46, 2.91)),
      # Default for k > 6
      list(lower = c(3.15, 2.48, 2.15), upper = c(4.36, 3.34, 2.82))
    )
    
    # t-test critical values (upper bound only relevant)
    cv_t <- switch(as.character(k),
      "1" = c(-3.96, -3.41, -3.13),
      "2" = c(-3.78, -3.22, -2.93),
      "3" = c(-3.66, -3.09, -2.81),
      "4" = c(-3.56, -3.00, -2.72),
      "5" = c(-3.49, -2.93, -2.65),
      "6" = c(-3.43, -2.87, -2.60),
      c(-3.38, -2.82, -2.55)  # default
    )
  } else if (case == 4) {
    # Case IV (with restricted trend)
    cv_F <- switch(as.character(k),
      "1" = list(lower = c(8.74, 6.56, 5.59), upper = c(9.63, 7.30, 6.26)),
      "2" = list(lower = c(6.34, 4.87, 4.19), upper = c(7.52, 5.85, 4.87)),
      "3" = list(lower = c(5.33, 4.09, 3.54), upper = c(6.50, 5.07, 4.31)),
      "4" = list(lower = c(4.79, 3.70, 3.21), upper = c(5.96, 4.66, 4.01)),
      list(lower = c(4.38, 3.42, 2.98), upper = c(5.52, 4.36, 3.78))
    )
    cv_t <- c(-4.40, -3.80, -3.50)
  } else {
    # Case V (unrestricted trend)
    cv_F <- switch(as.character(k),
      "1" = list(lower = c(8.17, 6.03, 5.09), upper = c(9.36, 6.76, 5.72)),
      "2" = list(lower = c(5.98, 4.57, 3.89), upper = c(7.16, 5.41, 4.59)),
      list(lower = c(5.15, 3.98, 3.43), upper = c(6.36, 4.88, 4.15))
    )
    cv_t <- c(-4.23, -3.64, -3.32)
  }
  
  return(list(
    F_lower = cv_F$lower,
    F_upper = cv_F$upper,
    t_upper = cv_t,
    significance = c("1%", "5%", "10%")
  ))
}


#' Bootstrap Bounds Test
#'
#' @description
#' Performs bootstrap-based bounds test for cointegration,
#' following McNown et al. (2018) methodology.
#'
#' @param y Dependent variable
#' @param X Independent variables
#' @param fourier Fourier terms
#' @param p Lag for y
#' @param q Lag for X
#' @param tau Quantiles
#' @param case Model case
#' @param n_boot Number of bootstrap replications
#' @param verbose Logical. Print progress messages (default: FALSE)
#'
#' @return List with bootstrap p-values
#'
#' @export
bootstrap_bounds_test <- function(y, X, fourier, p, q, tau, case, n_boot = 1000,
                                  verbose = FALSE) {
  
  n <- length(y)
  k <- ncol(X)
  
  # Estimate original model
  orig_result <- estimate_qardl(y, X, fourier, p, q, tau[1], case)
  
  # Original test statistics
  orig_t <- orig_result$t_statistics["y_lag1"]
  
  # Get level coefficients for F-test
  level_names <- grep("_lag1$", names(orig_result$coefficients), value = TRUE)
  orig_F <- mean(orig_result$t_statistics[level_names]^2)
  
  # Bootstrap under null of no cointegration
  # Generate data under I(1) process
  boot_F <- numeric(n_boot)
  boot_t <- numeric(n_boot)
  
  for (b in 1:n_boot) {
    # Generate random walk (I(1)) series
    y_boot <- cumsum(rnorm(n, sd = sd(diff(y))))
    
    # Estimate model on bootstrap sample
    tryCatch({
      boot_result <- estimate_qardl(y_boot, X, fourier, p, q, tau[1], case)
      
      boot_t[b] <- boot_result$t_statistics["y_lag1"]
      boot_F[b] <- mean(boot_result$t_statistics[level_names]^2, na.rm = TRUE)
    }, error = function(e) {
      boot_t[b] <- NA
      boot_F[b] <- NA
    })
    
    if (verbose && b %% 100 == 0) message(sprintf("   Bootstrap replication %d/%d", b, n_boot))
  }
  
  # Remove NAs
  boot_F <- boot_F[!is.na(boot_F)]
  boot_t <- boot_t[!is.na(boot_t)]
  
  # Calculate p-values
  p_value_F <- mean(boot_F >= orig_F)
  p_value_t <- mean(boot_t <= orig_t)  # t-stat is negative under alternative
  
  return(list(
    n_boot = n_boot,
    orig_F = orig_F,
    orig_t = orig_t,
    boot_F = boot_F,
    boot_t = boot_t,
    p_value_F = p_value_F,
    p_value_t = p_value_t,
    decision = if (p_value_F < 0.05 || p_value_t < 0.05) 
      "Reject null: Cointegration exists" 
    else 
      "Fail to reject null: No cointegration"
  ))
}


#' Quantile Wald Test for Coefficient Constancy
#'
#' @description
#' Tests whether coefficients are constant across quantiles.
#'
#' @param qardl_results List of QARDL results
#' @param coef_name Name of coefficient to test
#'
#' @return List with test results
#'
#' @export
quantile_wald_test <- function(qardl_results, coef_name) {
  
  n_tau <- length(qardl_results)
  
  # Extract coefficients and SEs
  coefs <- sapply(qardl_results, function(x) x$coefficients[coef_name])
  ses <- sapply(qardl_results, function(x) x$std_errors[coef_name])
  
  # Test: are all coefficients equal?
  # H0: beta_1 = beta_2 = ... = beta_q
  
  # Simple Wald test using chi-square
  mean_coef <- mean(coefs)
  var_coefs <- var(coefs)
  
  # Wald statistic
  W <- sum((coefs - mean_coef)^2 / ses^2)
  df <- n_tau - 1
  p_value <- 1 - pchisq(W, df)
  
  return(list(
    coefficient = coef_name,
    values = coefs,
    wald_stat = W,
    df = df,
    p_value = p_value,
    decision = if (p_value < 0.05) 
      "Reject constancy: Coefficients vary across quantiles" 
    else 
      "Fail to reject: Coefficients appear constant"
  ))
}

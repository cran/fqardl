#' =============================================================================
#' Multi-Threshold Nonlinear ARDL (MTNARDL)
#' Extension of NARDL with multiple threshold decomposition
#' Ported from Python: Dr. Merwan Roudane
#' R implementation: Muhammad Alkhalaf (Rufyq Elngeh)
#' =============================================================================

#' Multi-Threshold NARDL Estimation
#'
#' @description
#' Estimates NARDL model with multiple thresholds for more nuanced
#' asymmetric analysis. Allows decomposition into multiple regimes
#' (e.g., small positive, large positive, small negative, large negative).
#'
#' @param formula A formula of the form y ~ x1 + x2 + ...
#' @param data Data frame with time series
#' @param decompose Variables to decompose with thresholds
#' @param thresholds Named list of threshold values for each variable
#' @param max_p Maximum lag for dependent variable
#' @param max_q Maximum lag for independent variables
#' @param criterion Information criterion ("AIC", "BIC", "HQ")
#' @param case Model case (1-5)
#'
#' @return Object of class "mtnardl"
#'
#' @examples
#' \donttest{
#' result <- mtnardl(
#'   formula = gdp ~ oil_price,
#'   data = macro_data,
#'   decompose = "oil_price",
#'   thresholds = list(oil_price = c(-5, 0, 5)),
#'   max_p = 4, max_q = 4
#' )
#' summary(result)
#' }
#'
#' @param verbose Logical. Print progress messages (default: TRUE)
#' @export
mtnardl <- function(formula, data,
                    decompose = NULL,
                    thresholds = NULL,
                    max_p = 4,
                    max_q = 4,
                    criterion = c("BIC", "AIC", "HQ"),
                    case = 3,
                    verbose = TRUE) {
  
  criterion <- match.arg(criterion)
  
  # Extract variables
  vars <- all.vars(formula)
  y_name <- vars[1]
  x_names <- vars[-1]
  
  if (is.null(decompose)) {
    decompose <- x_names
  }
  
  # Default thresholds (0 only = standard NARDL)
  if (is.null(thresholds)) {
    thresholds <- lapply(decompose, function(v) 0)
    names(thresholds) <- decompose
  }
  
  # Extract data
  y <- data[[y_name]]
  n <- length(y)
  
  if (verbose) {
    message("=================================================================")
    message("   Multi-Threshold NARDL (MTNARDL) Estimation")
    message("   Multiple Regime Asymmetric Analysis")
    message("=================================================================\n")
  }
  
  # Step 1: Decompose with multiple thresholds
  if (verbose) message("Step 1: Multi-threshold decomposition...")
  decomposed <- list()
  regime_names <- list()
  
  for (var in decompose) {
    thresh <- thresholds[[var]]
    result <- decompose_multi_threshold(data[[var]], thresh)
    decomposed[[var]] <- result$components
    regime_names[[var]] <- result$names
    if (verbose) message(sprintf("   %s: %d regimes (thresholds: %s)", 
                var, length(result$names), 
                paste(thresh, collapse = ", ")))
  }
  if (verbose) message("")
  
  # Step 2: Build regressor matrix
  if (verbose) message("Step 2: Building regressor matrix...")
  X_full <- build_mtnardl_regressors(data, x_names, decompose, decomposed, regime_names)
  if (verbose) message(sprintf("   Total regressors: %d\n", ncol(X_full)))
  
  # Step 3: Select optimal lags
  if (verbose) message("Step 3: Selecting optimal lags...")
  lag_result <- select_mtnardl_lags(y, X_full, max_p, max_q, criterion)
  optimal_p <- lag_result$optimal_p
  optimal_q <- lag_result$optimal_q
  if (verbose) message(sprintf("   Optimal: p = %d, q = %d\n", optimal_p, optimal_q))
  
  # Step 4: Estimate model
  if (verbose) message("Step 4: Estimating MTNARDL model...")
  model_result <- estimate_mtnardl(y, X_full, optimal_p, optimal_q, case)
  if (verbose) message("   Model estimated.\n")
  
  # Step 5: Compute regime-specific multipliers
  if (verbose) message("Step 5: Computing regime-specific multipliers...")
  multipliers <- compute_regime_multipliers(model_result, decompose, regime_names)
  if (verbose) message("   Multipliers computed.\n")
  
  # Step 6: Test for regime asymmetry
  if (verbose) message("Step 6: Testing for regime asymmetry...")
  regime_tests <- test_regime_asymmetry(model_result, decompose, regime_names)
  
  if (verbose) {
    for (var in decompose) {
      message(sprintf("   %s:", var))
      for (test_name in names(regime_tests[[var]])) {
        test <- regime_tests[[var]][[test_name]]
        message(sprintf("      %s: Wald = %.3f (p = %.4f)",
                    test_name, test$wald_stat, test$p_value))
      }
    }
    message("")
  }
  
  # Step 7: Bounds test
  if (verbose) message("Step 7: Bounds test for cointegration...")
  bounds <- perform_mtnardl_bounds(model_result, n, ncol(X_full), case)
  if (verbose) {
    message(sprintf("   F-statistic: %.4f", bounds$F_stat))
    message(sprintf("   Decision: %s\n", bounds$decision))
    message("=================================================================")
    message("   Estimation complete!")
    message("=================================================================")
  }
  
  result <- list(
    call = match.call(),
    formula = formula,
    y_name = y_name,
    x_names = x_names,
    decompose = decompose,
    thresholds = thresholds,
    regime_names = regime_names,
    n = n,
    optimal_p = optimal_p,
    optimal_q = optimal_q,
    case = case,
    model = model_result,
    multipliers = multipliers,
    regime_tests = regime_tests,
    bounds_test = bounds
  )
  
  class(result) <- "mtnardl"
  return(result)
}


#' Multi-Threshold Decomposition
#'
#' @description
#' Decomposes a variable into multiple regimes based on thresholds.
#'
#' @param x Numeric vector
#' @param thresholds Threshold values (must include 0)
#'
#' @return List with regime components and names
#'
#' @export
decompose_multi_threshold <- function(x, thresholds) {
  
  # Ensure 0 is included and sort
  thresholds <- sort(unique(c(thresholds, 0)))
  
  n <- length(x)
  dx <- c(0, diff(x))
  
  # Create regime components
  n_regimes <- length(thresholds) + 1
  components <- matrix(0, nrow = n, ncol = n_regimes)
  
  # Name regimes
  regime_names <- character(n_regimes)
  for (i in 1:n_regimes) {
    if (i == 1) {
      regime_names[i] <- sprintf("lt_%g", thresholds[1])
    } else if (i == n_regimes) {
      regime_names[i] <- sprintf("gt_%g", thresholds[length(thresholds)])
    } else {
      regime_names[i] <- sprintf("%g_to_%g", thresholds[i-1], thresholds[i])
    }
  }
  
  # Allocate changes to regimes
  for (t in 1:n) {
    change <- dx[t]
    
    for (i in 1:n_regimes) {
      if (i == 1) {
        # Below lowest threshold
        components[t, i] <- min(change, thresholds[1])
        change <- change - components[t, i]
      } else if (i == n_regimes) {
        # Above highest threshold
        components[t, i] <- max(change - thresholds[length(thresholds)], 0) +
                            min(change, 0)
      } else {
        # Between thresholds
        lower <- thresholds[i-1]
        upper <- thresholds[i]
        if (change > 0) {
          components[t, i] <- max(0, min(change - lower, upper - lower))
        } else {
          components[t, i] <- min(0, max(change - upper, lower - upper))
        }
      }
    }
  }
  
  # Cumulative sums
  components <- apply(components, 2, cumsum)
  colnames(components) <- regime_names
  
  return(list(
    components = components,
    names = regime_names,
    thresholds = thresholds
  ))
}


#' Build MTNARDL Regressor Matrix
#'
#' @keywords internal
build_mtnardl_regressors <- function(data, x_names, decompose, decomposed, regime_names) {
  
  reg_list <- list()
  
  for (var in x_names) {
    if (var %in% decompose) {
      # Use regime components
      for (i in 1:ncol(decomposed[[var]])) {
        col_name <- paste0(var, "_", regime_names[[var]][i])
        reg_list[[col_name]] <- decomposed[[var]][, i]
      }
    } else {
      reg_list[[var]] <- data[[var]]
    }
  }
  
  X <- do.call(cbind, reg_list)
  return(X)
}


#' Select Optimal Lags for MTNARDL
#'
#' @keywords internal
select_mtnardl_lags <- function(y, X, max_p, max_q, criterion) {
  
  best_ic <- Inf
  best_p <- 1
  best_q <- 1
  
  for (p in 1:max_p) {
    for (q in 1:max_q) {
      tryCatch({
        result <- estimate_mtnardl(y, X, p, q, 3)
        ic <- if (criterion == "AIC") result$aic else if (criterion == "BIC") result$bic else result$hq
        if (ic < best_ic) {
          best_ic <- ic
          best_p <- p
          best_q <- q
        }
      }, error = function(e) NULL)
    }
  }
  
  return(list(optimal_p = best_p, optimal_q = best_q))
}


#' Estimate MTNARDL Model
#'
#' @keywords internal
estimate_mtnardl <- function(y, X, p, q, case) {
  
  n <- length(y)
  k <- ncol(X)
  max_lag <- max(p, q)
  n_eff <- n - max_lag
  
  # Build design matrix
  design_list <- list()
  col_names <- c()
  
  # Intercept
  design_list[["const"]] <- rep(1, n_eff)
  col_names <- c("const")
  
  # Trend (if case > 3)
  if (case >= 4) {
    design_list[["trend"]] <- 1:n_eff
    col_names <- c(col_names, "trend")
  }
  
  # Lagged y (ECT)
  y_lag1 <- y[max_lag:(n-1)]
  design_list[["y_lag1"]] <- y_lag1
  col_names <- c(col_names, "y_lag1")
  
  # Lagged differences of y
  dy <- diff(y)
  if (p > 1) {
    for (j in 1:(p-1)) {
      idx <- (max_lag - j):(n - 1 - j)
      if (length(idx) == n_eff) {
        design_list[[paste0("dy_lag", j)]] <- dy[idx]
        col_names <- c(col_names, paste0("dy_lag", j))
      }
    }
  }
  
  # X variables (lagged levels)
  for (i in 1:k) {
    x_lag <- X[max_lag:(n-1), i]
    var_name <- paste0(colnames(X)[i], "_lag1")
    design_list[[var_name]] <- x_lag
    col_names <- c(col_names, var_name)
  }
  
  # Differences of X
  dX <- rbind(0, diff(X))
  for (i in 1:k) {
    # Current
    dx_curr <- dX[(max_lag+1):n, i]
    var_name <- paste0("d_", colnames(X)[i])
    design_list[[var_name]] <- dx_curr
    col_names <- c(col_names, var_name)
    
    # Lagged
    if (q > 1) {
      for (j in 1:(q-1)) {
        idx <- (max_lag + 1 - j):(n - j)
        if (length(idx) == n_eff) {
          var_name <- paste0("d_", colnames(X)[i], "_lag", j)
          design_list[[var_name]] <- dX[idx, i]
          col_names <- c(col_names, var_name)
        }
      }
    }
  }
  
  # Dependent variable
  # CORRECTED in 1.0.3: the dependent variable of a PSS conditional error
  # correction model is Delta y, not y. Up to 1.0.2 the right-hand side was
  # the ECM design while the left-hand side was in levels, so the reported
  # coefficient on y_lag1 was 1 + rho instead of rho: positive, near unity,
  # with a large positive t-ratio at every quantile, and the long-run
  # multipliers -theta/phi carried the wrong sign.
  y_adj <- y[(max_lag+1):n] - y[max_lag:(n - 1)]
  
  # Combine
  X_design <- do.call(cbind, design_list)
  colnames(X_design) <- col_names
  
  # OLS
  model <- lm(y_adj ~ X_design - 1)
  summ <- summary(model)
  
  coefs <- coef(model)
  names(coefs) <- col_names
  
  se <- summ$coefficients[, 2]
  t_stats <- summ$coefficients[, 3]
  p_values <- summ$coefficients[, 4]
  names(se) <- names(t_stats) <- names(p_values) <- col_names
  
  # IC
  k_params <- length(coefs)
  ssr <- sum(residuals(model)^2)
  aic <- n_eff * log(ssr/n_eff) + 2 * k_params
  bic <- n_eff * log(ssr/n_eff) + k_params * log(n_eff)
  hq <- n_eff * log(ssr/n_eff) + 2 * k_params * log(log(n_eff))
  
  return(list(
    coefficients = coefs,
    std_errors = se,
    t_statistics = t_stats,
    p_values = p_values,
    fitted = fitted(model),
    residuals = residuals(model),
    r_squared = summ$r.squared,
    adj_r_squared = summ$adj.r.squared,
    aic = aic,
    bic = bic,
    hq = hq,
    n = n_eff,
    model = model,
    vcov = vcov(model)
  ))
}


#' Compute Regime-Specific Multipliers
#'
#' @keywords internal
compute_regime_multipliers <- function(model_result, decompose, regime_names) {
  
  coefs <- model_result$coefficients
  phi <- coefs["y_lag1"]
  
  multipliers <- list()
  
  for (var in decompose) {
    long_run <- list()
    short_run <- list()
    
    for (regime in regime_names[[var]]) {
      # Long-run
      lr_name <- paste0(var, "_", regime, "_lag1")
      if (lr_name %in% names(coefs)) {
        long_run[[regime]] <- -coefs[lr_name] / phi
      }
      
      # Short-run
      sr_name <- paste0("d_", var, "_", regime)
      if (sr_name %in% names(coefs)) {
        short_run[[regime]] <- coefs[sr_name]
      }
    }
    
    multipliers[[var]] <- list(
      long_run = unlist(long_run),
      short_run = unlist(short_run)
    )
  }
  
  multipliers$ect <- phi
  return(multipliers)
}


#' Test Regime Asymmetry
#'
#' @keywords internal
test_regime_asymmetry <- function(model_result, decompose, regime_names) {
  
  coefs <- model_result$coefficients
  vcov_mat <- model_result$vcov
  
  results <- list()
  
  for (var in decompose) {
    var_tests <- list()
    regimes <- regime_names[[var]]
    
    # Test pairs of regimes
    for (i in 1:(length(regimes)-1)) {
      for (j in (i+1):length(regimes)) {
        r1 <- regimes[i]
        r2 <- regimes[j]
        
        name1 <- paste0(var, "_", r1, "_lag1")
        name2 <- paste0(var, "_", r2, "_lag1")
        
        if (name1 %in% names(coefs) && name2 %in% names(coefs)) {
          idx1 <- which(names(coefs) == name1)
          idx2 <- which(names(coefs) == name2)
          
          diff_coef <- coefs[idx1] - coefs[idx2]
          var_diff <- vcov_mat[idx1, idx1] + vcov_mat[idx2, idx2] - 
                      2 * vcov_mat[idx1, idx2]
          
          wald <- (diff_coef^2) / var_diff
          p_val <- 1 - pchisq(wald, 1)
          
          test_name <- paste0(r1, "_vs_", r2)
          var_tests[[test_name]] <- list(
            wald_stat = as.numeric(wald),
            p_value = as.numeric(p_val),
            decision = if (p_val < 0.05) "Different" else "Similar"
          )
        }
      }
    }
    
    results[[var]] <- var_tests
  }
  
  return(results)
}


#' MTNARDL Bounds Test
#'
#' @keywords internal
perform_mtnardl_bounds <- function(model_result, n, k, case) {
  
  coefs <- model_result$coefficients
  t_stats <- model_result$t_statistics
  
  level_names <- grep("_lag1$", names(coefs), value = TRUE)
  t_levels <- t_stats[level_names]
  
  # CORRECTED in 1.0.3: a genuine Wald statistic, not mean(t^2), and the
  # critical values now come from the PSS table for the actual number of
  # level terms rather than a hard-coded pair.
  F_stat <- wald_bounds_F(coefs, model_result$vcov, level_names)
  t_phi  <- unname(t_stats["y_lag1"])
  cv <- get_pss_critical_values(max(1L, length(level_names) - 1L), case = 3)
  decision <- bounds_verdict(F_stat, cv$F_lower[2], cv$F_upper[2])

  return(list(
    F_stat   = F_stat,
    t_stat   = t_phi,          # was missing entirely in 1.0.2
    cv_5     = c(cv$F_lower[2], cv$F_upper[2]),
    t_cv_5   = c(cv$t_lower[2], cv$t_upper[2]),
    decision = decision
  ))
}


#' @export
print.mtnardl <- function(x, ...) {
  cat("\nMulti-Threshold NARDL Model\n")
  cat("===========================\n")
  cat(sprintf("Dependent: %s\n", x$y_name))
  cat(sprintf("Decomposed: %s\n", paste(x$decompose, collapse = ", ")))
  for (var in x$decompose) {
    cat(sprintf("  %s thresholds: %s\n", var, 
                paste(x$thresholds[[var]], collapse = ", ")))
    cat(sprintf("  %s regimes: %s\n", var,
                paste(x$regime_names[[var]], collapse = ", ")))
  }
  cat(sprintf("Lags: ARDL(%d, %d)\n", x$optimal_p, x$optimal_q))
  cat(sprintf("Bounds test: F = %.4f (%s)\n", 
              x$bounds_test$F_stat, x$bounds_test$decision))
  invisible(x)
}


#' @export
summary.mtnardl <- function(object, ...) {
  cat("\n")
  cat("=================================================================\n")
  cat("     Multi-Threshold NARDL (MTNARDL) - Summary\n")
  cat("=================================================================\n\n")
  
  cat("MODEL SPECIFICATION\n")
  cat("-------------------\n")
  cat(sprintf("Dependent: %s\n", object$y_name))
  cat(sprintf("Sample: %d observations\n", object$n))
  cat(sprintf("Lags: ARDL(%d, %d)\n\n", object$optimal_p, object$optimal_q))
  
  cat("REGIME STRUCTURE\n")
  cat("----------------\n")
  for (var in object$decompose) {
    cat(sprintf("%s:\n", var))
    cat(sprintf("  Thresholds: %s\n", 
                paste(object$thresholds[[var]], collapse = ", ")))
    cat(sprintf("  Regimes: %s\n\n", 
                paste(object$regime_names[[var]], collapse = ", ")))
  }
  
  cat("REGIME-SPECIFIC LONG-RUN MULTIPLIERS\n")
  cat("------------------------------------\n")
  for (var in object$decompose) {
    cat(sprintf("%s:\n", var))
    lr <- object$multipliers[[var]]$long_run
    for (regime in names(lr)) {
      cat(sprintf("  %s: %.4f\n", regime, lr[regime]))
    }
    cat("\n")
  }
  
  cat("REGIME ASYMMETRY TESTS\n")
  cat("----------------------\n")
  for (var in object$decompose) {
    cat(sprintf("%s:\n", var))
    tests <- object$regime_tests[[var]]
    for (test_name in names(tests)) {
      test <- tests[[test_name]]
      cat(sprintf("  %s: Wald = %.3f, p = %.4f -> %s\n",
                  test_name, test$wald_stat, test$p_value, test$decision))
    }
    cat("\n")
  }
  
  cat("BOUNDS TEST\n")
  cat("-----------\n")
  cat(sprintf("F-statistic: %.4f\n", object$bounds_test$F_stat))
  cat(sprintf("Decision: %s\n", object$bounds_test$decision))
  
  cat("\nERROR CORRECTION\n")
  cat("----------------\n")
  cat(sprintf("ECT (phi): %.4f\n", object$multipliers$ect))
  cat(sprintf("Half-life: %.2f periods\n", 
              log(0.5) / log(1 + object$multipliers$ect)))
  
  cat("\n=================================================================\n")
  invisible(object)
}

# =============================================================================
# Multi-Threshold Nonlinear ARDL (MTNARDL)
# Extension of NARDL with multiple threshold decomposition
# R implementation: Muhammad Alkhalaf (Rufyq Elngeh)
# =============================================================================

#' Multi-Threshold NARDL Estimation (deprecated)
#'
#' @description
#' \strong{Deprecated.} \code{fqardl::mtnardl()} is deprecated and will be
#' removed in fqardl 2.0.0; use \code{ardlverse::mtnardl()}. The first call
#' in a session gives a warning. The conventions differ: the lag order p of
#' \code{ardlverse::mtnardl()} equals p - 1 here (here p counts the lagged
#' level of y together with the p - 1 lagged differences), and the
#' thresholds are a single numeric vector there (a named list here).
#'
#' Estimates the NARDL model with the changes of each decomposed regressor
#' split at threshold 0 into positive and negative partial sums.
#'
#' @details
#' \strong{Changes in 1.1.0.}
#' \itemize{
#'   \item The bounds F statistic uses only the lagged level of y and the
#'     lagged levels of the regressors (partial sums and undecomposed
#'     variables), with k equal to the number of level regressors. Up to
#'     1.0.6 every coefficient whose name ended in \code{_lag1} entered the
#'     test, including \code{dy_lag1} and the lagged differences, and the
#'     table was read at the wrong k. When the selected p or q exceeds 1 the
#'     F statistic, the bounds and the verdict therefore differ from 1.0.6.
#'   \item Thresholds other than 0 give an error: the multi-threshold
#'     decomposition of earlier versions did not allocate the changes
#'     correctly (the regime pieces did not sum to the change of the
#'     variable).
#'   \item Only \code{case = 3} is supported; other cases give an error
#'     instead of being silently replaced by case 3.
#'   \item The returned object has class \code{c("fqardl_mtnardl",
#'     "mtnardl")} and the print and summary methods are registered for
#'     \code{"fqardl_mtnardl"} only, so that they no longer mask the methods
#'     of ardlverse.
#' }
#' The bounds are those of Pesaran, Shin and Smith (2001), Table CI(iii), at
#' k equal to the number of level regressors. This choice of k for partial
#' sums is a package choice, pending verification against Shin, Yu and
#' Greenwood-Nimmo (2014) and Pal and Mitra (2016).
#'
#' @param formula A formula of the form y ~ x1 + x2 + ...
#' @param data Data frame with time series
#' @param decompose Variables to decompose
#' @param thresholds Named list of threshold values for each variable; only
#'   0 is supported since 1.1.0.
#' @param max_p Maximum lag for dependent variable
#' @param max_q Maximum lag for independent variables
#' @param criterion Information criterion ("AIC", "BIC", "HQ")
#' @param case Model case; only 3 (unrestricted intercept, no trend) is
#'   supported.
#' @param verbose Logical. Print progress messages (default: TRUE)
#'
#' @references
#' Shin, Y., Yu, B. and Greenwood-Nimmo, M. (2014). Modelling asymmetric
#' cointegration and dynamic multipliers in a nonlinear ARDL framework. In
#' \emph{Festschrift in Honor of Peter Schmidt}, 281-314. Springer, New
#' York. \doi{10.1007/978-1-4899-8008-3_9}
#'
#' @return Object of class \code{c("fqardl_mtnardl", "mtnardl")}
#'
#' @examples
#' \donttest{
#' result <- suppressWarnings(mtnardl(
#'   formula = gdp ~ oil_price,
#'   data = oil_gdp_data,
#'   decompose = "oil_price",
#'   max_p = 2, max_q = 2,
#'   verbose = FALSE
#' ))
#' summary(result)
#' }
#'
#' @export
mtnardl <- function(formula, data,
                    decompose = NULL,
                    thresholds = NULL,
                    max_p = 4,
                    max_q = 4,
                    criterion = c("BIC", "AIC", "HQ"),
                    case = 3,
                    verbose = TRUE) {
  
  mtnardl_deprecation_warning()
  criterion <- match.arg(criterion)
  if (!identical(as.numeric(case), 3))
    stop("fqardl::mtnardl() supports only case = 3 (unrestricted intercept, ",
         "no trend). Up to 1.0.6 other cases were silently replaced by case 3 ",
         "in the lag selection and the bounds test. ardlverse::mtnardl() ",
         "supports cases 1 to 5.", call. = FALSE)
  
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
  for (var in decompose) check_mtnardl_thresholds(thresholds[[var]])
  
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
  
  class(result) <- c("fqardl_mtnardl", "mtnardl")
  return(result)
}


# Once-per-session deprecation warning of fqardl::mtnardl()
.fqardl_state <- new.env(parent = emptyenv())

mtnardl_deprecation_warning <- function() {
  if (isTRUE(.fqardl_state$mtnardl_warned)) return(invisible(FALSE))
  .fqardl_state$mtnardl_warned <- TRUE
  warning(
    "fqardl::mtnardl() is deprecated and will be removed in fqardl 2.0.0; ",
    "use ardlverse::mtnardl(). Conventions differ: p in ardlverse equals ",
    "p - 1 here, and thresholds are a single numeric vector there (a named ",
    "list here). Results changed in fqardl 1.1.0: the bounds F statistic was ",
    "computed incorrectly when p or q exceeded 1 and now differs; thresholds ",
    "other than 0 and cases other than 3 give an error. ",
    "This warning is shown once per session.", call. = FALSE)
  invisible(TRUE)
}

# Thresholds other than 0 are refused (the old decomposition was wrong)
check_mtnardl_thresholds <- function(thresholds) {
  th <- sort(unique(c(thresholds, 0)))
  if (length(th) > 1)
    stop("Thresholds other than 0 are not supported since fqardl 1.1.0 (got: ",
         paste(thresholds, collapse = ", "), "). The multi-threshold ",
         "decomposition of fqardl <= 1.0.6 was mathematically wrong: the regime ",
         "pieces did not sum to the change of the variable. Use ",
         "ardlverse::mtnardl() for several thresholds.", call. = FALSE)
  invisible(th)
}


#' Multi-Threshold Decomposition
#'
#' @description
#' Decomposes a variable into partial sums of its positive and negative
#' changes (threshold 0), starting at zero in the first observation.
#'
#' Since 1.1.0 thresholds other than 0 give an error: the multi-threshold
#' allocation of earlier versions was wrong (for example, with thresholds
#' -2, 0 and 2 a change of -1 was split into pieces that did not sum to -1).
#' Use \code{ardlverse::mtnardl()} for several thresholds.
#'
#' @param x Numeric vector
#' @param thresholds Threshold values; only 0 is supported.
#'
#' @return List with regime components and names
#'
#' @export
decompose_multi_threshold <- function(x, thresholds) {
  
  # Ensure 0 is included and sort; only 0 is supported since 1.1.0
  thresholds <- check_mtnardl_thresholds(thresholds)
  
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
    vcov = vcov(model),
    level_names = paste0(colnames(X), "_lag1")
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
#' Wald bounds F statistic on the lagged level of y and the lagged levels of
#' the regressors, with k equal to the number of level regressors (package
#' choice, pending verification against Shin, Yu and Greenwood-Nimmo 2014
#' and Pal and Mitra 2016). Corrected in 1.1.0: earlier versions included
#' every coefficient whose name ended in \code{_lag1} (also \code{dy_lag1}
#' and the lagged differences of the regressors).
#'
#' @param model_result Output of \code{estimate_mtnardl}.
#' @param n Sample size.
#' @param k Number of level regressors (used when the design does not record
#'   them).
#' @param case Model case; only 3 is supported.
#' @return A list with \code{F_stat}, \code{t_stat}, \code{k}, the 5
#'   percent bounds and the F-based \code{decision}.
#' @keywords internal
perform_mtnardl_bounds <- function(model_result, n, k, case) {
  
  if (!identical(as.numeric(case), 3))
    stop("Only case = 3 is supported.", call. = FALSE)
  coefs <- model_result$coefficients
  t_stats <- model_result$t_statistics
  
  # CORRECTED in 1.1.0: only y_lag1 and the lagged levels of the regressors.
  x_lev <- model_result$level_names
  if (is.null(x_lev))
    x_lev <- setdiff(grep("_lag1$", names(coefs), value = TRUE),
                     c("y_lag1", grep("^(dy_lag|d_)", names(coefs), value = TRUE)))
  level_names <- c("y_lag1", x_lev)
  
  # CORRECTED in 1.0.3: a genuine Wald statistic, not mean(t^2).
  F_stat <- wald_bounds_F(coefs, model_result$vcov, level_names)
  t_phi  <- unname(t_stats["y_lag1"])
  k_lev <- length(x_lev)
  cv <- get_pss_critical_values(k_lev, case = 3)
  decision <- bounds_verdict(F_stat, cv$F_lower[2], cv$F_upper[2])

  return(list(
    F_stat   = F_stat,
    t_stat   = t_phi,          # was missing entirely in 1.0.2
    k        = k_lev,
    cv_5     = c(cv$F_lower[2], cv$F_upper[2]),
    t_cv_5   = c(cv$t_lower[2], cv$t_upper[2]),
    decision = decision
  ))
}


#' Plot method for fqardl mtnardl objects
#'
#' fqardl provides no plot for its (deprecated) \code{mtnardl} objects. This
#' method stops with an informative error instead of letting \code{plot()}
#' dispatch, through the inherited class \code{"mtnardl"}, to a method of
#' another package that expects a different object.
#'
#' @param x An object of class \code{"fqardl_mtnardl"}.
#' @param ... Not used.
#' @return Does not return; always an error.
#' @export
plot.fqardl_mtnardl <- function(x, ...) {
  stop("No plot method for fqardl::mtnardl() results (deprecated). ",
       "Use ardlverse::mtnardl(), whose objects have a plot method.",
       call. = FALSE)
}


#' @export
print.fqardl_mtnardl <- function(x, ...) {
  cat("\nMulti-Threshold NARDL Model (fqardl, deprecated)\n")
  cat("================================================\n")
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
summary.fqardl_mtnardl <- function(object, ...) {
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

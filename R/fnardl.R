#' =============================================================================
#' Fourier Nonlinear ARDL (FNARDL) Implementation
#' Combines NARDL (Shin et al., 2014) with Fourier approximation
#' Ported from Stata: Dr. Merwan Roudane
#' R implementation: Muhammad Alkhalaf (Rufyq Elngeh)
#' =============================================================================

#' Fourier Nonlinear ARDL Estimation
#'
#' @description
#' Estimates the Fourier Nonlinear ARDL (FNARDL) model that combines:
#' - Nonlinear ARDL for asymmetric effects (Shin et al., 2014)
#' - Fourier approximation for smooth structural breaks
#'
#' @param formula A formula of the form y ~ x1 + x2 + ...
#' @param data A data frame containing the time series variables
#' @param decompose Vector of variable names to decompose into positive/negative
#' @param max_p Maximum lag for dependent variable
#' @param max_q Maximum lag for independent variables
#' @param max_k Maximum Fourier frequency
#' @param criterion Information criterion ("AIC", "BIC", "HQ")
#' @param case Model case (1-5)
#' @param bootstrap Perform bootstrap test
#' @param n_boot Number of bootstrap replications
#' @param verbose Logical. Print progress messages (default: TRUE)
#'
#' @return An object of class "fnardl"
#'
#' @examples
#' \donttest{
#' result <- fnardl(
#'   formula = gdp ~ oil_price + exchange_rate,
#'   data = macro_data,
#'   decompose = c("oil_price"),  # Decompose oil price into + and -
#'   max_k = 3,
#'   verbose = FALSE
#' )
#' summary(result)
#' plot(result, type = "asymmetry")
#' }
#'
#' @export
fnardl <- function(formula, data,
                   decompose = NULL,
                   max_p = 4,
                   max_q = 4,
                   max_k = 3,
                   criterion = c("BIC", "AIC", "HQ"),
                   case = 3,
                   bootstrap = FALSE,
                   n_boot = 1000,
                   verbose = TRUE) {
  
  criterion <- match.arg(criterion)
  
  # Extract variables
  vars <- all.vars(formula)
  y_name <- vars[1]
  x_names <- vars[-1]
  
  # Validate decompose variables

if (is.null(decompose)) {
    decompose <- x_names  # Decompose all by default
  }
  
  if (!all(decompose %in% x_names)) {
    stop("All 'decompose' variables must be in the formula")
  }
  
  # Extract data
  y <- data[[y_name]]
  n <- length(y)
  
  if (verbose) {
    message("=================================================================")
    message("   Fourier Nonlinear ARDL (FNARDL) Estimation")
    message("   Asymmetric Effects with Smooth Structural Breaks")
    message("=================================================================\n")
  }
  
  # Step 1: Decompose variables into positive and negative changes
  if (verbose) message("Step 1: Decomposing variables into positive/negative changes...")
  decomposed_data <- decompose_variables(data, decompose)
  if (verbose) message(sprintf("   Decomposed: %s\n", paste(decompose, collapse = ", ")))
  
  # Step 2: Build full regressor matrix
  X_full <- build_nardl_regressors(data, x_names, decompose, decomposed_data)
  x_names_full <- colnames(X_full)
  if (verbose) message(sprintf("   Total regressors: %d\n", ncol(X_full)))
  
  # Step 3: Select optimal Fourier frequency
  if (verbose) message("Step 2: Selecting optimal Fourier frequency...")
  k_results <- select_fourier_frequency(y, X_full, max_k, criterion)
  optimal_k <- k_results$optimal_k
  if (verbose) message(sprintf("   Optimal k = %d\n", optimal_k))
  
  # Step 4: Generate Fourier terms
  fourier_terms <- generate_fourier_terms(n, optimal_k)
  
  # Step 5: Select optimal lags
  if (verbose) message("Step 3: Selecting optimal lag structure...")
  lag_results <- select_optimal_lags(y, X_full, fourier_terms, max_p, max_q, criterion)
  optimal_p <- lag_results$optimal_p
  optimal_q <- lag_results$optimal_q
  if (verbose) message(sprintf("   Optimal lags: p = %d, q = %d\n", optimal_p, optimal_q))
  
  # Step 6: Estimate NARDL model
  if (verbose) message("Step 4: Estimating Fourier NARDL model...")
  nardl_result <- estimate_nardl(
    y = y,
    X = X_full,
    fourier = fourier_terms,
    p = optimal_p,
    q = optimal_q,
    case = case
  )
  if (verbose) message("   Model estimated.\n")
  
  # Step 7: Calculate asymmetric multipliers
  if (verbose) message("Step 5: Computing asymmetric multipliers...")
  multipliers <- compute_asymmetric_multipliers(
    nardl_result, 
    decompose, 
    x_names
  )
  if (verbose) message("   Long-run and short-run multipliers computed.\n")
  
  # Step 8: Wald tests for asymmetry
  if (verbose) message("Step 6: Testing for asymmetry...")
  asymmetry_tests <- test_asymmetry(nardl_result, decompose)
  if (verbose) {
    for (var in decompose) {
      message(sprintf("   %s: Wald = %.3f (p = %.4f) - %s", 
                  var,
                  asymmetry_tests[[var]]$wald_stat,
                  asymmetry_tests[[var]]$p_value,
                  asymmetry_tests[[var]]$decision))
    }
    message("")
  }
  
  # Step 9: Bounds test
  if (verbose) message("Step 7: Bounds test for cointegration...")
  bounds_result <- perform_nardl_bounds_test(nardl_result, n, length(x_names_full), case)
  if (verbose) {
    message(sprintf("   F-statistic: %.4f", bounds_result$F_stat))
    message(sprintf("   Decision: %s\n", bounds_result$decision))
  }
  
  # Step 10: Bootstrap (if requested)
  bootstrap_results <- NULL
  if (bootstrap) {
    if (verbose) message("Step 8: Bootstrap cointegration test...")
    bootstrap_results <- bootstrap_nardl(
      y, X_full, fourier_terms, optimal_p, optimal_q, case, n_boot
    )
    if (verbose) message(sprintf("   Bootstrap p-value: %.4f\n", bootstrap_results$p_value))
  }
  
  if (verbose) {
    message("=================================================================")
    message("   Estimation complete!")
    message("=================================================================")
  }
  
  # Compile results
  result <- list(
    call = match.call(),
    formula = formula,
    y_name = y_name,
    x_names = x_names,
    x_names_full = x_names_full,
    decompose = decompose,
    n = n,
    optimal_k = optimal_k,
    optimal_p = optimal_p,
    optimal_q = optimal_q,
    case = case,
    criterion = criterion,
    nardl_result = nardl_result,
    multipliers = multipliers,
    asymmetry_tests = asymmetry_tests,
    bounds_test = bounds_result,
    bootstrap = bootstrap_results,
    fourier_terms = fourier_terms,
    decomposed_data = decomposed_data
  )
  
  class(result) <- "fnardl"
  return(result)
}


#' Decompose Variables into Positive and Negative Changes
#'
#' @description
#' Decomposes time series into cumulative positive and negative partial sums.
#'
#' @param data Data frame
#' @param variables Variables to decompose
#'
#' @return List with positive and negative components
#'
#' @export
decompose_variables <- function(data, variables) {
  
  result <- list()
  
  for (var in variables) {
    x <- data[[var]]
    dx <- c(0, diff(x))  # First difference
    
    # Positive changes (partial sum)
    dx_pos <- pmax(dx, 0)
    x_pos <- cumsum(dx_pos)
    
    # Negative changes (partial sum)
    dx_neg <- pmin(dx, 0)
    x_neg <- cumsum(dx_neg)
    
    result[[var]] <- list(
      positive = x_pos,
      negative = x_neg,
      d_positive = dx_pos,
      d_negative = dx_neg
    )
  }
  
  return(result)
}


#' Build NARDL Regressor Matrix
#'
#' @param data Original data
#' @param x_names All independent variable names
#' @param decompose Variables that are decomposed
#' @param decomposed_data Output from decompose_variables
#'
#' @return Matrix of regressors
#'
#' @keywords internal
build_nardl_regressors <- function(data, x_names, decompose, decomposed_data) {
  
  n <- nrow(data)
  reg_list <- list()
  
  for (var in x_names) {
    if (var %in% decompose) {
      # Use decomposed versions
      reg_list[[paste0(var, "_pos")]] <- decomposed_data[[var]]$positive
      reg_list[[paste0(var, "_neg")]] <- decomposed_data[[var]]$negative
    } else {
      # Use original variable
      reg_list[[var]] <- data[[var]]
    }
  }
  
  X <- do.call(cbind, reg_list)
  return(X)
}


#' Estimate NARDL Model
#'
#' @param y Dependent variable
#' @param X Regressor matrix (with decomposed variables)
#' @param fourier Fourier terms
#' @param p Lag for y
#' @param q Lag for X
#' @param case Model case
#'
#' @return List with estimation results
#'
#' @keywords internal
estimate_nardl <- function(y, X, fourier, p, q, case = 3) {
  
  n <- length(y)
  k <- ncol(X)
  max_lag <- max(p, q)
  n_eff <- n - max_lag
  
  # Build design matrix
  design_list <- list()
  col_names <- c()
  
  # Intercept
  design_list[["intercept"]] <- rep(1, n_eff)
  col_names <- c(col_names, "intercept")
  
  # Lagged y (ECT component)
  y_lag1 <- y[max_lag:(n - 1)]
  design_list[["y_lag1"]] <- y_lag1
  col_names <- c(col_names, "y_lag1")
  
  # Lagged differences of y
  dy <- diff(y)
  if (p > 1) {
    for (j in 1:(p - 1)) {
      lag_idx <- (max_lag - j):(n - 1 - j)
      if (length(lag_idx) == n_eff) {
        design_list[[paste0("dy_lag", j)]] <- dy[lag_idx]
        col_names <- c(col_names, paste0("dy_lag", j))
      }
    }
  }
  
  # Independent variables in levels (lagged)
  for (i in 1:k) {
    x_lag1 <- X[max_lag:(n - 1), i]
    var_name <- paste0(colnames(X)[i], "_lag1")
    design_list[[var_name]] <- x_lag1
    col_names <- c(col_names, var_name)
  }
  
  # First differences of X (current and lagged)
  dX <- rbind(rep(0, k), diff(X))
  for (i in 1:k) {
    # Current difference
    dx_current <- dX[(max_lag + 1):n, i]
    var_name <- paste0("d_", colnames(X)[i])
    design_list[[var_name]] <- dx_current
    col_names <- c(col_names, var_name)
    
    # Lagged differences
    if (q > 1) {
      for (j in 1:(q - 1)) {
        lag_idx <- (max_lag + 1 - j):(n - j)
        if (length(lag_idx) == n_eff) {
          var_name <- paste0("d_", colnames(X)[i], "_lag", j)
          design_list[[var_name]] <- dX[lag_idx, i]
          col_names <- c(col_names, var_name)
        }
      }
    }
  }
  
  # Fourier terms
  fourier_adj <- fourier[(max_lag + 1):n, , drop = FALSE]
  for (j in 1:ncol(fourier_adj)) {
    design_list[[colnames(fourier_adj)[j]]] <- fourier_adj[, j]
    col_names <- c(col_names, colnames(fourier_adj)[j])
  }
  
  # Dependent variable (adjusted)
  y_adj <- y[(max_lag + 1):n]
  
  # Combine design matrix
  X_design <- do.call(cbind, design_list)
  colnames(X_design) <- col_names
  
  # Estimate OLS
  model <- lm(y_adj ~ X_design - 1)
  
  # Extract results
  coefs <- coef(model)
  names(coefs) <- col_names
  
  summary_model <- summary(model)
  se <- summary_model$coefficients[, 2]
  t_stats <- summary_model$coefficients[, 3]
  p_values <- summary_model$coefficients[, 4]
  
  names(se) <- names(t_stats) <- names(p_values) <- col_names
  
  return(list(
    coefficients = coefs,
    std_errors = se,
    t_statistics = t_stats,
    p_values = p_values,
    fitted = fitted(model),
    residuals = residuals(model),
    r_squared = summary_model$r.squared,
    adj_r_squared = summary_model$adj.r.squared,
    sigma = summary_model$sigma,
    df = summary_model$df,
    n = n_eff,
    model = model,
    design = list(y = y_adj, X = X_design)
  ))
}


#' Compute Asymmetric Multipliers
#'
#' @param nardl_result NARDL estimation results
#' @param decompose Decomposed variable names
#' @param x_names Original variable names
#'
#' @return List with long-run and short-run asymmetric multipliers
#'
#' @export
compute_asymmetric_multipliers <- function(nardl_result, decompose, x_names) {
  
  coefs <- nardl_result$coefficients
  phi <- coefs["y_lag1"]  # Speed of adjustment
  
  long_run_pos <- list()
  long_run_neg <- list()
  short_run_pos <- list()
  short_run_neg <- list()
  
  for (var in decompose) {
    # Long-run multipliers: -beta / phi
    pos_lag_name <- paste0(var, "_pos_lag1")
    neg_lag_name <- paste0(var, "_neg_lag1")
    
    if (pos_lag_name %in% names(coefs)) {
      long_run_pos[[var]] <- -coefs[pos_lag_name] / phi
    }
    if (neg_lag_name %in% names(coefs)) {
      long_run_neg[[var]] <- -coefs[neg_lag_name] / phi
    }
    
    # Short-run multipliers
    pos_diff_name <- paste0("d_", var, "_pos")
    neg_diff_name <- paste0("d_", var, "_neg")
    
    if (pos_diff_name %in% names(coefs)) {
      short_run_pos[[var]] <- coefs[pos_diff_name]
    }
    if (neg_diff_name %in% names(coefs)) {
      short_run_neg[[var]] <- coefs[neg_diff_name]
    }
  }
  
  return(list(
    long_run_positive = unlist(long_run_pos),
    long_run_negative = unlist(long_run_neg),
    short_run_positive = unlist(short_run_pos),
    short_run_negative = unlist(short_run_neg),
    ect = phi
  ))
}


#' Test for Asymmetry (Wald Test)
#'
#' @param nardl_result NARDL estimation results
#' @param decompose Decomposed variables
#'
#' @return List of Wald test results for each variable
#'
#' @export
test_asymmetry <- function(nardl_result, decompose) {
  
  coefs <- nardl_result$coefficients
  vcov_mat <- vcov(nardl_result$model)
  
  results <- list()
  
  for (var in decompose) {
    # Long-run asymmetry test: H0: L+ = L-
    pos_name <- paste0(var, "_pos_lag1")
    neg_name <- paste0(var, "_neg_lag1")
    
    if (pos_name %in% names(coefs) && neg_name %in% names(coefs)) {
      # Find indices
      idx_pos <- which(names(coefs) == pos_name)
      idx_neg <- which(names(coefs) == neg_name)
      
      # Coefficient difference
      diff_coef <- coefs[idx_pos] - coefs[idx_neg]
      
      # Variance of difference
      var_diff <- vcov_mat[idx_pos, idx_pos] + vcov_mat[idx_neg, idx_neg] - 
                  2 * vcov_mat[idx_pos, idx_neg]
      
      # Wald statistic
      wald_stat <- (diff_coef^2) / var_diff
      p_value <- 1 - pchisq(wald_stat, df = 1)
      
      results[[var]] <- list(
        wald_stat = as.numeric(wald_stat),
        p_value = as.numeric(p_value),
        coef_positive = coefs[idx_pos],
        coef_negative = coefs[idx_neg],
        difference = as.numeric(diff_coef),
        decision = if (p_value < 0.05) "Asymmetric" else "Symmetric"
      )
    }
  }
  
  return(results)
}


#' NARDL Bounds Test
#'
#' @keywords internal
perform_nardl_bounds_test <- function(nardl_result, n, k, case) {
  # Use same logic as regular bounds test
  coefs <- nardl_result$coefficients
  t_stats <- nardl_result$t_statistics
  
  # Get lagged level coefficients
  level_names <- grep("_lag1$", names(coefs), value = TRUE)
  t_stats_levels <- t_stats[level_names]
  
  F_stat <- mean(t_stats_levels^2, na.rm = TRUE)
  t_phi <- t_stats["y_lag1"]
  
  # Get critical values
  cv <- get_pss_critical_values(length(level_names) - 1, case)
  
  # Decision
  if (F_stat > cv$F_upper[2]) {
    decision <- "Cointegration exists"
  } else if (F_stat < cv$F_lower[2]) {
    decision <- "No cointegration"
  } else {
    decision <- "Inconclusive"
  }
  
  return(list(
    F_stat = F_stat,
    t_stat = t_phi,
    cv_5 = c(cv$F_lower[2], cv$F_upper[2]),
    decision = decision
  ))
}


#' Bootstrap NARDL Test
#'
#' @keywords internal
bootstrap_nardl <- function(y, X, fourier, p, q, case, n_boot) {
  
  n <- length(y)
  
  # Original test statistic
  orig <- estimate_nardl(y, X, fourier, p, q, case)
  orig_t <- orig$t_statistics["y_lag1"]
  
  boot_t <- numeric(n_boot)
  
  for (b in 1:n_boot) {
    # Generate I(1) process
    y_boot <- cumsum(rnorm(n, sd = sd(diff(y))))
    
    tryCatch({
      boot_res <- estimate_nardl(y_boot, X, fourier, p, q, case)
      boot_t[b] <- boot_res$t_statistics["y_lag1"]
    }, error = function(e) {
      boot_t[b] <- NA
    })
    
    # Progress is suppressed for CRAN compliance
  }
  
  boot_t <- boot_t[!is.na(boot_t)]
  p_value <- mean(boot_t <= orig_t)
  
  return(list(
    n_boot = n_boot,
    p_value = p_value,
    orig_t = orig_t
  ))
}


#' @export
print.fnardl <- function(x, ...) {
  cat("\nFourier Nonlinear ARDL (FNARDL) Model\n")
  cat("=====================================\n\n")
  cat("Call:\n")
  print(x$call)
  cat("\n")
  cat(sprintf("Sample size: %d\n", x$n))
  cat(sprintf("Decomposed variables: %s\n", paste(x$decompose, collapse = ", ")))
  cat(sprintf("Fourier frequency: k = %d\n", x$optimal_k))
  cat(sprintf("Lags: p = %d, q = %d\n", x$optimal_p, x$optimal_q))
  cat("\nAsymmetry Tests:\n")
  for (var in x$decompose) {
    cat(sprintf("  %s: %s (p = %.4f)\n", 
                var, 
                x$asymmetry_tests[[var]]$decision,
                x$asymmetry_tests[[var]]$p_value))
  }
  invisible(x)
}


#' @export
summary.fnardl <- function(object, ...) {
  cat("\n")
  cat("=================================================================\n")
  cat("       Fourier Nonlinear ARDL (FNARDL) - Summary\n")
  cat("=================================================================\n\n")
  
  cat("MODEL SPECIFICATION\n")
  cat("-------------------\n")
  cat(sprintf("Dependent: %s\n", object$y_name))
  cat(sprintf("Decomposed: %s\n", paste(object$decompose, collapse = ", ")))
  cat(sprintf("Fourier k: %d | Lags: ARDL(%d, %d)\n\n", 
              object$optimal_k, object$optimal_p, object$optimal_q))
  
  cat("ASYMMETRIC LONG-RUN MULTIPLIERS\n")
  cat("-------------------------------\n")
  for (var in object$decompose) {
    cat(sprintf("%s(+): %.4f\n", var, object$multipliers$long_run_positive[var]))
    cat(sprintf("%s(-): %.4f\n", var, object$multipliers$long_run_negative[var]))
    cat(sprintf("  Asymmetry test: Wald = %.3f, p = %.4f (%s)\n\n",
                object$asymmetry_tests[[var]]$wald_stat,
                object$asymmetry_tests[[var]]$p_value,
                object$asymmetry_tests[[var]]$decision))
  }
  
  cat("BOUNDS TEST\n")
  cat("-----------\n")
  cat(sprintf("F-stat: %.4f | Decision: %s\n\n", 
              object$bounds_test$F_stat, object$bounds_test$decision))
  
  cat("ERROR CORRECTION TERM\n")
  cat("---------------------\n")
  cat(sprintf("ECT (phi): %.4f\n", object$multipliers$ect))
  cat(sprintf("Half-life: %.2f periods\n", log(0.5) / log(1 + object$multipliers$ect)))
  
  cat("\n=================================================================\n")
  invisible(object)
}

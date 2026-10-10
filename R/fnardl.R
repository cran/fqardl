# =============================================================================
# Fourier Nonlinear ARDL (FNARDL) Implementation
# Combines NARDL (Shin et al., 2014) with Fourier approximation
# R implementation: Muhammad Alkhalaf (Rufyq Elngeh)
# =============================================================================

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
#' @param case Model case; only 3 (unrestricted intercept, no trend) is
#'   supported.
#' @param bootstrap Perform the recursive bootstrap bounds test
#'   (see Details).
#' @param n_boot Number of bootstrap replications
#' @param verbose Logical. Print progress messages (default: TRUE)
#' @param fourier_k Optional Fourier frequency fixed by the user in advance;
#'   \code{NULL} (default) selects it from \code{1:max_k} and the bootstrap
#'   re-selects it in every replication.
#' @param boot_scheme \code{"bvz"} (default) or \code{"mcnown"}; see
#'   \code{\link{bootstrap_bounds_test}}.
#' @param seed Optional seed for the bootstrap; the caller's random number
#'   state is restored afterwards.
#'
#' @details
#' The positive and negative partial sums start at zero in the first
#' observation (\code{decompose_variables}). The regression always contains
#' Fourier terms, for which the Pesaran, Shin and Smith (2001) bounds are not
#' valid: since 1.1.0 \code{bounds_test} reports the F, t and Find
#' statistics with the bounds labelled as the bounds for the model without
#' Fourier terms, not valid with Fourier terms, and no analytical verdict.
#' With \code{bootstrap = TRUE} the recursive bootstrap generates x* and
#' y*, rebuilds the partial sums from x*, holds the Fourier terms fixed,
#' re-selects the Fourier frequency (unless \code{fourier_k} is given) and
#' the lag orders in every replication (a package choice; Bertelli, Vacca
#' and Zoia 2022 do not discuss re-selection), and combines the three tests with
#' the AND rule (cointegration only if Fov, t and Find all reject).
#'
#' @references
#' Shin, Y., Yu, B. and Greenwood-Nimmo, M. (2014). Modelling asymmetric
#' cointegration and dynamic multipliers in a nonlinear ARDL framework. In
#' \emph{Festschrift in Honor of Peter Schmidt}, 281-314. Springer, New
#' York. \doi{10.1007/978-1-4899-8008-3_9}
#'
#' Bertelli, S., Vacca, G. and Zoia, M. (2022). Bootstrap cointegration tests
#' in ARDL models. \emph{Economic Modelling}, 116, 105987.
#' \doi{10.1016/j.econmod.2022.105987}
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
                   verbose = TRUE,
                   fourier_k = NULL,
                   boot_scheme = c("bvz", "mcnown"),
                   seed = NULL) {
  
  criterion <- match.arg(criterion)
  boot_scheme <- match.arg(boot_scheme)
  if (!is.null(fourier_k) &&
      (length(fourier_k) != 1 || fourier_k < 1 || fourier_k != round(fourier_k)))
    stop("'fourier_k' must be a single positive integer", call. = FALSE)
  
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
  if (is.null(fourier_k)) {
    k_results <- select_fourier_frequency(y, X_full, max_k, criterion)
    optimal_k <- k_results$optimal_k
    if (verbose) message(sprintf("   Optimal k = %d\n", optimal_k))
  } else {
    optimal_k <- as.integer(fourier_k)
    if (verbose) message(sprintf("   k = %d (fixed by the user)\n", optimal_k))
  }
  
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
  if (!bootstrap && isFALSE(bounds_result$bounds_valid))
    bounds_result$decision <- paste0(bounds_result$decision, fourier_boot_advice)
  if (verbose) {
    message(sprintf("   F-statistic: %.4f", bounds_result$F_stat))
    message(sprintf("   Decision: %s\n", bounds_result$decision))
  }
  
  # Step 10: Bootstrap (if requested)
  bootstrap_results <- NULL
  if (bootstrap) {
    if (verbose) message("Step 8: Bootstrap cointegration test...")
    bootstrap_results <- bootstrap_nardl(
      y, as.matrix(data[, x_names, drop = FALSE]), x_names, decompose,
      fourier_terms, optimal_p, optimal_q, case, n_boot,
      reselect = list(k = fourier_k, max_k = max_k, max_p = max_p,
                      max_q = max_q, criterion = criterion),
      scheme = boot_scheme, seed = seed
    )
    if (verbose) {
      message(sprintf("   Bootstrap p-values: Fov = %.4f, t = %.4f, Find = %.4f",
                      bootstrap_results$p_value_F, bootstrap_results$p_value_t,
                      bootstrap_results$p_value_Find))
      message(sprintf("   Decision: %s\n", bootstrap_results$label))
    }
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

  d <- build_nardl_design(y, X, fourier, p, q)
  y_adj <- d$y
  X_design <- d$X
  col_names <- colnames(X_design)
  n_eff <- d$n

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
    design = list(y = y_adj, X = X_design),
    level_names = d$level_names
  ))
}


#' Build the FNARDL Design Matrix
#'
#' Conditional ECM design of \code{estimate_nardl}: intercept, lagged level
#' of y, lagged differences of y, lagged levels of the regressors, current
#' and lagged differences of the regressors, Fourier terms.
#'
#' @inheritParams estimate_nardl
#' @return List with \code{y} (the differenced dependent variable),
#'   \code{X} (design matrix), \code{n} and \code{level_names} (the
#'   lagged-level columns of the regressors).
#' @keywords internal
build_nardl_design <- function(y, X, fourier, p, q) {
  
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
  # CORRECTED in 1.0.3: the dependent variable of a PSS conditional error
  # correction model is Delta y, not y. Up to 1.0.2 the right-hand side was
  # the ECM design while the left-hand side was in levels, so the reported
  # coefficient on y_lag1 was 1 + rho instead of rho: positive, near unity,
  # with a large positive t-ratio at every quantile, and the long-run
  # multipliers -theta/phi carried the wrong sign.
  y_adj <- y[(max_lag + 1):n] - y[max_lag:(n - 1)]
  
  # Combine design matrix
  X_design <- do.call(cbind, design_list)
  colnames(X_design) <- col_names
  
  list(y = y_adj, X = X_design, n = n_eff,
       level_names = paste0(colnames(X), "_lag1"))
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
#' Wald bounds statistics of the FNARDL model. The overall F statistic tests
#' the lagged level of y and the lagged levels of the regressors (partial
#' sums and undecomposed variables); k is the number of level regressors.
#' Up to 1.0.6 the F statistic also included the lagged differences
#' \code{dy_lag1} and \code{d_*_lag1} whenever p or q exceeded 1.
#'
#' With Fourier terms in the regression (always the case for models estimated
#' by \code{fnardl}) the PSS (2001) bounds do not apply (Monte Carlo
#' evidence): no analytical verdict is given
#' and the PSS bounds are returned only for reference, labelled as the bounds
#' for the model without Fourier terms, not valid with Fourier terms. Use the
#' recursive bootstrap (\code{fnardl(..., bootstrap = TRUE)}).
#'
#' @param nardl_result Output of \code{estimate_nardl}.
#' @param n Sample size.
#' @param k Number of level regressors (used when the design does not record
#'   them).
#' @param case Model case; only 3 is supported.
#' @return A list with \code{F_stat}, \code{t_stat}, \code{Find_stat},
#'   \code{k}, the 5 percent bounds \code{cv_5} and \code{t_cv_5},
#'   \code{bounds_valid}, \code{bounds_note} and \code{decision}.
#' @keywords internal
perform_nardl_bounds_test <- function(nardl_result, n, k, case) {
  coefs <- nardl_result$coefficients
  t_stats <- nardl_result$t_statistics
  nms <- names(coefs)

  # CORRECTED in 1.1.0: only the lagged levels of the regressors, not every
  # coefficient whose name ends in "_lag1".
  x_lev <- nardl_result$level_names
  if (is.null(x_lev))
    x_lev <- setdiff(grep("_lag1$", nms, value = TRUE),
                     c("y_lag1", grep("^(dy_lag|d_)", nms, value = TRUE)))
  level_names <- c("y_lag1", x_lev)
  
  V <- tryCatch(stats::vcov(nardl_result$model), error = function(e) NULL)
  if (!is.null(V)) dimnames(V) <- list(nms, nms)
  F_stat <- wald_bounds_F(coefs, V, level_names)
  Find_stat <- wald_bounds_F(coefs, V, x_lev)
  t_phi <- unname(t_stats["y_lag1"])
  
  k_lev <- length(x_lev)
  cv <- get_pss_critical_values(k_lev, case)
  fourier <- has_fourier_terms(nms)
  
  decision <- if (fourier) fourier_no_verdict else
    bounds_verdict(F_stat, cv$F_lower[2], cv$F_upper[2])
  
  return(list(
    F_stat = F_stat,
    t_stat = t_phi,
    Find_stat = Find_stat,
    k = k_lev,
    cv_5 = c(cv$F_lower[2], cv$F_upper[2]),
    t_cv_5 = c(cv$t_lower[2], cv$t_upper[2]),
    bounds_valid = !fourier,
    bounds_note = if (fourier) fourier_bounds_note else
                    "PSS (2001) Table CI(iii), 5 percent",
    decision = decision
  ))
}


# Partial-sum decomposition of fnardl as a function of the regressor matrix,
# with its one-step form for the bootstrap engine.
.fn_decomposer <- function(x_names, decompose) {
  full <- function(x) {
    x <- as.matrix(x)
    colnames(x) <- x_names
    df <- as.data.frame(x)
    X <- build_nardl_regressors(df, x_names, decompose,
                                decompose_variables(df, decompose))
    as.matrix(X)
  }
  isdec <- x_names %in% decompose
  step <- function(dx) {
    unlist(lapply(seq_along(dx), function(j)
      if (isdec[j]) c(max(dx[j], 0), min(dx[j], 0)) else dx[j]))
  }
  list(full = full, step = step)
}

# Statistics of fnardl's ECM (build_nardl_design), estimated by OLS
.fn_stats <- function(y, x, spec, dec) {
  Xf <- dec$full(x)
  d <- build_nardl_design(y, Xf, .fq_spec_fourier(spec, length(y)), spec$p, spec$q)
  .fq_ols_stats(d$y, d$X, "y_lag1", d$level_names)
}


#' Bootstrap NARDL Test
#'
#' Recursive bootstrap bounds test for the FNARDL model (rewritten in
#' 1.1.0). The regressors x are generated recursively and the partial sums
#' are rebuilt from x* with \code{decompose_variables}; the Fourier terms are
#' fixed deterministic regressors; the statistics (Fov, t, Find) are OLS
#' statistics of fnardl's own design; the Fourier frequency (unless fixed)
#' and the lag orders are re-selected in every replication when
#' \code{reselect} is given (a package choice); the decision is the AND rule. See
#' \code{\link{bootstrap_bounds_test}} for the schemes. Up to 1.0.6 this
#' function regressed a new Gaussian random walk on fixed regressors and
#' reported the t test only, without a warning.
#'
#' @param y Dependent variable.
#' @param x Matrix of the original regressors (columns \code{x_names}).
#' @param x_names Names of the regressors.
#' @param decompose Names of the decomposed regressors.
#' @param fourier Fourier terms of the observed model.
#' @param p,q Lag orders as in \code{estimate_nardl}.
#' @param case Model case; only 3 is supported.
#' @param n_boot Number of bootstrap replications.
#' @param reselect \code{NULL} or a list (\code{max_k}, \code{max_p},
#'   \code{max_q}, \code{criterion}, optional \code{k}); see
#'   \code{\link{bootstrap_bounds_test}}.
#' @param scheme \code{"bvz"} or \code{"mcnown"}.
#' @param level Significance level of the decision.
#' @param seed Optional seed; the caller's random number state is restored.
#' @return A list as returned by \code{\link{bootstrap_bounds_test}};
#'   \code{p_value} is the p-value of the t test (as in 1.0.6).
#' @keywords internal
bootstrap_nardl <- function(y, x, x_names, decompose, fourier, p, q, case = 3,
                            n_boot = 1000, reselect = NULL,
                            scheme = c("bvz", "mcnown"), level = 0.05,
                            seed = NULL) {
  if (!identical(as.numeric(case), 3))
    stop("Only case = 3 (unrestricted intercept, no trend) is supported.",
         call. = FALSE)
  scheme <- match.arg(scheme)
  x <- as.matrix(x)
  colnames(x) <- x_names
  dec <- .fn_decomposer(x_names, decompose)
  kw <- ncol(dec$full(x))
  sc <- .fq_scheme(scheme)
  spec <- list(fourier = as.matrix(fourier), p = p, q = q)
  eng <- .ardl_boot_engine(
    as.numeric(y), x, p = p - 1, q = rep(q - 1, kw), case = 3,
    t0 = max(p, q) + 1, det = as.matrix(fourier), decompose = dec,
    nulls = sc$nulls, xmodel = sc$xmodel, B = n_boot, init = sc$init,
    recentre = sc$recentre,
    stat_fun = function(yy, xx, sp) .fn_stats(yy, xx, sp, dec),
    select = .fq_make_select(reselect, transform = dec$full),
    spec = spec, use_find = TRUE, level = level, seed = seed)
  out <- .fq_boot_result(eng, scheme, n_boot, reselect,
                         model = "Fourier NARDL (conditional mean), fnardl design")
  out$orig_t <- c(y_lag1 = out$orig_t)   # named as in 1.0.6
  out
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
  if (!is.null(x$bounds_test$decision))
    cat(sprintf("\nBounds test: %s\n", x$bounds_test$decision))
  if (!is.null(x$bootstrap$label))
    cat(sprintf("Bootstrap bounds test: %s\n", x$bootstrap$label))
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
    # FIXED in 1.1.0 (printout only): the elements are named
    # "<var>.<var>_pos_lag1" by unlist(), so indexing by <var> printed NA.
    lrp <- object$multipliers$long_run_positive
    lrn <- object$multipliers$long_run_negative
    cat(sprintf("%s(+): %.4f\n", var,
                unname(lrp[paste0(var, ".", var, "_pos_lag1")])))
    cat(sprintf("%s(-): %.4f\n", var,
                unname(lrn[paste0(var, ".", var, "_neg_lag1")])))
    cat(sprintf("  Asymmetry test: Wald = %.3f, p = %.4f (%s)\n\n",
                object$asymmetry_tests[[var]]$wald_stat,
                object$asymmetry_tests[[var]]$p_value,
                object$asymmetry_tests[[var]]$decision))
  }
  
  cat("BOUNDS TEST\n")
  cat("-----------\n")
  cat(sprintf("F-stat: %.4f | t-stat: %.4f | Decision: %s\n",
              object$bounds_test$F_stat, object$bounds_test$t_stat %||% NA_real_,
              object$bounds_test$decision))
  if (isFALSE(object$bounds_test$bounds_valid))
    cat(sprintf("5%% PSS bounds: I(0) = %.3f, I(1) = %.3f. %s\n",
                object$bounds_test$cv_5[1], object$bounds_test$cv_5[2],
                fourier_bounds_note))
  cat("\n")
  if (!is.null(object$bootstrap)) {
    cat("BOOTSTRAP BOUNDS TEST\n")
    cat("---------------------\n")
    if (is.null(object$bootstrap$engine))
      cat(sprintf("Bootstrap p-value (t): %.4f\n", object$bootstrap$p_value))
    else
      .fq_print_boot(object$bootstrap)
    cat("\n")
  }
  
  cat("ERROR CORRECTION TERM\n")
  cat("---------------------\n")
  cat(sprintf("ECT (phi): %.4f\n", object$multipliers$ect))
  cat(sprintf("Half-life: %.2f periods\n", log(0.5) / log(1 + object$multipliers$ect)))
  
  cat("\n=================================================================\n")
  invisible(object)
}

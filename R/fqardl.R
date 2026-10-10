`%||%` <- function(a, b) if (is.null(a)) b else a

# =============================================================================
# Fourier Quantile ARDL (FQARDL) - Main Functions
# Ported from Stata to R
# R implementation: Muhammad Alkhalaf (Rufyq Elngeh)
# =============================================================================

#' Fourier Quantile ARDL Estimation
#'
#' @description
#' Estimates the Fourier Quantile Autoregressive Distributed Lag (FQARDL) model.
#' This methodology extends QARDL by incorporating Fourier trigonometric terms
#' to capture smooth structural breaks without prior knowledge of break timing.
#'
#' @param formula A formula of the form y ~ x1 + x2 + ...
#' @param data A data frame containing the time series variables
#' @param tau Numeric vector of quantiles to estimate (default: c(0.25, 0.5, 0.75))
#' @param max_p Maximum lag for dependent variable (default: 4)
#' @param max_q Maximum lag for independent variables (default: 4)
#' @param max_k Maximum Fourier frequency to test (default: 3)
#' @param criterion Information criterion for lag selection ("AIC", "BIC", "HQ")
#' @param case Model case; only 3 (unrestricted intercept, no trend) is
#'   supported.
#' @param bootstrap Logical, perform the recursive bootstrap bounds test
#'   (\code{\link{bootstrap_bounds_test}}).
#' @param n_boot Number of bootstrap replications (default: 1000)
#' @param seed Random seed, set with \code{set.seed()} at the start (the
#'   quantile regression standard errors and the bootstrap use random
#'   numbers).
#' @param verbose Logical. Print progress messages (default: TRUE).
#' @param fourier_k Optional Fourier frequency fixed by the user in advance.
#'   \code{NULL} (default) selects it from \code{1:max_k}; the bootstrap then
#'   re-selects it in every replication.
#' @param boot_scheme Bootstrap scheme, \code{"bvz"} (Bertelli, Vacca and
#'   Zoia 2022, default) or \code{"mcnown"} (the Fov null for all statistics,
#'   McNown, Sam and Goh 2018, applied to the conditional ECM); see
#'   \code{\link{bootstrap_bounds_test}}.
#'
#' @details
#' \strong{Inference on cointegration.} The regression always contains
#' Fourier terms, for which the Pesaran, Shin and Smith (2001) bounds are not
#' valid. Since 1.1.0 \code{bounds_test} reports the Wald F and t
#' statistics at each quantile with the PSS bounds labelled as the bounds
#' for the model without Fourier terms, and gives no analytical verdict.
#' Use \code{bootstrap = TRUE}: the recursive bootstrap of
#' \code{\link{bootstrap_bounds_test}} with OLS conditional ECM statistics,
#' the Fourier frequency (unless \code{fourier_k} is given) and the lag
#' orders re-selected in every replication (a package choice), and the AND
#' decision rule. The
#' quantile-specific statistics have no published null distribution (Cho,
#' Kim and Shin 2015) and are descriptive.
#'
#' @return An object of class "fqardl" containing:
#' \item{coefficients}{Estimated coefficients for each quantile}
#' \item{long_run}{Long-run multipliers}
#' \item{short_run}{Short-run multipliers}
#' \item{optimal_k}{Optimal Fourier frequency}
#' \item{optimal_lags}{Optimal lag structure}
#' \item{bounds_test}{Wald bounds statistics by quantile (no verdict with
#'   Fourier terms; see Details)}
#' \item{bootstrap}{Result of \code{\link{bootstrap_bounds_test}} when
#'   \code{bootstrap = TRUE}}
#' \item{diagnostics}{Model diagnostics}
#'
#' @references
#' Pesaran, M.H., Shin, Y. and Smith, R.J. (2001). Bounds testing approaches
#' to the analysis of level relationships. \emph{Journal of Applied
#' Econometrics}, 16(3), 289-326. \doi{10.1002/jae.616}
#'
#' Cho, J.S., Kim, T. and Shin, Y. (2015). Quantile cointegration in the
#' autoregressive distributed-lag modeling framework. \emph{Journal of
#' Econometrics}, 188(1), 281-300. \doi{10.1016/j.jeconom.2015.05.003}
#'
#' Bertelli, S., Vacca, G. and Zoia, M. (2022). Bootstrap cointegration tests
#' in ARDL models. \emph{Economic Modelling}, 116, 105987.
#' \doi{10.1016/j.econmod.2022.105987}
#'
#' @examples
#' \donttest{
#' data(macro_data)
#' result <- fqardl(gdp ~ inflation + interest_rate, 
#'                  data = macro_data,
#'                  tau = c(0.25, 0.5, 0.75),
#'                  max_k = 3)
#' summary(result)
#' plot(result)
#' }
#'
#' @export
fqardl <- function(formula, data, 
                   tau = c(0.25, 0.5, 0.75),
                   max_p = 4,
                   max_q = 4,
                   max_k = 3,
                   criterion = c("BIC", "AIC", "HQ"),
                   case = 3,
                   bootstrap = FALSE,
                   n_boot = 1000,
                   seed = NULL,
                   verbose = TRUE,
                   fourier_k = NULL,
                   boot_scheme = c("bvz", "mcnown")) {
  
  # Validate inputs
  criterion <- match.arg(criterion)
  boot_scheme <- match.arg(boot_scheme)
  if (!is.null(fourier_k) &&
      (length(fourier_k) != 1 || fourier_k < 1 || fourier_k != round(fourier_k)))
    stop("'fourier_k' must be a single positive integer", call. = FALSE)
  
  if (!inherits(formula, "formula")) {
    stop("'formula' must be a formula object")
  }
  
  if (!is.data.frame(data)) {
    stop("'data' must be a data frame")
  }
  
  if (any(tau <= 0 | tau >= 1)) {
    stop("'tau' must be between 0 and 1 (exclusive)")
  }
  
  if (!identical(as.numeric(case), 3)) {
    stop("Only case = 3 (unrestricted intercept, no trend) is supported.\n",
         "  In fqardl <= 1.0.2 `case` was accepted but never used: the design\n",
         "  matrix was always Case III while the critical values were switched.\n",
         "  See NEWS.md for version 1.0.3.", call. = FALSE)
  }

  if (!is.null(seed)) set.seed(seed)
  
  # Extract variables from formula
  vars <- all.vars(formula)
  y_name <- vars[1]
  x_names <- vars[-1]
  
  # Check variables exist
  if (!all(vars %in% names(data))) {
    missing <- vars[!vars %in% names(data)]
    stop(paste("Variables not found in data:", paste(missing, collapse = ", ")))
  }
  
  # Prepare data
  y <- data[[y_name]]
  X <- as.matrix(data[, x_names, drop = FALSE])
  n <- length(y)
  
  # Check for missing values
  if (any(is.na(y)) || any(is.na(X))) {
    warning("Missing values detected. Using complete cases only.")
    complete <- complete.cases(y, X)
    y <- y[complete]
    X <- X[complete, , drop = FALSE]
    n <- length(y)
  }
  
  if (verbose) {
    message("=================================================================")
    message("   Fourier Quantile ARDL (FQARDL) Estimation")
    message("=================================================================\n")
    message(sprintf("Dependent variable: %s", y_name))
    message(sprintf("Independent variables: %s", paste(x_names, collapse = ", ")))
    message(sprintf("Sample size: %d", n))
    message(sprintf("Quantiles: %s\n", paste(tau, collapse = ", ")))
  }
  
  # Step 1: Select optimal Fourier frequency
  if (verbose) message("Step 1: Selecting optimal Fourier frequency...")
  if (is.null(fourier_k)) {
    k_results <- select_fourier_frequency(y, X, max_k, criterion)
    optimal_k <- k_results$optimal_k
    if (verbose) message(sprintf("   Optimal k = %d (based on %s)\n", optimal_k, criterion))
  } else {
    k_results <- list(optimal_k = as.integer(fourier_k), ic_values = NA_real_,
                      criterion = "fixed by the user")
    optimal_k <- k_results$optimal_k
    if (verbose) message(sprintf("   k = %d (fixed by the user)\n", optimal_k))
  }
  
  # Step 2: Generate Fourier terms
  if (verbose) message("Step 2: Generating Fourier terms...")
  fourier_terms <- generate_fourier_terms(n, optimal_k)
  if (verbose) message(sprintf("   Generated sin and cos terms for k = %d\n", optimal_k))
  
  # Step 3: Select optimal lag structure
  if (verbose) message("Step 3: Selecting optimal lag structure...")
  lag_results <- select_optimal_lags(y, X, fourier_terms, max_p, max_q, criterion)
  optimal_p <- lag_results$optimal_p
  optimal_q <- lag_results$optimal_q
  if (verbose) message(sprintf("   Optimal lags: p = %d, q = %d (based on %s)\n", 
              optimal_p, optimal_q, criterion))
  
  # Step 4: Estimate QARDL for each quantile
  if (verbose) message("Step 4: Estimating Quantile ARDL for each tau...")
  qardl_results <- list()
  
  for (i in seq_along(tau)) {
    if (verbose) message(sprintf("   Estimating tau = %.2f...", tau[i]))
    qardl_results[[i]] <- estimate_qardl(
      y = y,
      X = X,
      fourier = fourier_terms,
      p = optimal_p,
      q = optimal_q,
      tau = tau[i],
      case = case
    )
  }
  names(qardl_results) <- paste0("tau_", tau)
  
  # Step 5: Calculate long-run and short-run multipliers
  if (verbose) message("Step 5: Computing multipliers...")
  multipliers <- compute_multipliers(qardl_results, x_names, tau)
  if (verbose) message("   Long-run and short-run multipliers computed.\n")
  
  # Step 6: Bounds test for cointegration
  if (verbose) message("Step 6: Performing bounds test for cointegration...")
  bounds_results <- perform_bounds_test(qardl_results, n, length(x_names), case)
  if (!bootstrap && isFALSE(bounds_results$bounds_valid))
    bounds_results$decision <- paste0(bounds_results$decision, fourier_boot_advice)
  if (verbose) {
    message(sprintf("   F-statistic: %.4f", bounds_results$F_stat))
    message(sprintf("   t-statistic: %.4f", bounds_results$t_stat))
    message(sprintf("   Decision: %s\n", bounds_results$decision))
  }
  
  # Step 7: Bootstrap cointegration test (if requested)
  bootstrap_results <- NULL
  if (bootstrap) {
    if (verbose) message("Step 7: Bootstrap cointegration testing...")
    bootstrap_results <- bootstrap_bounds_test(
      y, X, fourier_terms, optimal_p, optimal_q,
      tau, case, n_boot, verbose = verbose,
      reselect = list(k = fourier_k, max_k = max_k, max_p = max_p,
                      max_q = max_q, criterion = criterion),
      scheme = boot_scheme
    )
    if (verbose) {
      message(sprintf("   Bootstrap p-value (Fov): %.4f", bootstrap_results$p_value_F))
      message(sprintf("   Bootstrap p-value (t): %.4f", bootstrap_results$p_value_t))
      message(sprintf("   Bootstrap p-value (Find): %.4f", bootstrap_results$p_value_Find))
      message(sprintf("   Decision: %s\n", bootstrap_results$label))
    }
  }
  
  # Step 8: Diagnostics
  if (verbose) message("Step 8: Computing diagnostics...")
  diagnostics <- compute_diagnostics(qardl_results)
  if (verbose) message("   Diagnostics computed.\n")
  
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
    n = n,
    tau = tau,
    optimal_k = optimal_k,
    optimal_p = optimal_p,
    optimal_q = optimal_q,
    case = case,
    criterion = criterion,
    qardl_results = qardl_results,
    long_run = multipliers$long_run,
    short_run = multipliers$short_run,
    bounds_test = bounds_results,
    bootstrap = bootstrap_results,
    diagnostics = diagnostics,
    fourier_terms = fourier_terms,
    k_selection = k_results,
    lag_selection = lag_results
  )
  
  class(result) <- "fqardl"
  return(result)
}


#' @export
print.fqardl <- function(x, ...) {
  cat("\nFourier Quantile ARDL Model\n")
  cat("===========================\n\n")
  cat("Call:\n")
  print(x$call)
  cat("\n")
  cat(sprintf("Sample size: %d\n", x$n))
  cat(sprintf("Dependent variable: %s\n", x$y_name))
  cat(sprintf("Independent variables: %s\n", paste(x$x_names, collapse = ", ")))
  cat(sprintf("Optimal Fourier frequency: k = %d\n", x$optimal_k))
  cat(sprintf("Optimal lags: p = %d, q = %d\n", x$optimal_p, x$optimal_q))
  cat(sprintf("Model case: %d\n", x$case))
  cat(sprintf("Quantiles: %s\n", paste(x$tau, collapse = ", ")))
  cat("\n")
  cat("Bounds Test:\n")
  cat(sprintf("  F-statistic: %.4f\n", x$bounds_test$F_stat))
  cat(sprintf("  Decision: %s\n", x$bounds_test$decision))
  if (!is.null(x$bootstrap)) {
    cat("Bootstrap bounds test:\n")
    cat(sprintf("  p-values: Fov = %.4f, t = %.4f, Find = %.4f\n",
                x$bootstrap$p_value_F, x$bootstrap$p_value_t,
                x$bootstrap$p_value_Find))
    cat(sprintf("  Decision: %s\n", x$bootstrap$label))
  }
  invisible(x)
}


#' @export
summary.fqardl <- function(object, ...) {
  cat("\n")
  cat("=================================================================\n")
  cat("        Fourier Quantile ARDL (FQARDL) - Summary Results\n")
  cat("=================================================================\n\n")
  
  # Model specification
  cat("MODEL SPECIFICATION\n")
  cat("-------------------\n")
  cat(sprintf("Dependent variable: %s\n", object$y_name))
  cat(sprintf("Independent variables: %s\n", paste(object$x_names, collapse = ", ")))
  cat(sprintf("Sample size: %d\n", object$n))
  cat(sprintf("Fourier frequency (k): %d\n", object$optimal_k))
  cat(sprintf("Lag structure: ARDL(%d, %s)\n", 
              object$optimal_p, 
              paste(rep(object$optimal_q, length(object$x_names)), collapse = ", ")))
  cat(sprintf("Selection criterion: %s\n", object$criterion))
  cat(sprintf("Model case: %d\n\n", object$case))
  
  # Long-run coefficients
  cat("LONG-RUN COEFFICIENTS\n")
  cat("---------------------\n")
  print(round(object$long_run, 4))
  cat("\n")
  
  # Short-run coefficients  
  cat("SHORT-RUN COEFFICIENTS (ECT)\n")
  cat("----------------------------\n")
  print(round(object$short_run, 4))
  cat("\n")
  
  # Bounds test
  cat("BOUNDS TEST FOR COINTEGRATION\n")
  cat("-----------------------------\n")
  cat(sprintf("F-statistic: %.4f\n", object$bounds_test$F_stat))
  cat(sprintf("t-statistic: %.4f\n", object$bounds_test$t_stat))
  cat("\nCritical Values (Pesaran, Shin and Smith 2001, Table CI(iii)):\n")
  cat(sprintf("   5%%: I(0) = %.3f, I(1) = %.3f\n",
              object$bounds_test$cv_5[1], object$bounds_test$cv_5[2]))
  if (isFALSE(object$bounds_test$bounds_valid))
    cat(sprintf("   %s\n", fourier_bounds_note))
  cat("   1% and 10%: not supplied. Only the 5 percent column of the table has\n")
  cat("   been verified against the source; see ?get_pss_critical_values.\n")
  cat(sprintf("\nDecision: %s\n", object$bounds_test$decision))
  cat(sprintf("(reported at tau = %.2f; restrictions = %d)\n",
              object$bounds_test$reference_tau %||% NA_real_,
              object$bounds_test$n_restrictions %||% NA_integer_))
  if (!is.null(object$bounds_test$summary)) {
    cat("\nBounds test by quantile:\n")
    print(object$bounds_test$summary, row.names = FALSE)
  }
  cat("\nNote: the PSS critical values were simulated for the conditional mean\n")
  cat("of a model without Fourier terms. The quantile statistics are descriptive\n")
  cat("(no published null distribution); use bootstrap = TRUE for inference.\n\n")
  
  # Bootstrap results if available
  if (!is.null(object$bootstrap)) {
    cat("BOOTSTRAP COINTEGRATION TEST\n")
    cat("----------------------------\n")
    if (is.null(object$bootstrap$engine)) {
      cat(sprintf("Number of replications: %d\n", object$bootstrap$n_boot))
      cat(sprintf("Bootstrap p-value (F): %.4f\n", object$bootstrap$p_value_F))
      cat(sprintf("Bootstrap p-value (t): %.4f\n", object$bootstrap$p_value_t))
    } else {
      .fq_print_boot(object$bootstrap)
    }
    cat("\n")
  }
  
  # Diagnostics
  cat("DIAGNOSTICS\n")
  cat("-----------\n")
  for (i in seq_along(object$tau)) {
    tau_name <- names(object$diagnostics)[i]
    diag <- object$diagnostics[[i]]
    cat(sprintf("\nQuantile tau = %.2f:\n", object$tau[i]))
    cat(sprintf("  Pseudo R-squared: %.4f\n", diag$pseudo_r2))
    cat(sprintf("  Log-likelihood: %.4f\n", diag$loglik))
  }
  
  cat("\n=================================================================\n")
  invisible(object)
}

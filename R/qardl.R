# =============================================================================
# Quantile ARDL Estimation Functions
# Based on Cho et al. (2015) and extensions
# =============================================================================

#' Select Optimal Lag Structure
#'
#' @description
#' Selects optimal lag orders for ARDL model using grid search
#' over all combinations of p and q.
#'
#' @param y Dependent variable
#' @param X Matrix of independent variables
#' @param fourier Fourier terms matrix
#' @param max_p Maximum lag for dependent variable
#' @param max_q Maximum lag for independent variables
#' @param criterion Information criterion
#'
#' @return List with optimal lags
#'
#' @export
select_optimal_lags <- function(y, X, fourier, max_p, max_q, 
                                criterion = c("BIC", "AIC", "HQ")) {
  
  criterion <- match.arg(criterion)
  n <- length(y)
  k_vars <- ncol(X)
  
  best_ic <- Inf
  optimal_p <- 1
  optimal_q <- 1
  
  ic_matrix <- matrix(NA, nrow = max_p, ncol = max_q)
  
  for (p in 1:max_p) {
    for (q in 1:max_q) {
      
      # Calculate effective sample size
      max_lag <- max(p, q)
      n_eff <- n - max_lag
      
      if (n_eff < 20) next  # Skip if insufficient observations
      
      # Build design matrix
      design <- build_ardl_design(y, X, fourier, p, q)
      
      if (is.null(design)) next
      
      # Estimate OLS for IC comparison
      tryCatch({
        model <- lm(design$y ~ design$X - 1)
        
        loglik <- logLik(model)
        n_params <- ncol(design$X)
        
        ic <- switch(criterion,
          "AIC" = -2 * loglik + 2 * n_params,
          "BIC" = -2 * loglik + log(n_eff) * n_params,
          "HQ"  = -2 * loglik + 2 * log(log(n_eff)) * n_params
        )
        
        ic_matrix[p, q] <- ic
        
        if (ic < best_ic) {
          best_ic <- ic
          optimal_p <- p
          optimal_q <- q
        }
        
      }, error = function(e) NULL)
    }
  }
  
  return(list(
    optimal_p = optimal_p,
    optimal_q = optimal_q,
    best_ic = best_ic,
    ic_matrix = ic_matrix,
    criterion = criterion
  ))
}


#' Build ARDL Design Matrix
#'
#' @description
#' Constructs the design matrix for ARDL estimation including
#' lagged dependent and independent variables.
#'
#' @param y Dependent variable
#' @param X Independent variables matrix
#' @param fourier Fourier terms
#' @param p Lag order for y
#' @param q Lag order for X
#'
#' @return List with y and X design matrices
#'
#' @keywords internal
build_ardl_design <- function(y, X, fourier, p, q) {
  
  n <- length(y)
  k <- ncol(X)
  max_lag <- max(p, q)
  
  if (n <= max_lag + 1) return(NULL)
  
  # Effective sample
  n_eff <- n - max_lag
  
  # Dependent variable (adjusted)
  # CORRECTED in 1.0.3: the dependent variable of a PSS conditional error
  # correction model is Delta y, not y. Up to 1.0.2 the right-hand side was
  # the ECM design while the left-hand side was in levels, so the reported
  # coefficient on y_lag1 was 1 + rho instead of rho: positive, near unity,
  # with a large positive t-ratio at every quantile, and the long-run
  # multipliers -theta/phi carried the wrong sign.
  y_adj <- y[(max_lag + 1):n] - y[max_lag:(n - 1)]
  
  # Initialize design matrix list
  design_list <- list()
  col_names <- c()
  
  # Intercept
  design_list[["intercept"]] <- rep(1, n_eff)
  col_names <- c(col_names, "intercept")
  
  # Lagged y (levels and differences for ECM form)
  y_lag1 <- y[max_lag:(n - 1)]
  design_list[["y_lag1"]] <- y_lag1
  col_names <- c(col_names, "y_lag1")
  
  # Additional lags of first differences of y
  if (p > 1) {
    dy <- diff(y)
    for (j in 1:(p - 1)) {
      lag_idx <- (max_lag - j):(n - 1 - j)
      if (length(lag_idx) == n_eff) {
        design_list[[paste0("dy_lag", j)]] <- dy[lag_idx]
        col_names <- c(col_names, paste0("dy_lag", j))
      }
    }
  }
  
  # Independent variables in levels (for long-run)
  for (i in 1:k) {
    x_lag1 <- X[max_lag:(n - 1), i]
    design_list[[paste0("x", i, "_lag1")]] <- x_lag1
    col_names <- c(col_names, paste0("x", i, "_lag1"))
  }
  
  # First differences of X (current and lagged)
  dX <- diff(X)
  for (i in 1:k) {
    # Current difference
    dx_current <- dX[max_lag:(n - 1), i]
    design_list[[paste0("dx", i)]] <- dx_current
    col_names <- c(col_names, paste0("dx", i))
    
    # Lagged differences
    if (q > 1) {
      for (j in 1:(q - 1)) {
        lag_idx <- (max_lag - j):(n - 1 - j)
        if (length(lag_idx) == n_eff) {
          design_list[[paste0("dx", i, "_lag", j)]] <- dX[lag_idx, i]
          col_names <- c(col_names, paste0("dx", i, "_lag", j))
        }
      }
    }
  }
  
  # Fourier terms (adjusted)
  fourier_adj <- fourier[(max_lag + 1):n, , drop = FALSE]
  for (j in 1:ncol(fourier_adj)) {
    design_list[[colnames(fourier_adj)[j]]] <- fourier_adj[, j]
    col_names <- c(col_names, colnames(fourier_adj)[j])
  }
  
  # Combine into matrix
  X_design <- do.call(cbind, design_list)
  colnames(X_design) <- col_names
  
  return(list(
    y = y_adj,
    X = X_design,
    n = n_eff,
    max_lag = max_lag
  ))
}


#' Estimate Quantile ARDL Model
#'
#' @description
#' Estimates the Quantile ARDL model for a given quantile tau.
#'
#' @param y Dependent variable
#' @param X Independent variables
#' @param fourier Fourier terms
#' @param p Lag for dependent variable
#' @param q Lag for independent variables
#' @param tau Quantile (0 < tau < 1)
#' @param case Model case (1-5)
#'
#' @return List with estimation results
#'
#' @export
estimate_qardl <- function(y, X, fourier, p, q, tau, case = 3) {
  
  # Build design matrix
  design <- build_ardl_design(y, X, fourier, p, q)
  
  if (is.null(design)) {
    stop("Insufficient observations for the specified lag structure")
  }
  
  # Estimate quantile regression
  qr_model <- quantreg::rq(design$y ~ design$X - 1, tau = tau)
  
  # Extract coefficients
  coefs <- coef(qr_model)
  names(coefs) <- colnames(design$X)
  
  # Standard errors (using sandwich estimator)
  qr_summary <- summary(qr_model, se = "boot", R = 200, covariance = TRUE)
  se <- qr_summary$coefficients[, 2]
  names(se) <- colnames(design$X)   # FIXED in 1.0.3: in 1.0.2 only `coefs` was
                                    # renamed, so se["y_lag1"] silently returned NA
  # NEW in 1.0.3: retained so that perform_bounds_test() can form a genuine
  # Wald statistic instead of an average of squared marginal t-ratios.
  V <- tryCatch(qr_summary$cov, error = function(e) NULL)
  if (!is.null(V)) dimnames(V) <- list(colnames(design$X), colnames(design$X))
  
  # t-statistics
  t_stats <- coefs / se
  names(t_stats) <- colnames(design$X)
  
  # p-values
  p_values <- 2 * (1 - pnorm(abs(t_stats)))
  
  # Fitted values and residuals
  fitted <- as.numeric(design$X %*% coefs)
  residuals <- design$y - fitted
  
  # Pseudo R-squared
  rho <- function(u, tau) u * (tau - (u < 0))
  pseudo_r2 <- 1 - sum(rho(residuals, tau)) / sum(rho(design$y - quantile(design$y, tau), tau))
  
  # Log-likelihood (quantile regression)
  loglik <- -sum(rho(residuals, tau))
  
  return(list(
    tau = tau,
    coefficients = coefs,
    vcov = V,
    std_errors = se,
    t_statistics = t_stats,
    p_values = p_values,
    fitted = fitted,
    residuals = residuals,
    pseudo_r2 = pseudo_r2,
    loglik = loglik,
    n = design$n,
    design = design,
    model = qr_model
  ))
}


#' Compute Long-run and Short-run Multipliers
#'
#' @description
#' Calculates the long-run and short-run multipliers from QARDL estimates.
#'
#' @param qardl_results List of QARDL results for each quantile
#' @param x_names Names of independent variables
#' @param tau Vector of quantiles
#'
#' @return List with long-run and short-run multiplier matrices
#'
#' @export
compute_multipliers <- function(qardl_results, x_names, tau) {
  
  n_tau <- length(tau)
  n_x <- length(x_names)
  
  # Initialize matrices
  long_run <- matrix(NA, nrow = n_tau, ncol = n_x)
  short_run <- matrix(NA, nrow = n_tau, ncol = n_x)
  
  rownames(long_run) <- paste0("tau_", tau)
  rownames(short_run) <- paste0("tau_", tau)
  colnames(long_run) <- x_names
  colnames(short_run) <- x_names
  
  for (i in 1:n_tau) {
    coefs <- qardl_results[[i]]$coefficients
    
    # Get coefficient for y_lag1 (speed of adjustment)
    phi <- coefs["y_lag1"]
    
    for (j in 1:n_x) {
      # Long-run coefficient: -beta_x / phi
      x_lag_name <- paste0("x", j, "_lag1")
      if (x_lag_name %in% names(coefs)) {
        beta_x <- coefs[x_lag_name]
        long_run[i, j] <- -beta_x / phi
      }
      
      # Short-run coefficient (current difference)
      dx_name <- paste0("dx", j)
      if (dx_name %in% names(coefs)) {
        short_run[i, j] <- coefs[dx_name]
      }
    }
  }
  
  # ECT (Error Correction Term) = phi (speed of adjustment)
  ect <- sapply(qardl_results, function(x) x$coefficients["y_lag1"])
  names(ect) <- paste0("tau_", tau)
  
  return(list(
    long_run = long_run,
    short_run = short_run,
    ect = ect
  ))
}


#' Compute Model Diagnostics
#'
#' @description
#' Computes various diagnostic statistics for the QARDL models.
#'
#' @param qardl_results List of QARDL estimation results
#'
#' @return List of diagnostics for each quantile
#'
#' @export
compute_diagnostics <- function(qardl_results) {
  
  diagnostics <- lapply(qardl_results, function(res) {
    
    # Basic stats
    n <- res$n
    k <- length(res$coefficients)
    
    # Pseudo R-squared
    pseudo_r2 <- res$pseudo_r2
    
    # Log-likelihood
    loglik <- res$loglik
    
    # Information criteria
    aic <- -2 * loglik + 2 * k
    bic <- -2 * loglik + log(n) * k
    
    # Residual statistics
    resid_mean <- mean(res$residuals)
    resid_sd <- sd(res$residuals)
    
    list(
      n = n,
      k = k,
      pseudo_r2 = pseudo_r2,
      loglik = loglik,
      AIC = aic,
      BIC = bic,
      resid_mean = resid_mean,
      resid_sd = resid_sd
    )
  })
  
  return(diagnostics)
}

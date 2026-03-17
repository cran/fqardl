#' =============================================================================
#' Fourier Approximation Functions
#' For capturing smooth structural breaks
#' =============================================================================

#' Generate Fourier Trigonometric Terms
#'
#' @description
#' Creates sine and cosine terms for Fourier approximation of structural breaks.
#' Based on Enders & Lee (2012) methodology.
#'
#' @param n Sample size (number of observations)
#' @param k Fourier frequency (integer >= 1)
#' @param cumulative If TRUE, includes all frequencies from 1 to k
#'
#' @return A matrix with sine and cosine columns
#'
#' @details
#' The Fourier terms are computed as:
#' \deqn{sin(2\pi k t / T)}
#' \deqn{cos(2\pi k t / T)}
#' where t is the time index and T is the sample size.
#'
#' @export
generate_fourier_terms <- function(n, k, cumulative = FALSE) {
  
  if (k < 1 || k != round(k)) {
    stop("k must be a positive integer")
  }
  
  # Time index
  t <- 1:n
  
  if (cumulative) {
    # Include all frequencies from 1 to k
    fourier_mat <- matrix(NA, nrow = n, ncol = 2 * k)
    col_names <- character(2 * k)
    
    for (j in 1:k) {
      fourier_mat[, 2*j - 1] <- sin(2 * pi * j * t / n)
      fourier_mat[, 2*j] <- cos(2 * pi * j * t / n)
      col_names[2*j - 1] <- paste0("sin_", j)
      col_names[2*j] <- paste0("cos_", j)
    }
    colnames(fourier_mat) <- col_names
    
  } else {
    # Only include frequency k
    fourier_mat <- cbind(
      sin_k = sin(2 * pi * k * t / n),
      cos_k = cos(2 * pi * k * t / n)
    )
    colnames(fourier_mat) <- c(paste0("sin_", k), paste0("cos_", k))
  }
  
  return(fourier_mat)
}


#' Select Optimal Fourier Frequency
#'
#' @description
#' Selects the optimal Fourier frequency k based on information criteria.
#' Tests all frequencies from 1 to max_k and selects the one minimizing
#' the chosen criterion.
#'
#' @param y Dependent variable vector
#' @param X Matrix of independent variables
#' @param max_k Maximum Fourier frequency to test
#' @param criterion Information criterion ("AIC", "BIC", "HQ")
#'
#' @return A list containing:
#' \item{optimal_k}{The optimal Fourier frequency}
#' \item{ic_values}{Information criterion values for each k}
#' \item{criterion}{The criterion used}
#'
#' @export
select_fourier_frequency <- function(y, X, max_k = 3, 
                                     criterion = c("BIC", "AIC", "HQ")) {
  
  criterion <- match.arg(criterion)
  n <- length(y)
  
  # Storage for IC values
  ic_values <- numeric(max_k)
  
  for (k in 1:max_k) {
    # Generate Fourier terms
    fourier <- generate_fourier_terms(n, k)
    
    # Combine regressors
    X_full <- cbind(X, fourier)
    
    # Estimate simple OLS for IC comparison
    model <- lm(y ~ X_full)
    
    # Calculate information criterion
    loglik <- logLik(model)
    n_params <- length(coef(model))
    
    ic_values[k] <- switch(criterion,
      "AIC" = -2 * loglik + 2 * n_params,
      "BIC" = -2 * loglik + log(n) * n_params,
      "HQ"  = -2 * loglik + 2 * log(log(n)) * n_params
    )
  }
  
  # Select optimal k
  optimal_k <- which.min(ic_values)
  
  return(list(
    optimal_k = optimal_k,
    ic_values = ic_values,
    criterion = criterion
  ))
}


#' Fourier ADF Unit Root Test
#'
#' @description
#' Performs the Fourier Augmented Dickey-Fuller test for unit roots
#' with smooth structural breaks.
#'
#' @param y Time series vector
#' @param max_k Maximum Fourier frequency
#' @param max_lag Maximum number of lags for ADF
#' @param criterion Lag selection criterion
#'
#' @return A list with test results
#'
#' @export
fourier_adf <- function(y, max_k = 3, max_lag = 8, 
                        criterion = c("BIC", "AIC")) {
  
  criterion <- match.arg(criterion)
  n <- length(y)
  
  # First differences
  dy <- diff(y)
  y_lag <- y[-n]
  
  best_ic <- Inf
  best_result <- NULL
  
  for (k in 1:max_k) {
    # Generate Fourier terms
    fourier <- generate_fourier_terms(n - 1, k)
    
    for (p in 0:max_lag) {
      # Create lagged differences
      if (p > 0) {
        dy_lags <- embed(dy, p + 1)[, -1, drop = FALSE]
        colnames(dy_lags) <- paste0("dy_lag", 1:p)
        
        # Adjust vectors
        dy_adj <- dy[(p + 1):(n - 1)]
        y_lag_adj <- y_lag[(p + 1):(n - 1)]
        fourier_adj <- fourier[(p + 1):(n - 1), , drop = FALSE]
        
        X <- cbind(y_lag_adj, fourier_adj, dy_lags)
      } else {
        dy_adj <- dy
        X <- cbind(y_lag, fourier)
      }
      
      # Estimate model
      model <- lm(dy_adj ~ X)
      loglik <- logLik(model)
      n_obs <- length(dy_adj)
      n_params <- length(coef(model))
      
      # Information criterion
      ic <- switch(criterion,
        "AIC" = -2 * loglik + 2 * n_params,
        "BIC" = -2 * loglik + log(n_obs) * n_params
      )
      
      if (ic < best_ic) {
        best_ic <- ic
        
        # Extract t-statistic for y_lag coefficient
        coef_summary <- summary(model)$coefficients
        t_stat <- coef_summary[2, 3]  # t-stat for y_lag
        
        best_result <- list(
          k = k,
          lag = p,
          t_stat = t_stat,
          ic = ic,
          model = model
        )
      }
    }
  }
  
  # Critical values (Enders & Lee, 2012 - approximate)
  cv <- data.frame(
    significance = c("1%", "5%", "10%"),
    cv = c(-4.42, -3.81, -3.49)  # Approximate for k=1
  )
  
  # Decision
  if (best_result$t_stat < cv$cv[2]) {
    decision <- "Reject H0: Unit root (Series is stationary)"
  } else {
    decision <- "Fail to reject H0: Unit root exists"
  }
  
  result <- list(
    optimal_k = best_result$k,
    optimal_lag = best_result$lag,
    t_statistic = best_result$t_stat,
    critical_values = cv,
    decision = decision,
    ic = best_result$ic,
    criterion = criterion
  )
  
  class(result) <- "fadf"
  return(result)
}


#' @export
print.fadf <- function(x, ...) {
  cat("\nFourier ADF Unit Root Test\n")
  cat("==========================\n\n")
  cat(sprintf("Optimal Fourier frequency: k = %d\n", x$optimal_k))
  cat(sprintf("Optimal lag length: %d\n", x$optimal_lag))
  cat(sprintf("Test statistic: %.4f\n\n", x$t_statistic))
  cat("Critical values:\n")
  print(x$critical_values, row.names = FALSE)
  cat(sprintf("\nDecision: %s\n", x$decision))
  invisible(x)
}

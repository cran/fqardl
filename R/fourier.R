# =============================================================================
# Fourier Approximation Functions
# For capturing smooth structural breaks
# =============================================================================

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
#' @return The result of \code{\link{fourier_adf_test}} with a constant,
#'   with the additional fields \code{optimal_k} and \code{t_statistic}.
#'
#' @export
fourier_adf <- function(y, max_k = 3, max_lag = 8,
                        criterion = c("BIC", "AIC")) {
  criterion <- match.arg(criterion)
  # Same test as fourier_adf_test() with a constant; kept for compatibility
  res <- fourier_adf_test(y, model = "c", max_freq = max_k, max_lag = max_lag,
                          criterion = criterion, verbose = FALSE)
  res$optimal_k <- res$optimal_frequency
  res$t_statistic <- res$statistic
  res$criterion <- criterion
  res
}

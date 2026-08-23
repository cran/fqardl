
#' Genuine Wald bounds F statistic from an OLS conditional ECM
#'
#' Replaces the 1.0.2 construction \code{mean(t^2)}, which averaged squared
#' marginal t-ratios and so ignored the covariances among the level
#' coefficients. See PSS (2001) eq. (21).
#'
#' @keywords internal
wald_bounds_F <- function(coefs, V, level_names) {
  if (is.null(V)) return(NA_real_)
  pos <- match(level_names, names(coefs))
  if (anyNA(pos) || any(dim(V) != length(coefs))) return(NA_real_)
  R <- matrix(0, nrow = length(pos), ncol = length(coefs))
  R[cbind(seq_along(pos), pos)] <- 1
  Rb <- R %*% coefs
  W <- tryCatch(as.numeric(t(Rb) %*% solve(R %*% V %*% t(R)) %*% Rb),
                error = function(e) NA_real_)
  W / length(pos)
}

#' Three-way bounds verdict, reporting the inconclusive region
#' @keywords internal
bounds_verdict <- function(stat, lower, upper) {
  if (is.na(stat)) return("Statistic not available")
  if (stat > upper)  return("Cointegration (above the I(1) bound at 5%)")
  if (stat < lower)  return("No cointegration (below the I(0) bound at 5%)")
  "Inconclusive (inside the bounds at 5%)"
}

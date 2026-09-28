# =============================================================================
# Bounds Test for Cointegration
# Pesaran, Shin & Smith (2001), Journal of Applied Econometrics 16(3), 289-326.
#
# CONTAINMENT RELEASE 1.0.3 -- this file was rewritten. See NEWS.md.
# In 1.0.2 the "F-statistic" was mean(t^2) over the terms matching `_lag1$`.
# That is an average of squared marginal t-ratios: it ignores the covariances
# among the level coefficients, it has no distribution theory, and it was being
# compared against tables simulated for a genuine Wald statistic. It has been
# replaced by the Wald form of PSS eq. (21).
# =============================================================================

#' Wald Bounds Test for Cointegration
#'
#' @description
#' Computes the Pesaran, Shin and Smith (2001) bounds test from an estimated
#' conditional error correction model, at every quantile supplied.
#'
#' The F statistic is the Wald statistic for the joint null that all lagged
#' level coefficients are zero, divided by the number of restrictions:
#' \deqn{W = (Rb)' (R V R')^{-1} (Rb), \qquad F = W / m}
#' with \eqn{m = k+1} in Cases I, III and V and \eqn{m = k+2} in Cases II and IV.
#' The t statistic is the ratio on the lagged level of the dependent variable,
#' which is negative under cointegration.
#'
#' @param qardl_results List of quantile estimation results from
#'   \code{estimate_qardl}. The test is computed for every element.
#' @param n Sample size.
#' @param k Number of regressors (excluding the dependent variable).
#' @param case Model case. Only \code{case = 3} is supported; see Details.
#'
#' @return A list with one element per quantile plus a \code{summary} data frame.
#'   Each per-quantile element contains \code{F_stat}, \code{t_stat}, the
#'   critical value bounds, and a three-way \code{decision} that reports the
#'   inconclusive region rather than forcing a binary verdict.
#'
#' @details
#' \strong{Only Case III is supported.} In versions up to 1.0.2 the \code{case}
#' argument was accepted and never used: the design matrix was always that of
#' Case III (unrestricted intercept, no trend), while the critical values were
#' switched. Cases 4 and 5 therefore compared a Case III model against Case IV
#' and Case V tables. Rather than continue to report that, this release raises an
#' error for any case other than 3. Proper case handling requires changing the
#' design matrix, not the table, and is scheduled for 2.0.0.
#'
#' \strong{Quantile caveat.} The PSS critical values were simulated for the
#' conditional mean. Applying them quantile by quantile has no distribution
#' theory behind it. Treat quantile-specific verdicts as descriptive and prefer
#' \code{bootstrap_bounds_test} for inference.
#'
#' @references
#' Pesaran, M.H., Shin, Y. and Smith, R.J. (2001). Bounds testing approaches to
#' the analysis of level relationships. \emph{Journal of Applied Econometrics},
#' 16(3), 289-326. \doi{10.1002/jae.616}
#'
#' @export
perform_bounds_test <- function(qardl_results, n, k, case = 3) {

  if (!identical(as.numeric(case), 3)) {
    stop(
      "Only case = 3 (unrestricted intercept, no trend) is supported.\n",
      "  In fqardl <= 1.0.2 the `case` argument was accepted but never used: the\n",
      "  design matrix was always Case III while the critical values were switched,\n",
      "  so cases 4 and 5 compared a Case III model against Case IV/V tables.\n",
      "  Correct case handling changes the regressors, not the table, and is\n",
      "  scheduled for version 2.0.0.",
      call. = FALSE
    )
  }

  cv <- get_pss_critical_values(k, case = 3)
  m_restr <- k + 1L   # Cases I, III, V. Cases II and IV would use k + 2.

  per_tau <- lapply(qardl_results, function(res)
    bounds_one_quantile(res, k = k, m_restr = m_restr, cv = cv))

  summary_df <- data.frame(
    tau      = vapply(qardl_results, function(r) r$tau,          numeric(1)),
    F_stat   = vapply(per_tau,       function(r) r$F_stat,       numeric(1)),
    t_stat   = vapply(per_tau,       function(r) r$t_stat,       numeric(1)),
    ect      = vapply(per_tau,       function(r) r$ect,          numeric(1)),
    decision = vapply(per_tau,       function(r) r$decision,     character(1)),
    stringsAsFactors = FALSE,
    row.names = NULL
  )

  # The median (or, failing that, the first) quantile is promoted to the top
  # level so that code written against the 1.0.2 return shape keeps working.
  idx <- which(summary_df$tau == 0.5)
  if (length(idx) == 0) idx <- 1L

  out <- per_tau[[idx]]
  out$by_quantile   <- per_tau
  out$summary       <- summary_df
  out$reference_tau <- summary_df$tau[idx]
  out$case          <- 3
  out$k             <- k
  out$n             <- n
  out
}


#' Bounds test at a single quantile
#' @keywords internal
bounds_one_quantile <- function(res, k, m_restr, cv) {

  cf  <- res$coefficients
  nms <- names(cf)

  y_lev <- "y_lag1"
  x_lev <- grep("^x[0-9]+_lag1$", nms, value = TRUE)
  lev   <- c(y_lev, x_lev)

  if (!all(lev %in% nms))
    stop("Level terms not found in the design; cannot form the bounds test.",
         call. = FALSE)

  V <- res$vcov
  F_stat <- NA_real_

  if (!is.null(V) && all(dim(V) == length(cf))) {
    pos <- match(lev, nms)
    R <- matrix(0, nrow = length(pos), ncol = length(cf))
    R[cbind(seq_along(pos), pos)] <- 1
    Rb  <- R %*% cf
    RVR <- R %*% V %*% t(R)
    W <- tryCatch(as.numeric(t(Rb) %*% solve(RVR) %*% Rb),
                  error = function(e) NA_real_)
    F_stat <- W / m_restr
  }

  ect    <- unname(cf[y_lev])
  t_stat <- unname(ect / res$std_errors[y_lev])

  # Three-way verdict. The inconclusive region is part of the procedure and
  # must be reported; 1.0.2 collapsed it into a binary answer.
  dec_F <- if (is.na(F_stat)) "F not available"
           else if (F_stat > cv$F_upper[2]) "cointegration"
           else if (F_stat < cv$F_lower[2]) "no cointegration"
           else "inconclusive"

  dec_t <- if (is.na(t_stat)) "t not available"
           else if (t_stat < cv$t_upper[2]) "cointegration"
           else if (t_stat > cv$t_lower[2]) "no cointegration"
           else "inconclusive"

  decision <- if (dec_F == "cointegration" && dec_t == "cointegration")
                "Cointegration (both F and t reject at 5%)"
              else if (dec_F == "no cointegration" || dec_t == "no cointegration")
                "No cointegration at 5%"
              else if (dec_F == "inconclusive" || dec_t == "inconclusive")
                "Inconclusive (statistic inside the bounds)"
              else
                "Mixed evidence; use bootstrap_bounds_test"

  list(
    tau         = res$tau,
    F_stat      = F_stat,
    t_stat      = t_stat,
    ect         = ect,
    cv_1        = c(cv$F_lower[1], cv$F_upper[1]),
    cv_5        = c(cv$F_lower[2], cv$F_upper[2]),
    cv_10       = c(cv$F_lower[3], cv$F_upper[3]),
    t_cv_1      = c(cv$t_lower[1], cv$t_upper[1]),
    t_cv_5      = c(cv$t_lower[2], cv$t_upper[2]),
    t_cv_10     = c(cv$t_lower[3], cv$t_upper[3]),
    decision_F  = dec_F,
    decision_t  = dec_t,
    decision    = decision,
    n_restrictions = m_restr
  )
}


#' PSS (2001) Critical Values, Case III
#'
#' @description
#' Table CI(iii) for the F statistic and Table CII(iii) for the t statistic,
#' unrestricted intercept and no trend. Values are the 1 percent, 5 percent and
#' 10 percent points, each with an I(0) lower and an I(1) upper bound.
#'
#' In Table CII the I(0) bound does not depend on k.
#'
#' @param k Number of regressors, excluding the dependent variable.
#' @param case Deterministic case. Only 3 is supported.
#'
#' @return A list with \code{F_lower}, \code{F_upper}, \code{t_lower},
#'   \code{t_upper} and \code{significance}.
#'
#' @references
#' Pesaran, M.H., Shin, Y. and Smith, R.J. (2001), Tables CI(iii) and CII(iii).
#'
#' @keywords internal
get_pss_critical_values <- function(k, case = 3) {

  if (!identical(as.numeric(case), 3))
    stop("Only case = 3 is supported in this release.", call. = FALSE)

  kk <- max(1L, min(10L, as.integer(k)))

  # Table CI(iii), 5 percent points only, I(0) lower then I(1) upper, k = 1..10.
  # Verified line by line against the printed table and checked by a test.
  #
  # The 1 and 10 percent points are deliberately NOT supplied. An adversarial
  # audit of 1.0.3 found that the entries carried there had no traceable
  # provenance: they were neither the values shipped in 1.0.2 nor anything I
  # could tie to the paper. Shipping plausible-looking numbers of unknown origin
  # is precisely the fault this release exists to correct, so they are returned
  # as NA. The decision rule uses the 5 percent level only.
  #
  # For the record, 1.0.2's own Case III table was wrong for every k >= 2 and
  # collapsed to a single row for k > 6. At k = 2 it used an I(1) bound of 4.38
  # where the paper gives 4.85, which made rejection easier for a second reason
  # independent of the levels-versus-differences defect.
  F5_low <- c(4.94, 3.79, 3.23, 2.86, 2.62, 2.45, 2.32, 2.22, 2.14, 2.06)
  F5_upp <- c(5.73, 4.85, 4.35, 4.01, 3.79, 3.61, 3.50, 3.39, 3.30, 3.24)

  # Table CII(iii): the I(0) bound is constant in k.
  # VERIFIED against the paper: the 5 percent I(1) upper bounds, k = 0..10, are
  #   -2.86 -3.22 -3.53 -3.78 -3.99 -4.19 -4.38 -4.57 -4.72 -4.88 -5.03
  # and the I(0) column is -3.43 / -2.86 / -2.57 at 1 / 5 / 10 percent, equal to
  # the Dickey-Fuller values, as the paper notes it must be at k = 0.
  # The 1 and 10 percent I(1) upper bounds have NOT been verified against the
  # printed table and are therefore NA rather than guessed. Only the 5 percent
  # level is used for the t decision in this release.
  t_low <- c(-3.43, -2.86, -2.57)
  t_upp5 <- c(-3.22, -3.53, -3.78, -3.99, -4.19,
              -4.38, -4.57, -4.72, -4.88, -5.03)   # k = 1..10
  t_upp <- rbind(rep(NA_real_, 10), t_upp5, rep(NA_real_, 10))

  list(
    F_lower      = c(NA_real_, F5_low[kk], NA_real_),
    F_upper      = c(NA_real_, F5_upp[kk], NA_real_),
    t_lower      = unname(t_low),
    t_upper      = unname(t_upp[, kk]),
    significance = c("1%", "5%", "10%"),
    verified     = c(FALSE, TRUE, FALSE)
  )
}


#' Bootstrap Bounds Test
#'
#' @description
#' Bootstrap cointegration test in the spirit of McNown, Sam and Goh (2018).
#'
#' \strong{Known limitation, retained from 1.0.2 and scheduled for 2.0.0.}
#' The pseudo-samples generated here are not built by accumulating the
#' differenced series under the null; each replication draws a fresh random walk
#' whose scale is taken from the data. The three statistics also share a single
#' restricted null rather than each having its own. Treat the p-values as
#' indicative only.
#'
#' @param y Dependent variable.
#' @param X Independent variables.
#' @param fourier Fourier terms.
#' @param p Lag for y.
#' @param q Lag for X.
#' @param tau Quantiles.
#' @param case Model case; only 3 is supported.
#' @param n_boot Number of bootstrap replications.
#' @param verbose Print progress messages.
#'
#' @return A list with bootstrap p-values.
#'
#' @export
bootstrap_bounds_test <- function(y, X, fourier, p, q, tau, case, n_boot = 1000,
                                  verbose = FALSE) {

  warning("bootstrap_bounds_test(): the null DGP is not yet the McNown, Sam and ",
          "Goh (2018) construction. See ?bootstrap_bounds_test. p-values are ",
          "indicative only.", call. = FALSE)

  n <- length(y)
  k <- ncol(X)

  orig <- estimate_qardl(y, X, fourier, p, q, tau[1], case)
  nms  <- names(orig$coefficients)
  lev  <- c("y_lag1", grep("^x[0-9]+_lag1$", nms, value = TRUE))

  wald_F <- function(res) {
    V <- res$vcov
    if (is.null(V)) return(NA_real_)
    pos <- match(lev, names(res$coefficients))
    R <- matrix(0, nrow = length(pos), ncol = length(res$coefficients))
    R[cbind(seq_along(pos), pos)] <- 1
    Rb  <- R %*% res$coefficients
    out <- tryCatch(as.numeric(t(Rb) %*% solve(R %*% V %*% t(R)) %*% Rb),
                    error = function(e) NA_real_)
    out / (k + 1)
  }

  orig_F <- wald_F(orig)
  orig_t <- unname(orig$coefficients["y_lag1"] / orig$std_errors["y_lag1"])

  boot_F <- rep(NA_real_, n_boot)
  boot_t <- rep(NA_real_, n_boot)

  for (b in seq_len(n_boot)) {
    y_boot <- cumsum(stats::rnorm(n, sd = stats::sd(diff(y))))
    res <- tryCatch(estimate_qardl(y_boot, X, fourier, p, q, tau[1], case),
                    error = function(e) NULL)
    if (!is.null(res)) {
      boot_t[b] <- unname(res$coefficients["y_lag1"] / res$std_errors["y_lag1"])
      boot_F[b] <- wald_F(res)
    }
    if (verbose && b %% 100 == 0)
      message(sprintf("   Bootstrap replication %d/%d", b, n_boot))
  }

  boot_F <- boot_F[!is.na(boot_F)]
  boot_t <- boot_t[!is.na(boot_t)]

  list(
    n_boot    = n_boot,
    orig_F    = orig_F,
    orig_t    = orig_t,
    boot_F    = boot_F,
    boot_t    = boot_t,
    p_value_F = mean(boot_F >= orig_F),
    p_value_t = mean(boot_t <= orig_t),
    decision  = if (mean(boot_F >= orig_F) < 0.05 || mean(boot_t <= orig_t) < 0.05)
                  "Reject the null: evidence of cointegration"
                else
                  "Fail to reject the null: no evidence of cointegration"
  )
}


#' Quantile Wald Test for Coefficient Constancy
#'
#' @description
#' Tests whether a coefficient is constant across quantiles.
#'
#' \strong{Known limitation, scheduled for 2.0.0.} The statistic below treats the
#' estimates at different quantiles as independent. They are not: quantile
#' regression estimates from the same sample have covariance proportional to
#' \eqn{\min(\tau_i,\tau_j) - \tau_i\tau_j}, which is strictly positive. The
#' correct denominator is therefore smaller than the one used here, so this test
#' is conservative and under-rejects constancy. Fixing it requires the joint
#' cross-quantile covariance, which this release does not compute.
#'
#' @param qardl_results List of QARDL results.
#' @param coef_name Name of the coefficient to test.
#'
#' @return A list with test results.
#'
#' @export
quantile_wald_test <- function(qardl_results, coef_name) {

  warning("quantile_wald_test(): assumes independence across quantiles and is ",
          "therefore conservative. See ?quantile_wald_test.", call. = FALSE)

  n_tau <- length(qardl_results)
  coefs <- sapply(qardl_results, function(x) x$coefficients[coef_name])
  ses   <- sapply(qardl_results, function(x) x$std_errors[coef_name])

  mean_coef <- mean(coefs)
  W  <- sum((coefs - mean_coef)^2 / ses^2)
  df <- n_tau - 1
  p_value <- 1 - stats::pchisq(W, df)

  list(
    coefficient = coef_name,
    values      = coefs,
    wald_stat   = W,
    df          = df,
    p_value     = p_value,
    decision    = if (p_value < 0.05)
                    "Reject constancy: coefficients vary across quantiles"
                  else
                    "Fail to reject: coefficients appear constant"
  )
}

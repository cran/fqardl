# =============================================================================
# Bounds Test for Cointegration
# Pesaran, Shin and Smith (2001), Journal of Applied Econometrics 16(3), 289-326.
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
#'   critical value bounds, \code{bounds_valid}, \code{bounds_note} and
#'   \code{decision}. Without Fourier terms the decision is three-way and
#'   reports the inconclusive region; with Fourier terms it states that no
#'   analytical verdict is available.
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
#' \strong{Fourier caveat.} The PSS critical values are tabulated for the
#' model without Fourier terms. With Fourier terms in the regression the
#' null distributions change and the PSS (2001) bounds do not apply with
#' Fourier terms (Monte Carlo evidence; in a Monte Carlo
#' with n = 100 and independent random walks the 5 percent decision of
#' fqardl 1.0.6 rejected in 10 percent of samples, the fnardl decision in
#' 47.5 percent). Since 1.1.0, whenever the design contains Fourier terms
#' (columns named \code{sin_*} or \code{cos_*}, which is always the case
#' for models estimated by \code{fqardl} and \code{fnardl}), no analytical
#' verdict is given: \code{decision} says so, \code{bounds_valid} is
#' \code{FALSE}, and the bounds are returned only for reference, labelled as
#' the bounds for the model without Fourier terms, not valid with Fourier
#' terms. Use the recursive bootstrap (\code{bootstrap = TRUE} in
#' \code{\link{fqardl}}, see \code{\link{bootstrap_bounds_test}}); when
#' \code{fqardl} or \code{fnardl} is called without the bootstrap, their
#' decision text adds this advice.
#'
#' \strong{Quantile caveat.} The PSS critical values were simulated for the
#' conditional mean. Applying them quantile by quantile has no distribution
#' theory behind it (Cho, Kim and Shin 2015 give no bounds test for the null
#' of no cointegration), so the quantile statistics are descriptive.
#'
#' @references
#' Pesaran, M.H., Shin, Y. and Smith, R.J. (2001). Bounds testing approaches to
#' the analysis of level relationships. \emph{Journal of Applied Econometrics},
#' 16(3), 289-326. \doi{10.1002/jae.616}
#'
#' Cho, J.S., Kim, T. and Shin, Y. (2015). Quantile cointegration in the
#' autoregressive distributed-lag modeling framework. \emph{Journal of
#' Econometrics}, 188(1), 281-300. \doi{10.1016/j.jeconom.2015.05.003}
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
    bounds_one_quantile(res, k = k, m_restr = m_restr, cv = cv,
                        fourier = has_fourier_terms(names(res$coefficients))))

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


# TRUE when a design contains Fourier terms (columns sin_* / cos_*)
has_fourier_terms <- function(nms) any(grepl("^(sin|cos)_", nms))

# Label of PSS bounds shown for a model with Fourier terms
fourier_bounds_note <- paste(
  "Bounds for the model without Fourier terms; not valid with Fourier terms.")

# Decision text when the design contains Fourier terms
fourier_no_verdict <- paste(
  "No analytical verdict: the PSS bounds are for the model without Fourier",
  "terms and are not valid with Fourier terms")

# Advice appended by fqardl() and fnardl() when the bootstrap was not run
fourier_boot_advice <- "; use the bootstrap (bootstrap = TRUE)"



#' Bounds test at a single quantile
#' @keywords internal
bounds_one_quantile <- function(res, k, m_restr, cv, fourier = FALSE) {

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

  # NEW in 1.1.0: the PSS (2001) bounds do not apply with Fourier terms.
  if (fourier) {
    dec_F <- dec_t <- "not available (Fourier terms)"
    decision <- fourier_no_verdict
  }

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
    bounds_valid = !fourier,
    bounds_note = if (fourier) fourier_bounds_note else
                    "PSS (2001) Table CI(iii) and CII(iii), 5 percent",
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
#' Recursive bootstrap bounds test for the Fourier ARDL model estimated by
#' \code{\link{fqardl}}. Rewritten in 1.1.0.
#'
#' @details
#' The pseudo-samples are generated by the shared recursive bootstrap engine
#' (file \code{ardl_boot_engine.R}, engine version 1.1.0, copied verbatim
#' from the package 'ardlverse'). Under each null hypothesis the restricted
#' conditional error correction model (ECM) for \eqn{\Delta y} is estimated
#' by OLS together with a marginal model for \eqn{\Delta x}; \eqn{y^*} and
#' \eqn{x^*} are then generated recursively from resampled, recentred
#' residuals, with the Fourier terms held fixed as deterministic regressors.
#' \describe{
#'   \item{\code{scheme = "bvz"} (default)}{Bertelli, Vacca and Zoia (2022):
#'     a separate restricted model for each of the three statistics, a
#'     marginal VECM for \eqn{\Delta x}, residuals recentred after each
#'     draw, initial values from a random block of the data (the levels of
#'     the regressors and partial sums continue from the block).}
#'   \item{\code{scheme = "mcnown"}}{The Fov null for all statistics
#'     (McNown, Sam and Goh 2018, Steps 1-8) applied to the conditional ECM
#'     (package choice; MSG write the y equation in unconditional form): the
#'     null of the overall F test generates the samples for all statistics,
#'     the \eqn{\Delta x} equation is unrestricted (it includes
#'     \eqn{y_{t-1}}), residuals mean-centred once, the first observations
#'     as initial values.}
#' }
#' The three statistics are computed from fqardl's own design
#' (\code{build_ardl_design}) estimated by OLS: the overall F test on the
#' lagged levels of y and x (Fov), the t ratio on the lagged level of y (t)
#' and the F test on the lagged levels of x (Find). The same function is
#' applied to the data and to every bootstrap sample, so the bootstrap
#' distribution refers to the reported statistics. When \code{reselect} is
#' given, the Fourier frequency (unless fixed) and the lag orders are chosen
#' again on every bootstrap sample with \code{select_fourier_frequency} and
#' \code{select_optimal_lags}. Re-selection is a package choice: Bertelli,
#' Vacca and Zoia (2022, step 5(c)) re-estimate the unrestricted model on
#' each bootstrap sample but do not discuss re-selection. In the package's
#' Monte Carlo, holding a data-selected frequency fixed made the test
#' oversized.
#'
#' The decision combines the three tests with the AND rule: cointegration
#' only if Fov, t and Find all reject at level \code{level}; Fov rejecting
#' without t, or Fov and t without Find, are reported as the degenerate
#' cases of the first and second type. Critical values are order statistics
#' of the bootstrap distributions.
#'
#' The statistics are OLS (conditional mean) statistics. Cho, Kim and Shin
#' (2015) give no bounds test and no bootstrap for the null of no
#' cointegration in the quantile ARDL model, and no quantile-specific
#' bootstrap statistic is computed. Versions up to 1.0.6 regressed a new
#' Gaussian random walk on fixed regressors at the first quantile and
#' declared cointegration when F or t rejected; that procedure was oversized
#' and has been removed.
#'
#' @param y Dependent variable (no missing values).
#' @param X Matrix of regressors in levels.
#' @param fourier Matrix of Fourier terms (from \code{generate_fourier_terms})
#'   used for the observed data and as fixed deterministic terms of the
#'   bootstrap data generating process.
#' @param p Lag order for y, as in \code{build_ardl_design}.
#' @param q Lag order for X, as in \code{build_ardl_design}.
#' @param tau Not used since 1.1.0 (the statistics are OLS statistics); kept
#'   for compatibility.
#' @param case Model case; only 3 is supported.
#' @param n_boot Number of bootstrap replications.
#' @param verbose Print a progress message.
#' @param reselect \code{NULL} (hold the frequency and lags fixed) or a list
#'   with elements \code{max_k}, \code{max_p}, \code{max_q},
#'   \code{criterion} and optionally \code{k} (a frequency fixed by the
#'   user in advance) to re-run the selection on every bootstrap sample.
#' @param scheme \code{"bvz"} or \code{"mcnown"}; see Details.
#' @param level Significance level of the decision.
#' @param seed Optional seed for the bootstrap; the caller's random number
#'   state is restored afterwards.
#'
#' @return A list with the observed statistics (\code{orig_F},
#'   \code{orig_t}, \code{orig_Find}), the bootstrap distributions
#'   (\code{boot_F}, \code{boot_t}, \code{boot_Find}; failed replications are
#'   \code{NA}), the p-values (\code{p_value_F}, \code{p_value_t},
#'   \code{p_value_Find}), the critical values \code{cv} at 10, 5, 2.5 and 1
#'   percent, the decision code \code{decision} and its text \code{label},
#'   and the engine output \code{engine}.
#'
#' @references
#' Bertelli, S., Vacca, G. and Zoia, M. (2022). Bootstrap cointegration tests
#' in ARDL models. \emph{Economic Modelling}, 116, 105987.
#' \doi{10.1016/j.econmod.2022.105987}
#'
#' McNown, R., Sam, C.Y. and Goh, S.K. (2018). Bootstrapping the
#' autoregressive distributed lag test for cointegration. \emph{Applied
#' Economics}, 50(13), 1509-1521. \doi{10.1080/00036846.2017.1366643}
#'
#' Cho, J.S., Kim, T. and Shin, Y. (2015). Quantile cointegration in the
#' autoregressive distributed-lag modeling framework. \emph{Journal of
#' Econometrics}, 188(1), 281-300. \doi{10.1016/j.jeconom.2015.05.003}
#'
#' @export
bootstrap_bounds_test <- function(y, X, fourier, p, q, tau = NULL, case = 3,
                                  n_boot = 1000, verbose = FALSE,
                                  reselect = NULL, scheme = c("bvz", "mcnown"),
                                  level = 0.05, seed = NULL) {

  if (!identical(as.numeric(case), 3))
    stop("Only case = 3 (unrestricted intercept, no trend) is supported.",
         call. = FALSE)
  scheme <- match.arg(scheme)
  y <- as.numeric(y)
  X <- as.matrix(X)
  fourier <- as.matrix(fourier)
  if (nrow(fourier) != length(y))
    stop("'fourier' must have one row per observation", call. = FALSE)
  if (verbose)
    message(sprintf("   Recursive bootstrap (%s), %d replications", scheme, n_boot))

  sc <- .fq_scheme(scheme)
  spec <- list(fourier = fourier, p = p, q = q)
  eng <- .ardl_boot_engine(
    y, X, p = p - 1, q = rep(q - 1, ncol(X)), case = 3, t0 = max(p, q) + 1,
    det = fourier, nulls = sc$nulls, xmodel = sc$xmodel, B = n_boot,
    init = sc$init, recentre = sc$recentre,
    stat_fun = function(yy, xx, sp) .fq_stats_qardl(yy, xx, sp),
    select = .fq_make_select(reselect), spec = spec, use_find = TRUE,
    level = level, seed = seed)

  .fq_boot_result(eng, scheme, n_boot, reselect,
                  model = "Fourier ARDL (conditional mean), fqardl design")
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

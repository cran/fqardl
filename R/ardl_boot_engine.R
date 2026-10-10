# ---------------------------------------------------------------------------
# ardl_boot_engine.R: shared recursive bootstrap engine for ARDL bounds tests
# Engine version 1.1.0 (2026-10-03). Source of truth: package 'ardlverse'.
#
# This file is self-contained (base R and 'stats' only) and is meant to be
# copied verbatim into other packages. Do not edit a copy; edit the source in
# ardlverse and copy it again (packages check the copy with an md5 test).
#
# Schemes
#   nulls = "separate": Bertelli, Vacca and Zoia (2022, Section 3, eqs. 16-25):
#     one restricted model per statistic (Fov, t, Find).
#   nulls = "joint": the Fov null generates the samples for all statistics
#     (McNown, Sam and Goh, 2018, Steps 1-8). MSG write the y equation in
#     unconditional form (no contemporaneous dx, eq. 12); use j0 = 1 for that
#     form. With j0 = 0 (default) the joint null is applied to the
#     conditional ECM, a package choice.
#   xmodel = "vecm": marginal model for dx with x_{t-1} and lagged dy, dx
#     (BVZ eqs. 19-20); "var": unrestricted dx equation that also contains
#     y_{t-1} (MSG eq. 12); "rw": dx on the deterministic terms only.
#   recentre = "draw": residuals recentred after each resampling (BVZ eqs.
#     21-22); "once": recentred once before resampling (MSG eq. 13, which
#     subtracts the sum of the residuals divided by n - q - 1, the number of
#     residuals of their ARDL(1,1) model, i.e. the residual mean; no
#     rescaling factor is applied).
#   xlags: lags of dy and dx in the x equation; default p, the number of
#     lagged dy in the y equation (BVZ eq. 19 has p - 1 lags of dz for a VAR
#     of order p, whose conditional ECM has p - 1 lagged dy); may be 0.
#   j0: first lag of dw in the y equation (0 = conditional ECM with the
#     contemporaneous dw_t; 1 = unconditional form of MSG eq. 12).
#   init = "block": initial values are a random block of the data (BVZ step
#     5b); "observed": the first observations of the data. The initial levels
#     of w are the data's w in the block, so the generated w* continue the
#     levels of the estimated model.
# y* and x* are generated recursively; level regressors w = decompose(x) (for
# example partial sums) are rebuilt from x* with the caller's decomposition;
# deterministic columns 'det' (for example Fourier terms) are held fixed.
# ---------------------------------------------------------------------------

.abe_engine_version <- "1.1.0"

.abe_lag <- function(v, j) {
  n <- length(v)
  if (j == 0) return(v)
  c(rep(NA_real_, j), v[seq_len(n - j)])
}

# Canonical conditional ECM on rows t0..n:
# dy_t = const + trend + det_t + rho y_{t-1} + theta' w_{t-1}
#        + sum_{i=1}^{p} phi_i dy_{t-i} + sum_j sum_{l=0}^{q_j} pi_jl dw_{j,t-l} + u_t
.abe_design <- function(y, w, p, q, t0, case, det = NULL, j0 = 0) {
  n <- length(y)
  w <- as.matrix(w)
  kw <- ncol(w)
  rows <- t0:n
  dy <- c(NA_real_, diff(y))
  dw <- rbind(rep(NA_real_, kw), diff(w))
  cols <- list()
  nm <- character(0)
  grp <- character(0)
  lagv <- integer(0)
  colj <- integer(0)
  add <- function(v, name, g, l = NA_integer_, j = NA_integer_) {
    cols[[length(cols) + 1L]] <<- v[rows]
    nm <<- c(nm, name)
    grp <<- c(grp, g)
    lagv <<- c(lagv, l)
    colj <<- c(colj, j)
  }
  if (case >= 2) add(rep(1, n), "const", "const")
  if (case >= 4) add(as.numeric(seq_len(n)), "trend", "trend")
  if (!is.null(det) && ncol(det) > 0)
    for (j in seq_len(ncol(det))) add(det[, j], paste0("det", j), "det", NA, j)
  add(c(NA_real_, y[-n]), "ly", "ly")
  for (j in seq_len(kw)) add(c(NA_real_, w[-n, j]), paste0("lw", j), "lw", NA, j)
  for (i in seq_len(p)) add(.abe_lag(dy, i), paste0("dy", i), "dy", i)
  for (j in seq_len(kw)) {
    if (q[j] >= j0) for (i in j0:q[j])
      add(.abe_lag(dw[, j], i), paste0("dw", j, "_", i), "dw", i, j)
  }
  X <- do.call(cbind, cols)
  colnames(X) <- nm
  list(Y = dy[rows], X = X, grp = grp, lag = lagv, col = colj, rows = rows)
}

.abe_ols <- function(X, Y) {
  qx <- qr(X)
  if (qx$rank < ncol(X)) stop("rank-deficient design", call. = FALSE)
  b <- qr.coef(qx, Y)
  e <- qr.resid(qx, Y)
  list(b = b, e = e, qr = qx, df = NROW(X) - ncol(X))
}

# F (Wald / m with the OLS covariance) on columns idx, and t on column it
.abe_stats_fit <- function(f, iov, it, ind) {
  piv <- order(f$qr$pivot)
  XtXi <- chol2inv(qr.R(f$qr))[piv, piv, drop = FALSE]
  s2 <- sum(f$e^2) / f$df
  Fw <- function(i) {
    if (length(i) == 0) return(NA_real_)
    bi <- f$b[i]
    drop(crossprod(bi, solve(XtXi[i, i, drop = FALSE], bi))) / (length(i) * s2)
  }
  c(Fov = Fw(iov), t = unname(f$b[it] / sqrt(s2 * XtXi[it, it])), Find = Fw(ind))
}

.abe_null_cols <- function(grp, case, which) {
  switch(which,
         Fov = grp %in% c("ly", "lw") | (case == 2 & grp == "const") |
           (case == 4 & grp == "trend"),
         t = grp == "ly",
         Find = grp == "lw")
}

.abe_stats <- function(D, case) {
  f <- .abe_ols(D$X, D$Y)
  .abe_stats_fit(f, which(.abe_null_cols(D$grp, case, "Fov")),
                 which(D$grp == "ly"), which(D$grp == "lw"))
}

# Marginal model for dx on rows t0..n
.abe_xdesign <- function(y, x, xmodel, xlags, t0, case, det = NULL) {
  n <- length(y)
  x <- as.matrix(x)
  kx <- ncol(x)
  rows <- t0:n
  dy <- c(NA_real_, diff(y))
  dx <- rbind(rep(NA_real_, kx), diff(x))
  Z <- NULL
  if (case >= 2) Z <- cbind(Z, const = rep(1, length(rows)))
  if (case >= 4) Z <- cbind(Z, trend = rows)
  if (!is.null(det) && ncol(det) > 0) Z <- cbind(Z, det[rows, , drop = FALSE])
  if (xmodel == "var") Z <- cbind(Z, ly = y[rows - 1])
  if (xmodel %in% c("vecm", "var")) Z <- cbind(Z, x[rows - 1, , drop = FALSE])
  for (l in seq_len(xlags)) {
    Z <- cbind(Z, dy[rows - l], dx[rows - l, , drop = FALSE])
  }
  if (is.null(Z)) Z <- matrix(0, length(rows), 0)
  list(Y = dx[rows, , drop = FALSE], Z = Z)
}

# Recursive generation of (y*, x*). 'cy' holds the y-equation coefficients
# (canonical layout, zeros on restricted terms), 'cx' the x-equation ones.
.abe_generate <- function(n, t0, y0, x0, cy, cx, u, v, m, w0 = NULL) {
  kx <- m$kx
  kw <- m$kw
  ys <- numeric(n)
  xs <- matrix(0, n, kx)
  ws <- matrix(0, n, kw)
  dys <- rep(NA_real_, n)
  dxs <- matrix(NA_real_, n, kx)
  dws <- matrix(NA_real_, n, kw)
  h <- seq_len(t0 - 1)
  ys[h] <- y0
  xs[h, ] <- x0
  ws[h, ] <- if (!is.null(w0)) w0 else if (is.null(m$decompose)) x0 else
    as.matrix(m$decompose(xs[h, , drop = FALSE]))
  if (t0 > 2) {
    dys[2:(t0 - 1)] <- diff(ys[h])
    dxs[2:(t0 - 1), ] <- diff(xs[h, , drop = FALSE])
    dws[2:(t0 - 1), ] <- diff(ws[h, , drop = FALSE])
  }
  # deterministic contributions
  dety <- rep(0, n)
  if (m$case >= 2) dety <- dety + cy[m$iy$const]
  if (m$case >= 4) dety <- dety + cy[m$iy$trend] * seq_len(n)
  if (length(m$iy$det)) dety <- dety + drop(m$det %*% cy[m$iy$det])
  detx <- matrix(0, n, kx)
  if (length(m$ix$det0)) {
    D0 <- NULL
    if (m$case >= 2) D0 <- cbind(D0, rep(1, n))
    if (m$case >= 4) D0 <- cbind(D0, seq_len(n))
    if (!is.null(m$det) && ncol(m$det) > 0) D0 <- cbind(D0, m$det)
    detx <- D0 %*% cx[m$ix$det0, , drop = FALSE]
  }
  # stacked coefficients and linear-index offsets of the lagged terms
  by <- c(cy[m$iy$ly], cy[m$iy$lw], cy[m$iy$dy], cy[m$iy$dw])
  p <- length(m$iy$dy)
  offw <- (m$col_dw - 1L) * n - m$lag_dw
  ckw <- (seq_len(kw) - 1L) * n
  ckx <- (seq_len(kx) - 1L) * n
  xl <- m$xlags
  lagx <- seq_len(xl)
  offdx <- unlist(lapply(lagx, function(l) ckx - l))
  # x equation coefficients stacked as (ly, lx, dy lags, dx lags by lag)
  Cx <- rbind(cx[m$ix$ly, , drop = FALSE], cx[m$ix$lx, , drop = FALSE],
              cx[m$ix$dy, , drop = FALSE],
              do.call(rbind, lapply(lagx, function(l) cx[m$ix$dx[[l]], , drop = FALSE])))
  if (is.null(Cx) || nrow(Cx) == 0) Cx <- NULL
  has_ly <- length(m$ix$ly) > 0
  has_lx <- length(m$ix$lx) > 0
  lagp <- seq_len(p)
  dec <- m$decompose
  step <- m$step
  for (t in t0:n) {
    i <- t - t0 + 1L
    if (is.null(Cx)) {
      mx <- detx[t, ]
    } else {
      z <- c(if (has_ly) ys[t - 1], if (has_lx) xs[t - 1 + ckx],
             dys[t - lagx], dxs[t + offdx])
      mx <- detx[t, ] + drop(z %*% Cx)
    }
    dxs[t, ] <- mx + v[i, ]
    xs[t, ] <- xs[t - 1, ] + dxs[t, ]
    if (is.null(dec)) {
      dws[t, ] <- dxs[t, ]
    } else if (!is.null(step)) {
      dws[t, ] <- step(dxs[t, ])
    } else {
      W2 <- as.matrix(dec(xs[(t - 1):t, , drop = FALSE]))
      dws[t, ] <- W2[2, ] - W2[1, ]
    }
    ws[t, ] <- ws[t - 1, ] + dws[t, ]
    dys[t] <- dety[t] + sum(by * c(ys[t - 1], ws[t - 1 + ckw], dys[t - lagp],
                                   dws[t + offw])) + u[i]
    ys[t] <- ys[t - 1] + dys[t]
  }
  list(y = ys, x = xs, w = ws)
}

# Order-statistic critical values (MSG eqs. 15-16, BVZ eqs. 24-25)
.abe_cv <- function(d, alpha, lower = FALSE) {
  d <- sort(d[is.finite(d)])
  B <- length(d)
  if (B == 0) return(NA_real_)
  m <- floor(alpha * B + 1e-9)
  if (lower) d[min(B, m + 1)] else d[max(1, B - m)]
}

# Run 'expr' with a local seed and restore the caller's RNG state
.abe_with_seed <- function(seed, expr) {
  if (is.null(seed)) return(expr)
  env <- globalenv()
  had <- exists(".Random.seed", envir = env, inherits = FALSE)
  if (had) old <- get(".Random.seed", envir = env, inherits = FALSE)
  on.exit({
    if (had) assign(".Random.seed", old, envir = env) else
      if (exists(".Random.seed", envir = env, inherits = FALSE))
        rm(".Random.seed", envir = env)
  })
  set.seed(seed)
  expr
}

#' Shared recursive bootstrap engine for ARDL bounds tests (internal)
#'
#' @param y Numeric vector, dependent variable in levels (no missing values).
#' @param x Numeric matrix of the underlying regressors in levels.
#' @param p Number of lagged differences of y in the ECM.
#' @param q Integer vector, one entry per level regressor in
#'   \code{decompose(x)}: highest lag of its difference (0 = contemporaneous
#'   only, -1 = none).
#' @param case PSS case 1 to 5.
#' @param t0 First row of the estimation sample (default: the smallest
#'   admissible value).
#' @param det Optional n-row matrix of deterministic columns held fixed.
#' @param decompose Optional function mapping x (n x kx) to the level
#'   regressors w (n x kw); it must be causal with increments that depend
#'   only on the current change of x (identity and partial sums are). It may
#'   also be a list with elements \code{full} (that function) and
#'   \code{step}, a fast function mapping the vector of current changes of
#'   x to the changes of w; the reproduction check verifies both.
#' @param nulls "separate" (BVZ) or "joint" (MSG).
#' @param xmodel "vecm", "var" or "rw".
#' @param xlags Lags of dy and dx in the x equation (default \code{p}; 0
#'   allowed).
#' @param B Number of bootstrap replications.
#' @param init "block" or "observed".
#' @param recentre "draw" or "once".
#' @param stat_fun Function (y, x, spec) returning c(Fov, t, Find) computed
#'   with the caller's own design; NULL uses the canonical ECM.
#' @param select Optional function (y, x) returning a spec passed to
#'   \code{stat_fun}; re-run on every bootstrap sample.
#' @param spec Spec for the observed data (passed to \code{stat_fun}).
#' @param use_find Logical, include the Find test.
#' @param level Significance level of the combined decision.
#' @param seed Optional seed; the caller's RNG state is restored.
#' @param tol Tolerance of the reproduction and design checks.
#' @param j0 First lag of the differenced level regressors in the y equation:
#'   0 (conditional ECM, default) or 1 (unconditional form, MSG eq. 12).
#' @return A list (class "ardl_boot_engine").
#' @keywords internal
#' @noRd
.ardl_boot_engine <- function(y, x, p, q, case = 3, t0 = NULL, det = NULL,
                              decompose = NULL,
                              nulls = c("separate", "joint"),
                              xmodel = c("vecm", "var", "rw"), xlags = NULL,
                              B = 999, init = c("block", "observed"),
                              recentre = c("draw", "once"),
                              stat_fun = NULL, select = NULL, spec = NULL,
                              use_find = TRUE, level = 0.05, seed = NULL,
                              tol = 1e-6, j0 = 0) {
  nulls <- match.arg(nulls)
  xmodel <- match.arg(xmodel)
  init <- match.arg(init)
  recentre <- match.arg(recentre)
  y <- as.numeric(y)
  x <- as.matrix(x)
  storage.mode(x) <- "double"
  n <- length(y)
  if (nrow(x) != n) stop("'y' and 'x' must have the same length", call. = FALSE)
  if (anyNA(y) || anyNA(x)) stop("missing values in 'y' or 'x'", call. = FALSE)
  if (!is.null(det)) {
    det <- as.matrix(det)
    if (nrow(det) != n) stop("'det' must have one row per observation", call. = FALSE)
  }
  kx <- ncol(x)
  step <- NULL
  if (is.list(decompose)) {
    step <- decompose$step
    decompose <- decompose$full
  }
  w <- if (is.null(decompose)) x else as.matrix(decompose(x))
  kw <- ncol(w)
  if (length(q) == 1) q <- rep(q, kw)
  if (length(q) != kw) stop("'q' must have one entry per level regressor", call. = FALSE)
  if (!j0 %in% c(0, 1)) stop("'j0' must be 0 or 1", call. = FALSE)
  tmin <- max(p + 2, max(q) + 2, 2)
  if (is.null(t0)) t0 <- tmin
  if (t0 < tmin) stop("'t0' is too small for the lag orders", call. = FALSE)
  if (is.null(xlags)) xlags <- p
  xlags <- max(0, min(xlags, t0 - 2))
  if (xmodel == "rw") xlags <- 0

  # canonical design on the data
  D <- .abe_design(y, w, p, q, t0, case, det, j0)
  N <- length(D$rows)
  if (N - ncol(D$X) < 5) stop("too few observations for the bootstrap", call. = FALSE)
  stats_engine <- .abe_stats(D, case)
  if (is.null(stat_fun)) stat_fun <- function(y, x, spec) {
    ww <- if (is.null(decompose)) as.matrix(x) else as.matrix(decompose(as.matrix(x)))
    .abe_stats(.abe_design(y, ww, p, q, t0, case, det, j0), case)
  }
  obs <- stat_fun(y, x, spec)
  obs <- obs[c("Fov", "t", "Find")]
  if (!use_find) obs["Find"] <- NA_real_
  chk <- c(stats_engine[1:2], if (use_find) stats_engine[3])
  design_check <- max(abs(chk - obs[names(chk)]) / pmax(1, abs(chk)))
  if (!is.finite(design_check) || design_check > tol)
    warning("the statistics of 'stat_fun' differ from the canonical ECM of the ",
            "bootstrap engine (max relative difference ", signif(design_check, 3),
            "); the bootstrap model may not match the estimated model",
            call. = FALSE)

  # x equation
  XD <- .abe_xdesign(y, x, xmodel, xlags, t0, case, det)
  if (ncol(XD$Z) > 0) {
    fx <- qr(XD$Z)
    if (fx$rank < ncol(XD$Z)) stop("rank-deficient x equation", call. = FALSE)
    cx <- qr.coef(fx, XD$Y)
    ex <- qr.resid(fx, XD$Y)
  } else {
    cx <- matrix(0, 0, kx)
    ex <- XD$Y
  }
  cx <- as.matrix(cx)
  ex <- as.matrix(ex)
  # index map of the x equation
  nd0 <- (case >= 2) + (case >= 4) + (if (is.null(det)) 0 else ncol(det))
  pos <- nd0
  ix <- list(det0 = seq_len(nd0))
  ix$ly <- if (xmodel == "var") { pos <- pos + 1; pos } else integer(0)
  ix$lx <- if (xmodel %in% c("vecm", "var")) { r <- pos + seq_len(kx); pos <- pos + kx; r } else integer(0)
  ix$dy <- integer(0)
  ix$dx <- list()
  for (l in seq_len(xlags)) {
    ix$dy[l] <- pos + 1
    ix$dx[[l]] <- pos + 1 + seq_len(kx)
    pos <- pos + 1 + kx
  }
  iy <- list(const = which(D$grp == "const"), trend = which(D$grp == "trend"),
             det = which(D$grp == "det"), ly = which(D$grp == "ly"),
             lw = which(D$grp == "lw"), dy = which(D$grp == "dy"),
             dw = which(D$grp == "dw"))
  mdl <- list(kx = kx, kw = kw, case = case, det = det, iy = iy, ix = ix,
              lag_dw = D$lag[iy$dw], col_dw = D$col[iy$dw], xlags = xlags,
              decompose = decompose, step = step)

  # restricted y equations
  which_nulls <- if (nulls == "joint") "Fov" else c("Fov", "t", if (use_find) "Find")
  restr <- list()
  for (s in which_nulls) {
    keep <- !.abe_null_cols(D$grp, case, s)
    fr <- .abe_ols(D$X[, keep, drop = FALSE], D$Y)
    cy <- stats::setNames(numeric(ncol(D$X)), colnames(D$X))
    cy[keep] <- fr$b
    restr[[s]] <- list(cy = cy, e = fr$e)
  }

  # reproduction check: original residuals and initial values give the data
  dgpcheck <- sapply(which_nulls, function(s) {
    g <- .abe_generate(n, t0, y[seq_len(t0 - 1)], x[seq_len(t0 - 1), , drop = FALSE],
                       restr[[s]]$cy, cx, restr[[s]]$e, ex, mdl)
    max(abs(g$y - y), abs(g$x - x), abs(g$w - w))
  })
  if (any(!is.finite(dgpcheck)) || any(dgpcheck > tol * max(1, abs(y), abs(x))))
    warning("the bootstrap recursion does not reproduce the data (max abs ",
            "difference ", signif(max(dgpcheck), 3), ")", call. = FALSE)

  stat_names <- c("Fov", "t", "Find")
  boot <- matrix(NA_real_, B, 3, dimnames = list(NULL, stat_names))
  run <- function() {
    for (s in which_nulls) {
      e <- restr[[s]]$e
      cy <- restr[[s]]$cy
      ec <- e - mean(e)
      exc <- sweep(ex, 2, colMeans(ex))
      targets <- if (nulls == "joint") stat_names[c(TRUE, TRUE, use_find)] else s
      for (b in seq_len(B)) {
        j <- sample.int(N, N, replace = TRUE)
        if (recentre == "draw") {
          u <- e[j] - mean(e[j])
          v <- sweep(ex[j, , drop = FALSE], 2, colMeans(ex[j, , drop = FALSE]))
        } else {
          u <- ec[j]
          v <- exc[j, , drop = FALSE]
        }
        r <- if (init == "block") sample.int(n - t0 + 2, 1) else 1L
        hh <- r + seq_len(t0 - 1) - 1
        val <- tryCatch({
          g <- .abe_generate(n, t0, y[hh], x[hh, , drop = FALSE], cy, cx, u, v, mdl,
                             w0 = w[hh, , drop = FALSE])
          sp <- if (is.null(select)) spec else select(g$y, g$x)
          st <- stat_fun(g$y, g$x, sp)
          st[c("Fov", "t", "Find")]
        }, error = function(err) rep(NA_real_, 3))
        val[!is.finite(val)] <- NA_real_
        boot[b, targets] <<- val[targets]
      }
    }
    invisible(NULL)
  }
  .abe_with_seed(seed, run())

  used <- stat_names[c(TRUE, TRUE, use_find)]
  n_valid <- colSums(!is.na(boot))[used]
  n_fail <- B - n_valid
  if (any(n_fail > 0))
    warning(sprintf("%s of %d bootstrap replications failed and were set to NA",
                    paste(sprintf("%s: %d", used, n_fail), collapse = ", "), B),
            call. = FALSE)
  pval <- c(Fov = mean(boot[, "Fov"] >= obs["Fov"], na.rm = TRUE),
            t = mean(boot[, "t"] <= obs["t"], na.rm = TRUE),
            Find = if (use_find) mean(boot[, "Find"] >= obs["Find"], na.rm = TRUE) else NA_real_)
  levs <- c(0.10, 0.05, 0.025, 0.01)
  cv <- matrix(NA_real_, length(levs), 3,
               dimnames = list(paste0(levs * 100, "%"), stat_names))
  for (a in seq_along(levs)) {
    cv[a, "Fov"] <- .abe_cv(boot[, "Fov"], levs[a])
    cv[a, "t"] <- .abe_cv(boot[, "t"], levs[a], lower = TRUE)
    if (use_find) cv[a, "Find"] <- .abe_cv(boot[, "Find"], levs[a])
  }
  cvl <- c(Fov = .abe_cv(boot[, "Fov"], level),
           t = .abe_cv(boot[, "t"], level, lower = TRUE),
           Find = if (use_find) .abe_cv(boot[, "Find"], level) else NA_real_)
  rej <- c(Fov = unname(obs["Fov"] > cvl["Fov"]), t = unname(obs["t"] < cvl["t"]),
           Find = if (use_find) unname(obs["Find"] > cvl["Find"]) else NA)
  dec <- .abe_decision(rej, use_find)
  structure(list(statistic = obs, boot = boot, p_value = pval, cv = cv,
                 cv_level = cvl, reject = rej, decision = dec$decision,
                 label = dec$label, level = level, n_valid = n_valid,
                 n_fail = n_fail, dgpcheck = dgpcheck,
                 design_check = design_check,
                 settings = list(nulls = nulls, xmodel = xmodel, xlags = xlags, j0 = j0,
                                 init = init, recentre = recentre, B = B,
                                 case = case, p = p, q = q, t0 = t0,
                                 reselect = !is.null(select)),
                 engine_version = .abe_engine_version),
            class = "ardl_boot_engine")
}

# Combined decision: Fov AND t (AND Find); never an OR rule
.abe_decision <- function(rej, use_find = TRUE) {
  r <- rej[c("Fov", "t", if (use_find) "Find")]
  if (anyNA(r))
    return(list(decision = "UNDETERMINED",
                label = "Undetermined: a bootstrap statistic or critical value is missing"))
  tests <- if (use_find) "Fov, t and Find" else "Fov and t"
  if (!r["Fov"]) {
    list(decision = "NO_COINTEGRATION",
         label = "No cointegration: Fov does not reject")
  } else if (!r["t"]) {
    list(decision = "DEGENERATE_1",
         label = paste("No cointegration, degenerate case of the first type:",
                       "Fov rejects but t does not (lagged y not significant)"))
  } else if (use_find && !r["Find"]) {
    list(decision = "DEGENERATE_2",
         label = paste("No cointegration, degenerate case of the second type:",
                       "Fov and t reject but Find does not (lagged regressors",
                       "not significant)"))
  } else {
    list(decision = "COINTEGRATION",
         label = paste0("Cointegration: ", tests, " all reject"))
  }
}

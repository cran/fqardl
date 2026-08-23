# fqardl 1.0.3

**This is a correctness release. Results change. Anyone who has run version 1.0.2 should
re-run their analysis before relying on it.**

The defects below came to light after a question from Prof. Chi-Lu (Edward) Peng, National
Kaohsiung University of Science and Technology, who noticed that the bounds-test t statistic
was positive at every quantile when it should be negative. Our thanks to him for reporting it.

## Corrected

### The estimated equation had an error-correction right-hand side and a levels left-hand side

`build_ardl_design()`, `estimate_nardl()` and `estimate_mtnardl()` assembled a Pesaran, Shin
and Smith conditional error-correction design on the right-hand side, containing the lagged
levels `y_{t-1}` and `x_{t-1}` together with the differences `Delta y_{t-j}` and
`Delta x_{t-j}`, while setting the left-hand side to `y_t` in levels rather than `Delta y_t`.
The equation actually estimated was therefore the conditional levels ARDL.

Hence the coefficient reported as the error-correction term was `phi = 1 + rho` rather than
`rho`. Since `rho` lies in `(-1, 0)` under cointegration, the reported coefficient was
positive and close to unity, with a large positive t ratio at every quantile.

The defect propagated:

- **Long-run multipliers carried the wrong sign.** `compute_multipliers()` applies
  `-theta / phi`, which is valid only when the left-hand side is `Delta y`. In a Monte Carlo
  with a true long-run coefficient of `+0.80`, version 1.0.2 returned a median of `-0.471`
  and the wrong sign in 100 percent of replications. Version 1.0.3 returns `+0.791` with the
  correct sign in 100 percent of replications.
- **The bounds test manufactured cointegration** out of the unit root in `y`. On the package's
  own `macro_data`, version 1.0.2 reported `F = 65.57` against a 5 percent upper bound of 4.38
  and concluded "Evidence of cointegration". Version 1.0.3 reports `F = 2.34`, below the 5
  percent lower bound of 3.79, and concludes no cointegration.
- **`fnardl()` and `mtnardl()` were affected identically.** On `oil_gdp_data`, `fnardl()`
  previously reported an error-correction coefficient of `+0.894` with `t = +27.08` and
  `F = 246.07`; it now reports `-0.106` with `t = -3.21` and `F = 4.47`.

### The "F statistic" was not a Wald statistic

`perform_bounds_test()`, `perform_nardl_bounds_test()` and `perform_mtnardl_bounds()`
computed `mean(t^2)` over the terms matching `_lag1$`, that is, the arithmetic mean of the
squared marginal t ratios. That ignores the covariances among the level coefficients, has no
distribution theory, and was being compared against tables that Pesaran, Shin and Smith
simulated for a genuine Wald statistic.

All three now compute the Wald form of PSS equation (21),
`W = (Rb)' (R V R')^{-1} (Rb)`, with `F = W / (k+1)` for Case III.

### The inconclusive region was not reported

The bounds procedure has three outcomes, not two. When the statistic falls between the I(0)
and I(1) bounds, inference is inconclusive without knowing the cointegration rank of the
regressors. Version 1.0.2 collapsed this into a binary verdict. All bounds routines now
return a three-way decision.

### The bounds test is now computed at every quantile

`perform_bounds_test()` located `tau = 0.5`, computed one statistic there, and reported it as
though it were general. It now returns a per-quantile table in `$summary`, with the median
quantile promoted to the top level so that code written against the 1.0.2 return shape keeps
working.

### `case` was accepted and never used

`fqardl()` validated `case` and then ignored it: the design matrix was always Case III
(unrestricted intercept, no trend) while the critical values were switched. Cases 4 and 5
therefore compared a Case III model against Case IV and Case V tables, and `case = 3` and
`case = 5` returned identical coefficients. Any case other than 3 now raises an error.
Correct case handling changes the regressors, not the table, and is scheduled for 2.0.0.

### Standard errors were unnamed, so lookups returned NA silently

In `estimate_qardl()` the coefficient vector was renamed to the design column names while the
standard-error and t-statistic vectors were not. Any lookup such as `std_errors["y_lag1"]`
returned `NA` without warning. All three vectors are now named consistently.

### Critical-value table

The Case III tables are stated explicitly and checked by a test. The 5 percent column of
Tables CI(iii) and CII(iii) has been verified against the published table for every `k` from
1 to 10. The 1 and 10 percent columns of CI(iii) are inherited from earlier versions and have
not been independently verified; the 1 and 10 percent I(1) bounds of CII(iii) are returned as
`NA` rather than guessed. Only the 5 percent level is used in the decision rule.

### Smaller fixes

- `mtnardl()$bounds_test$t_stat` was `NULL`; it is now returned.
- `perform_mtnardl_bounds()` used a hard-coded pair of critical values regardless of the
  number of level terms; it now reads the PSS table for the actual `k`.

## Known limitations, documented rather than silently shipped

Each of these now raises a warning at the point of use. All are scheduled for 2.0.0.

1. **`bootstrap_bounds_test()` does not generate pseudo-samples under the null.** The McNown,
   Sam and Goh (2018) procedure requires each of the three statistics to have its own
   restricted residuals, and requires the level series to be built by accumulation,
   `y*_t = y*_{t-1} + Delta y*_t`. Neither is done. The p-values are indicative only.
2. **`quantile_wald_test()` assumes independence across quantiles.** Quantile-regression
   estimates from the same sample have covariance proportional to
   `min(tau_i, tau_j) - tau_i tau_j`, which is strictly positive, so the correct denominator
   is smaller than the one used. The test is conservative and under-rejects constancy. Fixing
   it requires the joint cross-quantile covariance, which this release does not compute.
3. **Quantile-specific bounds verdicts have no distribution theory.** The PSS critical values
   were simulated for the conditional mean. Treat the per-quantile table as descriptive.
4. **Deterministic cases I, II, IV and V are not available.**

## Verification

This release ships `verify_fqardl_103.R`, sixteen checks, all passing, including:

- the design's left-hand side equals `Delta y_t`;
- the package F statistic equals a hand-computed Wald statistic to machine precision;
- on a design with a true adjustment speed of `-0.35` and a true long-run coefficient of
  `+0.80`, the package recovers `-0.371` and `+0.791`, with a negative error-correction
  t ratio in 100 percent of replications;
- the shipped Case III critical values match Tables CI(iii) and CII(iii) at 5 percent for
  every `k` from 1 to 10.

---

---

# fqardl 1.0.0

## Initial CRAN Release (2026-02-25)

### New Features

#### FQARDL (Fourier Quantile ARDL)
* `fqardl()` - Main function for Fourier Quantile ARDL estimation
* Automatic selection of optimal Fourier frequency (k*)
* Automatic lag selection using BIC, AIC, or HQ criteria
* Estimation across multiple quantiles (τ)
* Long-run and short-run multiplier computation
* PSS bounds test for cointegration
* Bootstrap cointegration testing
* Publication-ready plots

#### FNARDL (Fourier Nonlinear ARDL)
* `fnardl()` - Main function for asymmetric ARDL with Fourier terms
* Variable decomposition into positive/negative partial sums
* Wald tests for long-run and short-run asymmetry
* Dynamic multiplier computation (positive vs negative shocks)
* Cumulative multiplier visualization
* Asymmetry comparison plots

#### Helper Functions
* `fourier_adf()` - Fourier ADF unit root test
* `generate_fourier_terms()` - Generate Fourier trigonometric terms
* `decompose_variables()` - Decompose variables for NARDL
* `perform_bounds_test()` - PSS bounds testing
* `bootstrap_bounds_test()` - Bootstrap cointegration test

#### Visualization
* Coefficient plots across quantiles
* Dynamic multiplier plots
* Cumulative multiplier plots
* Persistence profile plots
* 3D coefficient surfaces
* Heatmaps
* Asymmetry comparison plots

### References

* Pesaran, M. H., Shin, Y., & Smith, R. J. (2001). Bounds testing approaches 
  to the analysis of level relationships. Journal of Applied Econometrics, 
  16(3), 289-326.

* Shin, Y., Yu, B., & Greenwood-Nimmo, M. (2014). Modelling asymmetric 
  cointegration and dynamic multipliers in a nonlinear ARDL framework. 
  In Festschrift in Honor of Peter Schmidt (pp. 281-314). Springer.

* Enders, W., & Lee, J. (2012). A unit root test using a Fourier series to 
  approximate smooth breaks. Oxford Bulletin of Economics and Statistics, 
  74(4), 574-599.

* Cho, J. S., Kim, T., & Shin, Y. (2015). Quantile cointegration in the 
  autoregressive distributed-lag modeling framework. Journal of 
  Econometrics, 188(1), 281-300.

### Acknowledgments

* Dr. Merwan Roudane for the original Stata implementation
* Rufyq Elngeh (رفيق النجاح) for supporting this development

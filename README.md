# fqardl

Fourier ARDL methods for R: Fourier Quantile ARDL, Fourier Nonlinear ARDL, Multi-Threshold NARDL, and Fourier unit root tests.

CRAN: https://cran.r-project.org/package=fqardl

## Correctness notice for version 1.0.2 and earlier

**If you have produced results with fqardl 1.0.2 or earlier, please re-run them with 1.0.3 before relying on them.**

Versions up to 1.0.2 estimated the wrong equation. fqardl, fnardl and mtnardl assembled a Pesaran, Shin and Smith conditional error-correction right-hand side but kept the dependent variable in levels rather than in first differences. The coefficient reported as the error-correction term was therefore 1 + rho instead of rho: positive and close to unity, with a large positive t ratio at every quantile.

### What was affected

| Quantity | Status in 1.0.2 |
| --- | --- |
| Error-correction term and its t ratio | Wrong sign; positive at every quantile |
| Long-run multipliers | Wrong sign and wrong magnitude in essentially every case |
| Bounds test F and t | Rejected the null almost unconditionally; the F was an average of squared marginal t ratios, not a Wald statistic |
| Verdict on cointegration | Binary; the inconclusive region was never reported |
| case = 1, 2, 4, 5 | Accepted but ignored; the design was always Case III |
| Lookups such as std_errors of y_lag1 | Returned NA silently |

Not affected: the Fourier frequency selection, the lag selection, and the partial-sum decomposition in fnardl and mtnardl.

### Scale of the change

On the package's own macro_data, the bounds statistic falls from F = 65.57, reported as "Evidence of cointegration", to F = 2.34, which is below the 5 percent lower bound of 3.79. In a simulation with a true long-run coefficient of +0.80, version 1.0.2 returned a median of -0.471 with the wrong sign in 100 percent of replications; 1.0.3 returns +0.791 with the correct sign in 100 percent of replications. On oil_gdp_data, fnardl previously reported an error-correction coefficient of +0.894 with t = +27.08; it now reports -0.106 with t = -3.21.

### Still not fixed in 1.0.3

These now emit a warning at the point of use and are scheduled for 2.0.0. The bootstrap does not generate pseudo-samples under the null. The cross-quantile Wald test assumes independence across quantiles and is therefore conservative. Quantile-specific bounds verdicts use critical values that Pesaran, Shin and Smith simulated for the conditional mean, so they should be read as descriptive.

### Verification

Version 1.0.3 adds tests that assert the design's left-hand side equals Delta y, that the bounds F equals a hand-computed Wald statistic, that the estimator recovers a known long-run coefficient from simulated data, and that the shipped Case III critical values match Pesaran, Shin and Smith Tables CI(iii) and CII(iii) at 5 percent for every k from 1 to 10.

With thanks to Prof. Chi-Lu (Edward) Peng, National Kaohsiung University of Science and Technology, who reported the sign anomaly that led to this review.

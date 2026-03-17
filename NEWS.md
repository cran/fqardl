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

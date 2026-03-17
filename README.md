# fqardl: Fourier ARDL Methods for R

[![CRAN Status](https://www.r-pkg.org/badges/version/fqardl)](https://CRAN.R-project.org/package=fqardl)
[![R-CMD-check](https://github.com/muhammedalkhalaf/fqardl/workflows/R-CMD-check/badge.svg)](https://github.com/muhammedalkhalaf/fqardl/actions)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

**Comprehensive ARDL methodologies for cointegration analysis with structural breaks and asymmetric effects.**

## 🎯 Overview

`fqardl` provides a complete suite of advanced ARDL (Autoregressive Distributed Lag) methods for time series econometric analysis. Originally ported from Stata/Python implementations by Dr. Merwan Roudane.

## ✨ Features

| Method | Description | Reference |
|--------|-------------|-----------|
| **FQARDL** | Fourier Quantile ARDL | Quantile effects with smooth breaks |
| **FNARDL** | Fourier Nonlinear ARDL | Asymmetric effects (Shin et al., 2014) |
| **MTNARDL** | Multi-Threshold NARDL | Multiple regime asymmetry |
| **Fourier ADF** | Unit root test | Enders & Lee (2012) |
| **Fourier KPSS** | Stationarity test | Becker, Enders & Lee (2006) |

### Key Capabilities

- 📊 **Automatic selection** of lags and Fourier frequencies (AIC/BIC/HQ)
- 🔄 **Asymmetric decomposition** into positive/negative partial sums
- 📈 **PSS bounds testing** for cointegration (Pesaran et al., 2001)
- 🎲 **Bootstrap cointegration** tests for robust inference
- 📉 **Dynamic multipliers** with confidence intervals
- 📋 **Publication-ready** tables and visualizations

## 📦 Installation

### From CRAN (when available)
```r
install.packages("fqardl")
```

### From GitHub
```r
# install.packages("devtools")
devtools::install_github("muhammedalkhalaf/fqardl")
```

## 🚀 Quick Start

### Fourier Unit Root Test
```r
library(fqardl)

# Fourier ADF test for unit root with structural breaks
result <- fourier_adf_test(y, model = "c", max_freq = 3)
print(result)

# Complete analysis (ADF + KPSS)
analysis <- fourier_unit_root_analysis(y, name = "GDP")
```

### Fourier Nonlinear ARDL (Asymmetric Effects)
```r
# Test for asymmetric effects of oil price on GDP
result <- fnardl(
  formula = gdp ~ oil_price + exchange_rate,
  data = macro_data,
  decompose = c("oil_price"),  # Decompose into + and -
  max_k = 3,
  max_p = 4,
  max_q = 4
)

summary(result)
plot(result, type = "asymmetry")
plot(result, type = "dynamic", variable = "oil_price")
```

### Fourier Quantile ARDL (Quantile Effects)
```r
# Analyze effects across the distribution
result <- fqardl(
  formula = gdp ~ inflation + interest_rate,
  data = macro_data,
  tau = c(0.1, 0.25, 0.5, 0.75, 0.9),
  max_k = 3,
  bootstrap = TRUE
)

summary(result)
```

### Multi-Threshold NARDL (Multiple Regimes)
```r
# Multiple threshold analysis
result <- mtnardl(
  formula = stock_return ~ oil_change,
  data = market_data,
  decompose = "oil_change",
  thresholds = list(oil_change = c(-5, 0, 5)),  # 4 regimes
  max_p = 4,
  max_q = 4
)

summary(result)
```

## 📐 Theoretical Background

### NARDL (Shin et al., 2014)

Decomposes independent variable into positive and negative partial sums:

```
x⁺ₜ = Σ max(Δxⱼ, 0)
x⁻ₜ = Σ min(Δxⱼ, 0)
```

Error Correction Model:
```
Δyₜ = α + ρyₜ₋₁ + θ⁺x⁺ₜ₋₁ + θ⁻x⁻ₜ₋₁ + short-run dynamics + εₜ
```

Long-run multipliers: L⁺ = -θ⁺/ρ, L⁻ = -θ⁻/ρ

### Fourier Approximation (Enders & Lee, 2012)

Captures smooth structural breaks using trigonometric terms:

```
sin(2πkt/T) and cos(2πkt/T)
```

where k is the optimal frequency selected by minimizing information criterion.

## 📊 Output Example

```
=================================================================
   Fourier Nonlinear ARDL (FNARDL) - Summary
=================================================================

MODEL SPECIFICATION
-------------------
Dependent: gdp
Decomposed: oil_price
Fourier k: 2 | Lags: ARDL(3, 2)

ASYMMETRIC LONG-RUN MULTIPLIERS
-------------------------------
oil_price(+): -0.3421
oil_price(-):  0.5678
  Asymmetry test: Wald = 12.453, p = 0.0004 (Asymmetric)

BOUNDS TEST
-----------
F-stat: 8.2341 | Decision: Cointegration exists

ERROR CORRECTION TERM
---------------------
ECT (phi): -0.2156
Half-life: 2.87 periods
=================================================================
```

## 📚 References

- **Shin, Y., Yu, B., & Greenwood-Nimmo, M.** (2014). Modelling Asymmetric Cointegration and Dynamic Multipliers in a Nonlinear ARDL Framework. In: Festschrift in Honor of Peter Schmidt. Springer.

- **Enders, W., & Lee, J.** (2012). The flexible Fourier form and Dickey-Fuller type unit root tests. Economics Letters, 117(1), 196-199.

- **Becker, R., Enders, W., & Lee, J.** (2006). A stationarity test in the presence of an unknown number of smooth breaks. Journal of Time Series Analysis, 27(3), 381-409.

- **Pesaran, M. H., Shin, Y., & Smith, R. J.** (2001). Bounds testing approaches to the analysis of level relationships. Journal of Applied Econometrics, 16(3), 289-326.

## 👨‍💻 Authors

- **Muhammad Alkhalaf** - R implementation - [ORCID](https://orcid.org/0009-0002-2677-9246)
- **Dr. Merwan Roudane** - Original Stata/Python implementation

## 🏢 Affiliation

**Rufyq Elngeh** (رفيق النجاح) - Academic & Business Services  
🌐 [www.rufyqelngeh.com](https://www.rufyqelngeh.com)

## 📄 License

GPL-3 © Muhammad Alkhalaf

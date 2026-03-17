\name{fqardl}
\alias{fqardl}
\title{Fourier Quantile ARDL Estimation}
\description{
Estimates the Fourier Quantile Autoregressive Distributed Lag (FQARDL) model.
This methodology extends QARDL by incorporating Fourier trigonometric terms
to capture smooth structural breaks without prior knowledge of break timing.
}
\usage{
fqardl(formula, data, tau = c(0.25, 0.5, 0.75), max_p = 4, max_q = 4, 
       max_k = 3, criterion = c("BIC", "AIC", "HQ"), case = 3, 
       bootstrap = FALSE, n_boot = 1000, seed = NULL, verbose = TRUE)
}
\arguments{
  \item{formula}{A formula of the form y ~ x1 + x2 + ...}
  \item{data}{A data frame containing the time series variables.}
  \item{tau}{Numeric vector of quantiles to estimate.}
  \item{max_p}{Maximum lag for dependent variable.}
  \item{max_q}{Maximum lag for independent variables.}
  \item{max_k}{Maximum Fourier frequency to test.}
  \item{criterion}{Information criterion for lag selection ("BIC", "AIC", "HQ").}
  \item{case}{Model case (1-5) following Pesaran et al. (2001).}
  \item{bootstrap}{Logical, perform bootstrap cointegration test.}
  \item{n_boot}{Number of bootstrap replications.}
  \item{seed}{Random seed for reproducibility.}
  \item{verbose}{Logical. Print progress messages (default: TRUE).}
}
\value{
An object of class "fqardl" containing:
\item{coefficients}{Estimated coefficients for each quantile}
\item{long_run}{Long-run multipliers}
\item{short_run}{Short-run multipliers}
\item{optimal_k}{Optimal Fourier frequency}
\item{optimal_lags}{Optimal lag structure}
\item{bounds_test}{Results of bounds test for cointegration}
\item{diagnostics}{Model diagnostics}
}
\references{
Pesaran, M. H., Shin, Y., & Smith, R. J. (2001). Bounds testing approaches 
to the analysis of level relationships.
}
\author{Muhammad Alkhalaf}
\examples{
\donttest{
data(macro_data)
result <- fqardl(gdp ~ inflation + interest_rate, 
                 data = macro_data,
                 tau = c(0.25, 0.5, 0.75))
summary(result)
}
}

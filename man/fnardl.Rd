\name{fnardl}
\alias{fnardl}
\title{Fourier Nonlinear ARDL Estimation}
\description{
Estimates Fourier Nonlinear ARDL models with asymmetric cointegration 
following Shin, Yu & Greenwood-Nimmo (2014).
}
\usage{
fnardl(formula, data, decompose = NULL, max_p = 4, max_q = 4, max_k = 3,
       criterion = c("BIC", "AIC", "HQ"), case = 3, bootstrap = FALSE, 
       n_boot = 1000, verbose = TRUE)
}
\arguments{
  \item{formula}{A formula specifying the model.}
  \item{data}{A data frame containing the variables.}
  \item{decompose}{Character vector of variables to decompose into positive/negative. Default NULL.}
  \item{max_p}{Maximum lag order for dependent variable.}
  \item{max_q}{Maximum lag order for independent variables.}
  \item{max_k}{Maximum Fourier frequency.}
  \item{criterion}{Information criterion for selection ("BIC", "AIC", or "HQ").}
  \item{case}{PSS case (1-5), default is 3.}
  \item{bootstrap}{Logical. Whether to perform bootstrap bounds test.}
  \item{n_boot}{Number of bootstrap replications.}
  \item{verbose}{Logical. Print progress messages (default: TRUE).}
}
\value{
An object of class "fnardl" containing estimation results.
}
\references{
Shin, Y., Yu, B., & Greenwood-Nimmo, M. (2014). Modelling Asymmetric 
Cointegration and Dynamic Multipliers in a Nonlinear ARDL Framework.
}
\author{Muhammad Alkhalaf}
\examples{
\donttest{
  data(macro_data, package = "fqardl")
  result <- fnardl(gdp ~ oil_price, data = macro_data, decompose = "oil_price")
  summary(result)
}
}

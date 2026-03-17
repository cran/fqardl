\name{mtnardl}
\alias{mtnardl}
\title{Multi-Threshold Nonlinear ARDL}
\description{
Estimates Multi-Threshold NARDL models with multiple regime asymmetry.
}
\usage{
mtnardl(formula, data, decompose = NULL, thresholds = NULL,
        max_p = 4, max_q = 4, criterion = c("BIC", "AIC", "HQ"), case = 3, verbose = TRUE)
}
\arguments{
  \item{formula}{A formula specifying the model.}
  \item{data}{A data frame containing the variables.}
  \item{decompose}{Variable to decompose (default NULL).}
  \item{thresholds}{List of threshold values (default NULL).}
  \item{max_p}{Maximum lag for dependent variable.}
  \item{max_q}{Maximum lag for independent variables.}
  \item{criterion}{Information criterion for selection.}
  \item{case}{PSS bounds test case.}
  \item{verbose}{Logical. Print progress messages (default: TRUE).}
}
\value{
An object of class "mtnardl".
}
\author{Muhammad Alkhalaf}

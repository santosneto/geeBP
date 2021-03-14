\name{geeBP}
\alias{geeBP}
\title{
Fit Beta Prime Generalized Estimating Equations (GEEBP)
}
\description{
The geeBP function fits generalized estimating equations using the reparametrized beta prime distribution (Bourguignon et al., 2019).
}
\usage{
geeBP(formula, data, id, tol = 0.001, maxiter = 25,
      corstr = "independence", linkmu = "log")
}
\arguments{
  \item{formula}{
an object of class "formula" (or one that can be coerced to that class): a expression of the form "response ~ predictors". See the documentation of \code{\link{glm}} and \code{\link{formula}} for details.
}
  \item{data}{
an optional data frame, list or environment containing the variables in the formula. If not found in data, the variables are taken from environment(formula).
}
  \item{id}{
a vector which identifies the cluster of observations. The length of id should be the same as the number of observations.
}
  \item{tol}{
the tolerance used in the fitting algorithm.
}
  \item{maxiter}{
the maximum number of iterations.
}
  \item{corstr}{
a character string specifying the correlation structure. The following are permitted: "independence", "exchangeable", "AR-1", "unstructured", "one-dependent", "one-dependent-stat", "two-dependent", "two-dependent-stat".
}
  \item{linkmu}{
The link function for mean model. The following are permitted: "log" and "identity".
}
}
\details{
It is important that the input data be sorted by the argument called "id" so that the repeated measures structure is correctly specified.
}
\value{
An object of class "geeBP" representing the fit.
}
\references{
Bourguignon, M., Santos-Neto, M. and Castro, M. (2018+). A new regression model for positive data. ArXiv ~> https://arxiv.org/abs/1804.07734.

Freitas, J.V., Nobre, J., Bourguignon, M., and Santos-Neto, M. (2021) A new approach to modelling positive random variables with repeated measures.
}
\author{
João de Freitas, jvbfreitas@ime.unicamp.br.
}

\examples{
library(geeBP)

# modBP = geeBP(formula, data = data, id = data$subject, corstr = "unstructured", linkmu = "log")

# modBP
# summary(modBP)
# residuals(modBP, "standard")
# plot(modBP)
}

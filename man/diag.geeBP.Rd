\name{diag.geeBP}
\alias{diag.geeBP}
\title{
Compute diagnostics measures for geeBP model
}
\description{
This function can be used to compute some of the regression diagnostics for generalized estimating equations discussed in Venezuela et al. (2007) and Freitas et al.(2021).
}
\usage{
diag.geeBP(model, type = c("cooks distance", "hat"))
}
\arguments{
  \item{model}{
    An R object returned by geeBP.
}
  \item{type}{
    Type of measure: Cook's distance or the hat matrix.
}
}
\details{
}
\value{
}
\references{
Venezuela, M. K.; Botter, D. A.; Sandoval, M. C. Diagnostic techniques in generalized estimating equations. Journal of Statistical Computation and Simulation, Taylor & Francis, v. 77, n. 10, p. 879–888, 2007.

Freitas, J.V., Nobre, J., Bourguignon, M., and Santos-Neto, M. (2018+) A new approach to modelling positive random variables with repeated measures. ArXiv ~> link.

}
\author{
João de Freitas, jvbfreitas@ime.unicamp.br.
}
\note{
}

\seealso{
}
\examples{
library(geeBP)

#modBP = geeBP(formula, data = data, id = data$subject, corstr = "unstructured", linkmu = "log")

## Cook's distance chart
#CD = diag.geeBP(modBP, type = "cooks distance")
#cluster = rep(1:modBP$nclusters,each = modBP$max.id)
#qplot(cluster, CD/sum(CD))
}

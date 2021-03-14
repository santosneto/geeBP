\name{LI.geeBP}
\alias{LI.geeBP}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Compute local influence measures for geeBP model
}
\description{
This function can be used to compute the local influence measures under certain perturbations for beta prime generalized estimating equations based in Venezuela et al. (2011) and developed in Freitas et al.(2021).
}
\usage{
LI.geeBP(model, pert = c("case-weight", "response", "precision",
    "correlation", "covariate"), covariate = NULL,
    type = c("normal", "standard"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{
an R object returned by geeBP.
}
  \item{pert}{
a character string specifying the perturbation scheme. The following are permitted: "case-weight", "response", "precision", "correlation" and "covariate".
}
  \item{covariate}{
position of the covariate of interest in the coefficients vector if the covariance perturbation scheme is chosen.
}
  \item{type}{
type of measure returned vector: normal or standardized.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
Venezuela, M. K.; Sandoval, M. C.; Botter, D. A. Local influence in estimating equations. Computational Statistics & Data Analysis, v. 55, n. 4, p. 1867 – 1883, 2011.

Freitas, J.V., Nobre, J., Bourguignon, M., and Santos-Neto, M. (2019+) A new approach to modelling positive random variables with repeated measures. ArXiv ~> link.
}
\author{
João de Freitas, jvbfreitas@ime.unicamp.br.
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
library(nlme)
library(geeBP)

data(Oxboys)

cluster <- as.numeric(as.vector(Oxboys[,1]))
modBP = geeBP(height ~ age, data = Oxboys, id = cluster,
              corstr = "AR-1", linkmu = "identity")

qplot(1:modBP$nobs,LI.geeBP(modBP,pert = "case-weight", type = "standard"))
qplot(1:modBP$nobs,LI.geeBP(modBP,pert = "covariate", covariate = 2,
    type = "standard"))

}

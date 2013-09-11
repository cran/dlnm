\name{poly}
\alias{poly}

\title{ Generate a Basis Matrix of Polynomials }

\description{
The function generates a basis matrix of polynomial transformations. It is meant to be used internally by \code{\link{onebasis}} and \code{\link{crossbasis}} and not directly run by the users.
}

\usage{
poly(x, degree=1, scale, int=FALSE)

}

\arguments{
  \item{x }{ the predictor variable. Missing values are allowed.}
  \item{degree }{ Numerical scalar defining the degree of the polynomial.}
  \item{scale }{ Scaling factor. Default to the maximum of the absolute value of \code{x}.}
  \item{int }{ Logical. If \code{TRUE}, an intercept is included in the basis matrix. See Details below.}
}

\details{
The predictor vector is scaled by default to avoid numerical problem with powers of very high/low values.
}

\value{
A matrix object of class \code{"poly"}. It contains the attributes \code{degree}, \code{scale} and \code{int}, with values which can be different than the arguments provided due to internal reset.}

\author{Antonio Gasparrini, \email{antonio.gasparrini@lshtm.ac.uk}}

\note{
This function is mainly used internally thorugh \code{\link{onebasis}} and \code{\link{crossbasis}} to create basis and cross-basis matrices, respectively. It is not exported in the namespace to avoid conflicts with the function with the same name in the package \pkg{stats}, and can be accessed through the triple colon operator '\code{:::}' (see Examples below). 
}

\seealso{
\code{\link{onebasis}} to generate basis matrices and \code{\link{crossbasis}} to generate cross-basis matrices.

See \code{\link{dlnm-package}} for an introduction to the package and for links to package vignettes providing more detailed information.
}

\examples{
### simple use (accessing non-exported function through ':::')
dlnm:::poly(1:5, degree=3)
dlnm:::poly(1:5, degree=3, int=TRUE)

### use as an internal function in onebasis
b <- onebasis(chicagoNMMAPS$pm10, "poly", degree=3)
summary(b)
model <- glm(death ~ b, family=quasipoisson(), chicagoNMMAPS)
pred <- crosspred(b, model, at=0:60)
plot(pred, xlab="PM10", ylab="RR", main="RR for PM10")
}

\keyword{smooth}


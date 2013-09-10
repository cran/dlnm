\name{onebasis}
\alias{onebasis}
\alias{summary.onebasis}
\alias{mkbasis}
\alias{mklagbasis}

\title{ Generate a Basis Matrix for Different Functions }

\description{
The function generates the basis matrix for a predictor vector, optionally centering on a specific predictor value. The function operates as a wrapper to existing or user-defined functions. Amongst other options, default choices include splines, polynomials, strata and linear threshold functions.
}

\usage{
onebasis(x, fun="ns", cen, ...)

\method{summary}{onebasis}(object, ...)
}

\arguments{
  \item{x }{ the predictor variable. Missing values are allowed.}
  \item{fun }{ character scalar with the name of the function to be called. See Details below.}
  \item{cen }{ logical or a numeric scalar. It specifies the centering value, then used as a reference for predictions. Setting \code{cen=FALSE} is never recommended. See Note below.}
  \item{\dots }{ additional arguments to be passed to the function specified by \code{fun} or to \code{summary}.}
  \item{object }{ a object of class \code{"onebasis"}.}
}

\details{
The function \code{onebasis} is a wrapper to existing functions which are called internally to produce different types of basis matrices in a pre-defined format. Its main use is to be called by \code{\link{crossbasis}} to generate cross-basis matrices for modelling bi-dimensional exposure-lag-response associations thorugh distributed lag non-linear models. However, it can be used also for simplifying the modelling of uni-dimensional exposure-response relationship.

The function to be called is chosen through the argument \code{fun}. Standard choices are:

\itemize{
  \item \bold{\code{"ns"}} and \bold{\code{"bs"}}: natural cubic B-splines or B-splines of various degree. Performed through a call to functions \code{\link[splines]{ns}} or \code{\link[splines]{bs}} from package \pkg{splines}. Arguments passed through \code{\dots} may include \code{df}, \code{knots}, \code{intercept}, and \code{Boundary.knots}.
  
  \item \bold{\code{"poly"}}: polynomials functions. Performed through a call to the internal function \code{\link{poly}} (be aware that this is different from \code{\link[stats]{poly}} in the package \pkg{stats}). Arguments passed through \code{\dots} may include \code{degree} and \code{int}.
  
  \item \bold{\code{"strata"}}: indicator variables defining strata. Performed through a call to the function \code{\link{strata}}. Arguments passed through \code{\dots} may include \code{breaks}, \code{df} and \code{int}.
  
  \item \bold{\code{"thr"}}: high, low or double linear threshold functions. Performed through a call to the function \code{\link{thr}}. Arguments passed through \code{\dots} may include \code{thr.value}, \code{side} and \code{int}.
  
  \item \bold{\code{"integer"}}: indicator variables for each integer value. Performed through a call to the internal function \code{\link{integer}} (be aware that this is different from the function \code{\link[base]{integer}} in the package \pkg{base}). Arguments passed through \code{\dots} may include \code{int}.
  
  \item \bold{\code{"lin"}}: linear functions. Performed through a call to the internal function \code{\link{lin}}. Arguments passed through \code{\dots} may include \code{int}.
}

The help pages of the called functions provides additional information. In particular, the option \code{"lin"} and \code{"integer"} are usually applied for defining constrained and unconstrained distributed lag models.

In addition, any other existing or user-defined function can be potentially called through \code{onebasis}. The function should have the first argument \code{x} defining the vector to be transformed. It also should return a vector or matrix of transformed variables, with attributes including the arguments of the function itself which exactly define the transformations.

Results from models including basis functions are interpreted here relatively to a reference value of the predictor, determined automatically or through a centering point (see also Note below). In the latter case, the centering value is chosen by \code{cen} (if a numeric scalar), or fixed at the mean if \code{cen=TRUE}. The basis is uncentered for \code{cen=FALSE}. Use \code{\link{summary}} to check the result.
}

\value{
A matrix object of class \code{"onebasis"} which can be included in a model formula in order to estimate the association. It contains the attributes \code{fun}, \code{cen}, \code{range} (range of the original vector of observations) and additional attributes specific to the chosen function. The function \code{\link{summary.onebasis}} returns a summary of the basis matrix and the related attributes.
}

\author{Antonio Gasparrini, \email{antonio.gasparrini@lshtm.ac.uk}}

\note{
This function offers a wide range of options about modelling the shape of the exposure-response relationships, also extending the use of existing functions, which can be centered here. The function \code{\link{crosspred}} can be called on objects of class \code{"onebasis"} in order to obtain predictions and plotting of such unidimensional associations. If more than one variable is transformed through \code{onebasis} in the same model, different names must be specified. 

This function has replaced the two old functions \code{mkbasis} and \code{mklagbasis} since version 1.5.0.

Centering is always recommended for continuous functions such as splines or polynomials. Centering is prevented for functions \code{\link{strata}}, \code{\link{thr}} and \code{\link{integer}}. The inclusion of the intercept term nullifies the centering.
}

\section{Warnings}{
Meaningless combinations of arguments could lead to collinear variables, with identifiability problems in the model. The function \code{onebasis} does not perform many checks on the arguments provided. The user is expected to provide valid arguments.
}

\seealso{
\code{\link{crossbasis}} to generate cross-basis matrices. \code{\link{crosspred}} to obtain predictions after model fitting. The method function \code{\link[=plot.crosspred]{plot}} to plot several type of graphs.

See \code{\link{dlnm-package}} for an introduction to the package and for links to package vignettes providing more detailed information.
}

\examples{
### a polynomial transformation of a simple vector, centered and uncentered
onebasis(1:5, "poly", degree=3)
onebasis(1:5, "poly", degree=3, cen=FALSE)

### a low linear threshold parameterization, with and without intercept
onebasis(1:5, "thr", thr=3, side="l")
onebasis(1:5, "thr", thr=3, side="l", int=TRUE)

### relationship between PM10 and mortality estimated by a step function
b <- onebasis(chicagoNMMAPS$pm10, "strata", breaks=c(20,40))
summary(b)
model <- lm(death ~ b, family=quasipoisson(), chicagoNMMAPS)
pred <- crosspred(b, model, at=0:60)
plot(pred, xlab="PM10", ylab="RR", main="RR for PM10")

### relationship between temperature and mortality: double threshold
b <- onebasis(chicagoNMMAPS$temp, "thr", thr=c(10,25))
summary(b)
model <- lm(death ~ b, family=quasipoisson(), chicagoNMMAPS)
pred <- crosspred(b, model, by=1)
plot(pred, xlab="Temperature (C)", ylab="RR", main="RR for temperature")

### extending the example for the 'ns' function in package splines
b <- onebasis(women$height, df=5)
summary(b)
model <- lm(weight ~ b, data=women)
pred <- crosspred(b, model)
plot(pred, xlab="Height (in)", ylab="Weight (lb) difference",
  main="Association between weight and height")
  
### use with a user-defined function with proper attributes
mylog <- function(x,scale=min(x,na.rm=TRUE)+1) {
  basis <- log(x+scale)
  attributes(basis)$scale <- scale
  return(basis)
}
mylog(-2:5)
onebasis(0:5,"mylog",cen=FALSE)
onebasis(0:5,"mylog")
}

\keyword{smooth}


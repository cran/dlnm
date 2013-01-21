\name{onebasis}
\alias{onebasis}
\alias{summary.onebasis}
\alias{mkbasis}
\alias{mklagbasis}

\title{ Generate a Basis Matrix for Different Functions }

\description{
The function generates the basis matrix for a predictor vector, choosing among a set of possible basis functions. Specifically, different types of splines, polynomials, strata and linear threshold functions.
}

\usage{
onebasis(x, type="ns", df=1, degree=1, knots=NULL, bound, int=FALSE, cen)

\method{summary}{onebasis}(object, ...)
}

\arguments{
  \item{x }{ the predictor variable. Missing values are allowed.}
  \item{type }{ type of basis. See Details below for the list of possible choices.}
  \item{df }{ dimension of the basis, equivalent to number of degrees of freedom. They depend on \code{knots} if provided, or on \code{degree} for \code{type="poly"}.}
  \item{degree }{ degree of polynomial. Used only for \code{type} equal to \code{"bs"} (degree of the piecewise polynomial for the B-spline) or \code{"poly"} (degree of the polynomial).}
  \item{knots }{ knots location for the basis function. They specify the position of the internal knots for \code{"ns"} and \code{"bs"}, the internal cut-off points for \code{"strata"} (defining right-open intervals) and the threshold(s)/cut-off points for \code{"lthr"}, \code{"hthr"} and \code{"dthr"}. If provided, are automatically ordered and made unique, determining the value of \code{df}. If only \code{df} is provided, \code{knots} are placed at equally spaced quantiles.}
  \item{bound }{ boundary knots (sometimes called external knots) for \code{type} equal to \code{"ns"} and \code{"bs"}, or external boundaries for \code{type} equal \code{"strata"}. Default to the range of the variable.}
  \item{int }{ logical. If \code{TRUE}, an 'intercept' is included in the basis matrix, with different meanings depending on \code{type} above: see Details below.}
  \item{cen }{ logical or a numeric scalar. It specifies the value the basis functions for the space of predictor are centered at, then used as a reference for predictions. Setting \code{cen=FALSE} is never recommended. See Note below.}
  \item{object }{ a object of class \code{"onebasis"}.}
  \item{\dots }{ additional arguments to be passed to \code{summary}.}
}

\details{
The value in \code{type} defines the basis function. It must  be one of:

\bold{\code{"ns"}}: natural cubic B-splines (constrained to be linear beyond the boundary knots). Specified by \code{knots} (internal knots) and \code{bound} (boundary or external knots). If \code{knots} is provided, the dimension \code{df} is set to \code{length(knots)+1+int}. An intercept is included if \code{int=T} (corresponding to a vector of 1's if \code{df=1}, involving a more complex parameterization for \code{df>1}). See the functions \code{\link[splines]{ns}} for additional information. The transformed variables can be centered by \code{cen}.

\bold{\code{"bs"}}: B-splines characterized by \code{degree} (degree of the piecewise polynomial). Specified by \code{knots} (internal knots) and \code{bound} (boundary or external knots). If \code{knots} is provided, the dimension \code{df} is set to \code{length(knots)+degree+int}; if not, \code{df} must be higher than \code{degree+int}. An intercept is included if \code{int=T} (corresponding to a vector of 1's if \code{df=1}, involving a more complex parameterization for \code{df>1}). See the functions \code{\link[splines]{bs}} for additional information. The transformed variables can be centered by \code{cen}.

\bold{\code{"strata"}}: strata variables (dummy parameterization) determined by cut-off values specified in \code{knots} and \code{bound}, which represent the boundaries for right-open intervals. If \code{knots} is provided, the dimension \code{df} is set to \code{length(knots)+int}. A dummy variable for the reference stratum (the first one by default) is included if \code{int=T}. Never centered.

\bold{\code{"poly"}}: polynomial with power specified by \code{degree}. The dimension \code{df} is set to to \code{degree+int}. An intercept, corresponding to a vector of 1's (the power 0 of the polynomial) is included if \code{int=T}. The transformed variables can be centered by \code{cen}.

\bold{\code{"integer"}}: strata variables (dummy parameterization) for each integer values, expressly created to specify an unconstrained function in the space of lags. \code{df} is set automatically to the number of integer values minus 1 plus \code{int}. A dummy variable for the reference stratum (the first one by default) is included if \code{int=T}. Never centered.

\bold{\code{"hthr"}}, \bold{\code{"lthr"}}: high and low threshold parameterization, with a linear relationship above or below the threshold, respectively, and flat otherwise. The threshold is chosen by \code{knots}: if more than one is provided, a piecewise linear relationship is applied above the first knot or below the last one, respectively, with the slope changing at each further knot. \code{df} is automatically set to \code{length(knots)+int}. An intercept (corresponding to a vector of 1's) is included if \code{int=T}. Never centered.

\bold{\code{"dthr"}}: double threshold parameterization (2 independent linear relationships above the second and below the first threshold, flat between them). The thresholds are chosen by \code{knots}. If only one is provided, the threshold is unique (V-model). If more than 2 are provided, the first and the last ones are chosen. \code{df} is automatically set to \code{2+int}. An intercept (corresponding to a vector of 1's) is included if \code{int=T}. Never centered.

\bold{\code{"lin"}}: linear relationship (untransformed apart from optional centering). \code{df} is automatically set to \code{1+int}. An intercept (corresponding to a vector of 1's) is included if \code{int=T}. It can be centered by \code{cen}. 

Results from models including basis functions are interpreted here relatively to a reference value of the predictor, determined automatically or through a centering point (see also Note below). In the latter case, the centering value is chosen by \code{cen} (if a numeric scalar), or fixed at the mean if \code{cen=TRUE}. The basis is uncentered for \code{cen=FALSE}. Some other arguments can be automatically changed for not sensible combinations, or set to \code{NULL} if not required. Use \code{\link{summary.onebasis}} to check the result.

For a detailed illustration of the use of the function, see:

\code{vignette("dlnmOverview")}
}

\value{
A matrix object of class \code{"onebasis"} which can be included in a model formula in order to estimate the association. It contains the attribute \code{range} (range of the original vector of observations) and additional attributes corresponding to the arguments above. The function \code{\link{summary.onebasis}} returns a summary of the basis matrix and the related attributes, and can be used to check the options for the chosen basis function.
}

\author{Antonio Gasparrini, \email{antonio.gasparrini@lshtm.ac.uk}}

\note{
This function offers a wide range of options about modelling the shape of the exposure-response relationship, also extending the use of functions for splines or polynomials, which can be centered here. The function is called by \code{\link{crossbasis}} to generate cross-basis matrices. In addition, other functions in the package \pkg{dlnm} may be applied with these simple basis functions to obtain prediction (see \code{\link{crosspred}}) and to plot the results (see \code{\link{plot.crosspred}}).

This function has replaced the two old functions \code{mkbasis} and \code{mklagbasis} since version 1.5.0.

The name of the basis object will be used by \code{\link{crosspred}} in order to extract the related estimated parameters. If more than one variable is transformed through basis functions in the same model, different names must be specified. 

For continuous functions specified with \code{type} equal to \code{"ns"}, \code{"bs"}, \code{"poly"} or \code{"lin"}, the reference is set at the value specified by \code{cen}. For the other choices, the reference is automatic:  for \code{type} equal to \code{"strata"} and \code{"integer"}, the reference is the first interval, while for \code{vartype} equal to \code{"hthr"}, \code{"lthr"} and \code{"dthr"}, the reference is the region of null effect below, above or between the threshold(s), respectively. The inclusion of the intercept term nullifies the centering.
}

\section{Warnings}{
Meaningless combinations of arguments (for example the inclusion of knots lying outside the range for \code{type} equal to \code{"strata"} or \code{thr}-type) could lead to collinear variables, with identifiability problems in the model and the exclusion of some of them.
}

\seealso{
\code{\link{crossbasis}} to generate cross-basis matrices. \code{\link{crosspred}} to obtain predictions after model fitting. 
\code{\link{plot.crosspred}} to plot several type of graphs.

See \code{\link{dlnm-package}} for an overview of the package and type \code{'vignette(dlnmOverview)'} for a detailed description.
}

\examples{
### relationship between PM10 and mortality estimated by a step function
b <- onebasis(chicagoNMMAPS$pm10, "strata", knots=c(20,40))
summary(b)
model <- glm(death ~ b, family=quasipoisson(), chicagoNMMAPS)
pred <- crosspred(b, model, at=0:60)
plot(pred, xlab="PM10", ylab="RR", main="RR for PM10")

### relationship between temperature and mortality: double threshold
b <- onebasis(chicagoNMMAPS$temp, "dthr", knots=c(10,25))
summary(b)
model <- glm(death ~ b, family=quasipoisson(), chicagoNMMAPS)
pred <- crosspred(b, model, by=1)
plot(pred, xlab="Temperature (C)", ylab="RR", main="RR for temperature")

### extending the example for the 'ns' function in package splines
b <- onebasis(women$height, df=5)
summary(b)
model <- lm(weight ~ b, data=women)
pred <- crosspred(b, model)
plot(pred, xlab="Height (in)", ylab="Weight (lb) difference",
  main="Association between weight and height")
}

\keyword{smooth}


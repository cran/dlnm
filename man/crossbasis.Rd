\name{crossbasis}
\alias{crossbasis}
\alias{summary.crossbasis}

\title{ Generate a Cross-Basis Matrix for a DLNM }

\description{
Generate the basis matrices for the two dimensions of predictor and lags, choosing among a set of possible basis functions. Then, these functions are combined in order to create the related cross-basis matrix, which can be included in a model formula to fit a distributed lag non-linear model (DLNM).
}

\usage{
crossbasis(x, lag=c(0,0), argvar=list(), arglag=list(), group=NULL, ...)

\method{summary}{crossbasis}(object, ...)
}

\arguments{
  \item{x }{ the predictor variable, defined as a numeric vector representing a complete series of ordered observations.}
  \item{lag }{ either an integer scalar or vector of length 2, defining the the maximum lag or the lag range, respectively. If a scalar, the minimum is automatically set ot 0.}
  \item{argvar, arglag }{ lists of arguments to be passed to the function \code{\link{onebasis}} for generating the two basis matrices for predictor and lags, respectively. See Details below.}
 \item{group }{ a factor defining groups of observations, representing multiple series. Each series must be consecutive, complete and ordered.}
 \item{object }{ a object of class \code{"crossbasis"}.}
  \item{\dots }{ additional arguments. See Details below.}
}

\details{
Until version 1.5.0, the function adopted a completely different usage, with different arguments. The compatibility of the old code is retained by the additional arguments passed through \code{...}. The users are however suggested to adopt the current usage.

The arguments in \code{argvar} and \code{arglag} (optionally including \code{type}, \code{df}, \code{degree}, \code{knots}, \code{bound}, \code{int}, \code{cen}) define two set of basis functions for each dimension. The function \code{\link{onebasis}} is called internally, to build the related basis matrices. The \code{argvar} list is applied to \code{x}, in order to generate the matrix for the space of the predictor. The \code{arglag} list is applied to a new vector given by the sequence obtained by \code{lag}, in order to generate the matrix for the space of lags. Then, the two set of basis matrices are combined in order to create the related cross-basis matrix. See \code{\link{onebasis}} for additional information on how to specify each basis.

Results from DLNM are interpreted relatively to a reference value of the predictor, determined automatically or through a centering point. See \code{\link{onebasis}} for further details.

The basis functions for lags are defined with different default arguments than in \code{\link{onebasis}}: specifically, the knots are placed at equally spaced values on the log scale, an intercept is always included (see Warnings below), and the basis is never centered. Some arguments can be automatically changed for not sensible combinations, or set to \code{NULL} if not required. Use \code{\link{summary.crossbasis}} to check the result.

The argument \code{group} defines groups of observations representing independent series. Each series must be consecutive, complete and ordered. \code{crossbasis} is run on each of them applying the same cross-basis functions: default choices (knots position, range, etc.) are taken considering the pooled distribution.

For a detailed illustration of the use of the function, see:

\code{vignette("dlnmOverview")}
}

\value{
A matrix object of class \code{"crossbasis"} which can be included in a model formula in order to fit a DLNM. It contains the attributes \code{range} (range of the original vector of observations), \code{lag} (lag range), \code{argvar} and \code{arglag} (lists of arguments defining the basis functions in each space, which can be modified if compared to the arguments above). The function \code{\link{summary.crossbasis}} returns a summary of the cross-basis matrix and the related attributes, and can be used to check the options for the basis functions chosen for the two dimensions.
}

\references{
Gasparrini A. Distributed lag linear and non-linear models in R: the package dlnm. \emph{Journal of Statistical Software}. 2011; \bold{43}(8):1-20. URL http://www.jstatsoft.org/v43/i08/.

Gasparrini A., Armstrong, B.,Kenward M. G. Distributed lag non-linear models. \emph{Statistics in Medicine}. 2010; \bold{29}(21):2224-2234.
}

\author{Antonio Gasparrini, \email{antonio.gasparrini@lshtm.ac.uk}}

\note{
The values in \code{x} are expected to be equally-spaced (with the interval defining the lag unit) and ordered in time. The series must be complete. Each value in the series of transformed variables is computed also using previous observations included in the lag period considered: therefore, the first observations in the transformed variables up to the maximum lag defined in \code{lag} are set to \code{NA}. Missing values in \code{x} are allowed, but, for the same reason, the same and the next transformed values up to the maximum lag period will be set to \code{NA}. Although correct, this could generate computational problems for DLNMs with long lag periods in the presence of scattered missing observations. If \code{group} is defined, each groups is treated as a separate series (assumed ordered in time).

The name of the crossbasis object will be used by \code{\link{crosspred}} in order to extract the related estimated parameters. If more than one variable is transformed through cross-basis functions in the same model, different names must be specified. 
}

\section{Warnings}{
Meaningless combinations of arguments (for example the inclusion of knots lying outside the range for \code{type} equal to \code{"strata"} or \code{thr}-type) could lead to collinear variables, with identifiability problems in the model and the exclusion of some of them.

It is strongly recommended to avoid the inclusion of an intercept in the basis for \code{x} (\code{int} in \code{argvar} should be \code{FALSE}, as default), otherwise a rank-deficient cross-basis matrix will be specified, causing some of the cross-variables to be excluded in the regression model. Conversely, an intercept is included by default in the basis for the space of lags.
}

\seealso{
\code{\link{onebasis}} to generate one-dimensional basis matrices. \code{\link{crosspred}} to obtain predictions after model fitting. 
\code{\link{plot.crosspred}} to plot several type of graphs.

See \code{\link{dlnm-package}} for an overview of the package and type \code{'vignette(dlnmOverview)'} for a detailed description.
}

\examples{
### simple DLM
### space of predictor: linear relationship for PM10
### space of predictor: 5df natural cubic spline for temperature
### lag function: 4th degree polynomial for PM10 up to lag15
### lag function: strata intervals at lag 0 and 1-3 for temperature

# CREATE THE CROSS-BASIS FOR EACH PREDICTOR AND CHECK WITH SUMMARY
cb1.pm <- crossbasis(chicagoNMMAPS$pm10, lag=15, argvar=list(type="lin",
  cen=FALSE), arglag=list(type="poly",degree=4))
cb1.temp <- crossbasis(chicagoNMMAPS$temp, lag=3, argvar=list(df=5,cen=21),
  arglag=list(type="strata",knots=1))
summary(cb1.pm)
summary(cb1.temp)

# RUN THE MODEL AND GET THE PREDICTION FOR PM10
library(splines)
model1 <- glm(death ~ cb1.pm + cb1.temp + ns(time, 7*14) + dow,
	family=quasipoisson(), chicagoNMMAPS)
pred1.pm <- crosspred(cb1.pm, model1, at=0:20, cumul=TRUE)

# PLOT THE LINEAR ASSOCIATION OF PM10 ALONG LAGS
plot(pred1.pm, "slices", var=10, col=3, ylab="RR", ci.arg=list(density=15,lwd=2),
	main="Effects of a 10-unit increase in PM10 along lags")
plot(pred1.pm, "slices", var=10, cumul=TRUE, ylab="Cumulative RR",
	main="Cumulative effects of a 10-unit increase in PM10 along lags")

# GET THE FIGURES FOR THE OVERALL EFFECTS, WITH CI
pred1.pm$allRRfit["10"]
cbind(pred1.pm$allRRlow, pred1.pm$allRRhigh)["10",]

### See the vignette 'dlnmOverview' for a detailed explanation of this example
}

\keyword{smooth}
\keyword{ts}


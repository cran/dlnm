\name{crosspred}
\alias{crosspred}

\title{ Generate predicted effects for a DLNM }

\description{
Generate predicted effects from a distributed lag non-linear model (DLNM) for a set of values of the original predictor. It returns specific  effects for each combination of values and lags, plus overall and (optionally) cumulative effects (summed up along lags).
}

\usage{
crosspred(crossbasis, model, at=NULL, 
	from=NULL, to=NULL, by=NULL, cumul=FALSE)
}

\arguments{
  \item{crossbasis }{ an object of class \code{"crossbasis"}.}
  \item{model }{ a model object for which the prediction is desired. See Details below.}
  \item{at }{ vector of values used for prediction.}
  \item{from, to }{ range of values used for prediction.}
  \item{by }{ increment of the sequence.}
  \item{cumul }{ logical. If \code{TRUE}, cumulative effects are predicted. See details.}
}

\details{
\code{model} can be an object created by the following regression functions: \code{lm}, \code{glm}, \code{gam} (package \pkg{mgcv}), \code{glm.nb} (package \pkg{MASS}) or \code{gee} (package \pkg{geepack}), \code{clogit} and \code{coxph} (package \pkg{survival}). The model formula can contain more than one \code{crossbasis} objects.

The object \code{crossbasis} must be the same containing the cross-basis matrix included in \code{model}, including its attributes and class. The set of values for which the effects must be computed can be specified by \code{at} or alternatively  by \code{from}/\code{to}/\code{by}. If specified by \code{at}, the values are automatically ordered and made unique. By default, \code{from} and \code{to} correspond to the range of the original vector of observation stored in the \code{crossbasis} object (see \code{\link{crossbasis}}). If \code{by} is not provided, 30 equally spaced values are returned.

Matrices with cumulative effects summed upon lags for each values used for prediction are included if \code{cumul=TRUE}. For a long lag series (i.e. 1000 lags) the routine can be slow. These matrices are required by \code{\link{crossplot}} to graph the cumulative effects along lags.

For a detailed illustration of the use of the functions, see:

\code{vignette("dlnmOverview")}
}

\value{
A list object of class \code{"crosspred"} with the following components:
  \item{predvar }{ vector of observations used for prediction.}
  \item{maxlag }{ a positive value defining the maximum lag.}
  \item{coef, vcov }{ related coefficients and variance-covariance matrix from \code{model}.}
  \item{matfit, matse }{ matrices of effects and related standard errors for each value of \code{predvar} and lag.}
  \item{allfit, allse }{ vectors of total effects and related standard errors for each value of \code{predvar}.}
  \item{cumfit, cumse }{ matrices of cumulative effects (along lags) and related standard errors for each value of \code{predvar} and lag. Computed if \code{cumul=TRUE}.}
  \item{matRRfit }{ exponentiated effects from \code{matfit}.}
  \item{matRRlow, matRRhigh }{ matrices with low and high 95\% confidence intervals for \code{matRRfit}.}
  \item{allRRfit }{ exponentiated total effects from \code{allfit}.}
  \item{cumRRfit }{ exponentiated effects from \code{cumfit}. Computed if \code{cumul=TRUE}.}
  \item{cumRRlow, cumRRhigh }{ matrices with low and high 95\% confidence intervals for \code{cumRRfit}. Computed if \code{cumul=TRUE}.}
  \item{model.class }{ class of the model command used for estimation.}
  \item{model.link }{ a specification for the model link function.}
}

\references{ 
Armstrong, B. Models for the relationship between ambient temperature and daily mortality. \emph{Epidemiology}. 2006, \bold{17}(6):624-31.
}

\author{Antonio Gasparrini, \email{antonio.gasparrini@lshtm.ac.uk}}

\note{
All the effects are reported versus a reference value corresponding to the centering point for continuous functions or to the default values for the other options (see \code{\link{crossbasis}}). Exponentiated effects are included if \code{model.link} is equal to \code{log} or \code{logit}.
}

\section{Warnings}{
In case of collinear variables in the \code{crossbasis} object, some of them are discarded and the related parameters not included in \code{model}. \code{crosspred} will return an error. Check that the specification of the variables in \code{\link{crossbasis}} is meaningful.

The name of the object \code{crossbasis}  will be used by to extract the related estimated parameters from \code{model}. This name must not match the names of other predictors in the model formula. In addition, if more than one variable is transformed by cross-basis functions in the same model, different names must be specified.
}

\seealso{ \code{\link{crossbasis}}, \code{\link{crossplot}}}

\examples{
# Example 2. See crossbasis and crossplot for other examples

### seasonal analysis
### space of predictor: linear effect above 40.3 microgr/m3 for O3
### space of predictor: linear effects below 15C and above 25C for temperature
### lag function: integer lag parameterization (unconstrained) for O3 up to lag5
### lag function: strata intervals at lag 0-1, 2-5 and 6-10 for temperature

data(chicagoNMMAPS)
chicagoNMMAPSseas <- subset(chicagoNMMAPS, month>5 & month<10)
basis.o3 <- crossbasis(chicagoNMMAPSseas$o3, group=chicagoNMMAPSseas$year,
	vartype="hthr", varknots=40.3, lagtype="integer", maxlag=5)
basis.temp <- crossbasis(chicagoNMMAPSseas$temp, group=chicagoNMMAPSseas$year,
	vartype="dthr", varknots=c(15,25), lagtype="strata",
	lagknots=c(2,6), maxlag=10)
summary(basis.o3)
summary(basis.temp)
model <- glm(death ~  basis.o3 + basis.temp + ns(doy, 4) + dow,
	family=quasipoisson(), chicagoNMMAPSseas)
pred.o3 <- crosspred(basis.o3, model, at=c(0:65,40.3,50.3))

crossplot(pred.o3, "slices", var=50.3, ci="bars",
	title="Effects of a 10-unit increase in ozone along lags")
crossplot(pred.o3,"overall",label="Ozone", ylim=c(0.9,1.3),
	title="Overall effects of ozone over 5 days of lag")
# overall effect for a 10-unit increase in ozone over 15 days of lag, with CI
pred.o3$allRRfit["50.3"]
cbind(pred.o3$allRRlow, pred.o3$allRRhigh)["50.3",]
crossplot(pred.o3, label="Ozone", title="3D graph of ozone effect")

### See the vignette 'dlnmOverview' for a detailed explanation of this example
}

\keyword{smooth}
\keyword{ts}

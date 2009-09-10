\name{crosspred}
\alias{crosspred}

\title{ Generate predicted effects for a DLNM }

\description{
Generate predicted effects from a distributed lag non-linear model (DLNM) for a set of values of the original predictor. It reports specific effects for each combination of values and lags, and overall effects (summed up along lags).
}

\usage{
crosspred(crossbasis, model, at=NULL, from=NULL, to=NULL, by=NULL)
}

\arguments{
  \item{crossbasis }{ an object of class \code{"crossbasis"}.}
  \item{model }{ a model object for which the prediction is desired.}
  \item{at }{ vector of values used for prediction.}
  \item{from, to }{ range of values used for prediction.}
  \item{by }{ increment of the sequence.}
}

\details{
The object \code{crossbasis} must be the same containing the cross-basis matrix included in \code{model}, including its attributes and class. The set of values for which the effects must be computed can be specified by \code{at} or alternatively  by \code{from}/\code{to}/\code{by}. If specified by \code{at}, the values are automatically ordered and made unique. By default, \code{from} and \code{to} correspond to the range of the original vector of observation stored in the \code{crossbasis} object (see \code{\link{crossbasis}}). If \code{by} is not provided, 30 equally spaced values are returned.

For a detailed overview of the options, see:

\code{vignette("dlnmOverview")}
}

\value{
A list object of class \code{"crosspred"} with the following components:
  \item{predvar }{ vector of observations used for prediction.}
  \item{maxlag }{ a positive value defining the maximum lag.}
  \item{coef, vcov }{ related coefficients and variance-covariance matrix from \code{model}.}
  \item{matfit, matse }{ matrices of effects and related standard errors for each value of \code{predvar} and lag.}
  \item{allfit, allse }{ vectors of total effects and related standard errors for each value of \code{predvar}.}
  \item{matRRfit }{ exponentiated effects from \code{matfit}.}
  \item{matRRlow, matRRhigh }{ matrices with low and high 95\% confidence intervals for \code{matRRfit}.}
  \item{allRRfit }{ exponentiated total effects from \code{allfit}.}
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
The name of the object \code{crossbasis}  will be used by to extract the related estimated parameters from \code{model}. This name must not match the names of other predictors in the model formula. In addition, if more than one variable is transformed by cross-basis functions in the same model, different names must be specified.
}

\seealso{ \code{\link{crossbasis}}, \code{\link{crossplot}}}

\examples{
# Example 2. See crossbasis and crossplot for other examples

### DLM with threshold for the effect of O3 on mortality up to 5 days of lag
### space of predictor: linear effect above 40.3 microgr/m3 for O3
### space of predictor: linear effects below 10C and above 25C for temperature
### lag function: integer lag parameterization (unconstrained) for O3
### lag function: strata intervals at lag 0-1, 2-6 and 7-15 for temperature

data(chicagoNMMAPS)
basis.o3 <- crossbasis(chicagoNMMAPS$o3, vartype="hthr", varknots=40.3,
	lagtype="integer", maxlag=5)
basis.temp <- crossbasis(chicagoNMMAPS$temp, vartype="dthr", varknots=c(10,25),
	lagtype="strata", lagknots=c(2,7), maxlag=15)
summary(basis.o3)
summary(basis.temp)
model <- glm(death ~  basis.o3 + basis.temp, family=quasipoisson(), chicagoNMMAPS)
pred.o3 <- crosspred(basis.o3, model, at=c(0:65,40.3,50.3))

crossplot(pred.o3, "slices", var=50.3,
	title="Effect of a 10-unit increase in ozone along lags")
# overall effect for a 10-unit increase in ozone over 15 days of lag, with CI
pred.o3$allRRfit["50.3"]
cbind(pred.o3$allRRlow, pred.o3$allRRhigh)["50.3",]
crossplot(pred.o3, label="Ozone", title="3D graph of ozone effect")
crossplot(pred.o3, "contour", label="Ozone", 
	title="Contour graph of ozone effect")

### alternatively, a piecewise linear effect above 35 microgr/m3 for o3,
### 	with an additional change in slope at 50 microgr/m3
### same lag function as before

basis.o3 <- crossbasis(chicagoNMMAPS$o3, vartype="hthr", varknots=c(35,50),
	lagtype="integer", maxlag=5)
summary(basis.o3)
model <- update(model)
pred.o3 <- crosspred(basis.o3, model, at=c(0:65))

crossplot(pred.o3,"overall",label="Ozone",
	title="Overall effect of ozone over 5 days of lag")

### See the vignette 'dlnmOverview' for a detailed explanation of this example
}

\keyword{smooth}
\keyword{ts}

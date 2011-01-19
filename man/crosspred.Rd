\name{crosspred}
\alias{crosspred}
\alias{summary.crosspred}

\title{ Generate predicted effects for a DLNM }

\description{
Generate predicted effects from a distributed lag non-linear model (DLNM) for a set of values of the original predictor. It returns specific  effects for each combination of values and lags, plus overall and (optionally) cumulative effects (summed up along lags).
}

\usage{
crosspred(basis, model=NULL, model.link=NULL, coef=NULL,
	vcov=NULL, at=NULL, from=NULL, to=NULL, by=NULL,
	ci.level=0.95, cumul=FALSE)

\method{summary}{crosspred}(object, ...)
}

\arguments{
  \item{basis }{ an object of class \code{"crossbasis"}.}
  \item{model }{ a model object for which the prediction is desired. See Details below.}
  \item{coef, vcov, model.link }{ user-provided coefficients, (co)variance matrix and model link for the prediction. See Details below.}
  \item{at }{ vector of values used for prediction.}
  \item{from, to, by }{ range of values and increment of the sequence used for prediction.}
  \item{ci.level }{ confidence level for the computation of confidence intervals.}
  \item{cumul }{ logical. If \code{TRUE}, cumulative effects are predicted. See Details.}
  \item{object }{ a object of class \code{"crosspred"}.}
  \item{\dots }{ additional arguments to be passed to \code{summary}.}

}

\details{
\code{model} is the model object including \code{basis}, used to fit the DLNM. It must include methods for \code{\link{coef}} and \code{\link{vcov}}, applied to extract the parameters. For \code{model} classes without these methods, the user can manually extract the related parameters and include them in \code{coef}-\code{vcov}, also specifying \code{model.link}. In this case, the dimensions and order of the first two must match the variables included in \code{basis}.

The object \code{basis} must be the same containing the cross-basis matrix included in \code{model}, including its attributes and class. The set of values for which the effects must be computed can be specified by \code{at} or alternatively  by \code{from}/\code{to}/\code{by}. If specified by \code{at}, the values are automatically ordered and made unique. If \code{at} and \code{by} are not provided, approximately 50 equally-spaced rounded values are returned using \code{\link{pretty}}.

Exponentiated effects are included if \code{model.link} is equal to \code{log} or \code{logit}, together with confidence intervals computed using a normal approximation and a confidence level of \code{ci.level}. \code{model.link} is automatically selected from \code{model} for classes \code{"lm"}-\code{"glm"}-\code{"gam"}-\code{"clogit"}-\code{"coxph"}, but needs to be provided for different classes. Matrices with cumulative effects summed upon lags for each values used for prediction are included if \code{cumul=TRUE}. For a long lag series (i.e. 1000 lags) the routine can be slow.

For a detailed illustration of the use of the functions, see:

\code{vignette("dlnmOverview")}
}

\value{
A list object of class \code{"crosspred"} with the following components:
  \item{predvar }{ vector of observations used for prediction.}
  \item{maxlag }{ a positive value defining the maximum lag.}
  \item{coef, vcov }{ coefficients and their variance-covariance matrix.}
  \item{matfit, matse }{ matrices of effects and related standard errors for each value of \code{predvar} and lag.}
  \item{allfit, allse }{ vectors of total effects and related standard errors for each value of \code{predvar}.}
  \item{cumfit, cumse }{ matrices of cumulative effects (along lags) and related standard errors for each value of \code{predvar} and lag. Computed if \code{cumul=TRUE}.}
  \item{matRRfit }{ matrix of exponentiated effects from \code{matfit}.}
  \item{matRRlow, matRRhigh }{ matrices with low and high confidence intervals for \code{matRRfit}.}
  \item{allRRfit }{ vector of exponentiated total effects from \code{allfit}.}
  \item{allRRlow, allRRhigh }{ vectors with low and high confidence intervals for \code{allRRfit}.}
  \item{cumRRfit }{ matrix of exponentiated effects from \code{cumfit}. Computed if \code{cumul=TRUE}.}
  \item{cumRRlow, cumRRhigh }{ matrices with low and high confidence intervals for \code{cumRRfit}. Computed if \code{cumul=TRUE}.}
\item{ci.level }{ confidence level used for the computation of confidence intervals.}
  \item{model.class }{ class of the model command used for estimation.}
  \item{model.link }{ a specification for the model link function.}
The function \code{\link{summary.crosspred}} returns a summary of the list.
}

\references{
Gasparrini A., Armstrong, B.,Kenward M. G. Distributed lag non-linear models. \emph{Statistics in Medicine}. 2010; \bold{29}(21):2224-2234.
}

\author{Antonio Gasparrini, \email{antonio.gasparrini@lshtm.ac.uk}}

\note{
All the effects are reported versus a reference value corresponding to the centering point for continuous functions or to the default values for the other options (see \code{\link{crossbasis}}). Exponentiated effects are included if \code{model.link} (specified automatically by \code{model} or selected by the user) is equal to \code{log} or \code{logit}.
}

\section{Warnings}{
In case of collinear variables in the \code{basis} object, some of them are discarded and the related parameters not included in \code{model}. Then, \code{crosspred} will return an error. Check that the specification of the variables is meaningful through \code{\link{summary.crossbasis}} .

The name of the object \code{crossbasis}  will be used to extract the related estimated parameters from \code{model}. If more than one variable is transformed by cross-basis functions in the same model, different names must be specified.
}

\seealso{ \code{\link{crossbasis}}, \code{\link{plot.crosspred}}}

\examples{
# Example 2. See crossbasis and plot.crosspred for other examples

### seasonal analysis
### space of predictor: linear effect above 40.3 microgr/m3 for O3
### space of predictor: linear effects below 15C and above 25C for temperature
### lag function: integer lag parameterization (unconstrained) for O3 up to lag5
### lag function: strata intervals at lag 0-1, 2-5 and 6-10 for temperature

# SELECT SUMMER MONTHS OF THE SERIES
chicagoNMMAPSseas <- subset(chicagoNMMAPS, month>5 & month<10)

# CREATE THE CROSS-BASIS FOR EACH PREDICTOR, SPECIFYING THE GROUP VARIABLE
basis.o3 <- crossbasis(chicagoNMMAPSseas$o3, group=chicagoNMMAPSseas$year,
	vartype="hthr", varknots=40.3, lagtype="integer", maxlag=5)
basis.temp <- crossbasis(chicagoNMMAPSseas$temp, group=chicagoNMMAPSseas$year,
	vartype="dthr", varknots=c(15,25), lagtype="strata",
	lagknots=c(2,6), maxlag=10)
summary(basis.o3)
summary(basis.temp)

# RUN THE MODEL AND GET THE PREDICTION FOR 03
model <- glm(death ~  basis.o3 + basis.temp + ns(doy, 4) + dow,
	family=quasipoisson(), chicagoNMMAPSseas)
pred.o3 <- crosspred(basis.o3, model, at=c(0:65,40.3,50.3))

# PLOT THE LINEAR EFFECTS OF 03 ABOVE THE THRESHOLD (WITH 80%CI)
plot(pred.o3, "slices", var=50.3, ci="bars", type="p", pch=19, ci.level=0.80,
	main="Effects of 10-unit increase above the threshold (80CI)")
plot(pred.o3,"overall",xlab="Ozone", ci="lines", ylim=c(0.9,1.3), lwd=2,
	ci.level=0.80, main="Overall effect over 5 days of lag (80CI)")

# GET THE FIGURES FOR THE OVERALL EFFECT ABOVE, WITH CI
pred.o3$allRRfit["50.3"]
cbind(pred.o3$allRRlow, pred.o3$allRRhigh)["50.3",]

# 3-D PLOT WITH DEFAULT AND USER-DEFINED PERSPECTIVE
plot(pred.o3, xlab="Ozone", main="3D graph of ozone effect")
plot(pred.o3, xlab="Ozone", main="3D graph of ozone effect", theta=250, phi=40)

### See the vignette 'dlnmOverview' for a detailed explanation of this example
}

\keyword{smooth}
\keyword{ts}

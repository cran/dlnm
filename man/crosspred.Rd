\name{crosspred}
\alias{crosspred}
\alias{summary.crosspred}

\title{ Generate Predictions for a DLNM }

\description{
The function generates predictions from a distributed lag non-linear model (DLNM). These are interpreted as estimated associations defined on a grid of values of the original predictor and lags, computed versus a reference predictor value. Specifically, predictor-specific associations expressed along the lag space, lag-specific associations expressed along the predictor space, and the overall cumulative association on the whole lag period. The function can be used also for simple basis functions not including lag.
}

\usage{
crosspred(basis, model=NULL, coef=NULL, vcov=NULL, model.link=NULL, at=NULL,
  from=NULL, to=NULL, by=NULL, lag, bylag=1, ci.level=0.95, cumul=FALSE)

\method{summary}{crosspred}(object, ...)
}

\arguments{
  \item{basis }{ an object of class \code{"onebasis"} or \code{"crossbasis"}.}
  \item{model }{ a model object for which the prediction is desired. See Details below.}
  \item{coef, vcov, model.link }{ user-provided coefficients, (co)variance matrix and model link for the prediction. See Details below.}
  \item{at }{ vector of predictor values used for prediction.}
  \item{from, to}{ range of predictor values used for prediction.}
  \item{lag }{ either an integer scalar or vector of length 2, defining the lag range used for prediction. Defalut to values used for estimation.}
  \item{by, bylag}{ increment of the sequences of predictor and lag values used for prediction.}
  \item{ci.level }{ confidence level for the computation of confidence intervals.}
  \item{cumul }{ logical. If \code{TRUE}, incremental cumulative associations along lags are also predicted. See Details.}
  \item{object }{ an object of class \code{"crosspred"}.}
  \item{\dots }{ additional arguments to be passed to \code{summary}.}
}

\details{
\code{model} is the model object including \code{basis}. It must include methods for \code{\link{coef}} and \code{\link{vcov}}, applied to extract the parameters. For \code{model} classes without these methods, the user can manually extract the related parameters and include them in \code{coef}-\code{vcov}, also specifying \code{model.link}. In this case, the dimensions and order of the first two must match the variables included in \code{basis}.

The object \code{basis} must be the same containing the basis or cross-basis matrix included in \code{model}, preserving its attributes and class. The set of predictor values at which predictions must be computed can be specified by \code{at} or alternatively  by \code{from}/\code{to}/\code{by}. If specified by \code{at}, the values are automatically ordered and made unique. If \code{at} and \code{by} are not provided, approximately 50 equally-spaced rounded values are returned using \code{\link{pretty}}. \code{lag} defines the lag values used for prediction, while \code{bylag} determines the lag step (forced to 1 for incremental cumulative associations).

The function can be used to compute predictions for models with simple basis functions not including lag, computed with \code{\link{onebasis}}. In this case, only unlagged predicted associations are returned.

Exponentiated predictions are included if \code{model.link} is equal to \code{log} or \code{logit}, together with confidence intervals computed using a normal approximation and a confidence level of \code{ci.level}. \code{model.link} is automatically selected from \code{model} for classes \code{"lm"}-\code{"glm"}-\code{"gam"}-\code{"clogit"}-\code{"coxph"}, but needs to be provided for different classes. Matrices with incremental cumulative predicted associations along lags at each values used for prediction are included if \code{cumul=TRUE} and if \code{lag[1]} is equal to 0. For a long lag series (e.g. 1000 lags) the routine can be slow.

For a detailed illustration of the use of the functions, see:

\code{vignette("dlnmOverview")}
}

\value{
A list object of class \code{"crosspred"} with the following (optional) components:
  \item{predvar }{ vector of observations used for prediction.}
  \item{lag }{ integer vector defining the lag range used for prediction.}
  \item{bylag }{ increment of the sequence of lag values.}
  \item{coefficients, vcov }{ coefficients and their variance-covariance matrix.}
  \item{matfit, matse }{ matrices of specific predicted associations and related standard errors at each value of \code{predvar} and lag.}
  \item{allfit, allse }{ vectors of the overall cumulative predicted association and related standard errors at each value of \code{predvar} in the chosen lag period.}
  \item{cumfit, cumse }{ matrices of incremental cumulative predicted associations along lags and related standard errors at each value of \code{predvar} and lag. Computed if \code{cumul=TRUE}.}
  \item{matRRfit }{ matrix of exponentiated specific associations from \code{matfit}.}
  \item{matRRlow, matRRhigh }{ matrices with low and high confidence intervals for \code{matRRfit}.}
  \item{allRRfit }{ vector of exponentiated overall cumulative associations from \code{allfit}.}
  \item{allRRlow, allRRhigh }{ vectors with low and high confidence intervals for \code{allRRfit}.}
  \item{cumRRfit }{ matrix of exponentiated incremental cumulative associations from \code{cumfit}. Computed if \code{cumul=TRUE}.}
  \item{cumRRlow, cumRRhigh }{ matrices with low and high confidence intervals for \code{cumRRfit}. Computed if \code{cumul=TRUE} and if \code{lag[1]} is equal to 0.}
\item{ci.level }{ confidence level used for the computation of confidence intervals.}
  \item{model.class }{ class of the model command used for estimation.}
  \item{model.link }{ a specification for the model link function.}
The function \code{\link{summary.crosspred}} returns a summary of the list.
}

\references{
Gasparrini A. Distributed lag linear and non-linear models in R: the package dlnm. \emph{Journal of Statistical Software}. 2011; \bold{43}(8):1-20. [freely available \href{http://www.ag-myresearch.com/jss2011.html}{here}].
  
Gasparrini A., Armstrong, B.,Kenward M. G. Distributed lag non-linear models. \emph{Statistics in Medicine}. 2010; \bold{29}(21):2224-2234. [freely available \href{http://www.ag-myresearch.com/statmed2010}{here}]
}

\author{Antonio Gasparrini, \email{antonio.gasparrini@lshtm.ac.uk}}

\note{
All the predictions are generated using a reference value corresponding to the centering point for continuous functions or to the default values for the other options (see \code{\link{onebasis}} and \code{\link{crossbasis}}). Exponentiated predictions are included if \code{model.link} (specified automatically by \code{model} or selected by the user) is equal to \code{log} or \code{logit}.
}

\section{Warnings}{
In case of collinear variables in the \code{basis} object, some of them are discarded and the related parameters not included in \code{model}. Then, \code{crosspred} will return an error. Check that the specification of the variables is meaningful through \code{\link{summary.crossbasis}} or \code{\link{summary.onebasis}}.

The name of the object \code{basis}  will be used to extract the related estimated parameters from \code{model}. If more than one variable is transformed by cross-basis functions in the same model, different names must be specified.
}

\seealso{
\code{\link{onebasis}} to generate one-dimensional basis matrices. \code{\link{crossbasis}} to generate cross-basis matrices. \code{\link{crossreduce}} to reduce the fit to one dimension. 
The method function \code{\link[=plot.crosspred]{plot}} to plot several type of graphs.

See \code{\link{dlnm-package}} for an overview of the package and type \code{'vignette(dlnmOverview)'} for a detailed description.
}

\examples{
### seasonal analysis
### space of predictor: linear effect above 40.3 microgr/m3 for O3
### space of predictor: linear effects below 15C and above 25C for temperature
### lag function: integer lag parameterization (unconstrained) for O3 up to lag5
### lag function: strata intervals at lag 0-1, 2-5 and 6-10 for temperature

# SELECT SUMMER MONTHS OF THE SERIES
chicagoNMMAPSseas <- subset(chicagoNMMAPS, month>5 & month<10)

# CREATE THE CROSS-BASIS FOR EACH PREDICTOR, SPECIFYING THE GROUP VARIABLE
cb2.o3 <- crossbasis(chicagoNMMAPSseas$o3, lag=5, argvar=list(type="hthr",
  knots=40.3), arglag=list(type="integer"), group=chicagoNMMAPSseas$year)
cb2.temp <- crossbasis(chicagoNMMAPSseas$temp, lag=10,
  argvar=list(type="dthr",knots=c(15,25)), arglag=list(type="strata",
  knots=c(2,6)), group=chicagoNMMAPSseas$year)
summary(cb2.o3)
summary(cb2.temp)

# RUN THE MODEL AND GET THE PREDICTION FOR 03
library(splines)
model2 <- glm(death ~  cb2.o3 + cb2.temp + ns(doy, 4) + ns(year,3) + dow,
  family=quasipoisson(), chicagoNMMAPSseas)
pred2.o3 <- crosspred(cb2.o3, model2, at=c(0:65,40.3,50.3))

# PLOT THE LINEAR ASSOCIATION OF 03 ALONG LAGS (WITH 80%CI)
plot(pred2.o3, "slices", var=50.3, ci="bars", type="p", pch=19, ci.level=0.80,
  main="Association with a 10-unit increase above threshold (80CI)")
# PLOT THE OVERALL CUMULATIVE ASSOCIATION
plot(pred2.o3,"overall",xlab="Ozone", ci="lines", ylim=c(0.9,1.3), lwd=2,
  ci.arg=list(col=1,lty=3), main="Overall cumulative association for 5 lags")

# GET THE FIGURES FOR THE OVERALL CUMULATIVE ASSOCIATION ABOVE, WITH CI
pred2.o3$allRRfit["50.3"]
cbind(pred2.o3$allRRlow, pred2.o3$allRRhigh)["50.3",]

# 3-D PLOT WITH DEFAULT AND USER-DEFINED PERSPECTIVE
plot(pred2.o3, xlab="Ozone", main="3D: default perspective")
plot(pred2.o3, xlab="Ozone", main="3D: different perspective", theta=250, phi=40)
### See the vignette 'dlnmOverview' for a detailed explanation of this example
}

\keyword{smooth}
\keyword{ts}

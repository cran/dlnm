\name{crossreduce}
\alias{crossreduce}
\alias{summary.crossreduce}

\title{ Reduce the Fit of a DLNM to One-Dimensional summaries }

\description{
The function reduces the fit of a bi-dimensional DLNM to summaries defined in the the dimension of predictor or lags only, and re-expresses it in terms of modified parameters of the one-dimensional basis functions chosen for that space. Specifically, it summarizes lag-specific or overall cumulative predicted associations in the predictor dimension, or predictor-specific associations in the lag dimension. It returns the new parameters and predictions.
}

\usage{
crossreduce(basis, model=NULL, type="overall", value=NULL, coef=NULL, vcov=NULL,
  model.link=NULL, at=NULL, from=NULL, to=NULL, by=NULL, lag, bylag=1,
  ci.level=0.95)
  
\method{summary}{crossreduce}(object, ...)
}

\arguments{
  \item{basis }{ an object of class \code{"crossbasis"}.}
  \item{model }{ a model object for which the reduction and prediction are desired. See Details below.}
  \item{coef, vcov, model.link }{ user-provided coefficients, (co)variance matrix and model link for the reduction and then prediction. See Details below.}
  \item{type }{ type of reduction. Possible options are \code{"overall"} (default) for reduction to the overall cumulative association, \code{"lag"} for reduction to a lag-specific association, or \code{"var"} for reduction to a predictor-specific association. See Details below.}
  \item{value }{ the single value of predictor or lag at which predictor-specific or lag-specific associations must be defined, respectively. See Details below.}
  \item{at }{ vector of values used for prediction in the dimension of predictor.}
  \item{from, to}{ range of predictor values used for prediction.}
  \item{lag }{ either an integer scalar or vector of length 2, defining the lag range used for prediction. Defalut to values used for estimation.}
  \item{by, bylag}{ increment of the sequences of predictor and lag values used for prediction.}
  \item{ci.level }{ confidence level for the computation of confidence intervals.}
  \item{object }{  an object of class \code{"crossreduce"}.}
  \item{\dots }{ additional arguments to be passed to \code{summary}.}
}

\details{
\code{model} is the model object including \code{basis}. It must include methods for \code{\link{coef}} and \code{\link{vcov}}, applied to extract the parameters. For \code{model} classes without these methods, the user can manually extract the related parameters and include them in \code{coef}-\code{vcov}, also specifying \code{model.link}. In this case, the dimensions and order of the first two must match the variables included in \code{basis}.

The dimension to which the fit is reduced is chosen by \code{type}, computing summaries for overall cumulative or lag-specific associations expressed in the predictor space, or predictor-specific associations expressed in the lag space. For specific associations, the value at which the reduction is computed is chosen by \code{value}. The function then re-expresses the original fit of the model, defined by the parameters of the bi-dimensional cross-basis functions, in summaries defined by the one-dimensional basis for the related space and a (usually smaller) set of modified parameters.

Similarly to \code{\link{crosspred}}, the object \code{basis} must be the same containing the cross-basis matrix included in \code{object}, including its attributes and class. The optional arguments \code{at} and \code{from}/\code{to}/\code{by} provides the predictor values for predicted associations when the reduction is in the dimension of predictor. \code{lag} and \code{bylag} determine instead the lag values for predictor-specific associations. Exponentiated predictions and confidence intervals are also optionally returned. See \code{\link{crosspred}} for details.

For a detailed illustration of the use of the functions, see:

\code{vignette("dlnmOverview")}
}

\value{
A list object of class \code{"crossreduce"} with the following (optional) components:
  \item{coefficients, vcov }{ reduced parameters of the original fitted model for the chosen dimension.}
  \item{basis }{ basis matrix computed at \code{predvar} or for the sequence of lags defined by \code{lag}, depending on the chosen dimension.}
  \item{type, value }{ type of reduction and (optional) value, as arguments above.}
  \item{predvar }{ vector of observations used for prediction, if the reduction is in the dimension of predictor.}
  \item{lag }{ integer vector defining the lag range.}
  \item{bylag }{ increment of the sequence of lag values.}
  \item{fit, se }{ vectors of the predicted association and related standard errors.}
  \item{RRfit }{ vector of exponentiated predicted associations from \code{fit}.}
  \item{RRlow, RRhigh }{ vectors with low and high confidence intervals for \code{RRfit}.}
\item{ci.level }{ confidence level used for the computation of confidence intervals.}
  \item{model.class }{ class of the model command used for estimation.}
  \item{model.link }{ a specification for the model link function.}
}

\references{
Gasparrini A., Armstrong, B., Kenward M. G. Reducing and meta-analyzing estimates from distributed lag non-linear models.\emph{BMC Medical Research Methodology}. 2013; \bold{13}(1):1. [freely available \href{http://www.ag-myresearch.com/bmcmrm2013}{here}].

Gasparrini A. Distributed lag linear and non-linear models in R: the package dlnm. \emph{Journal of Statistical Software}. 2011; \bold{43}(8):1-20. [freely available \href{http://www.ag-myresearch.com/jss2011}{here}].
  
Gasparrini A., Armstrong, B.,Kenward M. G. Distributed lag non-linear models. \emph{Statistics in Medicine}. 2010; \bold{29}(21):2224-2234. [freely available \href{http://www.ag-myresearch.com/statmed2010}{here}]
}

\author{Antonio Gasparrini, \email{antonio.gasparrini@lshtm.ac.uk}}

\note{
All the predictions are generated using a reference value corresponding to the centering point for continuous functions or to the default values for the other options (see \code{\link{onebasis}} and \code{\link{crossbasis}}). Exponentiated predictions are included if \code{model.link} (specified automatically by \code{object} or selected by the user) is equal to \code{log} or \code{logit}.
}

\section{Warnings}{
In case of collinear variables in the \code{basis} object, some of them are discarded and the related parameters not included in \code{model}. Then, \code{crossreduce} will return an error. Check that the specification of the variables is meaningful through \code{\link{summary.crossbasis}} .

The name of the object \code{basis}  will be used to extract the related estimated parameters from \code{object}. If more than one variable is transformed by cross-basis functions in the same model, different names must be specified.
}

\seealso{
\code{\link{crossbasis}} to generate cross-basis matrices. \code{\link{crosspred}} to obtain predictions after model fitting. 
\code{\link{plot.crossreduce}} to plot several the graphs.

See \code{\link{dlnm-package}} for an overview of the package and type \code{'vignette(dlnmOverview)'} for a detailed description.
}

\examples{
# CREATE THE CROSS-BASIS: DOUBLE THRESHOLD AND NATURAL SPLINE
cb4 <- crossbasis(chicagoNMMAPS$temp, lag=30, argvar=list(type="dthr",
  knots=c(10,25)), arglag=list(df=5))

# RUN THE MODEL AND GET THE PREDICTION FOR TEMPERATURE
library(splines)
model4 <- glm(death ~  cb4 + ns(time, 7*14) + dow,
  family=quasipoisson(), chicagoNMMAPS)
pred4 <- crosspred(cb4, model4, by=1)

# REDUCE TO OVERALL CUMULATIVE ASSOCIATION
redall <- crossreduce(cb4, model4)
summary(redall)
# REDUCE TO LAG-SPECIFIC ASSOCIATION FOR LAG 5
redlag <- crossreduce(cb4, model4, type="lag", value=5)
# REDUCE TO PREDICTOR-SPECIFIC ASSOCIATION AT VALUE 33
redvar <- crossreduce(cb4, model4, type="var", value=33)

# NUMBER OF PARAMETERS OF THE ORIGINAL MODEL
length(coef(pred4))
# NUMBER OF PARAMETERS FOR OVERALL CUMULATIVE AND LAG-SPECIFIC SUMMARIES
length(coef(redall)) ; length(coef(redlag))
# NUMBER OF PARAMETERS FOR PREDICTOR-SPECIFIC SUMMARY
length(coef(redvar))

# TEST: IDENTICAL FIT BETWEEN ORIGINAL AND REDUCED FIT
plot(pred4, "overall", xlab="Temperature", ylab="RR",
  ylim=c(0.8,1.6), main="Overall cumulative association")
lines(redall, ci="lines",col=4,lty=2)
legend("top",c("Original","Reduced"),col=c(2,4),lty=1:2,ins=0.1)

# RECONSTRUCT THE FIT IN TERMS OF ONE-DIMENSIONAL BASIS
b4 <- onebasis(0:30,knots=attributes(cb4)$arglag$knots,int=TRUE,cen=FALSE)
pred4b <- crosspred(b4,coef=coef(redvar),vcov=vcov(redvar),model.link="log",by=1)

# TEST: IDENTICAL FIT BETWEEN ORIGINAL, REDUCED AND RE-CONSTRUCTED
plot(pred4, "slices", var=33, ylab="RR", ylim=c(0.9,1.2),
  main="Predictor-specific association at 33C")
lines(redvar, ci="lines", col=4, lty=2)
points(pred4b, col=1, pch=19, cex=0.6)
legend("top",c("Original","Reduced","Reconstructed"),col=c(2,4,1),lty=c(1:2,NA),
  pch=c(NA,NA,19),pt.cex=0.6,ins=0.1)
}

\keyword{smooth}
\keyword{ts}

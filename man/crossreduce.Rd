\name{crossreduce}
\alias{crossreduce}
\alias{summary.crossreduce}

\title{ Reduce the Fit of a DLNM to One-Dimensional Summaries }

\description{
The function reduces the fit of a bi-dimensional DLNM to summaries defined in the the dimension of predictor or lags only, and re-expresses it in terms of modified parameters of the one-dimensional basis functions chosen for that space.
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
The dimension to which the fit is reduced is chosen by \code{type}, computing summaries for overall cumulative or lag-specific associations defining an exposure-response relationship in the predictor space, or predictor-specific associations defining a lag-response relationship in the lag space. The function re-expresses the original fit of the model, defined by the parameters of the bi-dimensional cross-basis functions, in summaries defined by the one-dimensional basis for the related space and a (usually smaller) set of modified parameters.

Similarly to \code{\link{crosspred}}, the object \code{basis} must be the same containing the cross-basis matrix included in \code{model}, with its attributes and class. The optional arguments \code{at} and \code{from}/\code{to}/\code{by} provides the predictor values for predicted associations when the reduction is in the dimension of predictor. \code{lag} and \code{bylag} determine instead the lag values for predictor-specific associations. For specific associations, the value at which the reduction is computed is chosen by \code{value}. Exponentiated predictions and confidence intervals are also optionally returned. See \code{\link{crosspred}} for details.

The function automatically works with model objects from regression function \code{\link{lm}} and \code{\link{glm}}, \code{\link[mgcv]{gam}} (package \pkg{mgcv}), \code{\link[survival]{coxph}} and \code{\link[survival]{clogit}} (package \pkg{survival}), \code{\link[nlme]{lme}} and \code{\link[nlme]{nlme}} (package \pkg{nlme}), \code{\link[lme4]{lmer}} and \code{\link[lme4]{glmer}} and \code{\link[lme4]{nlmer}} (package \pkg{lme4}), \code{\link[gee]{gee}} (package \pkg{gee}), \code{\link[geepack]{geeglm}} (package \pkg{geepack}). The function also works with any regression function for which \code{\link{coef}} and \code{\link{vcov}} methods are available. Otherwise, the user needs to input the coefficients and associated (co)variance matrix related to the parameters of the crossbasis as arguments \code{coef} and \code{vcov}. In this case, their dimensions and order must match the variables included in \code{basis}.
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
  \item{low, high }{ vectors of confidence intervals for \code{fit}.}
  \item{RRfit }{ vector of exponentiated predicted associations from \code{fit}.}
  \item{RRlow, RRhigh }{ vectors of confidence intervals for \code{RRfit}.}
\item{ci.level }{ confidence level used for the computation of confidence intervals.}
  \item{model.class }{ class of the model command used for estimation.}
  \item{model.link }{ a specification for the model link function.}
}

\references{
Gasparrini A., Armstrong, B., Kenward M. G. Reducing and meta-analyzing estimates from distributed lag non-linear models.\emph{BMC Medical Research Methodology}. 2013; \bold{13}(1):1. [freely available \href{http://www.ag-myresearch.com/bmcmrm2013}{here}].
}

\author{Antonio Gasparrini, \email{antonio.gasparrini@lshtm.ac.uk}}

\note{
All the predictions are generated using a reference value corresponding to the centering point for continuous functions or to the default values for the other options (see \code{\link{onebasis}} and \code{\link{crossbasis}}). Exponentiated predictions are included if \code{model.link} (specified automatically by \code{object} or selected by the user) is equal to \code{log} or \code{logit}.
}

\section{Warnings}{
In case of collinear variables in the \code{basis} object, some of them are discarded and the related parameters not included in \code{model}. Then, \code{crossreduce} will return an error. Check that the specification of the variables is meaningful through \code{\link[=summary.crossbasis]{summary}}.

The name of the object \code{basis}  will be used to extract the related estimated parameters from \code{model}. If more than one variable is transformed by cross-basis functions in the same model, different names must be specified.
}

\seealso{
\code{\link{crossbasis}} to generate cross-basis matrices. \code{\link{crosspred}} to obtain predictions after model fitting. The method function \code{\link[=plot.crossreduce]{plot}} to plot the association.

See \code{\link{dlnm-package}} for an introduction to the package and for links to package vignettes providing more detailed information.
}

\examples{
# create the crossbasis object
lagnk <- 3
lagknots <- exp(((1+log(30))/(lagnk+1) * seq(lagnk))-1)
cb4 <- crossbasis(chicagoNMMAPS$temp, lag=30, argvar=list(fun="thr",
  thr=c(10,25)), arglag=list(knots=lagknots))

# # run the model and get the predictions
library(splines)
model4 <- glm(death ~  cb4 + ns(time, 7*14) + dow,
  family=quasipoisson(), chicagoNMMAPS)
pred4 <- crosspred(cb4, model4, by=1)

# reduce to overall cumulative association
redall <- crossreduce(cb4, model4)
summary(redall)
# reduce to exposure-response association for lag 5
redlag <- crossreduce(cb4, model4, type="lag", value=5)
# reduce to lag-response association for value 33
redvar <- crossreduce(cb4, model4, type="var", value=33)

# compare number of parameters
length(coef(pred4))
length(coef(redall))
length(coef(redlag))
length(coef(redvar))

# test
plot(pred4, "overall", xlab="Temperature", ylab="RR",
  ylim=c(0.8,1.6), main="Overall cumulative association")
lines(redall, ci="lines",col=4,lty=2)
legend("top",c("Original","Reduced"),col=c(2,4),lty=1:2,ins=0.1)

# reconstruct the fit in terms of uni-dimensional function
b4 <- onebasis(0:30,knots=attributes(cb4)$arglag$knots,int=TRUE,cen=FALSE)
pred4b <- crosspred(b4,coef=coef(redvar),vcov=vcov(redvar),model.link="log",by=1)

# test
plot(pred4, "slices", var=33, ylab="RR", ylim=c(0.9,1.2),
  main="Lag-response association at 33C")
lines(redvar, ci="lines", col=4, lty=2)
points(pred4b, col=1, pch=19, cex=0.6)
legend("top",c("Original","Reduced","Reconstructed"),col=c(2,4,1),lty=c(1:2,NA),
  pch=c(NA,NA,19),pt.cex=0.6,ins=0.1)
}

\keyword{smooth}
\keyword{ts}

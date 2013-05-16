\name{plot.crossreduce}
\alias{plot.crossreduce}
\alias{lines.crossreduce}
\alias{points.crossreduce}

\title{ Plot Predictions for a Reduced DLNM }

\description{
High and low-level method functions for graphs of predictions from reduced distributed lag non-linear models (DLNM).
}

\usage{
\method{plot}{crossreduce}(x, ci="area", ci.arg, ci.level=x$ci.level, exp=NULL, ...)

\method{lines}{crossreduce}(x, ci="n", ci.arg, ci.level=x$ci.level, exp=NULL, ...)

\method{points}{crossreduce}(x, ci="n", ci.arg, ci.level=x$ci.level, exp=NULL, ...)
}

\arguments{
  \item{x }{ an object of class \code{"crossreduce"}.}
  \item{ci }{ type of confidence intervals representation: one of \code{"area"}, \code{"bars"}, \code{"lines"} or \code{"n"}. Default to \code{"area"} in high level functions, \code{"n"} for low-level functions.}
  \item{ci.arg }{ list of arguments to be passed to low-level plotting functions to draw the confidence intervals. See Details.}
  \item{ci.level }{ confidence level for the computation of confidence intervals.}
  \item{exp }{ logical. It forces the choice about the exponentiation. See Details.}
  \item{\dots }{ optional graphical arguments. See Details.}
}

\details{
Differently than for plotting functions for \code{crosspred} objects (see the method function \code{\link[=plot.crosspred]{plot}} for objects of class \code{"crosspred"}), the type of the plot is automatically chosen by the dimension and value at which the model has been reduced. Namely, the lag-specific association at the chosen lag value, the predictor-specific association at the chosen predictor value, or the overall cumulative association.

These methods allow a great flexibility in the choice of graphical parameters, specified through arguments of the original plotting functions. See \code{\link[graphics]{plot.default}}, \code{\link[graphics]{lines}} and \code{\link[graphics]{points}} for information on additional graphical arguments. Some arguments, if not specified, are set to different default values than the original functions.

Confidence intervals are plotted for \code{ptype} equal to \code{"overall"} or \code{"slices"}. Their type is determined  by \code{ci}, with options \code{"area"} (default for \code{plot}), \code{"bars"}, \code{"lines"} or \code{"n"} (no confidence intervals, default for \code{points} and \code{lines}). The appearance may be modified through \code{ci.arg}, a list of arguments passed to to low-level plotting functions: \code{\link[graphics]{polygon}} for \code{"area"}, \code{\link[graphics]{segments}} for \code{"bars"} and \code{\link[graphics]{lines}} for \code{"lines"}. See the original functions for a complete list of the arguments. This option offers flexibility in the choice of confidence intervals display. As above, some unspecified arguments are set to different default values.

For a detailed illustration of the use of the functions, see:

\code{vignette("dlnmOverview")}
}

\references{
Gasparrini A., Armstrong, B., Kenward M. G. Reducing and meta-analyzing estimates from distributed lag non-linear models. 2012; (Submitted).
}

\author{Antonio Gasparrini, \email{antonio.gasparrini@lshtm.ac.uk}}

\note{
All the predictions are plotted using a reference value corresponding to the centering point for continuous functions or to the default values for the other options (see \code{\link{onebasis}} and \code{\link{crossbasis}}). Exponentiated predictions are returned by default if \code{x$model.link} is equal to \code{log} or \code{logit}.
}

\seealso{
\code{\link{onebasis}} to generate simple basis matrices. \code{\link{crosspred}} to obtain predictions after model fitting. \code{\link{crossreduce}} to reduce the fit ot one dimension.

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

\keyword{hplot}
\keyword{aplot}

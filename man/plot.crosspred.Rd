\name{plot.crosspred}
\alias{plot.crosspred}
\alias{lines.crosspred}
\alias{points.crosspred}
\alias{crossplot}

\title{ Plot Predictions for a DLNM }

\description{
High and low-level method functions for graphs (3d, contour, slices and overall) of predictions from distributed lag non-linear models (DLNM).
}

\usage{
\method{plot}{crosspred}(x, ptype, var=NULL, lag=NULL, ci="area", ci.arg,
  ci.level=x$ci.level, cumul=FALSE, exp=NULL, ...)

\method{lines}{crosspred}(x, ptype, var=NULL, lag=NULL, ci="n", ci.arg,
  ci.level=x$ci.level, cumul=FALSE, exp=NULL, ...)

\method{points}{crosspred}(x, ptype, var=NULL, lag=NULL, ci="n", ci.arg,
  ci.level=x$ci.level, cumul=FALSE, exp=NULL, ...)
}

\arguments{
  \item{x }{ an object of class \code{"crosspred"}.}
  \item{ptype }{ type of plot. Default to \code{"3d"} for lagged relationship, otherwise \code{"overall"}. See Details below.}
   \item{var, lag }{ vectors (for \code{plot}) or numeric scalars (for \code{lines}-\code{points}) of predictor or lag values for which specific effects must be plotted. Used only if \code{ptype="slices"}.}
  \item{ci }{ type of confidence intervals representation: one of \code{"area"}, \code{"bars"}, \code{"lines"} or \code{"n"}. Default to \code{"area"} in high level functions, \code{"n"} for low-level functions.}
  \item{ci.arg }{ list of arguments to be passed to low-level plotting functions to draw the confidence intervals. See Details.}
  \item{ci.level }{ confidence level for the computation of confidence intervals.}
  \item{cumul }{ logical. If \code{TRUE}, cumulative effects along lags are reported. Used only if \code{type="slices"}. See Details.}
  \item{exp }{ logical. It forces the choice about the exponentiation of effects. See Details.}
  \item{\dots }{ optional graphical arguments. See Details.}
}

\details{
Different plots can be obtained by choosing the following values for the argument \code{ptype}:

\bold{\code{"3d"}}: a 3-D plot of predicted effects on the grid of predictor-lag values. Additional graphical arguments can be included, such as \code{theta}-\code{phi} (perspective), \code{border}-\code{shade} (surface), \code{xlab}-\code{ylab}-\code{zlab} (axis labelling) or \code{col}. See \code{\link[graphics]{persp}} for additional information.

\bold{\code{"contour"}}: a contour/level plot of predicted effects on the grid of predictor-lag values.  Additional graphical arguments can be included, such as \code{plot.title}-\code{plot.axes}-\code{key.title} for titles and axis and key labelling. Arguments \code{x}-\code{y}-\code{z} and \code{col}-\code{level} are automatically set and cannot be specified by the user. See \code{\link[graphics]{filled.contour}} for additional information.

\bold{\code{"overall"}}: a plot of the overall effects (summing up all the single lag contributions). See \code{\link[graphics]{plot.default}}, \code{\link[graphics]{lines}} and \code{\link[graphics]{points}} for information on additional graphical arguments.

\bold{\code{"slices"}}: a (optionally multiple) plot of predictor-specific effects along the lag space, and/or lag-specific effects along the predictor space. Predictor and lag values are chosen by \code{var} and \code{lag}, respectively. See \code{\link[graphics]{plot.default}}, \code{\link[graphics]{lines}} and \code{\link[graphics]{points}} for information on additional graphical arguments.

The method function \code{plot} calls the high-level functions listed above for each \code{ptype}, while \code{lines}-\code{points} add lines or points for \code{ptype} equal to \code{"overall"} or \code{"slices"}. These methods allow a great flexibility in the choice of graphical parameters, specified through arguments of the original plotting functions. Some arguments, if not specified, are set to different default values than the original functions.

Confidence intervals are plotted for \code{ptype} equal to \code{"overall"} or \code{"slices"}. Their type is determined  by \code{ci}, with options \code{"area"} (default for \code{plot}), \code{"bars"}, \code{"lines"} or \code{"n"} (no confidence intervals, default for \code{points} and \code{lines}). The appearance may be modified through \code{ci.arg}, a list of arguments passed to to low-level plotting functions: \code{\link[graphics]{polygon}} for \code{"area"}, \code{\link[graphics]{segments}} for \code{"bars"} and \code{\link[graphics]{lines}} for \code{"lines"}. See the original functions for a complete list of the arguments. This option offers flexibility in the choice of confidence intervals display. As above, some unspecified arguments are set to different default values.

For \code{ptype="slices"}, up to 4 plots for each dimension of predictor and lags are allowed in \code{plot}, while for \code{lines}-\code{points} a single plot in one of the two dimension must be chosen. Cumulative effects along lags are reported if \code{cumul=TRUE}: in this case, the same option must have been set to obtain the prediction saved in the \code{crosspred} object (see \code{\link{crosspred}}). 

For a detailed illustration of the use of the functions, see:

\code{vignette("dlnmOverview")}
}

\references{
Gasparrini A. Distributed lag linear and non-linear models in R: the package dlnm. \emph{Journal of Statistical Software}. 2011; \bold{43}(8):1-20. [freely available \href{http://www.jstatsoft.org/v43/i08/}{here}].
  
Gasparrini A., Armstrong, B.,Kenward M. G. Distributed lag non-linear models. \emph{Statistics in Medicine}. 2010; \bold{29}(21):2224-2234. [freely available \href{http://www.ncbi.nlm.nih.gov/pubmed/20812303}{here}]
}

\author{Antonio Gasparrini, \email{antonio.gasparrini@lshtm.ac.uk}}

\note{
These methods for class \code{"crosspred"} have replaced the old function \code{crossplot} since version 1.3.0. The old function has been kept in the namespace of the package, but its use is discouraged.

All the predictions are plotted using a reference value corresponding to the centering point for continuous functions or to the default values for the other options (see \code{\link{onebasis}} and \code{\link{crossbasis}}). Exponentiated predicted outcomes are returned by default if \code{x$model.link} is equal to \code{log} or \code{logit}.
}

\section{Warnings}{
The values in \code{var} must match those specified in the object \code{crosspred} (see \code{\link{crosspred}}), while the values in \code{lag} must be included in the lag period specified by \code{\link{crossbasis}}.
}

\seealso{
\code{\link{crossbasis}} to generate cross-basis matrices. \code{\link{crosspred}} to obtain predictions after model fitting.

See \code{\link{dlnm-package}} for an overview of the package and type \code{'vignette(dlnmOverview)'} for a detailed description.
}

\examples{
### complex DLNM
### space of predictor: 5df quadratic spline for temperature
### space of predictor: linear effect for PM10
### lag function: 5df natural cubic spline for temperature up to lag30
### lag function: single strata at lag 0-1 for PM10

# CREATE THE CROSS-BASIS FOR EACH PREDICTOR AND CHECK WITH SUMMARY
cb3.pm <- crossbasis(chicagoNMMAPS$pm10, lag=1, argvar=list(type="lin",
  cen=FALSE), arglag=list(type="strata"))
cb3.temp <- crossbasis(chicagoNMMAPS$temp, lag=30, argvar=list(type="bs",
  df=5,degree=2,cen=21), arglag=list(df=5))
summary(cb3.pm)
summary(cb3.temp)

# RUN THE MODEL AND GET THE PREDICTION FOR TEMPERATURE
library(splines)
model3 <- glm(death ~  cb3.pm + cb3.temp + ns(time, 7*14) + dow,
  family=quasipoisson(), chicagoNMMAPS)
pred3.temp <- crosspred(cb3.temp, model3, by=1)

# 3-D AND CONTOUR PLOTS
plot(pred3.temp, xlab="Temperature", col="red", zlab="RR", shade=0.6,
  main="3D graph of temperature effect")
plot(pred3.temp, "contour", xlab="Temperature", key.title=title("RR"),
	plot.title=title("Contour plot",xlab="Temperature",ylab="Lag"))

# MULTIPLE SLICES
plot(pred3.temp, "slices", var=-20, ci="n", col=1, ylim=c(0.95,1.15), lwd=1.5,
	main="Lag-specific effects at different temperature, ref. 21C")
for(i in 1:3) lines(pred3.temp, "slices", var=c(0,27,33)[i], col=i+1, lwd=1.5)
legend("topright",paste("Temperature =",c(-20,0,27,33)), col=1:4, lwd=1.5)
plot(pred3.temp, "slices", var=c(-20,0,27,33), lag=c(0,5,15,28), col=4,
	ci.arg=list(density=40,col=grey(0.7)))

### See the vignette 'dlnmOverview' for a detailed explanation of this example
}

\keyword{hplot}
\keyword{aplot}
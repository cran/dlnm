\name{plot.crosspred}
\alias{plot.crosspred}
\alias{lines.crosspred}
\alias{points.crosspred}
\alias{crossplot}

\title{ Plot predicted effects for a DLNM }

\description{
High and low-level method functions for graphs (3d, contour, slices and overall effects) of predicted effects from distributed lag non-linear models (DLNM).
}

\usage{
\method{plot}{crosspred}(x, ptype="3d", var=NULL, lag=NULL, ci="area", ci.arg, 
	ci.level=x$ci.level, cumul=FALSE, exp=NULL, ...)

\method{lines}{crosspred}(x, ptype="overall", var=NULL, lag=NULL, ci="n", ci.arg, 
	ci.level=x$ci.level, cumul=FALSE, exp=NULL, ...)

\method{points}{crosspred}(x, ptype="overall", var=NULL, lag=NULL, ci="n", ci.arg, 
	ci.level=x$ci.level, cumul=FALSE, exp=NULL, ...)
}

\arguments{
  \item{x }{ an object of class \code{"crosspred"}.}
  \item{ptype }{ type of plot. See Details below.}
   \item{var, lag }{ vectors (for \code{plot}) or numeric scalars (for \code{lines}-\code{points}) of predictor values and lags for which specific effects must be plotted. Used only if \code{type="slices"}.}
  \item{ci }{ type of confidence intervals representation, one of \code{"area"}, \code{"bars"}, \code{"lines"} or \code{"n"}.}
  \item{ci.arg }{ list of arguments to be passed to low-level plotting functions to draw the confidence intervals. See Details.}
  \item{ci.level }{ confidence level for the computation of confidence intervals.}
  \item{cumul }{ logical. If \code{TRUE}, cumulative effects along lags are reported. Used only if \code{type="slices"}. See Details.}
   \item{exp }{ logical. It forces the choice about the exponentiation of effects. See Details.}
  \item{\dots }{ optional graphical arguments. See Details.}
}

\details{
Different plots can be obtained choosing the following values for the argument \code{ptype}:

\bold{\code{"3d"}}: a 3-D plot. Additional graphical arguments can be included, such as \code{theta}-\code{phi} (perspective), \code{border}-\code{shade} (surface), \code{xlab}-\code{ylab}-\code{zlab} (axis labelling) or \code{col}. See \code{\link[graphics]{persp}} for additional information.

\bold{\code{"contour"}}: a contour/level plot.  Additional graphical arguments can be included, such as \code{plot.title}-\code{plot.axes}-\code{key.title} for titles and axis and key labelling. Arguments \code{x}-\code{y}-\code{z} and \code{col}-\code{level} are automatically set and cannot be specified by the user. See \code{\link[graphics]{filled.contour}} for additional information.

\bold{\code{"overall"}}: a plot of the overall effects (summing up all the single lag contributions). See \code{\link[graphics]{plot.default}}, \code{\link[graphics]{lines}} and \code{\link[graphics]{points}} for information on additional graphical arguments.

\bold{\code{"slices"}}: a (optionally multiple) plot of effects along lags at specific predictor values and effects along the predictor at specific lags, chosen by \code{var} and \code{lag}, respectively. See \code{\link[graphics]{plot.default}}, \code{\link[graphics]{lines}} and \code{\link[graphics]{points}} for information on additional graphical arguments.

The method function \code{plot} calls the high-level functions listed above for each \code{ptype}, while \code{lines}-\code{points} add lines or points for \code{ptype} equal to \code{"overall"} or \code{"slices"}. These methods allow a great flexibility in the choice of graphical parameters, specified through arguments of the original plotting functions. Some arguments, if not specified, are set to different default values than the original functions.

Similarly, confidence intervals are plotted for \code{ptype} equal to \code{"overall"} or \code{"slices"} through \code{ci.arg}. This is a list of arguments to be passed to low-level plotting functions, determined by \code{ci}: \code{\link[graphics]{polygon}} for \code{"area"}, \code{\link[graphics]{segments}} for \code{"bars"} and \code{\link[graphics]{lines}} for \code{"lines"}. See the original functions for a complete list of the arguments. This option offers flexibility in the choice of confidence intervals display. As above some unspecified arguments are set to different default values.

For \code{ptype="slices"}, up to 4 plots for each dimension of predictor and lags are allowed in \code{plot}, while for \code{lines}-\code{points} a single plot in one of the two dimension must be chosen. Cumulative effects along lags are reported if \code{cumul=TRUE}: in this case, the same option must have been set to obtain the prediction saved in the \code{crosspred} object (see \code{\link{crosspred}}). 

For a detailed illustration of the use of the functions, see:

\code{vignette("dlnmOverview")}
}

\references{
Gasparrini A., Armstrong, B.,Kenward M. G. Distributed lag non-linear models. \emph{Statistics in Medicine}. 2010; \bold{29}(21):2224-2234.
}

\author{Antonio Gasparrini, \email{antonio.gasparrini@lshtm.ac.uk}}

\note{
These methods for class \code{"crosspred"} have replaced the old function \code{crossplot} since version 1.3.0. The old function has been kept in the namespace of the package, but its use is discouraged.

All the effects are reported versus a reference value. For continuous functions, this is specified by the centering point defined in the \code{crossbasis} object (see \code{\link{crossbasis}}). Exponentiated effects are automatically returned if \code{x} has \code{model.link} equal to \code{log} or \code{logit}, or forced with \code{exp=TRUE}.
}

\section{Warnings}{
The values in \code{var} must match those specified in the object \code{crosspred} (see \code{\link{crosspred}}), while the values in \code{lag} must be included in the lag period specified by \code{\link{crossbasis}}.
}

\seealso{ \code{\link{crossbasis}}, \code{\link{crosspred}}}

\examples{
# Example 3. See crossbasis and crosspred for other examples

### complex DLNM
### space of predictor: 5df quadratic spline for temperature
### space of predictor: linear effect for PM10
### lag function: 5df natural cubic spline for temperature up to lag30
### lag function: single strata at lag 0-1 for PM10

# CREATE THE CROSS-BASIS FOR EACH PREDICTOR AND CHECK WITH SUMMARY
basis.pm <- crossbasis(chicagoNMMAPS$pm10,vartype="lin", lagtype="strata",
	cen=FALSE, maxlag=1)
basis.temp <- crossbasis(chicagoNMMAPS$temp, vartype="bs", vardf=5, vardegree=2,
	lagdf=5, cenvalue=21, maxlag=30)
summary(basis.pm)
summary(basis.temp)

# RUN THE MODEL AND GET THE PREDICTION FOR TEMPERATURE
library(splines)
model <- glm(death ~  basis.pm + basis.temp + ns(time, 7*14) + dow,
	family=quasipoisson(), chicagoNMMAPS)
pred.temp <- crosspred(basis.temp, model, by=1)

# 3-D AND CONTOUR PLOTS
plot(pred.temp, xlab="Temperature", col="red", zlab="RR", shade=0.6,
	main="3D graph of temperature effect")
plot(pred.temp, "contour", xlab="Temperature", key.title=title("RR"),
	plot.title=title("Contour plot",xlab="Temperature",ylab="Lag"))

# MULTIPLE SLICES
plot(pred.temp, "slices", var=-20, ci="n", col=1, ylim=c(0.95,1.15), lwd=1.5,
	main="Lag-specific effects at different temperature, ref. 21C")
for(i in 1:3) lines(pred.temp, "slices", var=c(0,27,33)[i], col=i+1, lwd=1.5)
legend("topright",paste("Temperature =",c(-20,0,27,33)), col=1:4, lwd=1.5)
plot(pred.temp, "slices", var=c(-20,0,27,33), lag=c(0,5,15,28), col=4,
	ci.arg=list(density=40,col=grey(0.7)))

### See the vignette 'dlnmOverview' for a detailed explanation of this example
}

\keyword{hplot}
\keyword{aplot}
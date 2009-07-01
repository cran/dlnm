\name{crossbasis}
\alias{crossbasis}
\alias{summary.crossbasis}

\title{ Generate a cross-basis matrix for a DLNM }

\description{
Generate the basis functions for the two spaces of predictor and lags, choosing among a set of possible bases. Then, these functions are combined in order to create the related cross-basis matrix, which can be included in a model formula to fit a distributed lag non-linear model (DLNM).
}

\usage{
crossbasis(var, vartype="ns", vardf=1, varknots=NULL, varbound=range(var),
	varint=FALSE, cen=TRUE, cenvalue=mean(var), maxlag=0, lagtype="ns",
	lagdf=1, lagknots=NULL, lagbound=c(-1,maxlag), lagint=TRUE)

\method{summary}{crossbasis}(object, ...)
}

\arguments{
The arguments below define two set of basis functions. The first one is applied to \code{var}, in order to describe the relationship in the space of the predictor. The second one is applied to a new vector \code{0:maxlag}, in order to describe the relationship in the space of lags. Many arguments refer to the specific basis for each space (with stub \code{var-} or \code{lag-}). Then, the two set of basis functions are combined in order to create the related cross-basis functions.
  \item{var }{ the predictor variable, defined as a numeric vector of ordered observations.}
  \item{vartype, lagtype }{ type of basis. See Details below for the list of possible choices.}
  \item{vardf, lagdf }{ dimension of the basis (number of degrees of freedom spent to specify the relationship in each space). They depend on \code{varknots} and \code{lagknots} if provided.}
  \item{varknots, lagknots }{ knots location for the basis. They specify the position of the knots for \code{"ns"}, the cut-off points for \code{"strata"} (defining right-open intervals) and the threshold(s) for \code{"lthr"}, \code{"hthr"} and \code{"thr"}. They must be set within the range of \code{var} and \code{0:maxlag}, respectively.}
  \item{varbound, lagbound }{ boundary knots. Used only if \code{type="ns"}.}
  \item{varint, lagint }{ logical. If \code{TRUE}, an 'intercept' is included in the basis. The default values should not be changed: see Warnings below.}
  \item{cen }{ logical. If \code{TRUE}, the basis functions for the space of predictor are centered. See Note below.}
  \item{cenvalue }{ centering value.}
  \item{maxlag }{ a positive value defining the maximum lag.}
  \item{object }{ a object of class \code{"crossbasis"}.}
  \item{\dots }{ additional arguments to be passed to \code{summary}.}
}

\details{
The value in \code{type} defines the basis for each space (predictor and lags). It must  be one of:

\bold{\code{"ns"}}: natural cubic splines, specified by \code{knots} (internal knots) and \code{bound} (boundary knots). See the function \code{\link[splines::ns]{ns}} for additional information. An intercept is included if \code{int=T}. The transformed variables can be centered at \code{cenvalue}.

\bold{\code{"strata"}}: strata variables (dummy parameterization) determined by internal cut-off values specified in \code{knots}, which represent the lower boundaries for the right-open intervals. A dummy variable for the reference stratum (the first one by default) is included if \code{int=T}, generating a full rank basis. Never centered.

\bold{\code{"poly"}}: polynomial with degree equal to \code{df-int}. An intercept, corresponding to a vector of 1's (the power 0 of the polynomial) is included if \code{int=T}. The transformed variables can be centered at \code{cenvalue}.

\bold{\code{"integer"}}: strata variables (dummy parameterization) for each integer values, expressly created to specify an unconstrained function in the space of lags. \code{df} is set automatically to the number of integer values. A dummy variable for the reference stratum (the first one by default) is included if \code{int=T}, generating a full rank basis. Never centered.

\bold{\code{"hthr","lthr"}}: high and low threshold parameterization (linear relationship above or below the threshold, respectively, flat otherwise). The threshold is chosen by \code{knots} (if more than one are provided, the last and the first one is chosen, respectively). \code{df} is automatically set to 1. No intercept. Never centered.

\bold{\code{"thr"}}: double threshold parameterization (2 independent linear relationships above the second and below the first threshold, flat between them). The thresholds are chosen by \code{knots}. If only one is provided, the threshold is unique (V-model). If more than 2 are provided, the first 2 are chosen. \code{df} is automatically set to 2. No intercept. Never centered.

\bold{\code{"lin"}}: linear relationship (untransformed apart from optional centering). \code{df} is automatically set to 1. No intercept. It can be centered at \code{cenvalue}. 

The values in \code{knots}, if provided, automatically determine the value of \code{"df"} (equal to \code{length(knots)+int} for \code{"strata"} and \code{length(knots)+int+1} for \code{"ns"} and \code{"poly"}). If not provided, \code{varknots} are placed at equally spaced quantiles, and \code{lagknots} at equally spaced values on the log scale of lags. Some arguments can be automatically changed for not sensible combinations, or set to \code{NULL} if not required.

A vignette with a detailed overview of all the functionalities of the package will be added in the next releases.
}

\value{
A matrix object of class \code{"crossbasis"} which can be included in a model formula in order to fit a DLNM. It contains the attributes \code{crossdf} (global number of degrees of freedom) and \code{range} (range of the original vector of observations). Additional attributes are returned that correspond to the arguments to \code{crossbasis}, and explicitly give \code{type}, \code{df}, \code{knots}, \code{bound}, \code{cen}, \code{cenvalue} and \code{maxlag} related to the corresponding basis ( with stub \code{var-} or \code{lag-}) for use of \code{\link{crosspred}}. The function \code{\link{summary.crossbasis}} returns a summary of the cross-basis matrix and the related attributes.
}

\references{ 
Armstrong, B. Models for the relationship between ambient temperature and daily mortality. \emph{Epidemiology}. 2006, \bold{17}(6):624-31.
}

\author{Antonio Gasparrini, \email{antonio.gasparrini@lshtm.ac.uk}}

\note{
The values in \code{var} are expected to be equally spaced (with that space defining a lag unit) and ordered in time. \code{NA} values are allowed.

The name of the crossbasis object will be used by \code{\link{crosspred}} in order to extract the related estimated parameters. This name must not match the names of other predictors in the model formula. In addition, if more than one variable is transformed by cross-basis functions in the same model, different names must be specified. 

For continuous functions specified with \code{vartype} equal to \code{"ns"}, \code{"poly"} or \code{"lin"}, the reference for the effects predicted by \code{\link{crosspred}} is set at \code{cenvalue}. For the other choices, the reference is automatic:  for \code{vartype} equal to \code{"strata"} and \code{"integer"}, the reference is the first interval, while for \code{vartype} equal to \code{"hthr"}, \code{"lthr"} and \code{"thr"}, the reference is the region of null effect below, above or between the threshold(s), respectively.
}

\section{Warnings}{
It is strongly recommended to avoid the inclusion of an intercept in the basis for \code{var}, otherwise the presence of the additional intercept in the model used to fit the data will cause some of the cross-basis variables to be excluded. Conversely, an intercept should always be included in the basis for the space of lags when \code{lagtype} is equal to \code{"ns"}, \code{"strata"} or \code{"poly"}.
}

\seealso{ \code{\link{crosspred}}, \code{\link{crossplot}}}

\examples{
# Example 1. See crosspred and crossplot for other examples

# load and prepare the dataset
initDB()
data <- readCity("chic", collapseAge=TRUE)
data$temp <- (data$tmpd-32)*5/9
data$pm10 <- with(data, pm10mean+pm10trend)

### simple DLM for the effect of PM10 on mortality up to 15 days of lag
### space of predictor: linear effect for PM10
### space of predictor: 5df spline for temperature
### lag function: 4th degree polynomial for PM10
### lag function: strata at lag 0 and 1-3 for temperature

basis.pm <- crossbasis(data$pm10, vartype="lin", lagtype="poly",
	lagdf=5,cen=FALSE,maxlag=15)
basis.temp <- crossbasis(data$temp, vardf=5, lagtype="strata",
	lagknots=1,	cenvalue=21, maxlag=3)
summary(basis.pm)
summary(basis.temp)
model <- glm(death ~  basis.pm + basis.temp, family=quasipoisson(), data)
pred.pm <- crosspred(basis.pm, model, at=0:20)

crossplot(pred.pm,"slices",var=10,
	title="Effect of a 10-unit increase in PM10 along lags")
# overall effect for a 10-unit increase in PM over 15 days of lag, with CI
pred.pm$allRRfit["10"]
cbind(pred.pm$allRRlow, pred.pm$allRRhigh)["10",]
crossplot(pred.pm, "overall", label="PM10",
	title="Overall effect of PM10 over 15 days of lag")
}

\keyword{smooth}
\keyword{ts}


\name{dlnm-internal}
\alias{mkbasis}
\alias{mklagbasis}

\title{ Internal functions for package dlnm }

\description{
Generate the basis functions for the space of predictor and lags. These functions are built for internal use only.
}

\usage{
mkbasis(var, type="ns", df=1, degree=1, knots=NULL, bound=range(var),
	int=FALSE, cen=TRUE, cenvalue=mean(var))

mklagbasis(maxlag=0, type="ns", df=1, degree=1, knots=NULL,
	bound=c(0, maxlag), int=TRUE)
}

\arguments{
  \item{var }{ a numeric vector of ordered observations.}
  \item{type}{ type of basis.}
  \item{df}{ dimension of the basis. They depend on \code{knots} or \code{degree}  if provided.}
  \item{degree}{ degree of polynomial. Used only for \code{type} equal to \code{"bs"} of \code{"poly"}.}
  \item{knots}{ knots location for the basis.}
  \item{bound}{  boundary knots. Used only for \code{type} equal to \code{"ns"} of \code{"bs"}.}
  \item{int}{ logical. If \code{TRUE}, an intercept is included in the basis.}
  \item{cen }{ logical. If \code{TRUE}, the basis functions are centered.}
  \item{cenvalue }{ centering value.}
  \item{maxlag }{ maximum lag.}
}

\details{
These two internal functions are called by \code{\link{crossbasis}} to build the basis matrices for the two dimensions of predictor and lags. Even if they are not expected to be directly run by the users, they are included in the namespace of the package and therefore made accessible, with the intention to keep the process more transparent and give the opportunity to change or improve them.

\code{mkbasis} applies a set basis functions on the vector of the predictor specified in \code{var}. \code{mklagbasis} calls \code{mkbasis} to build the basis matrix for the space of lags. Basically, it creates a new vector \code{0:maxlag} (a vector of integers from 0 to the maximum lag allowed) and then apply the related basis functions. The functions for the 2 spaces are defined by the arguments above. See \code{\link{crossbasis}} for additional information.
}

\value{
  \item{basis }{ matrix of basis functions}
Additional values are returned that correspond to the arguments above, and explicitly give \code{type}, \code{df}, \code{ddegree}, \code{knots}, \code{bound}, \code{int}, \code{varcen}, \code{cenvalue} and \code{maxlag} related to the corresponding basis.
}

\references{ 
Armstrong, B. Models for the relationship between ambient temperature and daily mortality. \emph{Epidemiology}. 2006, \bold{17}(6):624-31.
}

\author{Antonio Gasparrini, \email{antonio.gasparrini@lshtm.ac.uk}}

\examples{
# natural cubic spline with a knot at 3
basis.var <- mkbasis(1:5, knots=3)
basis.var

# quadratic spline selected by df, automatic knots location
mkbasis(1:5, type="bs", df=4, degree=2)

# linear centered at 4
mkbasis(1:5, type="lin", cenvalue=4)

# polynomial with degree 3, with intercept
mklagbasis(maxlag=5, type="poly", degree=3)

# integer
mklagbasis(maxlag=5, type="integer")

# threshold-type: double and single threshold or piecewise
mkbasis(1:5, type="dthr", knots=c(2,3))
mkbasis(1:5, type="hthr", knots=3)
mkbasis(1:5, type="hthr", knots=c(2,3))

# the intercept: strata defined by 2 cut-off points
mklagbasis(maxlag=10, type="strata", knots=c(4,7))
mklagbasis(maxlag=10, type="strata", knots=c(4,7), int=FALSE)

# centering: polynomial
mkbasis(0:10, type="poly", degree=3)
mkbasis(0:10, type="poly", degree=3, cen=FALSE)

### See the vignette 'dlnmOverview' for a detailed explanation of this example
}

\keyword{internal}

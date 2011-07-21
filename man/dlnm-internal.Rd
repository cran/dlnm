\name{dlnm-internal}
\alias{mkbasis}
\alias{mklagbasis}

\title{Internal functions for package dlnm}

\description{
Generate different types of basis matrices for the space of predictor and lags. These functions are called internally by other functions and are not meant to be directly run by the users.
}

\usage{
mkbasis(x, type="ns", df=1, degree=1, knots=NULL,
	bound=range(x),int=FALSE, cen=TRUE, cenvalue=mean(x))

mklagbasis(maxlag=0, type="ns", df=1, degree=1, knots=NULL,
	bound=c(0, maxlag), int=TRUE)
}

\arguments{
  \item{x }{ the predictor variable. Missing values are allowed.}
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
\code{mkbasis} applies a set basis functions to the predictor vector \code{x}. \code{mklagbasis} calls \code{mkbasis} to build the basis matrix for the space of lags. Basically, it creates a new vector \code{0:maxlag} (a vector of integers from 0 to the maximum lag allowed) and then applies the related basis functions. The basis functions for each dimension are defined by the arguments above. See \code{\link{crossbasis}} for additional information on the specific bases.

These two functions are called by \code{\link{crossbasis}} to build the basis matrices for the two dimensions of predictor and lags in DLNMs. . Even if they are not expected to be directly run by the users, they are included in the namespace of the package and therefore made accessible, with the intention to keep the process more transparent and give the opportunity to change or improve them.

For a detailed illustration of the use of the functions, see:

\code{vignette("dlnmOverview")}
}

\value{
  \item{basis }{ matrix of basis functions}
Additional values are returned that correspond to the arguments above, and explicitly give \code{type}, \code{df}, \code{degree}, \code{knots}, \code{bound}, \code{int}, \code{varcen}, \code{cenvalue} and \code{maxlag} related to the corresponding basis.
}

\references{
Gasparrini A. Distributed lag linear and non-linear models in R: the package dlnm. \emph{Journal of Statistical Software}. 2011; \bold{43}(8):1-20. URL http://www.jstatsoft.org/v43/i08/.

Gasparrini A., Armstrong, B.,Kenward M. G. Distributed lag non-linear models. \emph{Statistics in Medicine}. 2010; \bold{29}(21):2224-2234.
}

\author{Antonio Gasparrini, \email{antonio.gasparrini@lshtm.ac.uk}}

\examples{
# NATURAL CUBIC SPLINE WITH A KNOT AT 3
basis.var <- mkbasis(1:5, knots=3)
basis.var

# QUADRATIC SPLINE SELECTED BY DF, AUTOMATIC KNOT LOCATION
mkbasis(1:5, type="bs", df=4, degree=2)

# LINEAR, CENTERED AT 4
mkbasis(1:5, type="lin", cenvalue=4)

# POLYNOMIAL WITH DEGREE 3, WITH INTERCEPT
mklagbasis(maxlag=5, type="poly", degree=3)

# INTEGER
mklagbasis(maxlag=5, type="integer")

# THRESHOLD-TYPE: DOUBLE AND SINGLE THRESHOLD OR PIECEWISE
mkbasis(1:5, type="dthr", knots=c(2,3))
mkbasis(1:5, type="hthr", knots=3)
mkbasis(1:5, type="hthr", knots=c(2,3))

# THE 'INTERCEPT': STRATA DEFINED BY 2 CUT-OFF POINTS
mklagbasis(maxlag=10, type="strata", knots=c(4,7))
mklagbasis(maxlag=10, type="strata", knots=c(4,7), int=FALSE)

# CENTERING: POLYNOMIAL
mkbasis(0:10, type="poly", degree=3)
mkbasis(0:10, type="poly", degree=3, cen=FALSE)

### See the vignette 'dlnmOverview' for a detailed explanation of this example
}

\keyword{internal}

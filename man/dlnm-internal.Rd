\name{dlnm-internal}
\alias{mkbasis}
\alias{mklagbasis}

\title{ Internal functions for package dlnm }

\description{
Generate the basis functions for the space of predictor and lags. These functions are built for internal use only.
}

\usage{
mkbasis(var, type="ns", df=1, knots=NULL, bound = range(var),
	int=FALSE, cen=TRUE, cenvalue=mean(var))

mklagbasis(maxlag=0, type="ns", df=1, knots=NULL,
	bound=c(0, maxlag), int=TRUE)
}

\arguments{
  \item{var }{ a numeric vector of ordered observations.}
  \item{type}{ type of basis.}
  \item{df}{ dimension of the basis. They depend on \code{varknots} and \code{lagknots} if provided.}
  \item{knots}{ knots location for the basis.}
  \item{bound}{  boundary knots. Used only if \code{type="ns"}.}
  \item{int}{ logical. If \code{TRUE}, an intercept is included in the basis.}
  \item{cen }{ logical. If \code{TRUE}, the basis functions are centered.}
  \item{cenvalue }{ centering value.}
  \item{maxlag }{ maximum lag.}
}

\details{
These functions are called by \code{crossbasis} in order to build the related cross-basis functions. See \code{\link{crossbasis}} for further details.
}

\value{
  \item{basis }{ matrix of basis functions}
Additional values are returned that correspond to the arguments above, and explicitly give \code{type}, \code{df}, \code{knots}, \code{bound}, \code{int}, \code{varcen}, \code{cenvalue} and \code{maxlag} related to the corresponding basis.
}

\references{ 
Armstrong, B. Models for the relationship between ambient temperature and daily mortality. \emph{Epidemiology}. 2006, \bold{17}(6):624-31.
}

\author{Antonio Gasparrini, \email{antonio.gasparrini@lshtm.ac.uk}}

\examples{
mkbasis(1:5, knots=3)
mkbasis(1:5, type="thr", knots=3)
mklagbasis(maxlag=5, type="poly", df=3)

# the intercept
mklagbasis(maxlag=10, type="strata", knots=c(4,7))
mklagbasis(maxlag=10, type="strata", knots=c(4,7), int=FALSE)

# centering
mkbasis(0:10, type="poly", df=3)
mkbasis(0:10, type="poly", df=3, cen=FALSE)

### See the vignette 'dlnmOverview' for a detailed explanation of this example
}

\keyword{internal}

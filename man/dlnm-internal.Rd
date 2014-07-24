\name{dlnm-internal}
\alias{mklag}
\alias{fci}
\alias{seqlag}
\alias{Lag2}
\alias{getcoef}
\alias{getvcov}
\alias{getlink}


\title{ Internal Functions for Package dlnm }

\description{
These functions are called internally by other functions and are not meant to be directly run by the users.
}

\usage{
mklag(lag)

fci(ci, x, high, low, ci.arg, plot.arg, noeff=NULL)

seqlag(lag, by=1)

Lag2(v, k, group=NULL)

getcoef(model, class, cond)

getvcov(model, class, cond)

getlink(model, class)
}

\arguments{
  \item{lag }{ either an integer scalar or vector of length 2, defining the the maximum lag or the lag range, respectively.}
  \item{ci }{ type of confidence intervals representation.}
  \item{x }{ the coordinates of the x axis.}
  \item{high, low }{ the coordinates of the y axis for the interval.}
  \item{ci.arg }{ list of arguments to draw the confidence intervals.}
  \item{plot.arg }{ list of arguments of the main plot.}
  \item{noeff }{ reference value of the null effect.}
  \item{lag }{ a numeric vector of length 2 specifying the range of sequence.}
  \item{by }{ increment of the sequence.}
  \item{v }{ a numeric vector for which lagged occurrences need to be computed.}
  \item{k }{ the lag values.}
  \item{group }{ a factor or a list of factors defining groups of observations.}
  \item{model }{ a regression model object.}
  \item{class }{ a regression model class.}
  \item{cond }{ a condition, usually defined through regular expressions (see \code{\link{regex}}).}
}

\details{
The function \code{mklag} is used internally for checking/defining the vector of length 2.

The function \code{fci} provides different options for representing confidence intervals, and it is called in plotting functions.

The function \code{seqlag} creates the vector with the sequence of lags given the range provided.

The function \code{Lag2} creates a matrix of lagged occurrences.

The functions \code{getcoef}, \code{getvcov} and \code{getlink} extract the part of coefficients and (co)variance matrix with names matching a specific condition, and the model link, depending on the class of the model
}

\author{Antonio Gasparrini, \email{antonio.gasparrini@lshtm.ac.uk}}

\note{
These and other undocumented internal functions are not exported in the namespace of the package \pkg{dlnm}. The user can visualize them through the use of the triple colon operator '\code{:::}' or through the function \code{\link{getAnywhere}} (see the examples below).  
}

\seealso{
See \code{\link{dlnm-package}} for an introduction to the package and for links to package vignettes providing more detailed information.
}

\examples{
### display internal functions
dlnm:::getcoef
getAnywhere(getcoef)

### display other undocumented functions
dlnm:::fci
getAnywhere(fci)
}

\keyword{internal}

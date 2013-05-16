\name{dlnm-internal}
\alias{.fci}
\alias{.seq}

\title{ Internal Functions for Package dlnm }

\description{
These functions are called internally by other functions and are not meant to be directly run by the users.
}

\usage{
.fci(ci, x, high, low, ci.arg, plot.arg, noeff=NULL)

.seq(lag, by=1)
}

\arguments{
  \item{ci }{ type of confidence intervals representation.}
  \item{x }{ the coordinates of the x axis.}
  \item{high, low }{ the coordinates of the y axis for the interval.}
  \item{ci.arg }{ list of arguments to draw the confidence intervals.}
  \item{plot.arg }{ list of arguments of the main plot.}
  \item{noeff }{ reference value of the null effect.}
  \item{lag }{ a numeric vector of length 2 specifying the range of sequence.}
  \item{by }{ increment of the sequence.}
}

\details{
The function \code{.fci} is called internally within plotting functions of the package \pkg{dlnm} and provides different options for representing confidence intervals.

The function \code{.seq} is called internally to create the vector with the sequence of lag given the range provided.
}

\author{Antonio Gasparrini, \email{antonio.gasparrini@lshtm.ac.uk}}

\seealso{
See \code{\link{dlnm-package}} for an overview of the package and type \code{'vignette(dlnmOverview)'} for a detailed description.
}

\keyword{internal}

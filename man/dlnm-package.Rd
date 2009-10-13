\name{dlnm-package}
\alias{dlnm-package}
\alias{dlnm}
\docType{package}

\title{
Distributed Lag Non-linear Models (DLNM)
}

\description{
The package \pkg{dlnm} contains functions to specify basis and cross-basis matrices in order to run distributed lag models (DLM) and their non-linear extension (DLNM), then to predict and plot the results for a fitted model.
}

\details{
Distributed lag non-linear models (DLNM) represent a modelling framework to describe simultaneously non-linear and delayed dependencies in time-series data. This methodology is based on the definition of a \emph{cross-basis}, a bi-dimensional space of functions specifying the dependency along the space of the predictor and along lags. The cross-basis functions are built combining the basis functions for the two dimensions, chosen among a set of possible bases. This family includes simple distributed lag models (DLM) as a special case. The DLNM framework is developed for time series data: the potential applications to other study designs are currently under evaluation and will be explored in next releases.

Given a series of observations ordered and equally spaced in time, \code{\link{crossbasis}} creates a matrix object of class \code{"crossbasis"} containing the transformed variables to be included in the model formula. The estimation is obtained by default model commands, usually \code{\link[stats]{glm}}. After the model fitting, \code{\link{crosspred}} predicts the results for a set of suitable values of the original predictor and stores them in a \code{"crosspred"} object. Finally, \code{\link{crossplot}} offers a set of choices to plot the results.

Use \code{citation("dlnm")} to cite this package.

A list of changes included in the current and previous versions can be found typing \code{file.show(system.file("ChangeLog", package = "dlnm"))}.

For further information on DLNM, see the references below. For a detailed description of the capabilities of the package, refer to:

\code{vignette("dlnmOverview")}
}

\author{
Antonio Gasparrini and  Ben Armstrong

Maintainer: Antonio Gasparrini, \email{antonio.gasparrini@lshtm.ac.uk}
}

\references{
Armstrong, B. Models for the relationship between ambient temperature and daily mortality. \emph{Epidemiology}. 2006, \bold{17}(6):624-31.
}

\keyword{package}

\seealso{
\code{\link{crossbasis}} to create the basis and cross-basis matrices. \code{\link{crosspred}} to predict the effects after model fitting. \code{\link{crossplot}} to plot several type of graphs.
}


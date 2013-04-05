\name{dlnm-package}
\alias{dlnm-package}
\alias{dlnm}
\docType{package}

\title{
Distributed Lag Non-linear Models (DLNM)
}

\description{
The package \pkg{dlnm} contains functions to specify and run distributed lag linear and non-linear models. These functions are used to specify basis and cross-basis matrices among a set of available options, then to predict and plot the results for a fitted model.
}

\details{
Distributed lag non-linear models (DLNM) represent a modelling framework to describe simultaneously non-linear and delayed dependencies in time-series data. This methodology is based on the definition of a \emph{cross-basis}, a bi-dimensional space of functions specifying the dependency along the space of the predictor and along lags. The cross-basis functions are built combining the basis functions for the two dimensions, chosen among a set of possible bases. This family includes simple distributed lag models (DLM) as a special case. A thorough methodological overview is given in the references below.

Given a series of observations ordered and equally spaced in time, \code{\link{crossbasis}} creates two set of basis functions to define the relationship in the two dimensions of predictor and lags, by calling the function \code{\link{onebasis}}, which provides a list of possible options (splines, strata, and threshold functions, among others). The two matrices are then combined in a matrix object of class \code{"crossbasis"}, containing the transformed variables to be included in the model formula. After the model fitting, \code{\link{crosspred}} generate predictions for a set of suitable values of the original predictor and stores them in a \code{"crosspred"} object. The fit of a DLNM can be reduced and re-expressed as the chosen function of one of the two dimensions through the function \code{\link{crossreduce}}. It returns a \code{"crossreduce"} object storing the new parameters and predictions.

Method functions are available for both objects \code{"crosspred"} and \code{"crossreduce"}. The plotting functions \code{\link[=plot.crosspred]{plot}}, \code{\link[=lines.crosspred]{lines}} and \code{\link[=points.crosspred]{points}}, offer a set of choices to plot the results, while \code{\link[=coef.crosspred]{coef}} and \code{\link[=vcov.crosspred]{vcov}} return the coefficients and associated (co)variance matrix for a (optionally reduced) DLNM.

The DLNM framework is developed for time series data: the data are assumed to be an equally-spaced, complete and ordered series of observations, where the interval defines the lag unit. The estimation can be carried out with the default regression functions, such as \code{\link{lm}}, \code{\link{glm}} or \code{gam} (package \pkg{mgcv}). The function \code{\link{crosspred}} extracts the relevant parameters directly from the model object if methods for \code{\link{coef}} and \code{\link{vcov}} are available for that model class. Alternatively, the user must extract the parameters manually.

Potential applications to other study desing and data structures are currently under evaluation. Preliminary results for matched case-control and cohort data are promising, and the functions included in the package already work with \code{clogit} and \code{coxph} (package \pkg{survival}).

Additonal information on the package \pkg{dlnm} are available in the vignette included in the installation. This document offers a detailed description of the capabilities of the package, and some examples of application to real data, with an extensive illustration of the use of the functions. The vignette is available by typing:

\code{vignette("dlnmOverview")}

Use \code{citation("dlnm")} to cite this package.

A list of changes included in the current and previous versions can be found by typing:

\code{file.show(system.file("ChangeLog",package="dlnm"))}

General information on the development and applications of the DLNM modelling framework, together with an updated version of the \R{} scripts for running the examples in published papers, can be found at \bold{\href{http://www.ag-myresearch.com}{www.ag-myresearch.com}}.
}

\author{
Antonio Gasparrini and  Ben Armstrong

Maintainer: Antonio Gasparrini, \email{antonio.gasparrini@lshtm.ac.uk}
}

\references{
Gasparrini A. Distributed lag linear and non-linear models in R: the package dlnm. \emph{Journal of Statistical Software}. 2011; \bold{43}(8):1-20. [freely available \href{http://www.ag-myresearch.com/jss2011.html}{here}].
  
Gasparrini A., Armstrong, B.,Kenward M. G. Distributed lag non-linear models. \emph{Statistics in Medicine}. 2010; \bold{29}(21):2224-2234. [freely available \href{http://www.ag-myresearch.com/statmed2010}{here}]

Gasparrini A., Armstrong, B., Kenward M. G. Reducing and meta-analyzing estimates from distributed lag non-linear models.\emph{BMC Medical Research Methodology}. 2013; \bold{13}(1):1. [freely available \href{http://www.ag-myresearch.com/bmcmrm2013}{here}].

Armstrong, B. Models for the relationship between ambient temperature and daily mortality. \emph{Epidemiology}. 2006, \bold{17}(6):624-31. [available \href{http://www.ncbi.nlm.nih.gov/pubmed/17028505}{here}]
}

\keyword{package}

\seealso{
\code{\link{onebasis}} to generate simple basis matrices.
\code{\link{crossbasis}} to generate cross-basis matrices. \code{\link{crosspred}} to obtain predictions after model fitting. \code{\link{crossreduce}} to reduce the fit to one dimension.
\code{\link{plot.crosspred}} and \code{\link{plot.crossreduce}} to plot several type of graphs.

Type \code{'vignette(dlnmOverview)'} for a detailed description.
}


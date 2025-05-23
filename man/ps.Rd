\name{ps}
\alias{ps}

\title{ Generate a Basis Matrix for P-Splines }

\description{
Generate the basis matrix for P-splines, namely a B-spline basis with difference penalties.
}

\usage{
ps(x, df=10, knots=NULL, degree=3, intercept=FALSE, fx= FALSE, S=NULL, diff=2)
}

\arguments{
  \item{x }{ the predictor variable. Missing values are allowed.}
  \item{df }{ degrees of freedom, basically the dimension of the basis matrix. If supplied in the absence of \code{knots}, it automatically selects \code{df+degree+2-intercept} equally-spaced knots (within and beyond the range of \code{x}). The minimum \code{df} allowed is \code{degree+1-intercept}.}
  \item{knots }{ breakpoints that define the spline. These are generally automatically selected, and not defined by the user. See Details below.}
  \item{degree }{ degree of the piecewise polynomial. Default is 3 for cubic splines.}
  \item{intercept }{ logical. If \code{TRUE}, an intercept is included in the basis matrix. See Details below.}
  \item{fx }{ logical. If \code{TRUE}, it removes the penalization. See Details below.}
  \item{S }{ penalty matrix, usually internally defined if \code{NULL} (default).}
  \item{diff }{ order difference of the penalty.}
}

\details{
The function has a usage similar to \code{\link{bs}} and \code{\link{ns}} in the \pkg{splines} package. It produces B-spline transformations through a call to \code{\link{splineDesign}}, plus a difference matrix to define penalties. The same results are returned by the related \code{\link[mgcv]{smooth.construct.ps.smooth.spec}} in the package \pkg{mgcv}.

The argument \code{knots} defines a vector of knots or (if of length 2) the lower and upper limits between which the splines can be evaluated. However, knots should be usually left automatically selected, and in particular these P-splines only have sense with equally-spaced knots, due to the nature of the penalization. It is important to highlight that, differently from \code{\link{bs}} where \emph{internal} and \emph{boundary} knots are defined, this function adopts a standard B-spline parameterization, including by default \code{2*(degree+1)} knots beyond the range of the variable.

The penalization is defined on the difference of adjacent coefficients during fitting procedure through a penalty matrix \code{S}. The argument \code{diff} selects the order difference (with the default 2 determining a second order difference, and 0 producing a ridge penalty), while setting \code{fx=TRUE} removes the penalization. 

Similarly to \code{\link{bs}} and \code{\link{ns}}, setting \code{intercept=FALSE} (default) determines the exclusion of the first transformed variables, and the corresponding first row and column in \code{S}, thus avoiding identifiability issues during the model fitting. Note how the procedure of imposing identifiability constraints is different from that adopted by \code{\link[mgcv]{smoothCon}} in the package \pkg{mgcv}, where a more complex reparameterization is produced. 
}

\value{
A matrix object of class \code{"ps"}. It contains the attributes \code{df}, \code{knots}, \code{degree}, \code{intercept}, \code{fx}, \code{S}, and \code{diff}, with values that can be different than the arguments provided due to internal reset.
}

\references{
Gasparrini A, Scheipl F, Armstrong B, Kenward MG. A penalized framework for distributed lag non-linear models. \emph{Biometrics}. 2017;\bold{73}(3):938-948. [freely available \href{http://www.ag-myresearch.com/2017_gasparrini_biomet.html}{here}]

Eilers P. H. C. and Marx B. D. Flexible smoothing with B-splines and penalties. \emph{Statistical Science}. 1996;\bold{11}(2):89-121.

Wood S. N. Generalized Additive Models: An Introduction with R. Chapman and Hall/CRC Press, 2006.
}

\author{Antonio Gasparrini <\email{antonio.gasparrini@lshtm.ac.uk}>, adapting code available from functions included in the package \pkg{mgcv} by Simon N. Wood.}

\note{
The function is primarily added here to specify penalized DLMs and DLNMs using the so-called \emph{external} method, \emph{i.e.} by including the penalty matrix in the argument \code{paraPen} of the \code{\link[mgcv]{gam}} regression function in \pkg{mgcv}  (see \code{\link{cbPen}}). However, this approach can be also used to fit standard uni-dimensional P-spline models as an alternative to the use of specific \code{\link[mgcv]{smooth.construct.ps.smooth.spec}}, as it takes advantage of the use of prediction and plotting functions in \pkg{dlnm}.
}

\seealso{
\code{\link{cr}} for penalized cubic regression splines. \code{\link{bs}} and \code{\link{ns}} for B-splines and natural cubic splines, respectively. \code{\link{cbPen}} for defining tensor-type bi-dimensional penalties in DLNMs. The related \code{\link[mgcv]{smooth.construct.ps.smooth.spec}} for P-spline smooths in \pkg{mgcv}. The \code{\link[=smooth.construct.cb.smooth.spec]{cb smooth constructor}} for cross-basis penalized spline smooths.

See \code{\link{dlnm-package}} for an introduction to the package and for links to package vignettes providing more detailed information.
}

\examples{
# to be added soon
}

\keyword{smooth}


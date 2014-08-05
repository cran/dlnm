###
### R routines for the R package dlnm (c) Antonio Gasparrini 2013-2014
#
integer <-
function(x, values, int=FALSE) {
#
################################################################################
#
  nx <- names(x)
  x <- as.vector(x)
#
  # DEFINE LEVELS AND TRANSFORM INTO A FACTOR
  levels <- if(!missing(values)) values else sort(unique(x))
  xfac <- factor(x,levels=levels)
#
  # TRANSFORMATION
  basis <- as.matrix(outer(xfac,levels,"==")+0L)
#
  # IF INTERCEPT IS NOT REQUIRED, DROP THE FIRST COLUMN
  if(ncol(basis)>1L) {
    if(!int) basis <- basis[,-1L,drop=FALSE]
  } else int <- TRUE
#
  # NAMES AND ATTRIBUTES
  dimnames(basis) <- list(nx,seq(ncol(basis)))
  attributes(basis) <- c(attributes(basis),list(values=levels,int=int))
#
  class(basis) <- c("integer","matrix")
#
  return(basis)
}

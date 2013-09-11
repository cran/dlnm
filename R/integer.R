###
### R routines for the R package dlnm (c) Antonio Gasparrini 2013
#
`integer` <-
function(x, values, int=FALSE) {
#
################################################################################
#
  nx <- names(x)
  x <- as.vector(x)
#
  # TRANSFORM INTO A FACTOR AND DEFINE LEVELS
  x <- round(x,0)
  levels <- if(!missing(values)) values else sort(unique(x))
  x <- factor(x,levels=levels)
#
  # IF INTERCEPT IS NOT REQUIRED, DROP THE FIRST LEVEL
  if(length(levels)>1L) {
    if(!int) levels <- levels[-1L]
  } else int <- TRUE
#
  # TRANSFORMATION
  basis <- as.matrix(outer(x,levels,"==")+0L)
#
  # NAMES AND ATTRIBUTES
  dimnames(basis) <- list(nx,seq(ncol(basis)))
  attributes(basis) <- c(attributes(basis),list(values=levels,int=int))
#
  class(basis) <- c("integer","matrix")
#
  return(basis)
}

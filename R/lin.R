###
### R routines for the R package dlnm (c) Antonio Gasparrini 2013
#
`lin` <-
function(x, int=FALSE) {
#
################################################################################
#
  nx <- names(x)
  x <- as.vector(x)
#
  # TRANSFORMATION
  basis <- as.matrix(x)
  if(int) basis <- cbind(1,basis)
#
  # NAMES AND ATTRIBUTES
  dimnames(basis) <- list(nx,seq(ncol(basis)))
  attributes(basis) <- c(attributes(basis),list(int=int))
#
  class(basis) <- c("lin","matrix")
#
  return(basis)
}

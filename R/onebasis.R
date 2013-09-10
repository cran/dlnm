###
### R routines for the R package dlnm (c) Antonio Gasparrini 2012-2013
#
`onebasis` <-
function(x, fun="ns", cen, ...) {
#
################################################################################
#
  nx <- names(x)
  x <- as.vector(x)
  range <- range(x,na.rm=TRUE)
  args <- list(...)
  args$x <- x
#
  # CHECK fun
  if(!is.character(fun)) stop("'fun' must be character")
  # CHECK OLD CODE
  .oldonebasis(fun,args)
  # CHECK fun HAS x ARGUMENT
  if(all(names(formals(fun))!="x")) stop("'fun' must contain argument 'x'")
#
###########################################################################
# TRANSFORMATION
#
  # CREATE THE BASIS
  basis <- do.call(fun,args)
  # FORCE TO BE A MATRIX (NOT WITH as.matrix AS IT DELETES ATTRIBUTES)
  if(is.null(dim(basis))) dim(basis) <- c(length(x),1)
  attr <- attributes(basis)
#
################################################################################
# CENTERING 
#
  if((!is.null(attr$int)&&attr$int) || 
    (fun%in%c("thr","strata","integer"))) cen <- FALSE
  if(missing(cen) || (is.logical(cen))&&cen) cen <- mean(x,na.rm=TRUE)
  ind <- match(names(formals(fun)),names(attr),nomatch=0)
  argscen <- if(all(ind==0)) list() else attr[ind]
  if(!is.logical(cen))
    basis <- sweep(basis,2,do.call(fun,args=c(list(x=cen),argscen)))
#
##########################################################################
#
  # NAMES AND ATTRIBUTES
  attributes(basis) <- c(list(fun=fun,cen=cen,range=range),attr)
  dimnames(basis) <- list(nx,paste("b",seq(ncol(basis)),sep=""))
#
  class(basis) <- c("onebasis","matrix")
#
  return(basis)
}

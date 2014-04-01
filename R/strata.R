###
### R routines for the R package dlnm (c) Antonio Gasparrini 2013-2014
#
strata <-
function(x, df=1, breaks=NULL, ref=1, int=FALSE) {
#
################################################################################
#
  nx <- names(x)
  x <- as.vector(x)
  range <- range(x,na.rm=TRUE)
#
  # DEFINE breaks AND df IF NEEDED
  if(!is.null(breaks)) {
    breaks <- sort(unique(breaks))
    df <- length(breaks)+int
  } else if(df-int>0) breaks <- quantile(x,1/(df-int+1)*1:((df-int)),na.rm=TRUE)
#
  # TRANSFORMATION
  xcat <- cut(x,c(range[1]-0.0001,breaks,range[2]+0.0001),right=FALSE)
  basis <- matrix(outer(xcat,levels(xcat),"==")+0,ncol=length(levels(xcat)))
#
  # DEFINE REFERENCE
  if(!ref%in%seq(ncol(basis)))
    stop("wrong value in 'ref' argument. See help('strata')")
  if(!is.null(breaks)) {
    basis <- basis[,-ref,drop=FALSE]
    if(int) basis <- cbind(1,basis)
  }
#
  # NAMES AND ATTRIBUTES
  dimnames(basis) <- list(nx,seq(ncol(basis)))
  attributes(basis) <- c(attributes(basis),list(df=df,breaks=breaks,ref=ref,
    int=int))
#
  class(basis) <- c("strata","matrix")
#
  return(basis)
}

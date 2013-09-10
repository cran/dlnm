###
### R routines for the R package dlnm (c) Antonio Gasparrini 2013
#
`strata` <-
function(x, df=1, breaks=NULL, int=FALSE) {
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
  x <- cut(x,c(range[1]-0.0001,breaks,range[2]+0.0001),right=FALSE)
  basis <- matrix(outer(x,levels(x),"==")+0,ncol=length(levels(x)))
  if(int==FALSE && !is.null(breaks)) basis <- basis[,-1,drop=FALSE]
#
  # NAMES AND ATTRIBUTES
  dimnames(basis) <- list(nx,seq(ncol(basis)))
  attributes(basis) <- c(attributes(basis),list(df=df,breaks=breaks,int=int))
#
  class(basis) <- c("strata","matrix")
#
  return(basis)
}

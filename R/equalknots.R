###
### R routines for the R package dlnm (c) Antonio Gasparrini 2013-2014
#
equalknots <-
function(x, nk=NULL, fun="ns", df=1, degree=3, int=FALSE) {
#
################################################################################
#
  x <- as.vector(x)
#
  range <- range(x,na.rm=TRUE) 
#
  # CHOOSE NUMBER OF KNOTS IF NOT PROVIDED
  if(is.null(nk)) {
    fun <- match.arg(fun,c("ns","bs","strata"))
    nk <- switch(fun,"ns"=df-1-int,"bs"=df-degree-int,"strata"=df-int)
  }
#
  # DEFINE KNOTS AT EQUALLY-SPACED VALUES ALONG range
  if(nk<1) stop("choice of arguments defines no knots")
  knots <- range[1] + (diff(range)/(nk+1))*seq(nk)
#
  return(knots)
}

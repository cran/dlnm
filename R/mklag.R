###
### R routines for the R package dlnm (c) Antonio Gasparrini 2013
#
`.mklag` <- 
function(lag) {
#
################################################################################
#
  #  lag MUST BE A POSITIVE INTEGER VECTOR 
  if(any(!is.numeric(lag))||length(lag)>2||any(lag<0)) 
    stop("'lag' must a positive integer vector or scalar")
  if(length(lag)==1L) lag <- c(0L,lag)
  if(diff(lag)<0L) stop("lag[1] must be <= lag[2]")
  return(round(lag[1L:2L]))
}


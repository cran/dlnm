###
### R routines for the R package dlnm (c) Antonio Gasparrini 2016
#
mkcen <-
function(cen, type, basis, range) {
#
################################################################################
#
  # IF GAM MODEL
  if(type=="gam") {
    # IF NULL OR TRUE, MID-RANGE
    if(is.null(cen) || (is.logical(cen)&&cen)) cen <- mean(range)
    # IF FALSE, NULL
    if(is.logical(cen)&&!cen) cen <- NULL
  }
#
  # ELSE, FOR ONEBASIS OR CROSSBASIS
  if(type%in%c("one","cb")) {
    # IF NULL, TRY TO EXTRACT IT FROM BASIS
    if(is.null(cen)) cen <- if(type=="one") attributes(basis)$cen else
      attributes(basis)$argvar$cen
#
    # DEPENDING ON FUNCTION
    fun <- if(type=="one") attributes(basis)$fun else 
      attributes(basis)$argvar$fun
    if(fun%in%c("thr","strata","integer","lin")) {
      if(is.logical(cen)) cen <- NULL
    } else {
      # IF NULL OR TRUE, MID-RANGE
      if(is.null(cen) || is.logical(cen)&&cen) cen <- mean(range)
      # IF FALSE, NULL
      if(is.logical(cen)&&!cen) cen <- NULL
    }
#
    # HOWEVER, IF INTERCEPT IS PRESENT, SET TO NULL
    int <- if(type=="one") attributes(basis)$intercept else 
      attributes(basis)$argvar$intercept
    if(is.logical(int)&&int) cen <- NULL
  }
#
  return(cen)
}

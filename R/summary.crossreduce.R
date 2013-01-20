###
### R routines for the R package dlnm (c) Antonio Gasparrini 2012-2013
#
`summary.crossreduce` <-
function(object, ...) {
#
################################################################################
#
  cat("REDUCED FIT\n")
  cat("dimension:",ifelse(is.null(object$var),"predictor","lag"),"\n")
  if(object$type!="overall") cat("value:",ifelse(is.null(object$var),
    paste("lag",object$lag),paste("var",object$var)),"\n")
  cat("reduced df:",length(coef(object)),"\n")
#
  cat("\nBASIS:\n")
  attr <- attributes(object$basis)
  cat("type:",attr$type)
  if(!is.null(attr$degree)) cat(" with degree",attr$degree)
  cat("\n")
  if(!is.null(attr$knots)) {
    cat("df:",attr$df,", knots at:",attr$knots,"\n")
  } else cat("df:",attr$df,"\n")
  if(!is.null(attr$bound)) cat("boundary knots at",attr$bound,"\n")
  if(!is.null(attr$cen)&&attr$cen==TRUE) cat("centered on",attr$cen,"\n")
  if(attr$int==TRUE) cat("with intercept\n")
#  
  cat("\nPREDICTIONS:\n")
  if(object$type!="var") cat("values:",length(object$predvar),"\n")
  if(object$type!="var") cat("range:",min(object$predvar),",",
    max(object$predvar),"\n")
  if(object$type=="var") cat("lag:",object$lag,"\n")
  cat("exponentiated:",ifelse(!is.null(object$RRfit),"yes","no"),"\n")
  cat("\n")
}

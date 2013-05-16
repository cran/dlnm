###
### R routines for the R package dlnm (c) Antonio Gasparrini 2012-2013
#
`summary.onebasis` <-
function(object, ...) {
#
################################################################################
#
  attr <- attributes(object)
  cat("BASIS FUNCTION\n")
  cat("observations:",nrow(object),"\n")
  cat("range:",attr$range[1],",",attr$range[2],"\n")
  cat("type:",attr$type)
  if(!is.null(attr$degree)) cat(" with degree",attr$degree)
  cat("\n")
  if(!is.null(attr$knots)) {
    cat("df:",attr$df,", knots at:",attr$knots,"\n")
  } else cat("df:",attr$df,"\n")
  if(!is.null(attr$bound)) cat("boundary knots at",attr$bound,"\n")
  if(is.logical(attr$cen)&&!attr$cen) {
    cat("not centered","\n")
  } else cat("centered at",attr$cen,"\n")
  cat(ifelse(attr$int==TRUE,"with","without"),"intercept\n")
}

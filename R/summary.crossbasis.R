`summary.crossbasis` <-
function(object, ...) {

  attr <- attributes(object)
  cat("CROSSBASIS FUNCTIONS\n")
  cat("observations:",nrow(object),"\n")
  if(!is.null(attr$group)) cat("groups:",attr$group,"\n")
  cat("range:",attr$range[1],",",attr$range[2],"\n")
  cat("total df:",attr$argvar$df*attr$arglag$df,"\n")
  cat("lag range:",attr$lag,"\n")

  cat("\nBASIS FOR VAR:\n")
  cat("type:",attr$argvar$type)
  if(!is.null(attr$argvar$degree)) cat(" with degree",attr$argvar$degree)
  cat("\n")
  if(!is.null(attr$argvar$knots)) {
    cat("df:",attr$argvar$df,", knots at:",attr$argvar$knots,"\n")
  } else cat("df:",attr$argvar$df,"\n")
  if(!is.null(attr$argvar$bound)) cat("boundary knots at",attr$argvar$bound,"\n")
  if(attr$argvar$cen==FALSE) {
    cat("not centered","\n")
  } else cat("centered on",attr$argvar$cen,"\n")
  cat(ifelse(attr$argvar$int==TRUE,"with","without"),"intercept\n")
  
  cat("\nBASIS FOR LAG:\n")
  cat("type:",attr$arglag$type)
  if(!is.null(attr$arglag$degree)) cat(" with degree",attr$arglag$degree)
  cat("\n")
  if(!is.null(attr$arglag$knots)) {
    cat("df:",attr$arglag$df,", knots at:",attr$arglag$knots,"\n")
  } else cat("df:",attr$arglag$df,"\n")
  if(!is.null(attr$arglag$bound)) cat("boundary knots at",attr$arglag$bound,"\n")
  cat(ifelse(attr$arglag$int==TRUE,"with","without"),"intercept\n")
  cat("\n")
}


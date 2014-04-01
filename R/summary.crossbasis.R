###
### R routines for the R package dlnm (c) Antonio Gasparrini 2012-2014
#
summary.crossbasis <-
function(object, ...) {
#
################################################################################
#
  attr <- attributes(object)
  cat("CROSSBASIS FUNCTIONS\n")
  cat("observations:",nrow(object),"\n")
  if(!is.null(attr$group)) cat("groups:",attr$group,"\n")
  cat("range:",attr$range[1],"to",attr$range[2],"\n")
  cat("lag period:",attr$lag,"\n")
  cat("total df: ",ncol(object),"\n")
#
  cat("\nBASIS FOR VAR:\n")
  for(i in seq(attr$argvar)[-2]) {
    cat(names(attr$argvar[i]),": ",sep="")
    cat(attr$argvar[[i]],"\n",sep=" ")
  }
  if(!is.logical(attr$argvar$cen) )
    cat("centered at",attr$argvar$cen,"\n") else cat("not centered","\n")
#  
  cat("\nBASIS FOR LAG:\n")
  for(i in seq(attr$arglag)[-2]) {
    cat(names(attr$arglag[i]),": ",sep="")
    cat(attr$arglag[[i]],"\n",sep=" ")
  }
  cat("\n")
}


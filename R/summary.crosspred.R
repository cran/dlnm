`summary.crosspred` <-
function(object, ...) {

	cat("PREDICTED VALUES:\n")
	cat("values:",length(object$predvar),"\n")
	cat("range:",min(object$predvar),",",max(object$predvar),"\n")
	cat("maxlag:",object$maxlag,"\n")
	cat("exponentiated:",ifelse(!is.null(object$allRRfit),"yes","no"),"\n")
	cat("cumulative:",ifelse(!is.null(object$cumfit),"yes","no"),"\n")

	cat("\nMODEL:\n")
	cat("parameters:",length(object$coef),"\n")
	cat("class:",object$model.class,"\n")
	cat("link:",object$model.link)
	cat("\n")
}


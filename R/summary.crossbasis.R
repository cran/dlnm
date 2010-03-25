`summary.crossbasis` <-
function(object, ...) {

	attr <- attributes(object)
	cat("CROSSBASIS FUNCTIONS\n")
	cat("observations:",nrow(object),"\n")
	if(attr$group>1) cat("groups:",attr$group,"\n")
	cat("range:",attr$range[1],",",attr$range[2],"\n")
	cat("total df:",attr$crossdf,"\n")
	cat("maxlag:",attr$maxlag,"\n")

	cat("\nBASIS FOR VAR:\n")
	cat("type:",attr$vartype)
	if(!is.null(attr$vardegree)) cat(" with degree",attr$vardegree)
	cat("\n")
	if(!is.null(attr$varknots)) {
		cat("df:",attr$vardf,", knots at:",attr$varknots,"\n")
	} else cat("df:",attr$vardf,"\n")
	if(!is.null(attr$varbound)) cat("boundary knots at",attr$varbound,"\n")
	if(attr$cen==TRUE) cat("centered on",attr$cenvalue,"\n")
	if(attr$varint==TRUE) cat("with intercept\n")

	cat("\nBASIS FOR LAG:\n")
	cat("type:",attr$lagtype)
	if(!is.null(attr$lagdegree)) cat(" with degree",attr$lagdegree)
	cat("\n")
	if(!is.null(attr$lagknots)) {
		cat("df:",attr$lagdf,", knots at:",attr$lagknots,"\n")
	} else cat("df:",attr$lagdf,"\n")
	if(!is.null(attr$lagbound)) cat("boundary knots at",attr$lagbound,"\n")
	if(attr$lagint==TRUE) cat("with intercept\n")

}


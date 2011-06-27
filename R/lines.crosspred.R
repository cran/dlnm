`lines.crosspred` <-
function(x, ptype="overall", var=NULL, lag=NULL, ci="n", ci.arg,
	ci.level=x$ci.level, cumul=FALSE, exp=NULL, ...) {

###########################################################################
# CHECKS

if(class(x)!="crosspred") stop("'x' must be of class 'crosspred'")
if(!ptype%in%c("slices","overall")) {
	stop("'ptype' must be one of 'slices' or 'overall'")
}
if(!xor(is.null(var),is.null(lag))&(ptype=="slices")) {
	stop("One (only) of 'var' or 'lag' must be provided when ptype='slices'")
}
if(!is.null(var)&&!is.numeric(var)&&length(var)>1&&ptype=="slices") {
	stop("'var' must be a numeric scalar")
}
if(!is.null(lag)&&!is.numeric(lag)&&length(lag)>1&&ptype=="slices") {
	stop("'lag' must be a numeric scalar")
}
if(!is.null(var)&&!var%in%x$predvar&&(ptype=="slices")) {
	stop("'var' must match values used for prediction")
}
if(!is.null(lag)&&!lag%in%0:x$maxlag&&(ptype=="slices")) {
	stop("'lag' must match values in 0:maxlag")
}
if(!ci%in%c("area","bars","lines","n")) {
	stop("'ci' must be one of 'area', 'bars', 'lines' or 'n'")
}
if(missing(ci.arg)) {
	ci.arg <- list()
} else if(!is.list(ci.arg)) stop("'ci.arg' must be a list")
if(!is.numeric(ci.level)||ci.level>=1||ci.level<=0) {
	stop("'ci.level' must be numeric and between 0 and 1")
}
if(cumul==TRUE&&is.null(x$cumfit)) {
	stop("Cumulative effects can be plotted if predicted in the 'crosspred'
	object. Set the argument 'cumul=TRUE' in the function crosspred()")
}
if(!is.null(exp)&&!is.logical(exp)) stop("'exp' must be logical")

##########################################################################
# COMPUTE EFFECTS

# CUMULATIVE IF CUMUL==T
if(cumul==TRUE) {
	x$matfit <- x$cumfit
	x$matse <- x$cumse
}

# SET THE Z LEVEL EQUAL TO THAT STORED IN OBJECT IF NOT PROVIDED
z <- qnorm(1-(1-ci.level)/2)
x$mathigh <- x$matfit+z*x$matse
x$matlow <- x$matfit-z*x$matse
x$allhigh <- x$allfit+z*x$allse
x$alllow <- x$allfit-z*x$allse

# EXPONENTIAL
if((is.null(exp)&&x$model.link%in%c("log","logit"))||
	(!is.null(exp)&&exp==TRUE)) {
	x$matfit <- exp(x$matfit)
	x$mathigh <- exp(x$mathigh)
	x$matlow <- exp(x$matlow)
	x$allfit <- exp(x$allfit)
	x$allhigh <- exp(x$allhigh)
	x$alllow <- exp(x$alllow)
}

##########################################################################
# GRAPHS

# FUNCTION TO CREATE CONFIDENCE INTERVALS
fci.lines <- function(ci,x,high,low,ci.arg,plot.arg) {
	if(ci=="area") {
		polygon.arg <- modifyList(list(col=grey(0.9),border=NA),ci.arg)
		polygon.arg <- modifyList(polygon.arg,
			list(x=c(x,rev(x)),y=c(high,rev(low))))
		do.call(polygon,polygon.arg)
	} else if(ci=="bars") {
		range <- diff(range(x))/300
		segments.arg <- modifyList(ci.arg,list(x0=x,y0=high,x1=x,y1=low))
		do.call(segments,segments.arg)
		segments.arg <- modifyList(segments.arg,list(x0=x-range,y0=high,
			x1=x+range,y1=high))
		do.call(segments,segments.arg)
		segments.arg <- modifyList(segments.arg,list(x0=x-range,y0=low,
			x1=x+range,y1=low))
		do.call(segments,segments.arg)
	} else if(ci=="lines") {
		lines.arg <- modifyList(list(lty=2,col=plot.arg$col),ci.arg)
		lines.arg <- modifyList(lines.arg,list(x=x,y=high))
		do.call(lines,lines.arg)
		lines.arg <- modifyList(lines.arg,list(x=x,y=low))
		do.call(lines,lines.arg)
	}
}

##########
# SLICES
##########

if(ptype=="slices") {

	# LAG
	if(!is.null(lag)) {

		xlag <- paste("lag",lag,sep="")

		# SET DEFAULT VALUES IF NOT INCLUDED BY THE USER
		plot.arg <- list(type="l",col=2)
		plot.arg <- modifyList(plot.arg,list(...))		

		# PLOT CONFIDENCE INTERVALS (IF ANY)
		fci.lines(ci=ci,x=x$predvar,high=x$mathigh[,xlag],
			low=x$matlow[,xlag],ci.arg,plot.arg)
		plot.arg <- modifyList(plot.arg,c(list(x=x$predvar,
			y=x$matfit[,xlag])))
		do.call(lines,plot.arg)
	}

	# VAR
	if(!is.null(var)) {

		xvar <- as.character(var)

		# SET DEFAULT VALUES IF NOT INCLUDED BY THE USER
		plot.arg <- list(type="l",col=2)
		plot.arg <- modifyList(plot.arg,list(...))		

		# PLOT CONFIDENCE INTERVALS (IF ANY)
		fci.lines(ci=ci,x=0:x$maxlag,high=x$mathigh[xvar,],
			low=x$matlow[xvar,],ci.arg,plot.arg)
		plot.arg <- modifyList(plot.arg,c(list(x=0:x$maxlag,
			y=x$matfit[xvar,])))
		do.call(lines,plot.arg)
	}

}

##########
# OVERALL
##########

if(ptype=="overall") {

	# SET DEFAULT VALUES IF NOT INCLUDED BY THE USER
	plot.arg <- list(type="l",col=2)
	plot.arg <- modifyList(plot.arg,list(...))
	
	# SET CONFIDENCE INTERVALS
	fci.lines(ci=ci,x=x$predvar,high=x$allhigh,
		low=x$alllow,ci.arg,plot.arg)
	plot.arg <- modifyList(plot.arg,c(list(x=x$predvar,y=x$allfit)))

	# PLOT
	do.call(lines,plot.arg)
}
}

#


`points.crosspred` <-
function(x, ptype="overall", var=NULL, lag=NULL, ci="n", 
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
if(cumul==T) {
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
fci.points <- function(ci,x,high,low,colarea=grey(0.9),collines=2) {
	if(ci=="area") {
		polygon(c(x,rev(x)),c(high,rev(low)),col=colarea,border="white")
	} else if(ci=="bars") {
		range <- diff(range(x))/300
		segments(x,high,x,low)
		segments(x-range,high,x+range,high)
		segments(x-range,low,x+range,low)
	} else if(ci=="lines") {
		lines(x,high,lty=2,col=collines)
		lines(x,low,lty=2,col=collines)
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
		argdef <- list(type="p",col=2)
		argdef <- modifyList(argdef,list(...))		

		# PLOT CONFIDENCE INTERVALS (IF ANY)
		fci.points(ci=ci,x=x$predvar,high=x$mathigh[,xlag],
			low=x$matlow[,xlag],colarea=grey(0.9),collines=argdef$col)
		argdef <- modifyList(argdef,c(list(x=x$predvar,
			y=x$matfit[,xlag])))
		do.call(points,argdef)
	}

	# VAR
	if(!is.null(var)) {

		xvar <- as.character(var)

		# SET DEFAULT VALUES IF NOT INCLUDED BY THE USER
		argdef <- list(type="p",col=2)
		argdef <- modifyList(argdef,list(...))		

		# PLOT CONFIDENCE INTERVALS (IF ANY)
		fci.points(ci=ci,x=0:x$maxlag,high=x$mathigh[xvar,],
			low=x$matlow[xvar,],colarea=grey(0.9),collines=argdef$col)
		argdef <- modifyList(argdef,c(list(x=0:x$maxlag,
			y=x$matfit[xvar,])))
		do.call(points,argdef)
	}

}

##########
# OVERALL
##########

if(ptype=="overall") {

	# SET DEFAULT VALUES IF NOT INCLUDED BY THE USER
	argdef <- list(type="p",col=2)
	argdef <- modifyList(argdef,list(...))

	# SET CONFIDENCE INTERVALS
	fci.points(ci=ci,x=x$predvar,high=x$allhigh,
		low=x$alllow,colarea=grey(0.9),collines=argdef$col)
	argdef <- modifyList(argdef,c(list(x=x$predvar,y=x$allfit)))

	# PLOT
	do.call(points,argdef)
}
}

#


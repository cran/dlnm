`plot.crosspred` <-
function(x, ptype="3d", var=NULL, lag=NULL, ci="area", 
	ci.level=x$ci.level, cumul=FALSE, exp=NULL, ...) {

###########################################################################
# CHECKS

if(class(x)!="crosspred") stop("'x' must be of class 'crosspred'")
if(!ptype%in%c("slices","3d","contour","overall")) {
	stop("'ptype' must be one of 'slices','3d','contour','overall'")
}
if(is.null(var)&&is.null(lag)&&(ptype=="slices")) {
	stop("at least 'var' or 'lag' must be provided when ptype='slices'")
}
if(!is.null(var)&&!is.numeric(var)&&length(var)>4&&ptype=="slices") {
	stop("'var' and 'lag' must be numeric and of length <=4")
}
if(!is.null(lag)&&!is.numeric(lag)&&length(lag)>4&&ptype=="slices") {
	stop("'var' and 'lag' must be numeric and of length <=4")
}
if((!is.null(var)&!is.null(lag))&&length(var)!=length(lag)&&ptype=="slices") {
	stop("if both are provided, length of 'var' and 'lag' must be the same")
}
if(!is.null(var)&&sum(var%in%x$predvar)!=length(var)&&(ptype=="slices")) {
	stop("'var' must match values used for prediction")
}
if(!is.null(lag)&&sum(lag%in%0:x$maxlag)!=length(lag)&&(ptype=="slices")) {
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
noeff <- 0

# EXPONENTIAL
if((is.null(exp)&&x$model.link%in%c("log","logit"))||
	(!is.null(exp)&&exp==TRUE)) {
	x$matfit <- exp(x$matfit)
	x$mathigh <- exp(x$mathigh)
	x$matlow <- exp(x$matlow)
	x$allfit <- exp(x$allfit)
	x$allhigh <- exp(x$allhigh)
	x$alllow <- exp(x$alllow)
	noeff <- 1
}

##########################################################################
# GRAPHS

# FUNCTION TO CREATE CONFIDENCE INTERVALS
fci <- function(ci,x,high,low,colarea=grey(0,9),collines=2,noeff) {
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
	abline(h=noeff)
}

##########
# SLICES
##########

if(ptype=="slices") {

	# SET FRAME AND GREY SCALE
	mar.old <- par()$mar
	mfrow.old <- par()$mfrow
	mgp.old <- par()$mgp
	grey <- grey(0.9)
	if(length(lag)+length(var)>1) {
		layout(matrix(1:(length(var)+length(lag)),
			ncol=sum(!is.null(var),!is.null(lag))))
		grey <- grey(0.8)
		par(mgp=c(2,0.7,0),mar=c(4.1,4.1,2.1,1.1))
	}

	# LAG
	if(!is.null(lag)) {

		# START LOOP FOR LAG
		xlag <- paste("lag",lag,sep="")
		for(i in seq(along=lag)) {

			# SET DEFAULT VALUES IF NOT INCLUDED BY THE USER
			argdef <- list(type="l",col=2,xlab="Var",ylab="Effect",
				ylim=c(min(x$matlow[,xlag]),
				max(x$mathigh[,xlag])),frame.plot=FALSE)
			if(length(lag)+length(var)>1)  argdef$cex.axis <- 0.7
			argdef <- modifyList(argdef,list(...))		

			# SET CONFIDENCE INTERVALS
			ci.list <- list(panel.first=call("fci",ci=ci,
				x=x$predvar,high=x$mathigh[,xlag[i]],
				low=x$matlow[,xlag[i]],colarea=grey,
				collines=argdef$col,noeff=noeff))
			argdef <- modifyList(argdef,c(ci.list,
				list(x=x$predvar,y=x$matfit[,xlag[i]])))
			if(length(lag)+length(var)>1) {
				argdef$main <- ""
				argdef$xlab <- "Var"
			}
			# PLOT
			do.call(plot,argdef)
			if(length(lag)>1) mtext(paste("Lag =",lag[i]),cex=0.8)
		}
	}

	# VAR
	if(!is.null(var)) {

		# START LOOP FOR VAR
		xvar <- as.character(var)
		for(i in seq(along=var)) {

			# SET DEFAULT VALUES IF NOT INCLUDED BY THE USER
			argdef <- list(type="l",col=2,xlab="Lag",ylab="Effect",
				ylim=c(min(x$matlow[xvar,]),
				max(x$mathigh[xvar,])),frame.plot=FALSE)
			if(length(lag)+length(var)>1)  argdef$cex.axis <- 0.7

			argdef <- modifyList(argdef,list(...))		

			# SET CONFIDENCE INTERVALS
			ci.list <- list(panel.first=call("fci",ci=ci,
				x=0:x$maxlag,high=x$mathigh[xvar[i],],
				low=x$matlow[xvar[i],],colarea=grey,
				collines=argdef$col,noeff=noeff))
			argdef <- modifyList(argdef,c(ci.list,
				list(x=0:x$maxlag,y=x$matfit[xvar[i],])))
			if(length(lag)+length(var)>1) {
				argdef$main <- ""
				argdef$xlab <- "Lag"
			}
			# PLOT
			do.call(plot,argdef)
			if(length(lag)>1) mtext(paste("Var =",var[i]),cex=0.8)
		}
	}
	if(length(lag)+length(var)>1) {
		par(mar=mar.old,mfrow=mfrow.old,mgp=mgp.old)
	}
}

##########
# OVERALL
##########
if(ptype=="overall") {

	# SET DEFAULT VALUES IF NOT INCLUDED BY THE USER
	argdef <- list(type="l",ylim=c(min(x$alllow),max(x$allhigh)),
		col=2,xlab="Var",ylab="Effect",frame.plot=FALSE)
	argdef <- modifyList(argdef,list(...))
	
	# SET CONFIDENCE INTERVALS
	ci.list <- list(panel.first=call("fci",ci=ci,x=x$predvar,
		high=x$allhigh,low=x$alllow,colarea=grey(0.9),
		collines=argdef$col,noeff=noeff))
	argdef <- modifyList(argdef,c(ci.list,
		list(x=x$predvar,y=x$allfit)))

	# PLOT
	do.call(plot,argdef)
}

##############
# CONTOURPLOT
##############

if(ptype=="contour") {

	# SET DEFAULT VALUES (NOT TO BE SPECIFIED BY THE USER)
	level <- pretty(x$matfit,20)
	col1 <- colorRampPalette(c("blue","white"))
	col2 <- colorRampPalette(c("white","red"))
	col <- c(col1(sum(level<noeff)),col2(sum(level>noeff)))

	filled.contour(x=x$predvar,y=0:x$maxlag,z=x$matfit,col=col,
		level=level,...)
}

#########
# 3-D
#########

if(ptype=="3d") {

	# SET DEFAULT VALUES IF NOT INCLUDED BY THE USER
	argdef <- list(ticktype="detailed",theta=210,phi=30,xlab="Var",
		ylab="Lag",	zlab="Effect",col="lightskyblue",
		zlim=c(min(x$matfit),max(x$matfit)),
		ltheta=290,shade=0.75,r=sqrt(3),d=5)
	argdef <- modifyList(argdef,list(...))
	argdef <- modifyList(argdef,list(x=x$predvar,y=0:x$maxlag,z=x$matfit))

	# PLOT
	do.call(persp,argdef)
}
}


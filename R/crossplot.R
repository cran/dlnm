`crossplot` <-
function(crosspred, type="3d",
	var=NULL, lag=NULL, ylim=NULL, title=NULL, label="var") {

if(!type%in%c("slices","3d","contour","overall")) {
	stop("Type must be one of 'slices','3d','contour','overall'")
}

if(is.null(var)&is.null(lag)&(type=="slices")) {
	stop("at least 'var' or 'lag' must be provided when type='slices'")
}
if((length(var)>4|length(lag)>4)&type=="slices") {
	stop("length of 'var' and 'lag' must be <=4")
}
if((length(var)!=length(lag))&(!is.null(var)&!is.null(lag))&type=="slices") {
	stop("if both are provided, length of 'var' and 'lag' must be the same")
}
if(!is.null(var)&sum(var%in%crosspred$predvar)!=length(var)&
	(type=="slices"|type=="slice")) {
	stop("'var' must match values used for prediction")
}
if(!is.null(lag)&sum(lag%in%0:crosspred$maxlag)!=length(lag)&
	(type=="slices"|type=="slice")) {
	stop("'lag' must match values in 0:maxlag")
}

##########################################################################
# GRAPHS
###########

if(crosspred$model.link %in% c("log","logit")) {
	crosspred$matfit <- crosspred$matRRfit
	crosspred$mathigh <- crosspred$matRRhigh
	crosspred$matlow <- crosspred$matRRlow
	crosspred$allfit <- crosspred$allRRfit
	crosspred$allhigh <- crosspred$allRRhigh
	crosspred$alllow <- crosspred$allRRlow
	plotlab <- "RR"
	noeff <- 1
} else {
	crosspred$mathigh <- crosspred$matfit+1.96*crosspred$matse
	crosspred$matlow <- crosspred$matfit-1.96*crosspred$matse
	crosspred$allhigh <- crosspred$allfit+1.96*crosspred$allse
	names(crosspred$allhigh) <- names(crosspred$allfit)
	crosspred$alllow <- crosspred$allfit-1.96*crosspred$allse
	names(crosspred$alllow) <- names(crosspred$allfit)
	plotlab <- "Effect"
	noeff <- 0
}

##########
# SLICES
##########


if(type=="slices") {
	mar.old <- par()$mar
	mfrow.old <- par()$mfrow
	grey <- grey(0.9)
	if(length(lag)+length(var)>1) {
		layout(matrix(1:(length(var)+length(lag)),
			ncol=sum(!is.null(var),!is.null(lag))))
		title <- ""
		grey <- grey(0.8)
		par(mar=c(3.1,4.1,2.1,1.1))
	}
	if(!is.null(lag)) {
		predvar <- crosspred$predvar
		x <- paste("lag",lag,sep="")
		if(is.null(ylim)) {
			ylimlag <- c(min(crosspred$matlow[,x]),
				max(crosspred$mathigh[,x]))
		} else ylimlag <- ylim
		for(i in 1:length(lag)) {
			plot(predvar,crosspred$matfit[,x[i]],type="l",
				frame.plot=FALSE,xlab="",ylab="",ylim=ylimlag,axes=FALSE)
			axis(1,cex.axis=0.7,mgp=c(3,0.7,0))
			axis(2,cex.axis=0.7,mgp=c(3,0.7,0))
			polygon(c(predvar,rev(predvar)),c(crosspred$mathigh[,x[i]],
				rev(crosspred$matlow[,x[i]])),col=grey,
				border="white")
			abline(h=noeff)
			lines(predvar,crosspred$matfit[,x[i]],col="red")
			if(length(lag)>1) mtext(paste("Lag =",lag[i]),cex=0.8)
			title(main=title,ylab=plotlab,xlab=label,mgp=c(2,0.7,0),cex=0.7)
			if(length(lag)==1) title(title,mgp=c(3,0.7,0))
		}
	}
	if(!is.null(var)) {
		predlag <- 0:crosspred$maxlag
		x <- as.character(var)
		if(is.null(ylim)) {
			ylimvar <- c(min(crosspred$matlow[x,]),
				max(crosspred$mathigh[x,]))
		} else ylimvar <- ylim
		for(i in 1:length(var)) {
			plot(predlag,crosspred$matfit[x[i],],type="n",frame.plot=FALSE,
				xlab="",ylab="",ylim=ylimvar,axes=FALSE)
			axis(1,cex.axis=0.7,mgp=c(3,0.7,0))
			axis(2,cex.axis=0.7,mgp=c(3,0.7,0))
			polygon(c(predlag,rev(predlag)),c(crosspred$mathigh[x[i],],
				rev(crosspred$matlow[x[i],])),col=grey,
				border="white")
			abline(h=noeff)
			lines(0:crosspred$maxlag,crosspred$matfit[x[i],],col="red")
			if(length(var)>1) mtext(paste(label,"=",var[i]),cex=0.8)
			title(main=title,ylab=plotlab,xlab="Lag",mgp=c(2,0.7,0),cex=0.7)
			if(length(lag)==1) title(title,mgp=c(3,0.7,0))
		}
	}
	if(length(lag)+length(var)>1) {
		par(mar=mar.old)
		par(mfrow=mfrow.old)
	}
}

##########
# OVERALL
##########

if(type=="overall") {
	predvar <- crosspred$predvar
	if(is.null(ylim)) {
		ylim <- c(min(crosspred$alllow),max(crosspred$allhigh))
	}
	plot(predvar,crosspred$allfit,type="n",frame.plot=FALSE,xlab="",
		ylab=plotlab,ylim=ylim,axes=FALSE)
	axis(1,xlab=predvar)
	axis(2,xlab=crosspred$allfit)
	polygon(c(predvar,rev(predvar)),
		c(crosspred$allhigh,rev(crosspred$alllow)),
			col=grey(0.9), border="white")
	abline(h=noeff)
	lines(predvar,crosspred$allfit,col="red")
	title(main=title,xlab=label)

}

##############
# CONTOURPLOT
##############

if(type=="contour") {
	col1 <- colorRampPalette(c("blue","white"))
	col2 <- colorRampPalette(c("white","red"))
	level <- pretty(crosspred$matfit,20)
	col <- c(col1(sum(level<noeff)),col2(sum(level>noeff)))
	filled.contour(crosspred$predvar,0:crosspred$maxlag,
		crosspred$matfit,col= col,
		plot.title = title(main = title,
			xlab = label, ylab = "Lags"),
	    	key.title = title(main=plotlab))
}

#########
# 3-D
#########

if(type=="3d") {
	mar.old <- par()$mar
	par(mar=c(2.1,1.1,2.1,1.1))
	if(is.null(ylim)) {
		ylim <- c(min(crosspred$matfit),max(crosspred$matfit))
	}
	persp(crosspred$predvar,0:crosspred$maxlag,crosspred$matfit,
		ticktype="detailed",theta = 210,phi = 30,xlab=label, ylab="Lag",
		zlab=plotlab,main=title,col="lightskyblue",zlim=ylim,
		ltheta=290,shade=0.75,r=sqrt(3),d=5)
	par(mar=mar.old)
}
}


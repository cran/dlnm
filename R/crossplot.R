`crossplot` <-
function(crosspred, type="3d", cumul=FALSE, ci="area",
  var=NULL, lag=NULL, ylim=NULL, title=NULL, label="var") {
  
  
warning("The function 'crossplot' has been replaced by method functions for class
'crosspred' since version 1.3.0, and kept only for compatibility reasons.
The users are strongly suggested to use the new functions\n")

if(!type%in%c("slices","3d","contour","overall")) {
  stop("type must be one of 'slices','3d','contour','overall'")
}
if(!ci%in%c("area","bars","lines")) {
  stop("'ci' must be one of 'area', 'bars' or 'lines'")
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
  (type=="slices")) {
  stop("'var' must match values used for prediction")
}
if(!is.null(lag)&sum(lag%in%.seq(crosspred$lag))!=length(lag)&
  (type=="slices")) {
  stop("'lag' must match lag values used for prediction")
}
if(cumul==TRUE&is.null(crosspred$cumfit)) {
  stop("Cumulative effects can be plotted if predicted in the 'crosspred'
object. Set the argument 'cumul=TRUE' in the function crosspred()")
}

##########################################################################
# GRAPHS
###########

if(cumul==T) {
  crosspred$matfit <- crosspred$cumfit
  crosspred$matse <- crosspred$cumse
  crosspred$matRRfit <- crosspred$cumRRfit
  crosspred$matRRhigh <- crosspred$cumRRhigh
  crosspred$matRRlow <- crosspred$cumRRlow
}
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
  if(cumul==TRUE) plotlab <- paste("Cumulative",plotlab)
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
      ylimlag <- c(min(crosspred$matlow[,x]),max(crosspred$mathigh[,x]))
    } else ylimlag <- ylim
    for(i in 1:length(lag)) {
      plot(predvar,crosspred$matfit[,x[i]],type="n",
        frame.plot=FALSE,xlab="",ylab="",ylim=ylimlag,axes=FALSE)
      axis(1,cex.axis=0.7,mgp=c(3,0.7,0))
      axis(2,cex.axis=0.7,mgp=c(3,0.7,0))
      if(ci=="area") {
        polygon(c(predvar,rev(predvar)),c(crosspred$mathigh[,x[i]],
          rev(crosspred$matlow[,x[i]])),col=grey,border="white")
        abline(h=noeff)
        lines(predvar,crosspred$matfit[,x[i]],col="red")
      } else if(ci=="bars") {
        segments(predvar,crosspred$mathigh[,x[i]],
          predvar,crosspred$matlow[,x[i]])
        segments(predvar-0.1,crosspred$mathigh[,x[i]],
          predvar+0.1,crosspred$mathigh[,x[i]])
        segments(predvar-0.1,crosspred$matlow[,x[i]],
          predvar+0.1,crosspred$matlow[,x[i]])
        abline(h=noeff)
        points(predvar,crosspred$matfit[,x[i]],col="red",pch=19)				
      } else if(ci=="lines") {
        lines(predvar,crosspred$mathigh[,x[i]],lty=2)
        lines(predvar,crosspred$matlow[,x[i]],lty=2)
        abline(h=noeff)
        lines(predvar,crosspred$matfit[,x[i]],col="red")
      }
      lines(predvar,crosspred$matfit[,x[i]],col="red")
      if(length(lag)>1) mtext(paste("Lag =",lag[i]),cex=0.8)
      title(main=title,ylab=plotlab,xlab=label,mgp=c(2,0.7,0),cex=0.7)
      if(length(lag)==1) title(title,mgp=c(3,0.7,0))
    }
  }
  if(!is.null(var)) {
    predlag <- .seq(crosspred$lag)
    x <- as.character(var)
    if(is.null(ylim)) {
      ylimvar <- c(min(crosspred$matlow[x,]),max(crosspred$mathigh[x,]))
    } else ylimvar <- ylim
    for(i in 1:length(var)) {
    plot(predlag,crosspred$matfit[x[i],],type="n",frame.plot=FALSE,
      xlab="",ylab="",ylim=ylimvar,axes=FALSE)
    axis(1,cex.axis=0.7,mgp=c(3,0.7,0))
    axis(2,cex.axis=0.7,mgp=c(3,0.7,0))
    if(ci=="area") {
      polygon(c(predlag,rev(predlag)),c(crosspred$mathigh[x[i],],
        rev(crosspred$matlow[x[i],])),col=grey,border="white")
      abline(h=noeff)
      lines(predlag,crosspred$matfit[x[i],],col="red")
    } else if(ci=="bars") {
      segments(predlag,crosspred$mathigh[x[i],],
        predlag,crosspred$matlow[x[i],])
      segments(predlag-0.1,crosspred$mathigh[x[i],],
        predlag+0.1,crosspred$mathigh[x[i],])
      segments(predlag-0.1,crosspred$matlow[x[i],],
        predlag+0.1,crosspred$matlow[x[i],])
      abline(h=noeff)
      points(.seq(crosspred$lag),crosspred$matfit[x[i],],col="red",pch=19)			
      } else if(ci=="lines") {
        lines(predlag,crosspred$mathigh[x[i],],lty=2)
        lines(predlag,crosspred$matlow[x[i],],lty=2)
        abline(h=noeff)
        lines(.seq(crosspred$lag),crosspred$matfit[x[i],],col="red")
      }
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
  if(ci=="area") {
    polygon(c(predvar,rev(predvar)),c(crosspred$allhigh,
      rev(crosspred$alllow)),col=grey(0.9), border="white")
    abline(h=noeff)
    lines(predvar,crosspred$allfit,col="red")
  } else if(ci=="bars") {
    segments(predvar,crosspred$allhigh,predvar,crosspred$alllow)
    segments(predvar-0.1,crosspred$allhigh,predvar+0.1,crosspred$allhigh)
    segments(predvar-0.1,crosspred$alllow,predvar+0.1,crosspred$alllow)
    abline(h=noeff)
    points(predvar,crosspred$allfit,col="red",pch=19)			
  } else if(ci=="lines") {
    lines(predvar,crosspred$allhigh,lty=2)
    lines(predvar,crosspred$alllow,lty=2)
    abline(h=noeff)
    lines(predvar,crosspred$allfit,col="red")
  }
  abline(h=noeff)
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
  filled.contour(crosspred$predvar,.seq(crosspred$lag),
    crosspred$matfit,col= col,
    plot.title = title(main = title,xlab = label, ylab = "Lags"),
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
  persp(crosspred$predvar,.seq(crosspred$lag),crosspred$matfit,
    ticktype="detailed",theta = 210,phi = 30,xlab=label, ylab="Lag",
    zlab=plotlab,main=title,col="lightskyblue",zlim=ylim,
    ltheta=290,shade=0.75,r=sqrt(3),d=5)
  par(mar=mar.old)
}

}


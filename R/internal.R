###
### R routines for the R package dlnm (c) Antonio Gasparrini 2012-2013
#
`.onAttach` <- 
function(lib, pkg) {
#
################################################################################
#
  meta <- packageDescription("dlnm")
  attachmsg <- paste("This is dlnm ",meta$Version,
    ". For details: help(dlnm) and vignette('dlnmOverview').",sep="")
  #attachmsg <- paste("Important changes since version 1.5.1\nSee:",
  # "'file.show(system.file('Changesince151',package='dlnm'))'")
  packageStartupMessage(attachmsg, domain = NULL, appendLF = TRUE)
}
###
### R routines for the R package dlnm (c) Antonio Gasparrini 2012-2013
#
`.fci` <-
function(ci, x, high, low, ci.arg, plot.arg, noeff=NULL) {
#
################################################################################
#
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
  if(!is.null(noeff)) abline(h=noeff)
}
###
### R routines for the R package dlnm (c) Antonio Gasparrini 2012-2013
#
`.seq` <-
function(lag,by=1) seq(from=lag[1],to=lag[2],by=by)
#
################################################################################
#
###
### Code from Roger Peng, included in the function Lag() in the package tsModel
#
`.Lag` <-
function(v, k, group = NULL) 
{
  stopifnot(length(k) > 0)
  v <- as.numeric(v)
  if (max(abs(k)) >= length(v)) 
    stop("largest lag in 'k' must be less than 'length(v)'")
  lag.f <- function(x) {
    lagmat <- matrix(nrow = length(x), ncol = length(k))
    n <- length(x)
    for (i in seq(along = k)) {
      lag <- k[i]
      if (lag > 0) 
        lagmat[, i] <- c(rep(NA, lag), x[1:(n - lag)])
      else if (lag < 0) 
        lagmat[, i] <- c(x[(-lag + 1):n], rep(NA, -lag))
      else lagmat[, i] <- x
    }
    lagmat
  }
  lagmat <- if (!is.null(group)) {
    groupLag <- tapply(v, group, lag.f)
    lagmat <- matrix(nrow = length(v), ncol = length(k))
    split(lagmat, group) <- groupLag
    lagmat
  }
  else lag.f(v)
  colnames(lagmat) <- as.character(k)
  drop(lagmat)
}

#

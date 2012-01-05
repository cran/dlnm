`.onAttach` <- function(lib, pkg) {
  meta <- packageDescription("dlnm")
  #attachmsg <- paste("This is dlnm ",meta$Version,
  # ". For details: help(dlnm) and vignette('dlnmOverview').",sep="")
  attachmsg <- paste("Important changes since version 1.5.1\nSee:",
    "'file.show(system.file('Changesince151',package='dlnm'))'")
  packageStartupMessage(attachmsg, domain = NULL, appendLF = TRUE)
}


`.fci` <-
function(ci, x, high, low, ci.arg, plot.arg, noeff=NULL) {
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

`.seq` <-
function(lag) lag[1]:lag[2]

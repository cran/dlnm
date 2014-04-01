###
### R routines for the R package dlnm (c) Antonio Gasparrini 2013-2014
#
checkgroup <- 
function(group,x,basisvar,lag) {
#
################################################################################
#
  if(NCOL(x)>1L) stop("'group' allowed only for time series data")
  if(min(tapply(x,group,length))<=diff(lag))
    stop("each group must have length > diff(lag)")
}
#
#
checkoldonebasis <- 
function(fun,args) {
#
################################################################################
#
  # ARGUMENT bound FOR SPLINE FUNCTIONS
  if(fun%in%c("ns","bs")&&!is.null(args$bound)) {
    names(args)[names(args)=="bound"] <- "Boundary.knots"
    warning("use the default argument 'Boundary.knots' for fun 'ns'-'bs'")
  }  
  # OLD THRESHOLD FUNCTIONS
  if(fun%in%c("hthr","lthr","dthr")) {
    args$side <- switch(fun,hthr="h",lthr="l",dthr="d")
    fun <- "thr"
    warning("function 'hthr'-'lthr'-'dthr' replaced by 'thr'. See ?thr")
  }
  if(fun=="thr"&&!is.null(args$knots)) {
    names(args)[names(args)=="knots"] <- "thr.value"
    warning("argument 'knots' replaced by 'thr.value' in function thr. See ?thr")
  }
  # OLD STRATA FUNCTION
  if(fun=="strata"&&!is.null(args$knots)) {
    names(args)[names(args)=="knots"] <- "breaks"
    warning("argument 'knots' replaced by 'breaks' in function strata. See ?strata")
  }
  assign("fun",fun,parent.frame())
  assign("args",args,parent.frame())
}
#
#
checkoldcrossbasis <- 
function(argvar,arglag,addarg) {
#
################################################################################
#
  # OLD ARGUMENT type
  if(is.null(argvar$fun)&&!is.null(argvar$type)) {
    names(argvar)[names(argvar)=="type"] <- "fun"
    assign("argvar",argvar,parent.frame())
    warning("argument 'type' replaced by 'fun'. See ?onebasis")
  }
  if(is.null(arglag$fun)&&!is.null(arglag$type)) {
    names(arglag)[names(arglag)=="type"] <- "fun"
    assign("arglag",arglag,parent.frame())
    warning("argument 'type' replaced by 'fun'. See ?onebasis")
  }
  # OLD DEFAULT KNOTS PLACEMENT FOR LAG SPACE
  checklag <- function(fun=NULL,df=NULL,knots=NULL,...) {
    #browser()
    if((is.null(fun)||fun%in%c("ns","bs","strata")) && 
      is.null(knots) && (!is.null(df)&&df>1))
      warning("default knots placement along lags has changed since version 2.0.0.",
        "\n","See 'file.show(system.file('Changesince200',package='dlnm'))'.",
        "\n","See also help(logknots) for setting the knots",
        "\n","consistently with the previous versions")
  }
  do.call(checklag,arglag)
  # 'VERY' OLD USAGE
  if(any(c("vartype","vardf","vardegree","varknots","varbound","varint",
    "cen","cenvalue","maxlag","lagtype","lagdf","lagdegree","lagknots",
    "lagbound","lagint") %in% addarg))
    stop("old usage not allowed any more. See ?crossbasis and ?onebasis")
}


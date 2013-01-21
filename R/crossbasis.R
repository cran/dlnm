###
### R routines for the R package dlnm (c) Antonio Gasparrini 2012-2013
#
`crossbasis` <-
  function(x, lag=c(0,0), argvar=list(), arglag=list(), group=NULL, ...) {         
#
################################################################################
#
  #  lag MUST BE A POSITIVE INTEGER VECTOR 
  if(!is.numeric(lag)||length(lag)>2||any(lag<0)) {
    stop("'lag' must a positive integer vector or scalar")
  }
  if(length(lag)==1L) lag <- c(0L,lag)
  if(diff(lag)<0L) stop("lag[1] must be <= lag[2]")
  lag <- round(lag[1L:2L])
#  
  if(!is.list(argvar)) stop("'var' must be a list")
  if(!is.list(arglag)) stop("'arglag' must be a list")
#
############################################################################
# SET DEFAULT VALUES (NEEDED FOR BASIS FOR THE LAG, UNCENTERED LATER)
#############
#
  argvar <- modifyList(list(type="ns",df=1,degree=1,int=FALSE),argvar)
  arglag  <- modifyList(list(type="ns",df=1,degree=1,int=TRUE),arglag)
#
############################################################################
# ASSURE PORTABILITY OF OLD USAGE
############# 
#
  oldarg <- list(...)
  oldlab <-   c("vartype","vardf","vardegree","varknots","varbound","varint",
    "cen","cenvalue","maxlag","lagtype","lagdf","lagdegree","lagknots",
      "lagbound","lagint")
  if(length(oldarg)>0) {
    # STOP IF NOT IN OLD ARGUMENTS
    if(!all(names(oldarg)%in%oldlab)) stop("unused argument(s)")
    # SET LAG
    if(!is.null(oldarg$maxlag)) lag <- round(c(0,oldarg$maxlag))
    # TRANSFORM CENVALUE
    if(!is.null(oldarg$cen) && oldarg$cen==TRUE) {
      oldarg$cen <- oldarg$cenvalue
    }
    # REPLACE IN argvar AND arglag LISTS
    oldvar <- list(type=oldarg$vartype,df=oldarg$vardf,
      degree=oldarg$vardegree,knots=oldarg$varknots,bound=oldarg$varbound,
      int=oldarg$varint,cen=oldarg$cen)
    argvar <- modifyList(argvar,oldvar[sapply(oldvar,is.null)==FALSE])
    oldlag <- list(type=oldarg$lagtype,df=oldarg$lagdf,
      degree=oldarg$lagdegree,knots=oldarg$lagknots,bound=oldarg$lagbound,
      int=oldarg$lagint)
    arglag <- modifyList(arglag,oldlag[sapply(oldlag,is.null)==FALSE])
    # WARNING
    warning("the function 'crossbasis' has changed arguments since version 1.5.0
  The old arguments are only kept for compatibility reasons.
  The users are strongly suggested to adopt the new usage: see '?crossbasis'")
  }
#
############################################################################
# CREATE THE BASIS FOR THE PREDICTOR SPACE
#############
#
  # x MUST BE A VECTOR OR MATRIX WITH NUMBER OF COLUMNS COMPATIBLE WITH lag
  # IF A VECTOR, x  IS TREATED AS A TIME SERIES
  # OTHERWISE, x IS TREATED AS A MATRIX OF LAGGED OCCURRENCES
  dim <- dim(as.matrix(x))
  if(!dim[2]%in%c(1L,diff(lag)+1L)) {
    stop("ncol(x) must be equal to 1 (if x is a time series vector),
  otherwise to the lag period (for x as a matrix of lagged occurrences)")
  }
#
  # THE BASIS TRANSFORMATION CREATES DIFFERENT MATRICES DEPENDING THE DATA :
  #   IF TIME SERIES, EACH COLUMN CONTAINS THE UNLAGGED TRANSFORMATION
  #   IF NOT, EACH COLUMN CONTAINS THE TRANFORMATION FOR ALL THE LAGGED 
  basisvar <- do.call("onebasis",modifyList(argvar,list(x=x)))
#
############################################################################
# CREATE THE BASIS FOR THE LAG SPACE
#############
#
  # IF LAG=0, THAN TYPE="STRATA" AND DF=1
  if(diff(lag)==0L) {
    arglag <- list(type="strata",df=1,knots=NULL,int=FALSE)
  }
  # SET THE KNOTS (DEFAULT AT EQUALLY-SPACED VALUES IN THE LOG SCALE)
  if(is.null(arglag$knots)&&arglag$df>1+arglag$int&&(arglag$type=="ns")) {
    arglag$knots <- lag[1] + exp(((1+log(diff(lag)))/(arglag$df-arglag$int))*
      1:(arglag$df-arglag$int-1)-1)
  }
  if(is.null(arglag$knots)&&(arglag$type=="bs")) {
    if(arglag$df>arglag$degree+arglag$int) {
      arglag$knots <- lag[1] + exp(((1+log(diff(lag)))/(arglag$df-arglag$int-
        arglag$degree+1))*1:(arglag$df-arglag$int-arglag$degree)-1)
    }
  }
  if(is.null(arglag$knots)&&arglag$df>1+arglag$int&
    arglag$type%in%c("strata","dthr","hthr","lthr")) {
      arglag$knots <- lag[1] + exp(((1+log(diff(lag)))/(arglag$df+1-arglag$int))*
        1:(arglag$df-arglag$int)-1)
  }
#
  basislag <- do.call("onebasis",modifyList(arglag,list(x=.seq(lag),cen=FALSE))) 
#
############################################################################
# CROSSBASIS COMPUTATION
#############
#
  # REDEFINE ARGUMENTS FOR BASES, THEY MIGHT HAVE BEEN CHANGED BY onebasis
  argvar <- attributes(basisvar)[names(attributes(basisvar))%in%
    c("type","df","degree","knots","bound","int","cen")]
  arglag <- attributes(basislag)[names(attributes(basislag))%in%
    c("type","df","degree","knots","bound","int","cen")]
#
  # COMPUTE CROSS-BASIS:
  #   FOR TIME SERIES DATA, COMPUTE THE MATRIX OF LAGGED OCCURRENCES FIRST
  #   IF x WAS ALREADY A MATRIX, JUST RECOMPUTE THE APPROPRIATE DIMENSIONS
  crossbasis <- matrix(0,nrow=dim[1],ncol=argvar$df*arglag$df)
  for(v in seq(length=argvar$df)) {
    if(dim[2]==1L) {
      mat <- as.matrix(Lag(basisvar[, v],.seq(lag)))
    } else mat <- matrix(basisvar[,v],ncol=diff(lag)+1)
    for(l in seq(length=arglag$df)) {
      crossbasis[,argvar$df*(l-1)+v] <- mat%*%(basislag[,l])
    }
  }
#
  # NAMES TO THE NEW CROSS-VARIABLES
  colnames(crossbasis) <- as.vector(outer(paste("v",seq(length=argvar$df),
    sep=""),paste("l",seq(length=arglag$df),sep=""),paste,sep="."))
#
############################################################################
# GROUP
#############
#
  # IF GROUP, JUST SET TO NA ALL THE FIRST  
  if(!is.null(group)) {
    if(dim[2]>1L) stop("'group' allowed only for time series data")
    # FORCE GROUP TO FACTOR, ELIMINATE UNUSED LEVELS
    group <- factor(group)
    # FIRST COHERENCE CHECKS
    if(any(is.na(group))) stop("missing values in 'group' are not allowed")
    if(length(group)!=length(x)) {
      stop("length(group) must be equal to length(x)")
    }
    if(sum(diff(as.numeric(group))!=0)!=length(levels(group))-1) {
      stop("series defined by 'group' must be consecutive")
    }
    if(min(tapply(x,group,length))<=argvar$df) {
      stop("each group must have length > vardf")
    }
    if(min(tapply(x,group,length))<=diff(lag)) {
      stop("each group must have length > diff(lag)")
    }
    # SET TO NA ALL THE FIRST MAXLAG OBS FOR EACH GROUP
    crossbasis[sequence(tapply(seq(length=nrow(crossbasis)),
      group,length))%in%seq(length=lag[2]),] <- NA
  }
#
############################################################################
# ATTRIBUTES
############
#
  attributes(crossbasis) <- c(attributes(crossbasis),list(
    range = attributes(basisvar)$range,
    lag = lag,
    argvar = argvar,
    arglag = arglag
  ))
  if(!is.null(group)) attributes(crossbasis)$group <- length(unique(group))
#
  class(crossbasis) <- c("crossbasis","matrix")
#
  return(crossbasis)
}


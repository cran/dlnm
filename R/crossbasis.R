###
### R routines for the R package dlnm (c) Antonio Gasparrini 2012-2014
#
crossbasis <-
  function(x, lag, argvar=list(), arglag=list(), group=NULL, ...) {         
#
################################################################################
# COHERENCE CHECKS
#
  # OLD USAGE
  checkoldcrossbasis(argvar,arglag,list(...))
#
  #  lag MUST BE A POSITIVE INTEGER VECTOR 
  lag <- if(missing(lag)) c(0,NCOL(x)-1) else mklag(lag)
#  
  if(!is.list(argvar)) stop("'var' must be a list")
  if(!is.list(arglag)) stop("'arglag' must be a list")
#
############################################################################
# CREATE THE BASIS FOR THE PREDICTOR SPACE
#
  # x MUST BE A VECTOR OR MATRIX WITH NUMBER OF COLUMNS COMPATIBLE WITH lag
  # IF A VECTOR, x  IS TREATED AS A TIME SERIES
  # OTHERWISE, x IS TREATED AS A MATRIX OF LAGGED OCCURRENCES
  x <- as.matrix(x)
  dim <- dim(x)
  if(!dim[2]%in%c(1L,diff(lag)+1L)) stop("NCOL(x) must be equal to 1 (if x is ",
    "a time series vector), otherwise to the lag period (for x as a matrix of ",
    "lagged occurrences)")
#
  # THE BASIS TRANSFORMATION CREATES DIFFERENT MATRICES DEPENDING THE DATA :
  #   IF TIME SERIES, EACH COLUMN CONTAINS THE UNLAGGED TRANSFORMATION
  #   IF NOT, EACH COLUMN CONTAINS THE TRANFORMATION FOR ALL THE LAGGED 
  basisvar <- do.call("onebasis",modifyList(argvar,list(x=x)))
#
############################################################################
# CREATE THE BASIS FOR THE LAG SPACE
#
  # SET FUN="STRATA" AND DF=1 UNDER SPECIFIC CIRCUMSTANCES
  if(length(arglag)==0L || diff(lag)==0L) 
    arglag <- list(fun="strata",df=1,int=TRUE)
#
  # IF NOT SPECIFIED AND AN ARGUMENT, INCLUDE AN INTERCEPT BY DEFAULT
  if(is.null(arglag$int)) arglag$int <- TRUE
#
  basislag <- do.call("onebasis",modifyList(arglag,list(x=seqlag(lag),cen=FALSE))) 
#
############################################################################
# CROSSBASIS COMPUTATION
#
  # GROUP
  if(!is.null(group)) checkgroup(group,x,basisvar,lag)
#
  # COMPUTE CROSS-BASIS:
  #   FOR TIME SERIES DATA, COMPUTE THE MATRIX OF LAGGED OCCURRENCES FIRST
  #   IF x WAS ALREADY A MATRIX, JUST RECOMPUTE THE APPROPRIATE DIMENSIONS
  crossbasis <- matrix(0,nrow=dim[1],ncol=ncol(basisvar)*ncol(basislag))
  for(v in seq(length=ncol(basisvar))) {
    if(dim[2]==1L) {
      mat <- as.matrix(Lag2(basisvar[, v],seqlag(lag),group=group))
    } else mat <- matrix(basisvar[,v],ncol=diff(lag)+1)
    for(l in seq(length=ncol(basislag))) {
      crossbasis[,ncol(basisvar)*(l-1)+v] <- mat%*%(basislag[,l])
    }
  }
#
############################################################################
# ATTRIBUTES AND NAMES
#
  # NAMES
  cn <- as.vector(outer(paste("v",seq(length=ncol(basisvar)),
    sep=""),paste("l",seq(length=ncol(basislag)),sep=""),paste,sep="."))
  dimnames(crossbasis) <- list(rownames(x),cn)
#
  # REDEFINE ARGUMENTS FOR BASES, THEY MIGHT HAVE BEEN CHANGED BY onebasis
  ind <- match(names(formals(attributes(basisvar)$fun)),
    names(attributes(basisvar)),nomatch=0)
  argvar <- c(attributes(basisvar)[c("fun","cen")],attributes(basisvar)[ind])
  ind <- match(names(formals(attributes(basislag)$fun)),
    names(attributes(basislag)),nomatch=0)
  arglag <- c(attributes(basislag)[c("fun","cen")],attributes(basislag)[ind])
#
#
  attributes(crossbasis) <- c(attributes(crossbasis),
    list(df=c(ncol(basisvar),ncol(basislag)),range=range(x,na.rm=T),lag=lag,
      argvar=argvar,arglag=arglag))
  if(!is.null(group)) attributes(crossbasis)$group <- length(unique(group))
#
  class(crossbasis) <- c("crossbasis","matrix")
#
  return(crossbasis)
}


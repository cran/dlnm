###
### R routines for the R package dlnm (c) Antonio Gasparrini 2012-2014
#
crosspred <-
function(basis, model=NULL, coef=NULL, vcov=NULL, model.link=NULL, at=NULL,
  from=NULL, to=NULL, by=NULL, lag, bylag=1, ci.level=0.95, cumul=FALSE) {
#
################################################################################
#
  list <- vector("list",0)
  name <- deparse(substitute(basis))
#
###########################################################################
# CHECK CROSSBASIS AND WRITE CONDITION (REGULAR EXPRESSION) TO EXTRACT COEF-VCOV
#
  # IF SIMPLE BASIS, THEN RE-CREATE ATTRIBUTES
  attr <- attributes(basis)
  if(any(class(basis)=="onebasis")) {
    ind <- match(names(formals(attr$fun)),names(attr),nomatch=0)
    attr <- list(range=attr$range,lag=c(0,0),argvar=c(attr[c("fun","cen")],
      attr[ind]),arglag=list(fun="strata",df=1,int=TRUE))
    cond <- paste(name,"[[:print:]]*b[0-9]{1,2}",sep="")
  } else cond <- paste(name,"[[:print:]]*v[0-9]{1,2}\\.l[0-9]{1,2}",sep="")
  if(ncol(basis)==1L) cond <- name
#
  if(!any(class(basis)%in%c("crossbasis","onebasis")))
    stop("the first argument must be an object of class 'crossbasis' or 'onebasis'")
#
###########################################################################
# OTHER COHERENCE CHECKS
#
  if(is.null(model)&&(is.null(coef)||is.null(vcov)))
    stop("At least 'model' or 'coef'-'vcov' must be provided")
  if(!is.numeric(ci.level)||ci.level>=1||ci.level<=0)
    stop("'ci.level' must be numeric and between 0 and 1")
#
  #  lag MUST BE A POSITIVE INTEGER VECTOR, BY DEFAULT THAT USED FOR ESTIMATION
  lag <- if(missing(lag)) attr$lag else mklag(lag)
  if(lag!=attr$lag && attr$arglag$fun=="integer")
      stop("prediction for lag sub-period not allowed for type 'integer'")
#
  # CUMULATIVE EFFECTS ONLY WITH LAGGED EFFECTS AND lag[1]==0
  if(cumul==TRUE && (diff(lag)==0L || lag[1]!=0L)) {
    cumul <- FALSE
    warning("Cumulative predictions only computed if diff(lag)>0 and lag[1]=0")
  }
#
###########################################################################
# SET COEF, VCOV CLASS AND LINK FOR EVERY TYPE OF MODELS
#
  # IF MODEL PROVIDED, EXTRACT FROM HERE, OTHERWISE DIRECTLY FROM COEF AND VCOV
  if(!is.null(model)) {
    model.class <- class(model)
    coef <- getcoef(model,model.class,cond)
    vcov <- getvcov(model,model.class,cond)
    model.link <- getlink(model,model.class)
  } else model.class <- NA
#
  # CHECK COEF AND VCOV
  if(length(coef)!=ncol(basis) || length(coef)!=dim(vcov)[1] ||
    any(is.na(coef))|| any(is.na(vcov))) {
    stop("number of estimated parameters does not match number of cross-functions.
  Possible reasons:
  1) 'model' and 'basis' objects do not match
  2) wrong 'coef' or 'vcov' arguments or methods
  3) model dropped some cross-functions because of collinearity (set to NA)
  4) name of basis matrix matches other parameters in the model formula
  Unlikely, but in this case change the name of the onebasis-crossbasis object")
  }
#
##########################################################################
# PREDVAR
#
  # SET PREDVAR FROM AT, FROM/TO/BY OR AUTOMATICALLY
  # IF at IS A MATRIX, PREDVAR IS ITS ROWNAMES (IF ANY) OR DEFAULT
  if(is.null(at)) {
    if(is.null(from)) from <- attr$range[1]
    if(is.null(to)) to <- attr$range[2]
    nobs <- ifelse(is.null(by),50,max(1,diff(attr$range)/by))
    pretty <- pretty(c(from,to),n=nobs)
    pretty <- pretty[pretty>=from&pretty<=to]	
    predvar <- if(is.null(by)) pretty else seq(from=min(pretty),
      to=to,by=by)
  } else if(is.matrix(at)) {
    if(dim(at)[2]!=diff(lag)+1L)
      stop("matrix in 'at' must have ncol=diff(lag)+1")
    if(!is.null(rownames(at))) {
      predvar <- as.numeric(rownames(at))
      if(any(is.na(predvar))) stop("rownames(at) must represent numurical values")
    } else predvar <- seq(nrow(at))
  } else predvar <- sort(unique(at))
#
##########################################################################
# PREDICTION OF LAG-SPECIFIC EFFECTS
#####################################
#
  # PREDLAG
  predlag <- seqlag(lag,bylag)
  # EXPAND at IF A VECTOR HAS BEEN PROVIDED
  if(is.matrix(at)) {
    matlag <- at
  } else matlag <- matrix(rep(predvar,length(predlag)),nrow=length(predvar))
  # CREATE THE BASIS FOR VAR AND LAG
  basisvar <- do.call("onebasis",c(list(x=matlag),attr$argvar))
  basislag <- do.call("onebasis",c(list(x=predlag),attr$arglag))
#
  # CREATE LAG-SPECIFIC EFFECTS AND SE + (OPTIONAL) CUMULATIVE
  matfit <- matse <- matrix(0,length(predvar),length(predlag))
  for(i in seq(length(predlag))) {
    ind <- seq(length(predvar))+length(predvar)*(i-1)
    ZM <- basislag[i,,drop=FALSE] %x% basisvar[ind,,drop=FALSE]
    matfit[,i] <- ZM %*% coef
    matse[,i] <- sqrt(diag(ZM %*% vcov %*% t(ZM)))
  }
#
  # NAMES
  rownames(matfit) <- rownames(matse) <- predvar
  colnames(matfit) <- colnames(matse) <- outer("lag",seqlag(lag,bylag),
    paste,sep="")
#
##########################################################################
# PREDICTION OF OVERALL+CUMULATIVE EFFECTS
#####################################
#
  # PREDLAG
  predlag <- seqlag(lag)
  # EXPAND at IF A VECTOR HAS BEEN PROVIDED
  if(is.matrix(at)) {
    matlag <- at
  } else matlag <- t(rep(1,length(predlag)))%x%predvar 
  # CREATE THE BASIS FOR VAR AND LAG
  basisvar <- do.call("onebasis",c(list(x=matlag),attr$argvar))
  basislag <- do.call("onebasis",c(list(x=seqlag(lag)),attr$arglag))
#
  # CREATE OVERALL AND (OPTIONAL) CUMULATIVE EFFECTS AND SE
  ZMall <- 0
  if(cumul) {
    cumfit <- cumse <- matrix(0,length(predvar),diff(lag)+1)
  }
  for (i in 1:(diff(lag)+1)) {
    ind <- seq(length(predvar))+length(predvar)*(i-1)
    ZMall <- ZMall + basislag[i,,drop=FALSE] %x% basisvar[ind,,drop=FALSE]
    if(cumul) {
      cumfit[, i] <- ZMall %*% coef
      cumse[, i] <- sqrt(diag(ZMall %*% vcov %*% t(ZMall)))
    }
  }
  allfit <- as.vector(ZMall %*% coef)
  allse <- sqrt(diag(ZMall %*% vcov %*% t(ZMall)))
#
  # NAMES
  names(allfit) <- names(allse) <- predvar
  if(cumul) {
    rownames(cumfit) <- rownames(cumse) <- predvar
    colnames(cumfit) <- colnames(cumse) <- outer("lag",seqlag(lag),paste,sep="")
  }
#
###########################################################################
#
  list$predvar <- predvar
  list$lag <- lag
  list$bylag <- bylag
  list$coefficients <- coef
  list$vcov <- vcov
  list$matfit <- matfit
  list$matse <- matse
  list$allfit <- allfit
  list$allse <- allse
  if(cumul) {
    list$cumfit <- cumfit
    list$cumse <- cumse
  }
#
  # MATRICES AND VECTORS WITH EXPONENTIATED EFFECTS AND CONFIDENCE INTERVALS
  z <- qnorm(1-(1-ci.level)/2)
  if(!is.null(model.link) && model.link %in% c("log","logit")) {
    list$matRRfit <- exp(matfit)
    list$matRRlow <- exp(matfit-z*matse)
    list$matRRhigh <- exp(matfit+z*matse)
    list$allRRfit <- exp(allfit)
    list$allRRlow <- exp(allfit-z*allse)
    names(list$allRRlow) <- names(allfit)
    list$allRRhigh <- exp(allfit+z*allse)
    names(list$allRRhigh) <- names(allfit)
    if(cumul) {
      list$cumRRfit <- exp(cumfit)
      list$cumRRlow <- exp(cumfit-z*cumse)
      list$cumRRhigh <- exp(cumfit+z*cumse)
    }
  } else {
    list$matlow <- matfit-z*matse
    list$mathigh <- matfit+z*matse
    list$alllow <- allfit-z*allse
    names(list$alllow) <- names(allfit)
    list$allhigh <- allfit+z*allse
    names(list$allhigh) <- names(allfit)
    if(cumul) {
      list$cumlow <- cumfit-z*cumse
      list$cumhigh <- cumfit+z*cumse
    }
  }
#
  list$ci.level <- ci.level
  list$model.class <- model.class
  list$model.link <- model.link
#
  class(list) <- "crosspred"
#
  return(list)
}


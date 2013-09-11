###
### R routines for the R package dlnm (c) Antonio Gasparrini 2012-2013
#
`crossreduce` <-
function(basis, model=NULL, type="overall", value=NULL, coef=NULL, vcov=NULL,
  model.link=NULL, at=NULL, from=NULL, to=NULL, by=NULL, lag, bylag=1,
  ci.level=0.95) {
#
################################################################################
#
  list <- vector("list",0)
  name <- deparse(substitute(basis))
#
###########################################################################
# CHECK BASIS AND WRITE CONDITION (REGULAR EXPRESSION) TO EXTRACT COEF-VCOV
#
  if(all(class(basis)!="crossbasis")) {
    stop("the first argument must be an object of class 'crossbasis'")
  }
  attr <- attributes(basis)
  cond <- paste(name,"[[:print:]]*v[0-9]{1,2}\\.l[0-9]{1,2}",sep="")
  if(ncol(basis)==1) cond <- name
#
###########################################################################
# COHERENCE CHECKS (SEE CROSSPRED)
#
  if(is.null(model)&&(is.null(coef)||is.null(vcov))) {
    stop("At least 'model' or 'coef'-'vcov' must be provided")
  }
  type <- match.arg(type,c("overall","var","lag"))
  if(type!="overall") {
    if(is.null(value)) stop("'value' must be provided for type 'var' or 'lag'")
    else if(!is.numeric(value)||length(value)>1) {
      stop("'value' must be a numeric scalar")
    }
    if(type=="lag" && (any(value<attr$lag[1]) ||any(value>attr$lag[2]))) {
      stop("'value' of lag-specific effects must be within the lag range")
    }
  } else value <- NULL
#
  #  lag MUST BE A POSITIVE INTEGER VECTOR, BY DEFAULT THAT USED FOR ESTIMATION
  lag <- if(missing(lag)) attr$lag else .mklag(lag)
  if(lag!=attr$lag && attr$arglag$fun=="integer")
      stop("prediction for lag sub-period not allowed for type 'integer'")
#
  if(!is.numeric(ci.level)||ci.level>=1||ci.level<=0) {
    stop("'ci.level' must be numeric and between 0 and 1")
  }
#
###########################################################################
# SET COEF, VCOV CLASS AND LINK FOR EVERY TYPE OF MODELS
#
  # IF MODEL PROVIDED, EXTRACT FROM HERE, OTHERWISE DIRECTLY FROM COEF AND VCOV
  if(!is.null(model)) {
    model.class <- class(model)
    coef <- .getcoef(model,model.class,cond)
    vcov <- .getvcov(model,model.class,cond)
    model.link <- .getlink(model,model.class)
  } else {
    model.class <- NA
    model.link <- NA
  }
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
  if(is.null(at)) {
    if(is.null(from)) from <- attr$range[1]
    if(is.null(to)) to <- attr$range[2]
    pretty <- pretty(c(from,to),n=50,min.n=30)
    pretty <- pretty[pretty>=from&pretty<=to]  
    if(is.null(by)) {
      predvar <- pretty
    } else predvar <- seq(from=min(pretty),to=max(pretty),by=by)
  } else predvar <- sort(unique(at))
#
##########################################################################
# REDUCTION
#
  # CREATE TRANSFORMATION MATRIX AND BASIS
  if(type=="overall") {
    lagbasis <- do.call("onebasis",c(list(x=.seq(lag)),attr$arglag))
    M <- (t(rep(1,diff(lag)+1)) %*% lagbasis) %x% 
      diag(ncol(basis)/ncol(lagbasis))
    newbasis <- do.call("onebasis",c(list(x=predvar),attr$argvar))
  }else if(type=="lag") {
    lagbasis <- do.call("onebasis",c(list(x=value),attr$arglag))
    M <- lagbasis %x% diag(ncol(basis)/ncol(lagbasis))
    newbasis <- do.call("onebasis",c(list(x=predvar),attr$argvar))
  } else if(type=="var") {
    varbasis <- do.call("onebasis",c(list(x=value),attr$argvar))
    M <- diag(ncol(basis)/ncol(varbasis)) %x% varbasis
    newbasis <- do.call("onebasis",c(list(x=.seq(lag,bylag)),attr$arglag))
  }
#
  # CREATE NEW SET OF COEF AND VCOV
  newcoef <- as.vector(M%*%coef)
  names(newcoef) <- colnames(newbasis)
  newvcov <- M%*%vcov%*%t(M)
  dimnames(newvcov) <- list(colnames(newbasis),colnames(newbasis))
#
##########################################################################
# PREDICTION
#
  fit <- as.vector(newbasis%*%newcoef)
  se <- sqrt(diag(newbasis%*%newvcov%*%t(newbasis)))
  if(type=="var") {
    names(fit) <- names(se) <- outer("lag",.seq(lag,bylag),paste,sep="")
  }else names(fit) <- names(se) <- predvar
#
###########################################################################
#
  list$coefficients <- newcoef
  list$vcov <- newvcov
  list$basis <- newbasis
  list$type <- type
  list$value <- value
  if(type!="var") list$predvar <- predvar
  list$lag <- lag
  list$bylag <- bylag
  list$fit <- fit
  list$se <- se
#
  # VECTORS WITH EXPONENTIATED EFFECTS AND CONFIDENCE INTERVALS
  z <- qnorm(1-(1-ci.level)/2)
  if(model.link %in% c("log","logit")) {
    list$RRfit <- exp(fit)
    list$RRlow <- exp(fit-z*se)
    names(list$RRlow) <- names(fit)
    list$RRhigh <- exp(fit+z*se)
    names(list$RRhigh) <- names(fit)
  } else {
    list$low <- fit-z*se
    names(list$low) <- names(fit)
    list$high <- fit+z*se
    names(list$high) <- names(fit)
  }
#
  list$ci.level <- ci.level
  list$model.class <- model.class
  list$model.link <- model.link
#
  class(list) <- "crossreduce"
#
  return(list)
}

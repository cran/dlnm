`crossreduce` <-
function(basis, model=NULL, type="overall", value=NULL, coef=NULL, vcov=NULL,
  model.link=NULL, at=NULL, from=NULL, to=NULL, by=NULL, lag, bylag=1,
  ci.level=0.95) {

list <- vector("list",0)
name <- deparse(substitute(basis))

###########################################################################
# CHECK BASIS AND WRITE CONDITION (REGULAR EXPRESSION) TO EXTRACT COEF-VCOV

if(all(class(basis)!="crossbasis")) {
  stop("the first argument must be an object of class 'crossbasis'")
}
attr <- attributes(basis)
cond <- paste(name,"[[:print:]]*v[0-9]{1,2}\\.l[0-9]{1,2}",sep="")
if(ncol(basis)==1) cond <- name

###########################################################################
# COHERENCE CHECKS (SEE CROSSPRED)

if(is.null(model)&&(is.null(coef)||is.null(vcov))) {
  stop("At least 'model' or 'coef'-'vcov' must be provided")
}
if(!is.null(vcov)&&!is.numeric(vcov)) stop("'vcov' must be numeric")
if(!is.null(coef)&&!is.numeric(coef)) stop("'coef' must be numeric")

if(!type%in%c("overall","var","lag")) {
  stop("'type' must be one of 'overall', 'var' or lag'")
}
if(type!="overall") {
  if(is.null(value)) stop("'value' must be provided for type 'var' or 'lag'")
  else if(!is.numeric(value)||length(value)>1) {
    stop("'value' must be a numeric scalar")
  }
  if(type=="lag" && (any(value<attr$lag[1]) ||any(value>attr$lag[2]))) {
    stop("'value' of lag-specific effects must be within the lag range")
  }
} else value <- NULL

#  lag MUST BE A POSITIVE INTEGER VECTOR, BY DEFAULT THAT USED FOR ESTIMATION
if(missing(lag)) {
  lag <- attr$lag
} else {
  if(lag!=attr$lag && attr$arglag$type=="integer") {
    stop("prediction for lag sub-period not allowed for type 'integer'")
  }
  if(!is.numeric(lag)||length(lag)>2||any(lag<0)) {
    stop("'lag' must a positive integer vector or scalar")
  }
  if(length(lag)==1L) lag <- c(0L,lag)
  if(diff(lag)<0L) stop("lag[1] must be <= lag[2]")
  lag <- round(lag[1L:2L])
}

if(!is.null(at)&&!is.numeric(at)) stop("'at' must be numeric")
if(!is.null(from)&&!is.numeric(from)) stop("'from' must be numeric")
if(!is.null(to)&&!is.numeric(to)) stop("'to' must be numeric")
if(!is.null(by)&&!is.numeric(by)) stop("'by' must be numeric")
if(!is.numeric(ci.level)||ci.level>=1||ci.level<=0) {
  stop("'ci.level' must be numeric and between 0 and 1")
}

###########################################################################
# EXTRACT ORIGINAL PARAMETERS FROM FITTED MODEL

# IF MODEL PROVIDED, EXTRACT FROM HERE, OTHERWISE DIRECTLY FROM COEF AND VCOV
model.class <- NA
if(!is.null(model)) {

  # EXTRACT COEF AND VCOV FROM THE MODEL, IF METHODS EXIST
  model.class <- class(model)
  modelcoef <- tryCatch(coef(model),error=function(w) "error")
  modelvcov <- tryCatch(vcov(model),error=function(w) "error")
  if(identical(modelcoef,"error")||identical(modelvcov,"error")) {
    stop("methods for coef() and vcov() must exist for the class
of object 'model'. If not, extract them manually and use 
the argumetns 'coef' and 'vcov'")
  }

  # SELECT COEF-VCOV
  indcoef <- grep(cond,names(modelcoef))
  coef <- modelcoef[indcoef]
  indvcov <- grep(cond,dimnames(modelvcov)[[1]])
  vcov <- modelvcov[indvcov,indvcov,drop=FALSE]
}

# CHECK COEF AND VCOV
if(length(coef)!=attr$argvar$df*attr$arglag$df || length(coef)!=dim(vcov)[1] ||
  any(is.na(coef))|| any(is.na(vcov))) {
  stop("number of estimated parameters does not match number of cross-functions.
Possible reasons:
1) 'model' and 'basis' objects do not match
2) wrong 'coef' or 'vcov' arguments
3) model dropped some cross-functions because of collinearity (set to NA)
4) name of crossbasis matrix matches other parameters in the model formula
Unlikely, but in this case change the name of the crossbasis object")
}

###########################################################################
# MODEL LINK

if(all(model.class %in% c("lm"))) {
  model.link <- "identity"
} else if(any(model.class %in% c("clogit"))) {
  model.link <- "logit"
} else if(all(model.class %in% c("coxph"))) {
  model.link <- "log"
} else if(any(model.class %in% c("glm"))) {
  model.link <- model$family$link
} else if(is.null(model.link)) model.link <- NA


##########################################################################
# PREDVAR

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

##########################################################################
# REDUCTION

# CREATE TRANSFORMATION MATRIX AND BASIS
if(type=="overall") {
  lagbasis <- do.call("onebasis",c(list(x=.seq(lag)),attr$arglag))
  M <- (t(rep(1,diff(attr$lag)+1)) %*% lagbasis) %x% diag(attr$argvar$df)
  newbasis <- do.call("onebasis",c(list(x=predvar),attr$argvar))
}else if(type=="lag") {
  lagbasis <- do.call("onebasis",c(list(x=value),attr$arglag))
  M <- lagbasis %x% diag(attr$argvar$df)
  newbasis <- do.call("onebasis",c(list(x=predvar),attr$argvar))
} else if(type=="var") {
  varbasis <- do.call("onebasis",c(list(x=value),attr$argvar))
  M <- diag(attr$arglag$df) %x% varbasis
  newbasis <- do.call("onebasis",c(list(x=.seq(lag,bylag)),attr$arglag))
}

# CREATE NEW SET OF COEF AND VCOV
newcoef <- as.vector(M%*%coef)
names(newcoef) <- colnames(newbasis)
newvcov <- M%*%vcov%*%t(M)
dimnames(newvcov) <- list(colnames(newbasis),colnames(newbasis))
      
##########################################################################
# PREDICTION

fit <- as.vector(newbasis%*%newcoef)
se <- sqrt(diag(newbasis%*%newvcov%*%t(newbasis)))
if(type=="var") {
  names(fit) <- names(se) <- outer("lag",.seq(lag,bylag),paste,sep="")
}else names(fit) <- names(se) <- predvar

###########################################################################

list$newcoef <- newcoef
list$newvcov <- newvcov
list$newbasis <- newbasis
list$type <- type
list$value <- value
if(type!="var") list$predvar <- predvar
list$lag <- lag
list$bylag <- bylag
list$fit <- fit
list$se <- se

# VECTORS WITH EXPONENTIATED EFFECTS AND CONFIDENCE INTERVALS
z <- qnorm(1-(1-ci.level)/2)
if(model.link %in% c("log","logit")) {
  list$RRfit <- exp(fit)
  list$RRlow <- exp(fit-z*se)
  names(list$RRlow) <- names(fit)
  list$RRhigh <- exp(fit+z*se)
  names(list$RRhigh) <- names(fit)
}

list$ci.level <- ci.level
list$model.class <- model.class
list$model.link <- model.link

class(list) <- "crossreduce"
return(list)
}

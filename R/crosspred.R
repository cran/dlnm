`crosspred` <-
function(basis, model=NULL, coef=NULL, vcov=NULL, model.link=NULL,
  at=NULL, from=NULL, to=NULL, by=NULL,  ci.level=0.95, cumul=FALSE) {

list <- vector("list",0)
name <- deparse(substitute(basis))

###########################################################################
# CHECK CROSSBASIS AND WRITE CONDITION (REGULAR EXPRESSION) TO EXTRACT COEF-VCOV

# IF SIMPLE BASIS, THEN RE-CREATE ATTRIBUTES
if(any(class(basis)=="onebasis")) {
  attr <- attributes(basis)
  attr <- list(range=attr$range,lag=c(0,0),argvar=list(type=attr$type,df=attr$df,
    degree=attr$degree,knots=attr$knots,bound=attr$bound,int=attr$int,
    cen=attr$cen),arglag=list(type="strata",df=1,int=FALSE))
  cond <- paste(name,"[[:print:]]*b[0-9]{1,2}",sep="")
} else {
  attr <- attributes(basis)
  cond <- paste(name,"[[:print:]]*v[0-9]{1,2}\\.l[0-9]{1,2}",sep="")
}
if(ncol(basis)==1) cond <- name

if(!any(class(basis)%in%c("crossbasis","onebasis"))) {
  stop("the first argument must be an object of class 'crossbasis' or 'onebasis'")
}

###########################################################################
# OTHER COHERENCE CHECKS

if(is.null(model)&&(is.null(coef)||is.null(vcov))) {
  stop("At least 'model' or 'coef'-'vcov' must be provided")
}
if(!is.null(vcov)&&!is.numeric(vcov)) stop("'vcov' must be numeric")
if(!is.null(coef)&&!is.numeric(coef)) stop("'coef' must be numeric")

if(!is.null(at)&&!is.numeric(at)) stop("'at' must be numeric")
if(!is.null(from)&&!is.numeric(from)) stop("'from' must be numeric")
if(!is.null(to)&&!is.numeric(to)) stop("'to' must be numeric")
if(!is.null(by)&&!is.numeric(by)) stop("'by' must be numeric")
if(!is.numeric(ci.level)||ci.level>=1||ci.level<=0) {
  stop("'ci.level' must be numeric and between 0 and 1")
}
# CUMULATIVE EFFECTS ONLY WITH LAGGED EFFECTS
if(diff(attr$lag)==0L) cumul <- FALSE

###########################################################################
# SET COEF, VCOV AND LINK FOR EVERY TYPE OF MODELS

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
4) name of basis matrix matches other parameters in the model formula
Unlikely, but in this case change the name of the onebasis-crossbasis object")
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
# IF at IS A MATRIX, PREDVAR IS ITS ROWNAMES (IF ANY) OR DEFAULT
if(is.null(at)) {
  if(is.null(from)) from <- attr$range[1]
  if(is.null(to)) to <- attr$range[2]
  pretty <- pretty(c(from,to),n=50,min.n=30)
  pretty <- pretty[pretty>=from&pretty<=to]	
  if(is.null(by)) {
    predvar <- pretty
  } else predvar <- seq(from=min(pretty),to=max(pretty),by=by)
} else if(is.matrix(at)) {
  if(dim(at)[2]!=diff(attr$lag)+1L) stop("matrix in 'at' must have ncol=diff(lag)+1")
  if(!is.null(rownames(at))||is.numeric(rownames(at))) {
    predvar <- rownames(at)
  } else predvar <- seq(nrow(at))
} else predvar <- sort(unique(at))

##########################################################################
# PREDICTION
#############

# EXPAND predlag IF A VECTOR HAS BEEN PROVIDED
if(is.matrix(at)) {
  matlag <- at
} else matlag <- t(rep(1,diff(attr$lag)+1))%x%predvar 
# CREATE THE BASIS FOR VAR AND LAG
basisvar <- do.call("onebasis",c(list(x=matlag),attr$argvar))
basislag <- do.call("onebasis",c(list(x=.seq(attr$lag)),attr$arglag))

# CREATE OVERAL, LAG-SPECIFIC AND (OPTIONAL) CUMULATIVE EFFECTS AND SE
matfit <- matse <- matrix(0,length(predvar),diff(attr$lag)+1)
ZMall <- 0
if(cumul) {
  cumfit <- cumse <- matrix(0,length(predvar),diff(attr$lag)+1)
  ZMcum <- 0
}
for (i in 1:(diff(attr$lag)+1)) {
  ind <- seq(length(predvar))+length(predvar)*(i-1)
  ZM <- basislag[i,,drop=FALSE] %x% basisvar[ind,,drop=FALSE]
  matfit[,i] <- ZM %*% coef
  matse[,i] <- sqrt(diag(ZM %*% vcov %*% t(ZM)))
  ZMall <- ZM+ZMall
  if(cumul) {
    cumfit[, i] <- ZMall %*% coef
    cumse[, i] <- sqrt(diag(ZMall %*% vcov %*% t(ZMall)))
  }
}
allfit <- as.vector(ZMall %*% coef)
allse <- sqrt(diag(ZMall %*% vcov %*% t(ZMall)))

# NAMES
rownames(matfit) <- rownames(matse) <- predvar
colnames(matfit) <- colnames(matse) <- outer("lag",.seq(attr$lag),paste,sep="")
names(allfit) <- names(allse) <- predvar
if(cumul) {
  rownames(cumfit) <- rownames(cumse) <- predvar
  colnames(cumfit) <- colnames(cumse) <- outer("lag",.seq(attr$lag),paste,sep="")
}

###########################################################################

list$predvar <- predvar
list$lag <- attr$lag
list$coef <- coef
list$vcov <- vcov
list$matfit <- matfit
list$matse <- matse
list$allfit <- allfit
list$allse <- allse
if(cumul) {
  list$cumfit <- cumfit
  list$cumse <- cumse
}

# MATRICES AND VECTORS WITH EXPONENTIATED EFFECTS AND CONFIDENCE INTERVALS
z <- qnorm(1-(1-ci.level)/2)
if(model.link %in% c("log","logit")) {
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
}

list$ci.level <- ci.level
list$model.class <- model.class
list$model.link <- model.link

class(list) <- "crosspred"
return(list)
}


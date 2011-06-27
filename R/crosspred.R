`crosspred` <-
function(basis, model=NULL, model.link=NULL, coef=NULL,
	vcov=NULL, at=NULL, from=NULL, to=NULL, by=NULL,
	ci.level=0.95, cumul=FALSE) {

list <- vector("list",0)

###########################################################################
# COHERENCE CHECKS

if(all(class(basis)!="crossbasis")) {
	stop("the first argument must be an object of class 'crossbasis'")
}
attr <- attributes(basis)

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
if(attr$maxlag==0) cumul <- FALSE

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

	# WRITE THE CONDITION (WITH REGULAR EXPRESSIONS) AND SELECT
	if(ncol(basis)==1) {
		cond <- deparse(substitute(basis))
	} else {
		cond <- paste(deparse(substitute(basis)),
			"[[:print:]]*v[0-9]{1,2}\\.l[0-9]{1,2}",sep="")
	}
	indcoef <- grep(cond,names(modelcoef))
	coef <- modelcoef[indcoef]
	indvcov <- grep(cond,dimnames(modelvcov)[[1]])
	vcov <- modelvcov[indvcov,indvcov,drop=FALSE]
}
if(length(coef)!=attr$crossdf || length(coef)!=dim(vcov)[1] ||
	any(is.na(coef))|| any(is.na(vcov))) {
	stop("number of estimated parameters does not match number of
cross-functions. Possible reasons:
1) model dropped some cross-functions because of collinearity (set to NA)
It may happens when knots specified beyond the range of var/lag
2) wrong 'coef' or 'vcov' arguments
3) name of crossbasis matrix matches other parameters in the model formula
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
if(is.null(from)) from <- attr$range[1]
if(is.null(to)) to <- attr$range[2]
if(is.null(at)) {
	pretty <- pretty(c(from,to),n=50,min.n=30)
	pretty <- pretty[pretty>=from&pretty<=to]	
	if(is.null(by)) {
		predvar <- pretty
	} else predvar <- seq(from=min(pretty),to=max(pretty),by=by)
} else predvar <- sort(unique(at))

if(length(predvar)<attr$vardf+attr$varint) {
	stop("number of predicted values must be >= vardf+varint")
}

##########################################################################
# PREDICTION
#############

maxlag <- attr$maxlag

# CREATE VARBASIS AND LAGBASIS
predvarbasis <- mkbasis(predvar,type=attr$vartype,df=attr$vardf,
	degree=attr$vardegree,knots=attr$varknots,int=attr$varint,
	bound=attr$varbound,cen=attr$cen,cenvalue=attr$cenvalue)$basis
rownames(predvarbasis) <- predvar
lagbasis <- mklagbasis(maxlag=attr$maxlag,type=attr$lagtype,df=attr$lagdf,
	degree=attr$lagdegree,knots=attr$lagknots,
	int=attr$lagint,bound=attr$lagbound)$basis

# CREATE PREDARRAY: DIFFERENTLY FROM CROSSBASIS, VARBASIS ALWAYS THE SAME
predarray <- array(0,dim=c(length(predvar),attr$crossdf,maxlag+1))
for(i in 1:(maxlag+1)) {
	predarray[,,i] <- matrix(outer(predvarbasis,lagbasis[i,],"*"),
		nrow=length(predvar))
}
dimnames(predarray) <- with(attr, list(predvar,colnames(basis),
	rownames(lagbasis)))
predcrossbasis <- apply(predarray,c(1,2),sum)

# CREATE THE MATRIX OF LAG-SPECIFIC EFFECTS AND SE
matfit <- matse <- matrix(0,dim(predarray)[1],dim(predarray)[3])
for (i in 1:(maxlag + 1)) {
	matfit[, i] <- as.matrix(predarray[, , i]) %*% coef
	matse[, i] <- sqrt(diag(as.matrix(predarray[, , i]) %*% vcov %*% 
		t(as.matrix(predarray[, , i]))))
}
rownames(matfit) <- rownames(matse) <- predvar
colnames(matfit) <- colnames(matse) <- dimnames(predarray)[[3]]

# CREATE THE VECTOR OF OVERALL EFFECTS AND SE
allfit <- as.vector(predcrossbasis%*%coef)
allse <- sqrt(diag(predcrossbasis%*%vcov%*%t(predcrossbasis)))
names(allfit) <- names(allse) <- predvar

# MATRICES AND VECTORS FOR CUMULATIVE EFFECTS AND SE
if(cumul==TRUE) {
	cumfit <- cumse <- matrix(0,dim(predarray)[1],dim(predarray)[3])
	for (i in 1:(maxlag + 1)) {
		# THIS WAY, OTHERWISE ARRAY LOSES DIM IF 1 COL
		predcumarray <- array(predarray[,,1:i],c(dim(predarray)[1:2],i))
		predcumbasis <- apply(predcumarray,c(1,2),sum)
		cumfit[,i] <- predcumbasis%*%coef
		cumse[,i] <- sqrt(diag(predcumbasis%*%vcov%*%t(predcumbasis)))
	}
	rownames(cumfit) <- rownames(cumse) <- predvar
	colnames(cumfit) <- colnames(cumse) <- dimnames(predarray)[[3]] 
}

###########################################################################

list$predvar <- predvar
list$maxlag <- maxlag
list$coef <- coef
list$vcov <- vcov
list$matfit <- matfit
list$matse <- matse
list$allfit <- allfit
list$allse <- allse
if(cumul==TRUE) {
	list$cumfit <- cumfit
	list$cumse <- cumse
}

# MATRICES AND VECTORS WITH EXPONENTIATED EFFECTS AND CONFIDENCE INTERVALS
z <- qnorm(1-(1-ci.level)/2)
if(model.link %in% c("log","logit")) {
	list$matRRfit <- exp(matfit)
	list$matRRhigh <- exp(matfit+z*matse)
	list$matRRlow <- exp(matfit-z*matse)
	list$allRRfit <- exp(allfit)
	list$allRRhigh <- exp(allfit+z*allse)
	names(list$allRRhigh) <- names(allfit)
	list$allRRlow <- exp(allfit-z*allse)
	names(list$allRRlow) <- names(allfit)
	if(cumul==TRUE) {
		list$cumRRfit <- exp(cumfit)
		list$cumRRhigh <- exp(cumfit+z*cumse)
		list$cumRRlow <- exp(cumfit-z*cumse)
	}
}

list$ci.level <- ci.level
list$model.class <- model.class
list$model.link <- model.link
class(list) <- "crosspred"

return(list)
}


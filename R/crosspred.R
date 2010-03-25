`crosspred` <-
function(crossbasis, model, at=NULL,
	from=NULL, to=NULL, by=NULL, cumul=FALSE) {
list <- vector("list",0)

if(class(crossbasis)!="crossbasis") {
	stop("the first argument must be an object of class 'crossbasis'")
}
attr <- attributes(crossbasis)
model.class <- class(model)
if(!any(model.class %in% c("glm","lm","clogit","coxph","gam"))) {
	stop("model class must be one of 'glm','lm','clogit','coxph','gam'
crosspred() needs to be modified in order to include other model functions")
}
index <- grep(deparse(substitute(crossbasis)),names(coef(model)),fixed=T)
coef <- coef(model)[index]
vcov <- vcov(model)[index,index]

if(any(model.class %in% c("glm","gam"))) {
	model.link <- model$family$link
} else if(all(model.class %in% "lm")) {
	model.link <- "identity"
} else model.link <- "log"

if(length(coef)!=attr$crossdf | any(is.na(coef))) {
	stop("number of estimated parameters does not match number of cross-functions
Possible reasons:
1) model dropped some cross-functions because of collinearity
It may happens when knots specify all-0 variables for var/lag
2) name of crossbasis matrix matches other parameters in the model formula
In this case change the name of the crossbasis object")
}

##########################################################################
# PREDVAR
#############

if(is.null(from)) from <- attributes(crossbasis)$range[1]
if(is.null(to)) to <- attributes(crossbasis)$range[2]

if(is.null(at)) {
	if(is.null(by)) {
		predvar <- seq(from=from,to=to,length.out=30)
	} else predvar <- seq(from=from,to=to,by=by)
} else predvar <- sort(unique(at))

predvar <- c(attr$range[1],predvar,attr$range[2])

if(length(predvar)<attr$vardf+attr$varint+2) {
	stop("number of predicted values must be > vardf+varint")
}

##########################################################################
# PREDICTION
#############

maxlag <- attr$maxlag
predvarbasis <- mkbasis(predvar,type=attr$vartype,df=attr$vardf,
	degree=attr$vardegree,knots=attr$varknots,int=attr$varint,
	bound=attr$varbound,cen=attr$cen,cenvalue=attr$cenvalue)$basis
rownames(predvarbasis) <- predvar
lagbasis <- mklagbasis(maxlag=attr$maxlag,type=attr$lagtype,df=attr$lagdf,
	degree=attr$lagdegree,knots=attr$lagknots,
	int=attr$lagint,bound=attr$lagbound)$basis
predarray <- array(0,dim=c(length(predvar),attr$crossdf,maxlag+1))
for(i in 1:(maxlag+1)) {
	predarray[,,i] <- matrix(outer(predvarbasis,lagbasis[i,],"*"),
		nrow=length(predvar))
}
dimnames(predarray) <- with(attr, list(predvar,colnames(crossbasis),
	rownames(lagbasis)))
predcrossbasis <- apply(predarray,c(1,2),sum)

matfit <- matse <- matrix(0,dim(predarray)[1],dim(predarray)[3])
for (i in 1:(maxlag + 1)) {
	matfit[, i] <- as.matrix(predarray[, , i]) %*% coef
	matse[, i] <- sqrt(diag(as.matrix(predarray[, , i]) %*% vcov %*% 
		t(as.matrix(predarray[, , i]))))
}
rownames(matfit) <- rownames(matse) <- predvar
colnames(matfit) <- colnames(matse) <- dimnames(predarray)[[3]]

allfit <- predcrossbasis%*%coef
allse <- sqrt(diag(predcrossbasis%*%vcov%*%t(predcrossbasis)))
names(allfit) <- names(allse) <- predvar

if(cumul==TRUE) {
	cumfit <- cumse <- matrix(0,dim(predarray)[1],dim(predarray)[3])
	for (i in 1:(maxlag + 1)) {
		predcumbasis <- apply(predarray[,,1:i],c(1,2),sum)
		cumfit[,i] <- predcumbasis%*%coef
		cumse[,i] <- sqrt(diag(predcumbasis%*%vcov%*%t(predcumbasis)))
	}
	rownames(cumfit) <- rownames(cumse) <- predvar
	colnames(cumfit) <- colnames(cumse) <- dimnames(predarray)[[3]] 
}

###########################################################################

matfit <- matfit[-c(1,length(predvar)),]
matse <- matse[-c(1,length(predvar)),]
allfit <- allfit[-c(1,length(predvar))]
allse <- allse[-c(1,length(predvar))]
if(cumul==TRUE) {
	cumfit <- cumfit[-c(1,length(predvar)),]
	cumse <- cumse[-c(1,length(predvar)),]
}
predvar <- predvar[-c(1,length(predvar))]

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

if(model.link %in% c("log","logit")) {
	list$matRRfit <- exp(matfit)
	list$matRRhigh <- exp(matfit+1.96*matse)
	list$matRRlow <- exp(matfit-1.96*matse)
	list$allRRfit <- exp(allfit)
	list$allRRhigh <- exp(allfit+1.96*allse)
	names(list$allRRhigh) <- names(allfit)
	list$allRRlow <- exp(allfit-1.96*allse)
	names(list$allRRlow) <- names(allfit)
	if(cumul==TRUE) {
		list$cumRRfit <- exp(cumfit)
		list$cumRRhigh <- exp(cumfit+1.96*cumse)
		list$cumRRlow <- exp(cumfit-1.96*cumse)
	}
}
list$model.class <- model.class
list$model.link <- model.link
class(list) <- "crosspred"
list
}


`crosspred` <-
function(crossbasis,model,at=NULL,from=NULL,to=NULL,by=NULL) {
list <- vector("list",0)

attr <- attributes(crossbasis)
index <- grep(deparse(substitute(crossbasis)),
	rownames(summary(model)$coeff),fixed=T)
coef <- summary(model)$coeff[index,1]
vcov <- vcov(model)[index,index]

if(length(coef)!=attr$crossdf) {
	stop("number of estimated parameters does not match number of cross-functions
Possible reasons:
1) model dropped some cross-functions because of collinearity
2) name of crossbasis matrix matches other parameters in the model formula
In this last case change the name of crossbasis object")
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

if(!is.null(attr$varknots)) {
	if(min(predvar)>min(attr$varknots)|max(predvar)<max(attr$varknots)) {
	stop("predicted range should contains original knots")
}}

##########################################################################
# PREDICTION
#############

maxlag <- attr$maxlag
predvarbasis <- mkbasis(predvar,type=attr$vartype,
	df=attr$vardf,knots=attr$varknots,int=attr$varint,bound=attr$varbound,
	cen=attr$cen,cenvalue=attr$cenvalue)$basis
rownames(predvarbasis) <- predvar
lagbasis <- mklagbasis(maxlag=attr$maxlag,type=attr$lagtype,df=attr$lagdf,
	knots=attr$lagknots,int=attr$lagint,bound=attr$lagbound)$basis
predarray <- array(0,dim=c(length(predvar),attr$crossdf,maxlag+1))
for(i in 1:(maxlag+1)) {
	predarray[,,i] <- matrix(outer(predvarbasis,lagbasis[i,],"*"),
		nrow=length(predvar))
}
dimnames(predarray) <- with(attr, list(predvar,colnames(crossbasis),
	rownames(lagbasis)))
predcrossbasis <- apply(predarray,c(1,2),sum)

matfit <- matrix(0,dim(predarray)[1],dim(predarray)[3])
matse <- matrix(0,dim(predarray)[1],dim(predarray)[3])
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

###########################################################################

list$predvar <- predvar
list$maxlag <- maxlag

list$coef <- coef
list$vcov <- vcov

list$matfit <- matfit
list$matse <- matse
list$allfit <- allfit
list$allse <- allse

if(model$family$link %in% c("log","logit")) {
	list$matRRfit <- exp(matfit)
	list$matRRhigh <- exp(matfit+1.96*matse)
	list$matRRlow <- exp(matfit-1.96*matse)
	list$allRRfit <- exp(allfit)
	list$allRRhigh <- exp(allfit+1.96*allse)
	names(list$allRRhigh) <- names(allfit)
	list$allRRlow <- exp(allfit-1.96*allse)
	names(list$allRRlow) <- names(allfit)
}

class(list) <- "crosspred"
list
}


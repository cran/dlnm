`crossbasis` <-
function(x, vartype="ns", vardf=1, vardegree=1, varknots=NULL,
	varbound=range(x), varint=FALSE, cen=TRUE, cenvalue=mean(x),
	maxlag=0, lagtype="ns", lagdf=1, lagdegree=1, lagknots=NULL,
	lagbound=c(0,maxlag), lagint=TRUE, group=NULL) {

############################################################################
# CROSSBASIS 
#############

# CREATE THE BASIS FOR THE PREDICTOR SPACE
basisvar <- mkbasis(x=x,type=vartype,df=vardf,degree=vardegree,
	knots=varknots,int=varint,bound=varbound,cen=cen,cenvalue=cenvalue)

# CREATE THE BASIS FOR THE LAG SPACE
basislag <- mklagbasis(maxlag=maxlag,type=lagtype,df=lagdf,degree=lagdegree,
	knots=lagknots,int=lagint,bound=lagbound)

# CROSSBASIS COMPUTATION
basis <- matrix(0,nrow=length(x),ncol=basisvar$df*basislag$df)
for(v in seq(length=basisvar$df)) {
	mat <- as.matrix(Lag(basisvar$basis[, v],0:maxlag))
	for(l in seq(length=basislag$df)) {
		basis[,basisvar$df*(l-1)+v] <- mat%*%(basislag$basis[,l])
	}
}

# NAMES TO THE NEW CROSS-VARIABLES
colnames(basis) <- outer(paste("v",seq(length=basisvar$df),sep=""),
	paste("l",seq(length=basislag$df),sep=""), function(x,y) paste(x,y,sep="."))

# IF GROUP, JUST SET TO NA ALL THE FIRST  
if(!is.null(group)) {
	# FIRST COHERENCE CHECKS
	if(any(is.na(group))) stop("missing values in 'group' are not allowed")
	if(length(group)!=length(x)) {
		stop("length(group) must be equal to length(x)")
	}
	if(min(tapply(x,group,length))<=basisvar$df) {
		stop("each group must have length > vardf")
	}
	if(min(tapply(x,group,length))<=maxlag) {
		stop("each group must have length > maxlag")
	}
	# SET TO NA ALL THE FIRST MAXLAG OBS FOR EACH GROUP
	basis[sequence(tapply(seq(length=nrow(basis)),
		group,length))%in%seq(length=maxlag),] <- NA
}

############################################################################
# ATTRIBUTES
############

attributes(basis) <- c(attributes(basis),list(
	range = range(x,na.rm=TRUE),
	crossdf = basisvar$df*basislag$df,

	vartype = basisvar$type,
	vardf = basisvar$df,
	vardegree=basisvar$degree,
	varknots = basisvar$knots,
	varbound = basisvar$bound,
	varint = basisvar$int,
	cen = basisvar$cen,
	cenvalue = basisvar$cenvalue,

	maxlag = basislag$maxlag,
	lagtype = basislag$type,
	lagdf = basislag$df,
	lagdegree = basislag$degree,
	lagknots = basislag$knots,
	lagbound = basislag$bound,
	lagint = basislag$int
))
if(!is.null(group)) attributes(basis)$group <- length(unique(group))

class(basis) <- c("crossbasis","matrix")

return(basis)
}


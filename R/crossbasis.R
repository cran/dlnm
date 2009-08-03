`crossbasis` <-
function(var, vartype="ns", vardf=1, vardegree=1, varknots=NULL,
	varbound=range(var), varint=FALSE, cen=TRUE, cenvalue=mean(var),
	maxlag=0, lagtype="ns", lagdf=1, lagdegree=1, lagknots=NULL,
	lagbound=c(0,maxlag), lagint=TRUE) {

############################################################################
# CROSSBASIS 
#############

basisvar <- mkbasis(var=var,type=vartype,df=vardf,degree=vardegree,
	knots=varknots,int=varint,bound=varbound,cen=cen,cenvalue=cenvalue)

basislag <- mklagbasis(maxlag=maxlag,type=lagtype,df=lagdf,degree=lagdegree,
	knots=lagknots,int=lagint,bound=lagbound)

# CROSSBASIS COMPUTATION
basis <- matrix(0,nrow=length(var),ncol=basisvar$df*basislag$df)
for(v in 1:basisvar$df) {
	mat <- as.matrix(Lag(basisvar$basis[,v],0:maxlag))
	for(l in 1:basislag$df) {
		basis[,basisvar$df*(l-1)+v] <- mat%*%(basislag$basis[,l])
	}
}
colnames(basis) <- outer(paste("v",1:basisvar$df,sep=""),
	paste("l",1:basislag$df,sep=""),	function(x,y) paste(x,y,sep="."))

############################################################################
# ATTRIBUTES
############

attributes(basis) <- c(attributes(basis),list(
	range = range(var,na.rm=TRUE),
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
	lagdegree=basislag$degree,
	lagknots = basislag$knots,
	lagbound = basislag$bound,
	lagint = basislag$int
))
class(basis) <- "crossbasis"
basis
}


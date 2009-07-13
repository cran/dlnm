`mklagbasis` <-
function(maxlag=0,type="ns",df=1,knots=NULL,
	bound=c(0,maxlag),int=TRUE) {

# MAXLAG MUST BE >=0
if(maxlag<0) {
	stop("maxlag must be >= 0")
}

# IF MAXLAG=0, THAN TYPE="STRATA" AND DF=1
if(maxlag==0) {
	type <- "strata"
	df=1
}

lag <- 0:maxlag

# DF MUST BE <= LENGTH(LAG)
if(df>1&df>(length(lag)+1)) {
	stop("for df>1,  df for lag must be <= maxlag+1")
}

if(is.null(knots)&df>1+int&(type=="ns"|type=="poly")) {
	knots <- exp(((1+log(maxlag))/(df-int))*1:(df-int-1)-1)
}
if(is.null(knots)&df>1+int&type=="strata") {
	knots <- exp(((1+log(maxlag))/(df+1-int))*1:(df-int)-1)
}

list <- mkbasis(lag,type=type,df=df,knots=knots,bound=bound,int=int,cen=FALSE)
rownames(list$basis) <- outer("lag",lag, function(x,y) paste(x,y,sep=""))
list$maxlag <- maxlag
list$cen <- NULL
list

}


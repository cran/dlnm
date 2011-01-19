`mklagbasis` <-
function(maxlag=0, type="ns", df=1, degree=1, knots=NULL,
	bound=c(0,maxlag), int=TRUE) {

# MAXLAG MUST BE >=0
if(!is.numeric(maxlag)||maxlag<0) {
	stop("maxlag must be numeric and >= 0")
}

# IF MAXLAG=0, THAN TYPE="STRATA" AND DF=1
if(maxlag==0) {
	type <- "strata"
	df=1
	knots=NULL
	int=F
}

lag <- 0:maxlag

# SET THE KNOTS (DEFAULT AT EQUALLY-SPACED VALUES IN TH LOG SCALE)
if(is.null(knots)&&df>1+int&&(type=="ns")) {
	knots <- exp(((1+log(maxlag))/(df-int))*1:(df-int-1)-1)
}
if(is.null(knots)&&(type=="bs")) {
	if(df>degree+int) {
		knots <- exp(((1+log(maxlag))/(df-int-degree+1))*
			1:(df-int-degree)-1)
	}
}
if(is.null(knots)&&df>1+int&type%in%c("strata","dthr","hthr","lthr")) {
	knots <- exp(((1+log(maxlag))/(df+1-int))*1:(df-int)-1)
}

list <- mkbasis(lag,type=type,df=df,degree=degree,knots=knots,
	bound=bound,int=int,cen=FALSE)
rownames(list$basis) <- outer("lag",lag,paste,sep="")

list$maxlag <- maxlag
list$cen <- NULL

return(list)
}


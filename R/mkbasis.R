`mkbasis` <-
function(var,type="ns",df=1,knots=NULL,bound=range(var),
	int=FALSE,cen=TRUE,cenvalue=mean(var)) {
list <- vector("list",0)

if(!is.null(bound)) {
	if(any(is.na(bound))) bound <- range(var,na.rm=TRUE)
}
if(!is.null(cenvalue)) {
	if(is.na(cenvalue)) cenvalue <- mean(var,na.rm=TRUE)
}

###########################################################################
# COHERENCE CHECKS
##################

# DEFINE THE POSSIBLE TYPES
if(!type%in%c("ns","strata","poly","integer","thr","lthr","hthr","lin")) {
	stop("type must be one of ns,strata,poly,integer,thr,lthr,hthr,lin")
}
# ONE OF DF OR KNOTS MUST BE GIVEN FOR NS, STRATA OR POLY
if(is.null(df)&is.null(knots) & type%in%c("ns","poly","strata")) {
	stop("for type 'ns', 'poly' or 'strata' at least df or knots must be specified")
}
# KNOTS FORCED WITHIN THE RANGE OF THE VARIABLE
if(!is.null(knots)) {
	if(min(knots)<=min(var,na.rm=TRUE)|max(knots)>=max(var,na.rm=TRUE)) {
	stop("knots must be within the range of var")
}}
# KNOTS ORDERED AND MADE UNIQUE, THEY OVERCOMES DF
if(!is.null(knots)&type%in%c("ns","strata","poly")) {
	knots <- sort(unique(knots))
	df <- length(knots)+1+int
	if(type=="strata") df <- length(knots)+int
}
# DF MUST BE <= N. OBS
if(df>length(var)) {
	stop("df must be <= length(var)")
}
# MINIMUM DF = 1
if((is.null(df)|df<1)&is.null(knots)&type%in%c("ns","strata","poly")) { 
	df <- 1
	stop("df must be >1") 
}
# CENTERING
if(!type%in%c("ns","poly","lin")) cen <- F
if(cen==FALSE) cenvalue <- NULL

# IF KNOTS NOT PROVIDED FOR THRESHOLD MODELS, FIXED AT MEAN
if(is.null(knots)&type%in%c("thr","lthr","hthr")) {
		knots <- mean(var,na.rm=TRUE)
}


###########################################################################

#######
# NS
#######
# IF NOT PROVIDED, KNOTS SET TO EQUALLY SPACED QUANTILES
if(type=="ns")	{
	if(df==1) {
		int <- F
	}
	if(is.null(knots)&df>(1+int)) {
		knots <- quantile(var,1/(df-int)*1:((df-int)-1),na.rm=TRUE)
	}
	if(cen==TRUE) {
	list$basis <- as.matrix(ns(var,df=df,knots=knots,Bound=bound,
		int=int)[,] - matrix(ns(cenvalue,df=df,knots=knots,
		Bound=bound,int=int)[,], 
		nrow=length(var),ncol=length(knots)+1+int,byrow=TRUE))
	} else {
		list$basis <- as.matrix(ns(var,df=df,knots=knots,
		Bound=bound,int=int)[,])
	}
}

##########
# STRATA
##########
# NEVER CENTERED, INTERNAL KNOTS SPECIFY STRATA LOWER BOUNDARIES
# IF ONLY DF PROVIDED, KNOTS PLACED AT EQUALLY SPACED QUANTILES
# KNOTS SPECIFYING EMPTY COLUMNS (AND THE COLUMNS THEMSELVES) ARE EXCLUDED
if(type=="strata")	{
	range <- range(var,na.rm=TRUE)
	if(is.null(knots)&df>1) {
		knots <- quantile(var,1/(df+1-int)*1:(df-int),na.rm=TRUE)
	}
	x <- cut(var,c(range[1],knots,range[2]+0.1),right=FALSE)
	temp <- as.matrix(outer(x,levels(x),"==")+0)
	list$basis <- as.matrix(temp[,apply(temp,2,sum)!=0])
	knots <- knots[(apply(temp,2,sum)!=0)[1:length(knots)]]
	if(df==1) {
		int <- F
	}
	if(int==FALSE&df>1) list$basis <- as.matrix(list$basis[,-1])
	df <- ncol(list$basis)
	bound <- NULL
	cen <- FALSE
	cenvalue <- NULL
}

########
# POLY
########
# THE NUMBER OF KNOTS CAN BE USED TO SPECIFY THE DF (NKNOTS+1+INT)
if(type=="poly")	{
	if(df==1) {
		int <- F
	}
	if(!is.null(knots)) {
		df <- length(knots)+1+int
	} 
	if(cen==TRUE) list$basis <- outer(var-(cenvalue),(1-int):(df-int),"^")
	if(cen==FALSE) list$basis <- outer(var,(1-int):(df-int),"^")
	knots <- NULL
	bound <- NULL
}

###########
# INTEGER
###########
# NEVER CENTERED, IT DOESN'T TAKE INTO ACCOUNT DF OR KNOTS
if(type=="integer")	{
	x <- factor(round(var,0))
	list$basis <- as.matrix(outer(x,levels(x),"==")+0)
	if(int==FALSE&length(levels(x))>1) list$basis <- list$basis[,-1]
	df <- ncol(list$basis)
	knots <- NULL
	bound <- NULL
	cen <- FALSE
	cenvalue <- NULL
}

#######
# THR
#######
# NEVER CENTERED, IT DOESN'T TAKE INTO ACCOUNT DF, ONLY KNOTS
# ONLY THE FIRST 2 KNOTS CONSIDERED, KNOT REPLIED IF SINGLE (VMODEL)
if(type=="thr")	{
	if(length(knots)==1) { 
		knots[2] <- knots[1]
	}
	list$basis <- as.matrix(cbind(-pmin(var-knots[1],0),
		pmax(var-knots[2],0)))
	df <- 2
	knots <- knots[1:2]
	bound <- NULL
	cen <- FALSE
	cenvalue <- NULL
	int <- FALSE
}

#######
# LTHR
#######
# NEVER CENTERED, IT DOESN'T TAKE INTO ACCOUNT DF, ONLY KNOTS
# ONLY THE FIRST KNOT CONSIDERED
if(type=="lthr")	{
	list$basis <- as.matrix(-pmin(var-knots[1],0))
	df <- 1
	knots <- knots[1]
	bound <- NULL
	cen <- FALSE
	cenvalue <- NULL
	int <- FALSE
}

#######
# HTHR
#######
# NEVER CENTERED, IT DOESN'T TAKE INTO ACCOUNT DF, ONLY KNOTS
# ONLY THE LAST KNOT CONSIDERED
if(type=="hthr")	{
	list$basis <- as.matrix(pmax(var-knots[length(knots)],0))
	df <- 1
	knots <- knots[1]
	bound <- NULL
	cen <- FALSE
	cenvalue <- NULL
	int <- FALSE
}

#######
# LIN
#######
# IT DOESN'T TAKE INTO ACCOUNT DF AND KNOTS
if(type=="lin")	{
	if(cen==TRUE) list$basis <- as.matrix(var-cenvalue)
	if(cen==FALSE) list$basis <- as.matrix(var)
	df <- 1
	knots <- NULL
	bound <- NULL
	int <- FALSE
}

##########################################################################

colnames(list$basis) <- paste("b",1:df,sep="")
list$type <- type
list$df <- df
list$knots <- knots
list$bound <- bound
list$int <- int
list$cen <- cen
list$cenvalue <- cenvalue
list
}


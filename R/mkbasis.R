`mkbasis` <-
function(var, type="ns", df=1, degree=1, knots=NULL, bound=range(var),
	int=FALSE, cen=TRUE, cenvalue=mean(var)) {
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
# CHANGE IN DOUBLE THRESHOLD
if(type=="thr") {
	stop("Type for double threshold changed from 'thr' to 'dthr' since version 0.4.0.
See help(crossbasis) for further information")
}
# DEFINE THE POSSIBLE TYPES
if(!type%in%c("ns","bs","strata","poly","integer","dthr","lthr","hthr","lin")) {
	stop("type must be one of ns,bs,strata,poly,integer,dthr,lthr,hthr,lin")
}
# CHECK ARGUMENT TYPE
if(!is.numeric(var)) stop("'var' must be a numeric vector")
if(!is.null(df)&!is.numeric(df)) stop("'df' must be numeric")
if(!is.null(knots)&!is.numeric(knots)) stop("'knots' must be numeric")
if(!is.null(bound)&!is.numeric(bound)) stop("'bound' must be numeric")
if(!is.logical(int)) stop("'int' must be logical")
if(!is.logical(cen)) stop("'cen' must be logical")
if(!is.null(cenvalue)&!is.numeric(cenvalue)) stop("'cenvalue' must be numeric")
if(!is.null(degree)&!is.numeric(degree)) stop("'df' must be numeric")

# ONE OF DF OR KNOTS MUST BE GIVEN
if(is.null(df)&is.null(knots) & !type%in%c("poly","integer","dthr","lin")) {
	stop("at least df or knots must be specified")
}
# DEGREE MUST BE GIVEN FOR POLY, AND OVERCOMES DF
if(type=="poly") {
	if(degree<1) stop("for type 'ns' and 'poly' degree must be >=1")
	df <- degree+int
	knots <- NULL
}
# KNOTS FORCED WITHIN THE RANGE OF THE VARIABLE
if(!is.null(knots)) {
	if(min(knots)<=min(var,na.rm=TRUE)|max(knots)>=max(var,na.rm=TRUE)) {
	stop("knots must be within the range of var")
}}
# DF FIXED FOR TYPES INTEGER (SOLVED LATER), DTHR AND LIN
if(type=="integer") df <- 1
if(type=="dthr") df <- 2+int
if(type=="lin") df <- 1+int

# KNOTS ORDERED AND MADE UNIQUE, THEY OVERCOMES DF
if(!is.null(knots)) {
	knots <- sort(unique(knots))
	if(type=="ns") df <- length(knots)+1+int
	if(type=="bs") df <- length(knots)+degree+int
	if(type%in%c("strata","hthr","lthr")) df <- length(knots)+int
}
# MINIMUM DF=1 AND INT ONLY IF DF>1
if(df<1) stop("df must be >=1") 
if(df==1 & !type%in%c("integer","dthr","lin")) int <- F

# DF MUST BE >=DEGREE+INT FOR BS
if(type=="bs") {
	if(df==1) {
		degree <- 1
	}
	if(df<degree+int) stop("df must be >=degree+int for type='bs'")
}
# DF MUST BE <= N. OBS
if(df>length(var)) {
	stop("df+int must be <= length(var)")
}
# CENVALUE ONLY WITH CEN=T
if(cen==F) cenvalue <- NULL

###########################################################################

#######
# NS
#######
# IF NOT PROVIDED, KNOTS SET TO EQUALLY SPACED QUANTILES
if(type=="ns")	{
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
	degree <- NULL
}

#######
# BS
#######
# IF NOT PROVIDED, KNOTS SET TO EQUALLY SPACED QUANTILES
if(type=="bs")	{
	if(is.null(knots)&df>(degree+int)) {
		knots <- quantile(var,1/(df-int-degree+1)*1:((df-int)-degree),
			na.rm=TRUE)
	}
	if(cen==TRUE) {
	list$basis <- as.matrix(bs(var,df=df,knots=knots,Bound=bound,
		int=int,degree=degree)[,] - matrix(bs(cenvalue,df=df,knots=knots,
		Bound=bound,int=int,degree=degree)[,], 
		nrow=length(var),ncol=length(knots)+degree+int,byrow=TRUE))
	} else {
		list$basis <- as.matrix(bs(var,df=df,knots=knots,
		Bound=bound,int=int,degree=degree)[,])
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
	list$basis <- as.matrix(temp[,apply(temp,2,sum,na.rm=T)!=0])
	knots <- knots[(apply(temp,2,sum,na.rm=T)!=0)[1:length(knots)]]
	if(int==FALSE & !is.null(knots)) {
		list$basis <- as.matrix(list$basis[,-1])
	}
	if(!is.null(knots)) df <- length(knots)+int
	bound <- NULL
	cen <- FALSE
	cenvalue <- NULL
	degree <- NULL
}

########
# POLY
########
# THE NUMBER OF KNOTS CAN BE USED TO SPECIFY THE DF (NKNOTS+1+INT)
if(type=="poly")	{
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
	degree <- NULL
}

#######
# DTHR
#######
# NEVER CENTERED, IT DOESN'T TAKE INTO ACCOUNT DF, ONLY KNOTS
# ONLY THE FIRST AND LAST KNOTS CONSIDERED, KNOT REPLIED IF SINGLE (VMODEL)
if(type=="dthr")	{
	if(is.null(knots)) {
		knots <- quantile(var,(1/3)*1:2,na.rm=TRUE)
	}
	knots <- c(knots[1],knots[length(knots)])
	list$basis <- as.matrix(cbind(-pmin(var-knots[1],0),
		pmax(var-knots[2],0)))
	if(int==TRUE) list$basis <- cbind(1,list$basis)
	bound <- NULL
	cen <- FALSE
	cenvalue <- NULL
	degree <- NULL
}

#######
# LTHR
#######
# NEVER CENTERED, KNOTS SPECIFY THRESHOLDS OR CUT-OFF POINTS
# IF ONLY DF PROVIDED, KNOTS PLACED AT EQUALLY SPACED QUANTILES
if(type=="lthr")	{
	if(is.null(knots)) {
		knots <- quantile(var,1/(df+1-int)*1:(df-int),na.rm=TRUE)
	}
	list$basis <- as.matrix(outer(var,knots, function(a,b) -pmin(a-b,0)))
	if(int==TRUE) list$basis <- cbind(1,list$basis)
	bound <- NULL
	cen <- FALSE
	cenvalue <- NULL
	degree <- NULL
}

#######
# HTHR
#######
# NEVER CENTERED, KNOTS SPECIFY THRESHOLDS OR CUT-OFF POINTS
# IF ONLY DF PROVIDED, KNOTS PLACED AT EQUALLY SPACED QUANTILES
if(type=="hthr")	{
	if(is.null(knots)) {
		knots <- quantile(var,1/(df+1-int)*1:(df-int),na.rm=TRUE)
	}
	list$basis <- as.matrix(outer(var,knots, function(a,b) pmax(a-b,0)))
	if(int==TRUE) list$basis <- cbind(1,list$basis)
	bound <- NULL
	cen <- FALSE
	cenvalue <- NULL
	degree <- NULL
}

#######
# LIN
#######
# IT DOESN'T TAKE INTO ACCOUNT DF AND KNOTS
if(type=="lin")	{
	if(cen==TRUE) list$basis <- as.matrix(var-cenvalue)
	if(cen==FALSE) list$basis <- as.matrix(var)
	if(int==TRUE) list$basis <- cbind(1,list$basis)
	knots <- NULL
	bound <- NULL
	degree <- NULL
}

##########################################################################

colnames(list$basis) <- paste("b",1:df,sep="")
list$type <- type
list$df <- df
list$degree <- degree
list$knots <- knots
list$bound <- bound
list$int <- int
list$cen <- cen
list$cenvalue <- cenvalue
list
}


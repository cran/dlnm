`onebasis` <-
function(x, type="ns", df=1, degree=1, knots=NULL, bound, int=FALSE, cen) {

x <- as.vector(x)
# SET DEFAULT TO BOUND AND CEN IF NA OR NULL
range <- range(x,na.rm=TRUE)
if(missing(bound)||is.null(bound)||any(is.na(bound))) bound <- range
if(missing(cen)||is.null(cen)||is.na(cen)) cen <- mean(x,na.rm=TRUE)

###########################################################################
# COHERENCE CHECKS
##################

# DEFINE THE POSSIBLE TYPES
if(!type%in%c("ns","bs","strata","poly","integer","dthr","lthr","hthr","lin")) {
  stop("type must be one of ns,bs,strata,poly,integer,dthr,lthr,hthr,lin")
}
# CHECK OTHER ARGUMENTS
if(!is.numeric(x)) stop("'x' must be a numeric vector or matrix")
if(!is.numeric(df)||df<1) stop("'df' must be numeric and >=1")
if(!is.null(knots)&&!is.numeric(knots)) stop("'knots' must be numeric")
if(!is.numeric(bound)) stop("'bound' must be numeric")
if(!is.logical(int)) stop("'int' must be logical")
if(!is.logical(cen)&&!is.numeric(cen)) stop("'cen' must be logical or numeric")
if(!is.null(degree)&&!is.numeric(degree)&&degree<1)  stop("'degree' must be numeric")

# DF FIXED FOR TYPES INTEGER (SOLVED LATER), DTHR AND LIN
if(type=="poly") df <- degree+int
if(type=="integer") df <- 1+int
if(type=="dthr") df <- 2+int
if(type=="lin") df <- 1+int

# KNOTS ORDERED AND MADE UNIQUE, THEY OVERCOMES DF
if(!is.null(knots)) {
  knots <- sort(unique(knots))
  if(type=="ns") df <- length(knots)+1+int
  if(type=="bs") df <- length(knots)+degree+int
  if(type%in%c("strata","hthr","lthr")) df <- length(knots)+int
}

# CENTERING VALUE
if(is.logical(cen)) cen <- ifelse(cen,mean(x,na.rm=TRUE),0)

###########################################################################

#######
# NS
#######
if(type=="ns")	{
  # IF NOT PROVIDED, KNOTS SET TO EQUALLY SPACED QUANTILES
  if(is.null(knots)&df>(1+int)) {
    knots <- quantile(x,1/(df-int)*1:((df-int)-1),na.rm=TRUE)
  }
  # IF NOT INTERCEPT-ONLY, CREATE THE BASIS, OTHERWISE A COLUMN OF 1'S
  if(df>=(1+int)) {
    basis <- matrix(ns(x,df=df,knots=knots,Boundary.knots=bound,
      intercept=int)[,],ncol=df)
    if(cen) basis <- sweep(basis,2,ns(cen,df=df,knots=knots,
      Boundary.knots=bound,intercept=int)[,],"-")
  } else basis <- as.matrix(rep(1,length(x)))
  degree <- 3
  if(cen==0) cen <- FALSE
}

#######
# BS
#######
# IF NOT PROVIDED, KNOTS SET TO EQUALLY SPACED QUANTILES
if(type=="bs")	{
  if(df<degree+int) df <- degree+int
  if(is.null(knots)&df>(degree+int)) {
    knots <- quantile(x,1/(df-int-degree+1)*1:((df-int)-degree),na.rm=TRUE)
  }
	# IF NOT INTERCEPT-ONLY, CREATE THE BASIS, OTHERWISE A COLUMN OF 1'S
  if(df>=(1+int)) {
    basis <- matrix(bs(x,df=df,degree=degree,knots=knots,
      Boundary.knots=bound,intercept=int)[,],ncol=df)
    if(cen) basis <- sweep(basis,2,bs(cen,df=df,degree=degree,knots=knots,
      Boundary.knots=bound,intercept=int)[,],"-")
  } else basis <- as.matrix(rep(1,length(x)))
  if(cen==0) cen <- FALSE
}

##########
# STRATA
##########
# NEVER CENTERED, INTERNAL KNOTS SPECIFY STRATA LOWER BOUNDARIES
# IF ONLY DF PROVIDED, KNOTS PLACED AT EQUALLY SPACED QUANTILES
# KNOTS SPECIFYING EMPTY COLUMNS (AND THE COLUMNS THEMSELVES) ARE EXCLUDED
if(type=="strata")	{
  # IF NOT PROVIDED, KNOTS SET TO EQUALLY SPACED QUANTILES
  if(is.null(knots)&&df>1) {
    knots <- quantile(x,1/(df+1-int)*1:(df-int),na.rm=TRUE)
  }
  # CREATE A DESIGN MATRIX WITH DUMMY VARIABLES
  x <- cut(x,c(range[1],knots,range[2]+0.0001),right=FALSE)
  basis <- matrix(outer(x,levels(x),"==")+0,ncol=length(levels(x)))
  # IF WITH INTERCEPT AND MORE THAN 1 COLUMN, ELIMINATE THE FIRST COLUMN
  if(int==FALSE && !is.null(knots)) {
  basis <- basis[,-1,drop=FALSE]
  }
  if(!is.null(knots)) df <- length(knots)+int
  bound <- NULL
  cen <- FALSE
  degree <- NULL
}

########
# POLY
########
# THE NUMBER OF KNOTS CAN BE USED TO SPECIFY THE DF (NKNOTS+1+INT)
if(type=="poly")	{
  basis <- outer(x-(cen),(1-int):(df-int),"^")
  knots <- NULL
  bound <- NULL
  if(cen==0) cen <- FALSE
}

###########
# INTEGER
###########
# NEVER CENTERED, IT DOESN'T TAKE INTO ACCOUNT DF OR KNOTS
if(type=="integer")	{
  x <- factor(round(x,0))
  basis <- as.matrix(outer(x,levels(x),"==")+0)
  if(int==FALSE&length(levels(x))>1) basis <- basis[,-1]
  df <- ncol(basis)
  knots <- NULL
  bound <- NULL
  cen <- FALSE
  degree <- NULL
}

#######
# DTHR
#######
# NEVER CENTERED, IT DOESN'T TAKE INTO ACCOUNT DF, ONLY KNOTS
# ONLY THE FIRST AND LAST KNOTS CONSIDERED, KNOT REPLIED IF SINGLE (VMODEL)
if(type=="dthr")	{
  if(is.null(knots)) {
    knots <- quantile(x,(1/3)*1:2,na.rm=TRUE)
  }
  knots <- c(knots[1],knots[length(knots)])
  basis <- as.matrix(cbind(-pmin(x-knots[1],0),
    pmax(x-knots[2],0)))
  if(int==TRUE) basis <- cbind(1,basis)
  bound <- NULL
  cen <- FALSE
  degree <- NULL
}

#######
# LTHR
#######
# NEVER CENTERED, KNOTS SPECIFY THRESHOLDS OR CUT-OFF POINTS
# IF ONLY DF PROVIDED, KNOTS PLACED AT EQUALLY SPACED QUANTILES
if(type=="lthr")	{
  if(df-int>0) {
    if(is.null(knots)) {
      knots <- quantile(x,1/(df+1-int)*1:(df-int),na.rm=TRUE)
    }
    basis <- as.matrix(outer(x,knots, function(a,b) -pmin(a-b,0)))
    if(int==TRUE) basis <- cbind(1,basis)
  } else basis <- as.matrix(rep(1,length(x)))
  bound <- NULL
  cen <- FALSE
  degree <- NULL
}

#######
# HTHR
#######
# NEVER CENTERED, KNOTS SPECIFY THRESHOLDS OR CUT-OFF POINTS
# IF ONLY DF PROVIDED, KNOTS PLACED AT EQUALLY SPACED QUANTILES
if(type=="hthr")	{
  if(df-int>0) {
    if(is.null(knots)) {
      knots <- quantile(x,1/(df+1-int)*1:(df-int),na.rm=TRUE)
    }
    basis <- as.matrix(outer(x,knots, function(a,b) pmax(a-b,0)))
    if(int==TRUE) basis <- cbind(1,basis)
  } else basis <- as.matrix(rep(1,length(x)))
  bound <- NULL
  cen <- FALSE
  degree <- NULL
}

#######
# LIN
#######
# IT DOESN'T TAKE INTO ACCOUNT DF AND KNOTS
if(type=="lin")	{
  basis <- as.matrix(x-cen)
  if(int==TRUE) basis <- cbind(1,basis)
  knots <- NULL
  bound <- NULL
  degree <- NULL
  if(cen==0) cen <- FALSE
}

##########################################################################

dimnames(basis) <- list(names(x),paste("b",1:df,sep=""))
attributes(basis) <- c(attributes(basis),list(range=range,type=type,df=df,
  degree=degree,knots=knots,bound=bound,int=int,cen=cen))

class(basis) <- c("onebasis","matrix")
return(basis)
}


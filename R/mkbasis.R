###
### R routines for the R package dlnm (c) Antonio Gasparrini 2012-2013
#
`mkbasis` <-
function(x, type="ns", df=1, degree=1, knots=NULL, bound=range(x),
  int=FALSE, cen=TRUE, cenvalue=mean(x)) {
#
################################################################################
#
warning("the functions 'mkbasis' and 'mklagbasis' has been replaced since version
1.5.0 by the new function 'onebasis', and only kept for compatibility reasons.
The users are strongly suggested to apply the new function: see '?onebasis'")
  


# CREATE AN EMPTY LIST
list <- vector("list",0)

# SET DEFAULT TO BOUND AND CENVALUE IF NA OR NULL
if(!is.null(bound)&&any(is.na(bound))) bound <- range(x,na.rm=TRUE)
if(!is.null(cenvalue)&&is.na(cenvalue)) cenvalue <- mean(x,na.rm=TRUE)

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
if(!is.numeric(x)) stop("'x' must be a numeric vector")
if(!is.numeric(df)||df<1) stop("'df' must be numeric and >=1")
if(!is.null(knots)&&!is.numeric(knots)) stop("'knots' must be numeric")
if(!is.null(bound)&&!is.numeric(bound)) stop("'bound' must be numeric")
if(!is.logical(int)) stop("'int' must be logical")
if(!is.null(cen)&&!is.logical(cen)) stop("'cen' must be logical")
if(!is.null(cenvalue)&&!is.numeric(cenvalue)) stop("'cenvalue' must be numeric")
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

# DF MUST BE <= N. OBS
if(df+int>length(x)) {
  stop("df+int must be <= length(x)")
}
# CENVALUE ONLY WITH CEN=T
if(cen==FALSE) cenvalue <- NULL

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
    list$basis <- as.matrix(ns(x,df=df,knots=knots,
      Boundary.knots=bound,intercept=int)[,])
    if(cen==TRUE) {
      list$basis <- sweep(list$basis,2,ns(cenvalue,df=df,
        knots=knots,Boundary.knots=bound,intercept=int)[,],"-")
    }
  } else list$basis <- as.matrix(rep(1,length(x)))
  degree <- NULL
}

#######
# BS
#######
# IF NOT PROVIDED, KNOTS SET TO EQUALLY SPACED QUANTILES
if(type=="bs")	{
  if(df<degree+int) stop("df must be >=degree+int for type='bs'")
  if(is.null(knots)&df>(degree+int)) {
    knots <- quantile(x,1/(df-int-degree+1)*1:((df-int)-degree),na.rm=TRUE)
  }
	# IF NOT INTERCEPT-ONLY, CREATE THE BASIS, OTHERWISE A COLUMN OF 1'S
  if(df>=(1+int)) {
    list$basis <- as.matrix(bs(x,df=df,degree=degree,knots=knots,
      Boundary.knots=bound,intercept=int)[,])
    if(cen==TRUE) {
      list$basis <- sweep(list$basis,2,bs(cenvalue,df=df,
        degree=degree,knots=knots,Boundary.knots=bound,intercept=int)[,],"-")
    }
  } else list$basis <- as.matrix(rep(1,length(x)))
}

##########
# STRATA
##########
# NEVER CENTERED, INTERNAL KNOTS SPECIFY STRATA LOWER BOUNDARIES
# IF ONLY DF PROVIDED, KNOTS PLACED AT EQUALLY SPACED QUANTILES
# KNOTS SPECIFYING EMPTY COLUMNS (AND THE COLUMNS THEMSELVES) ARE EXCLUDED
if(type=="strata")	{
  range <- range(x,na.rm=TRUE)
  # IF NOT PROVIDED, KNOTS SET TO EQUALLY SPACED QUANTILES
  if(is.null(knots)&df>1) {
    knots <- quantile(x,1/(df+1-int)*1:(df-int),na.rm=TRUE)
  }
  # CREATE A DESIGN MATRIX WITH DUMMY VARIABLES
  x <- cut(x,c(range[1],knots,range[2]+0.1),right=FALSE)
  list$basis <- as.matrix(outer(x,levels(x),"==")+0)
  # IF WITH INTERCEPT AND MORE THAN 1 COLUMN, ELIMINATE THE FIRST COLUMN
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
  if(cen==TRUE) list$basis <- outer(x-(cenvalue),(1-int):(df-int),"^")
  if(cen==FALSE) list$basis <- outer(x,(1-int):(df-int),"^")
  knots <- NULL
  bound <- NULL
}

###########
# INTEGER
###########
# NEVER CENTERED, IT DOESN'T TAKE INTO ACCOUNT DF OR KNOTS
if(type=="integer")	{
  x <- factor(round(x,0))
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
    knots <- quantile(x,(1/3)*1:2,na.rm=TRUE)
  }
  knots <- c(knots[1],knots[length(knots)])
  list$basis <- as.matrix(cbind(-pmin(x-knots[1],0),
    pmax(x-knots[2],0)))
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
    knots <- quantile(x,1/(df+1-int)*1:(df-int),na.rm=TRUE)
  }
  list$basis <- as.matrix(outer(x,knots, function(a,b) -pmin(a-b,0)))
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
    knots <- quantile(x,1/(df+1-int)*1:(df-int),na.rm=TRUE)
  }
  list$basis <- as.matrix(outer(x,knots, function(a,b) pmax(a-b,0)))
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
  if(cen==TRUE) list$basis <- as.matrix(x-cenvalue)
  if(cen==FALSE) list$basis <- as.matrix(x)
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

return(list)
}


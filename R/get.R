###
### R routines for the R package dlnm (c) Antonio Gasparrini 2013
#
`.getcoef` <-
function(model, class, cond) {
#
################################################################################
#
  # EXTRACT COEF
  # NB: gam, gee AND geeglm HAVE CLASS glm AS WELL
  coef <- if(any(class%in%c("glm","gam","coxph"))) coef(model) else
    if(any(class%in%c("lme","mer"))) fixef(model) else
    tryCatch(coef(model),error=function(w) "error")
  if(identical(coef,"error")) stop("methods for coef() and vcov() must",
    "exist for the class of object 'model'. If not, extract them manually and",
    "use the arguments 'coef' and 'vcov'")
  # SELECT COEF
  indcoef <- grep(cond,names(coef))
  coef <- coef[indcoef]
#
  return(coef)
}
#
#
`.getvcov` <-
function(model, class, cond) {
#
################################################################################
#
  # EXTRACT VCOV
  # NB: gam, gee AND geeglm HAVE CLASS glm AS WELL
  vcov <- if(any(class%in%c("lm","glm","lme","coxph")) &&
    !any(class%in%c("gee"))) vcov(model) else if(identical(class,c("gee","glm")))
    model$robust.variance else if(any(class%in%c("geeglm")))
    summary(model)$cov.scaled else if(any(class%in%c("mer")))
    as.matrix(vcov(model)) else tryCatch(vcov(model),error=function(w) "error")
  if(identical(vcov,"error")) stop("methods for coef() and vcov() must",
    "exist for the class of object 'model'. If not, extract them manually and",
    "use the arguments 'coef' and 'vcov'")
  # SELECT VCOV
  indvcov <- grep(cond,dimnames(vcov)[[1]])
  vcov <- vcov[indvcov,indvcov,drop=FALSE]
#
  return(vcov)
}
#
#
`.getlink` <-
function(model, class) {
#
################################################################################
#
  # EXTRACT MODEL LINK
  link <- if(all(class%in%c("lm")) || all(class%in%c("lme")) ||
    any(class%in%"nlme")) "identity" else if(any(class %in% c("clogit")))
    "logit" else if(all(class %in% c("coxph"))) "log" else
    if(any(class %in% c("glm")) || any(class %in% c("glmmPQL")))
    model$family$link else NA
#
  return(link)
}


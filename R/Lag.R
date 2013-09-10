###
### Code from Roger Peng, included in the function Lag() in the package tsModel
#
`.Lag` <-
function(v, k, group = NULL) {
  stopifnot(length(k) > 0)
  v <- as.numeric(v)
  if (max(abs(k)) >= length(v)) 
    stop("largest lag in 'k' must be less than 'length(v)'")
  lag.f <- function(x) {
    lagmat <- matrix(nrow = length(x), ncol = length(k))
    n <- length(x)
    for (i in seq(along = k)) {
      lag <- k[i]
      if (lag > 0) 
        lagmat[, i] <- c(rep(NA, lag), x[1:(n - lag)])
      else if (lag < 0) 
        lagmat[, i] <- c(x[(-lag + 1):n], rep(NA, -lag))
      else lagmat[, i] <- x
    }
    lagmat
  }
  lagmat <- if (!is.null(group)) {
    groupLag <- tapply(v, group, lag.f)
    lagmat <- matrix(nrow = length(v), ncol = length(k))
    split(lagmat, group) <- groupLag
    lagmat
  }
  else lag.f(v)
  colnames(lagmat) <- as.character(k)
  drop(lagmat)
}


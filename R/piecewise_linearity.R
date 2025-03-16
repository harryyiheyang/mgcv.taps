#' Generate a piecewise linear basis matrix with specified changepoints
#'
#' This function creates a design matrix representing a piecewise linear function with one or multiple changepoints.
#' Each changepoint introduces a knot at which the slope of the function may change.
#'
#' @param x A numeric vector representing the covariate values at which to evaluate the piecewise linear function.
#' @param para A numeric scalar or numeric vector indicating the location(s) of the changepoint(s) (knot positions).
#'
#' @return A numeric matrix with columns corresponding to the intercept, the original covariate `x`, and additional columns representing the positive parts of `(x - knot)` for each specified changepoint.
#'
#' @examples
#' x <- seq(0, 10, length.out=100)
#' piecewise_linearity(x, para=c(3,7))
#'
#' @export
piecewise_linearity <- function(x, para){
  softthres <- function(x,b){
    y <- x - b
    y[y < 0] <- 0
    return(y)
  }
  knot <- para
  p <- length(knot)
  G <- matrix(0, length(x), p)
  for(i in 1:p){
    G[,i] <- softthres(x, knot[i])
  }
  return(cbind(1, x, G))
}

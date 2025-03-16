#' Generate basis functions for Regression Discontinuity Design with multiple breakpoints
#'
#' This function creates a design matrix for a regression discontinuity model with multiple cutoff points.
#' Each cutoff introduces a discontinuity and a slope change at the corresponding threshold.
#'
#' @param x A numeric vector representing the covariate values.
#' @param para A numeric scalar or vector indicating the location(s) of the changepoint(s).
#'
#' @return A numeric matrix where columns represent the intercept, the covariate `x`,
#'         jump indicators `I(x > para[i])`, and slope change terms `(x - para[i]) I(x > para[i])`.
#'
#' @examples
#' x <- seq(0, 10, length.out=100)
#' linearity_discontinuity(x, para=c(3, 7))
#'
#' @export
linearity_discontinuity <- function(x, para) {
  if (!is.numeric(x)) stop("x must be a numeric vector.")
  if (!is.numeric(para)) stop("para must be a numeric vector or scalar.")

  p <- length(para)
  n <- length(x)

  # Initialize the basis matrix with intercept and x
  B <- matrix(0, n, 2 + 2 * p)
  B[,1] <- 1  # Intercept
  B[,2] <- x  # Linear term

  for (i in 1:p) {
    knot <- para[i]
    B[, 2 + 2*i - 1] <- as.numeric(x > knot)  # Jump at the breakpoint
    B[, 2 + 2*i] <- (x - knot) * (x > knot)   # Slope change after the breakpoint
  }

  colnames(B) <- c("Intercept", "x",
                   as.vector(t(outer(para, c("_jump", "_slope"), paste0))))

  return(B)
}

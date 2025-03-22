#' Generate a polynomial basis matrix
#'
#' This function constructs a polynomial basis matrix for a given covariate `x` up to a specified order.
#' The resulting matrix includes columns for each power of `x` from 0 (intercept) to the specified polynomial order.
#'
#' @param x A numeric vector representing the covariate values at which to evaluate the polynomial basis.
#' @param order A non-negative integer specifying the highest polynomial order to include in the basis matrix.
#'
#' @return A numeric matrix with `length(x)` rows and `order + 1` columns.
#'         The first column corresponds to the intercept (constant term, `x^0 = 1`),
#'         the second column to `x^1`, and so on up to `x^order`.
#'
#' @examples
#' x <- seq(-1, 1, length.out = 10)
#' polynomial(x, para = 3)
#'
#' @export
polynomial <- function(x, para) {
order=para
if (!is.numeric(x)) stop("x must be a numeric vector.")
if (!is.numeric(order) || length(order) != 1 || order < 0 || order %% 1 != 0) {
stop("order must be a non-negative integer.")
}

n <- length(x)
B <- matrix(0, n, order + 1)
for (i in 1:(order + 1)) {
B[, i] <- x^(i - 1)
}
return(B)
}

#' @title Construct a Fixed-Effect-Only Smooth Term
#'
#' @description
#' This function constructs a smooth term based solely on a fixed effect structure defined by `getA`.
#' The fixed effect structure is user-defined via `xt$getA`, with corresponding parameters in `xt$para`.
#' The function ensures that only fixed effects are modeled by checking that `fx=TRUE`. If `fx=FALSE`, it
#' stops with an error message. The penalty matrix `S` is set to a zero matrix with dimensions matching
#' the number of columns in `A`.
#'
#' @usage
#' \method{smooth.construct}{Fixed.smooth.spec}(object, data, knots)
#'
#' @param object A smooth specification object created by `s()`, containing user-defined smoothing parameters.
#' @param data A data frame containing the covariate for the smooth term.
#' @param knots A list of knots (not used in this fixed-effect-only implementation).
#'
#' @return A smooth term object of class `"Fixed.smooth"`, `"mgcv.smooth"`, containing:
#' \itemize{
#'   \item `X`: The design matrix consisting solely of the fixed effect components from `A`.
#'   \item `S`: A zero penalty matrix with dimensions matching the number of columns in `A`.
#'   \item `rank`: The rank of the penalty matrix (set to 0 as it is a zero matrix).
#'   \item `null.space.dim`: The dimension of the fixed effect (null space), equal to the number of columns in `A`.
#'   \item `df`: The degrees of freedom, equal to the number of columns in `A`.
#' }
#'
#' @importFrom mgcv smooth.construct
#' @importFrom Matrix bdiag
#'
#' @export
smooth.construct.Fixed.smooth.spec <- function(object, data, knots) {

# Check if fx is TRUE; stop if FALSE
if (!object$fixed) {
stop("This function only fits fixed parameter functions; please set fx=TRUE.")
}

x <- data[[object$term]]
n <- length(x)

# Obtain or define default getA
if (is.null(object$xt$getA)) {
getA <- function(x, para) cbind(1, x)  # Default to using x itself
para <- 0
A <- getA(x, para)
} else {
getA <- object$xt$getA
para <- object$xt$para
A <- getA(x, para)
}

# Check dimensions
if (!is.matrix(A)) {
stop("'getA(x)' must return a matrix.")
}
if (nrow(A) != n) {
stop("The number of rows in matrix 'A' returned by getA(x) must match length of x.")
}

# Store getA and para for later use
object$getA <- getA
object$para <- para
m <- ncol(A)

# Set design matrix X to A
object$X <- A

# Set penalty matrix S as a zero matrix with dimensions m x m
object$S[[1]] <- matrix(0, m, m)

# Set object properties
object$rank <- 0  # Rank of a zero penalty matrix
object$null.space.dim <- m
object$df <- m

# Assign class
class(object) <- c("Fixed.smooth", "mgcv.smooth")

return(object)
}

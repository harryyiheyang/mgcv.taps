#' Generate basis functions for a Smoothed Regression Discontinuity Design
#'
#' This function constructs a design matrix for a regression discontinuity model with multiple breakpoints,
#' incorporating a smoothed transition around each breakpoint instead of an abrupt jump.
#'
#' @param x A numeric vector representing the covariate values.
#' @param para A list containing:
#'   \itemize{
#'     \item `para[[1]]`: A numeric vector indicating the location(s) of the breakpoint(s).
#'     \item `para[[2]]`: A numeric scalar controlling the smoothing transition width `c`. Larger `c` results in a smoother connection.
#'   }
#'
#' @return A numeric matrix where columns represent the intercept, the covariate `x`,
#'         smoothed jump indicators, and smoothed slope change terms.
#'
#' @details
#' Instead of a sharp discontinuity, this function replaces the traditional jump indicator `I(x > v)`
#' with a **steep but continuous linear transition** in the neighborhood `[v - c, v + c]`. The smoothed function is:
#' \deqn{
#' f(x) = 4x + \frac{2}{c} (x - v + c)_+ - \frac{2}{c} (x - v - c)_+
#' }
#' where `(x - v + c)_+` and `(x - v - c)_+` are hinge functions (ReLU functions).
#'
#' @examples
#' x <- seq(0, 10, length.out=100)
#' smoothed_linearity_discontinuity(x, para=list(c(3, 7), 0.2))
#'
#' @export
smoothed_linearity_discontinuity <- function(x, para) {
if (!is.numeric(x)) stop("x must be a numeric vector.")
if (!is.list(para) || length(para) != 2) stop("para must be a list with two elements: breakpoints and smoothing parameter.")

knots <- para[[1]]  # Breakpoints (numeric vector)
c_val <- para[[2]]  # Smoothing width (scalar)

if (!is.numeric(knots) || !is.numeric(c_val)) stop("Both elements of para must be numeric.")
if (c_val <= 0) stop("Smoothing parameter c must be positive.")

p <- length(knots)
n <- length(x)

# Initialize the design matrix
B <- matrix(0, n, 2 + 2 * p)
B[,1] <- 1  # Intercept
B[,2] <- x  # Linear term

# Smooth transition function using hinge-like basis
smoothed_jump <- function(x, v, c) {
  return((x - v + c) * (x >= v - c & x <= v + c) / (2 * c) + (x > v + c))
}

smoothed_slope <- function(x, v, c) {
  return(((x - v + c)^2 / (4 * c)) * (x >= v - c & x <= v + c) + (x - v) * (x > v + c))
}

# Populate the matrix
for (i in 1:p) {
  knot <- knots[i]
  B[, 2 + 2*i - 1] <- smoothed_jump(x, knot, c_val)   # Smoothed jump term
  B[, 2 + 2*i] <- smoothed_slope(x, knot, c_val)      # Smoothed slope term
}

colnames(B) <- c("Intercept", "x",
                 as.vector(t(outer(knots, c("_jump", "_slope"), paste0))))

return(B)
}

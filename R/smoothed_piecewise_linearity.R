#' Generate a smoothed piecewise linear basis matrix with specified changepoints
#'
#' This function creates a design matrix representing a piecewise linear function with one or multiple changepoints,
#' incorporating a smooth approximation near each changepoint using a quadratic smoothing function.
#'
#' @param x A numeric vector representing the covariate values at which to evaluate the smoothed piecewise linear function.
#' @param para A list containing:
#'   \itemize{
#'     \item `para[[1]]`: A numeric vector indicating the location(s) of the changepoint(s) (knot positions).
#'     \item `para[[2]]`: A numeric scalar representing the smoothing parameter `b`, which controls the transition smoothness around the changepoints.
#'   }
#'
#' @return A numeric matrix with columns corresponding to the intercept, the original covariate `x`, and additional columns representing the smoothed version of `(x - knot)_+` for each specified changepoint.
#'
#' @details
#' This function extends the standard piecewise linear basis by introducing a **quadratic smoothing** around each changepoint.
#' The smoothing is achieved by approximating `(x - v)_+` with a quadratic spline in the interval `[v - b, ν + b]`:
#' \deqn{s(x, v, b) = \frac{(x - v + b)^2}{4b} 1_{[v - b, ν + b]} + (x - v) 1_{[x \geq v + b]}}
#' where `b` is the user-specified smoothing parameter.
#'
#' @examples
#' x <- seq(0, 10, length.out=100)
#' smoothed_piecewise_linearity(x, para=list(c(3,7), 0.5))
#'
#' @export
smoothed_piecewise_linearity <- function(x, para){
  # Extract changepoints and smoothing parameter
  knots <- para[[1]]  # Vector of changepoints
  b <- para[[2]]      # Smoothing parameter

  # Quadratic smoothing function
  smooth_softthres <- function(x, v, b) {
    smoothed_part <- ((x - v + b)^2 / (4 * b)) * (x >= v - b & x <= v + b)
    linear_part <- (x - v) * (x >= v + b)
    return(smoothed_part + linear_part)
  }

  # Construct the basis matrix
  p <- length(knots)
  G <- matrix(0, length(x), p)
  for (i in seq_len(p)) {
    G[, i] <- smooth_softthres(x, knots[i], b)
  }

  # Create output matrix
  basis_matrix <- cbind(1, x, G)

  # Ensure correct column names
  colnames(basis_matrix) <- c("Intercept", "x", paste0("Smooth_", seq_len(p)))

  return(basis_matrix)
}

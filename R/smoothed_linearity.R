#' Generate a Smoothed Linearity Basis Matrix
#'
#' This function constructs a design matrix for testing linearity, incorporating a smooth deviation
#' using a Gaussian-like function centered at `0.5` after rescaling `x` to `[0,1]`.
#'
#' @param x A numeric vector representing the covariate values.
#' @param a A numeric scalar controlling the magnitude of the deviation from linearity.
#'
#' @return A numeric matrix where columns represent:
#'   \itemize{
#'     \item `Intercept`: A column of ones.
#'     \item `x`: The original covariate.
#'     \item `Nonlinear term`: The added nonlinear Gaussian-like bump.
#'   }
#'
#' @details
#' This function rescales `x` to `[0,1]`, applies a **Gaussian perturbation** centered at `0.5`,
#' and then rescales the transformation back to the original range.
#' The transformed function is:
#' \deqn{ f(x) = 4x + a \exp(-16(x - 0.5)^2) }
#' where `a` controls the magnitude of the deviation.
#'
#' @examples
#' x <- seq(0, 10, length.out=100)
#' smoothed_linearity(x, a=2)
#'
#' @export
#'
smoothed_linearity <- function(x,a) {
  if (!is.numeric(x) || length(x) < 2) stop("x must be a numeric vector with at least two elements.")
  if (!is.numeric(a)) stop("a must be a numeric scalar.")
  # Scale x to [0,1]
  x_min <- min(x)
  x_max <- max(x)
  x_scaled <- (x - x_min) / (x_max - x_min)

  # Apply the smoothed deviation
  nonlinear_component <-  a*exp(-16 * (x_scaled - 0.5)^2)

  # Rescale back to original range
  nonlinear_component <- nonlinear_component * (x_max - x_min) + x_min

  return(x+nonlinear_component)
}

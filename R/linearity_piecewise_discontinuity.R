#' Generate a basis matrix for piecewise linear and regression discontinuity effects
#'
#' This function constructs a basis matrix that captures both piecewise linear changes at
#' changepoints and discontinuous jumps at breakpoints.
#'
#' @param x A numeric vector representing the covariate values.
#' @param para A list with two numeric vectors:
#'   \itemize{
#'     \item `para[[1]]`: Numeric vector of changepoints where slope changes occur.
#'     \item `para[[2]]`: Numeric vector of discontinuity points where jumps occur.
#'   }
#'
#' @return A numeric matrix where columns represent:
#'   - An intercept term (1)
#'   - The original covariate `x`
#'   - Additional columns for each changepoint capturing slope changes
#'   - Additional columns for each discontinuity capturing jumps and post-breakpoint slope changes
#'
#' @export
linearity_piecewise_discontinuity <- function(x, para) {
if (!is.list(para) || length(para) != 2) stop("para must be a list of two numeric vectors.")
if (!is.numeric(x)) stop("x must be a numeric vector.")

changepoints <- para[[1]]
breakpoints <- para[[2]]

n <- length(x)
p1 <- length(changepoints)
p2 <- length(breakpoints)

# Initialize the basis matrix with intercept and x
B <- matrix(0, n, 2 + p1 + 2 * p2)
B[,1] <- 1  # Intercept
B[,2] <- x  # Linear term

# Add changepoint basis (piecewise linear terms)
if (p1 > 0) {
  for (i in 1:p1) {
    B[, 2 + i] <- pmax(x - changepoints[i], 0)  # (x - changepoint)_{+}
  }
}

# Add breakpoint basis (discontinuities + post-jump slopes)
if (p2 > 0) {
  for (i in 1:p2) {
    B[, 2 + p1 + 2*i - 1] <- as.numeric(x > breakpoints[i])  # Jump term I(x > ν)
    B[, 2 + p1 + 2*i] <- (x - breakpoints[i]) * (x > breakpoints[i])  # Slope change term
  }
}

colnames(B) <- c("Intercept", "x",
                 paste0("Change_", changepoints),
                 as.vector(t(outer(breakpoints, c("_jump", "_slope"), paste0))))

return(B)
}

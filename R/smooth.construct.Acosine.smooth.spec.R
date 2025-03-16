#' Construct a Cosine Basis Smooth Term with Fixed and Random Effects
#'
#' This function constructs a smooth term using cosine basis functions (Fourier series representation)
#' with support in `[0,1]`. The transformation to `[0,1]` is achieved by mapping `x` through its
#' empirical cumulative distribution function (`ecdf`). The smooth term includes both fixed effects,
#' defined by a user-specified function `getA`, and random effects modeled using cosine basis functions.
#'
#' The fixed effects structure is defined via `xt$getA`, which must be a function taking `x` (data)
#' and `para` (a list of parameters) as inputs. The user-specified parameters must be stored in `xt$para`.
#'
#' The final smooth term dimension (`k`) consists of both the fixed effect (null space) dimension and the
#' random effect (penalized) dimension. The cosine basis expansion generates `k - m + 1` basis functions,
#' where `m` is the null space dimension (fixed effect dimension).
#'
#' @param object A smooth specification object created by `s()`, containing user-defined smoothing parameters.
#' @param data A data frame containing the covariate for the smooth term.
#' @param knots A list of knots supplied by the user or automatically generated from `x` (not used in this method).
#'
#' @return A smooth term object of class `"Acosine.smooth"`, `"mgcv.smooth"`, containing:
#' \itemize{
#'   \item `X`: The final design matrix combining fixed and random effect components.
#'   \item `S`: The smoothing penalty matrix.
#'   \item `rank`: The effective rank of the penalty matrix.
#'   \item `null.space.dim`: The dimension of the fixed effect (null space).
#'   \item `null.project`: Projection matrix onto the null space.
#'   \item `df`: The effective degrees of freedom of the smooth term.
#'   \item `delta`: The power exponent controlling basis function scaling.
#' }
#'
#' @details
#' This function constructs a smooth term using cosine basis functions derived from the Fourier series with
#' support in `[0,1]`. The covariate `x` is first mapped to `[0,1]` using its empirical cumulative distribution
#' function (`ecdf`). The cosine basis functions are defined as:
#' \deqn{ B_k(x) = \cos(k \pi x) }
#' for `k = 1, ..., (k - m + 1)`, where `m` is the null space dimension determined by the fixed effects.
#'
#' The penalty matrix is constructed using a diagonal weight matrix based on the sequence `k^delta`, where
#' `delta` is specified by the user through `s(...,m=delta,...)`.
#'
#' @importFrom CppMatrix matrixMultiply
#' @importFrom Matrix bdiag
#' @importFrom mgcv smooth.construct Predict.matrix
#'
#' @examples
#' \dontrun{
#' library(mgcv)
#' set.seed(42)
#' dat <- data.frame(x = runif(100), y = rnorm(100))
#' fit <- gam(y ~ s(x, bs="Acosine", k=12, m = 6, xt=list(getA=function(x, para) cbind(1, x), para=list())),
#'            data=dat, method="REML")
#' }
#'
#' @export

smooth.construct.Acosine.smooth.spec <- function(object, data, knots) {
  delta <- object$p.order[1]
  x <- data[[object$term]]
  n <- length(x)

  # Obtain or define default getA
  if (is.null(object$xt$getA)) {
    getA <- function(x,para) cbind(1,x)  # Default to using x itself
    para = 0
    A <- getA(x)
  } else {
    getA <- object$xt$getA
    para =  object$xt$para
    A <- getA(x,para)
  }

  # Check dimensions
  if (!is.matrix(A)) {
    stop("'getA(x)' must return a matrix.")
  }
  if (nrow(A) != n) {
    stop("The number of rows in matrix 'A' returned by getA(x) must match length of x.")
  }

  object$getA <- getA  # store function for later use
  object$para = para
  m <- ncol(A)

  # Set default bs.dim
  if (object$bs.dim < 0) object$bs.dim <- 10
  if (object$bs.dim <= m) stop("The smoothing dimension 'k' is too small given the fixed effects.")

  object$m <- m
  if (is.na(delta)) delta <- 4
  object$delta <- delta

  # Compute matrix B
  funi <- ecdf(x)
  xi <- funi(x)
  B <- outer(xi, 1:(object$bs.dim-m+1), function(x, k) cos(k * pi * x))
  svec <- (1:(object$bs.dim-m+1))^delta

  # Combine matrices A and B
  C <- cbind(A, B)

  # QR decomposition
  G <- matrixMultiply(t(C), A)/n
  fitqr <- qr(G)
  Q <- qr.Q(fitqr, complete = TRUE)[, (m + 1):nrow(G), drop = FALSE]
  X <- cbind(A, matrixMultiply(C,Q))
  svec <- c(rep(0, m), svec)
  Omega <- matrixMultiply(t(Q),Q*svec)
  Omega <- Matrix::bdiag(matrix(0, m, m), Omega)
  Omega <- as.matrix(Omega)

  object$X <- X
  if (!object$fixed) {
    object$S[[1]] <- (Omega + t(Omega)) / 2
  }
  object$rank <- object$bs.dim - m + 1
  object$null.space.dim <- m
  object$null.project <- Q
  object$df <- ncol(X)

  class(object) <- c("Acosine.smooth", "mgcv.smooth")
  object
}

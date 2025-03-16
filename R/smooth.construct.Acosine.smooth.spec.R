#' @importFrom CppMatrix matrixMultiply
#' @importFrom Matrix bdiag
#' @importFrom mgcv smooth.construct Predict.matrix
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

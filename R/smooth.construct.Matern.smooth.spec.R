#' @title Construct a Low-Rank Matern Basis Smooth Term
#'
#' @description
#' This function constructs a smooth term using a low-rank Matern basis representation.
#' The Matern basis is constructed based on a set of selected knots and a smoothing
#' scale parameter. Principal component analysis (PCA) is applied to reduce the rank,
#' ensuring that the number of retained components matches the specified smoothing
#' dimension (`k`).
#'
#' The `m` parameter determines:
#'   - The number of knots (`nk`) used for generating the Matern basis (selected as
#'     quantiles of `x`).
#'   - The scale parameter (`lambda_matern`), which controls the smoothness of the Matern kernel.
#'
#' The function returns a smooth term object containing the final design matrix, penalty matrix,
#' and metadata such as the selected Matern quantiles and scale parameters.
#'
#' @usage
#' \method{smooth.construct}{Matern.smooth.spec}(object, data, knots)
#'
#' @param object A smooth specification object created by `s()`, containing user-defined smoothing parameters.
#' @param data A data frame containing the covariate for the smooth term.
#' @param knots A list of knots supplied by the user or automatically generated from `x`.
#'
#' @return A smooth term object of class `"Matern.smooth"`, `"mgcv.smooth"`, containing:
#' \itemize{
#'   \item `X`: The final design matrix based on a low-rank Matern basis.
#'   \item `S`: The smoothing penalty matrix.
#'   \item `rank`: The effective rank of the penalty matrix.
#'   \item `null.space.dim`: The dimension of the unpenalized null space.
#'   \item `null.project`: Projection matrix onto the null space.
#'   \item `smoothfun`: The function used to generate the Matern basis.
#'   \item `kappa_matern`: The selected quantiles used as knots for Matern splines.
#'   \item `lambda_matern`: The scale parameter for the Matern spline.
#' }
#'
#' @importFrom mgcv smooth.construct Predict.matrix
#' @importFrom CppMatrix matrixMultiply matrixEigen
#' @importFrom Matrix bdiag
#'
#' @examples
#' \dontrun{
#' library(mgcv)
#' set.seed(42)
#' dat <- data.frame(x = runif(100), y = rnorm(100))
#' fit <- gam(y ~ s(x, bs="Matern", k=12, xt=list(para=list(scale=2))),
#'            data=dat, method="REML")
#' }
#'
#' @export

smooth.construct.Matern.smooth.spec <- function(object, data, knots) {

  nk <- object$p.order[1]
  lambda_matern <- object$p.order[2]
  x <- data[[object$term]]
  n <- length(x)
  A = cbind(1,x)
  m <- 2

  # Set default bs.dim
  if (object$bs.dim < 0) object$bs.dim <- 10
  if (object$bs.dim <= 2) stop("The smoothing dimension 'k' is too small given the fixed effects.")

  object$m <- 2
  if (is.na(nk)) nk <- min(round(n/3),300)
  object$nk <- nk
  if(nk<object$bs.dim) stop("The number of knots 'm' is too small.")
  if (is.na(lambda_matern)) lambda_matern <- 10*sd(x)
  object$lambda_matern <- lambda_matern

  kappa_matern=kappa_quantile(x,nk)
  object$kappa_matern=kappa_matern

  fit_design=matern_basis(x,kappa=kappa_matern,lambda=lambda_matern)
  B=fit_design$DX
  Omega=fit_design$DK
  C <- cbind(A, B)
  Omega=as.matrix(bdiag(diag(2)*0,Omega))
  # QR decomposition
  G <- matrixMultiply(t(C), A)/n
  fitqr <- qr(G)
  Q <- qr.Q(fitqr, complete = TRUE)[, (m + 1):nrow(G), drop = FALSE]
  B = matrixMultiply(C,Q)
  fiteigen=matrixEigen(matrixMultiply(t(B),B))
  v=fiteigen$vectors[,1:(object$bs.dim-1)]
  B=matrixMultiply(B,v)
  X=cbind(A,B)
  Omega <- matrixListProduct(list(t(v),t(Q),Omega,Q,v))
  Omega <- Matrix::bdiag(diag(2)*0, Omega)
  Omega <- as.matrix(Omega)

  object$X <- X
  object$v = v
  if (!object$fixed) {
    object$S[[1]] <- (Omega + t(Omega)) / 2
  }
  object$rank <- object$bs.dim-1
  object$null.space.dim <- 2
  object$null.project <- matrixMultiply(Q,v)
  object$df <- ncol(X)
  object$smoothfun=matern_basis
  class(object) <- c("Matern.smooth", "mgcv.smooth")
  object
}

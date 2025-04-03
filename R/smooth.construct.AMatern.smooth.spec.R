#' @title Construct a Matern-based Smooth Term with Fixed and Random Effects
#'
#' @description
#' This function constructs a smooth term based on Matern splines, incorporating both fixed and random effects.
#' The fixed effects structure is user-defined via `getA`, which is stored in `xt$getA`, while the second
#' element of `xt` (i.e., `xt$para`) contains the corresponding parameters for `getA`. Users can customize
#' the fixed effect structure, but the associated parameters must be provided as a list, ensuring that the
#' function only has `x` (data) and `para` as inputs.
#'
#' The final smooth term dimension (`k`) includes both the fixed effect (null space) dimension and the
#' random effect dimension. The `m` parameter allows two values: the number of knots (`nk`) for generating
#' the Matern spline (which is always based on quantiles of `x`), and the scale parameter used for the
#' Matern spline construction. The resulting Matern spline undergoes PCA, extracting the top principal
#' components, ensuring that the number of retained PCs is `k - fixed effect dimension`.
#'
#' This function constructs a smooth term using a Matern spline approach. The fixed effect structure is
#' defined through `getA(x, para)`, which is extracted from `xt$getA` and takes two inputs: `x` (data) and
#' `para` (a list of parameters). The number of knots (`nk`) is determined based on quantiles of `x`. The
#' smoothing term undergoes PCA, extracting the top components such that the number of retained principal
#' components equals `k - fixed effect dimension`.
#'
#' The function returns a smooth term object containing the final design matrix, penalty matrix, and
#' additional metadata such as the selected Matern quantiles and scale parameters.
#'
#' @usage
#' \method{smooth.construct}{AMatern.smooth.spec}(object, data, knots)
#'
#' @param object A smooth specification object created by `s()`, containing user-defined smoothing parameters.
#' @param data A data frame containing the covariate for the smooth term.
#' @param knots A list of knots supplied by the user or automatically generated from `x`.
#'
#' @return A smooth term object of class `"AMatern.smooth"`, `"mgcv.smooth"`, containing:
#' \itemize{
#'   \item `X`: The final design matrix combining fixed and random effect components.
#'   \item `S`: The smoothing penalty matrix.
#'   \item `rank`: The effective rank of the penalty matrix.
#'   \item `null.space.dim`: The dimension of the fixed effect (null space).
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
#' piecewise_linearity <- function(x, para){
#'  soft_thresholding <- function(x,b){
#'  y <- x - b
#'  y[y < 0] <- 0
#'  return(y)
#'  }
#'  knot <- para
#'  p <- length(knot)
#'  G <- matrix(0, length(x), p)
#'  for(i in 1:p){
#'  G[,i] <- soft_thresholding(x, knot[i])
#'  }
#'  basis_matrix <- cbind(1, x, G)
#'  colnames(basis_matrix) <- c("Intercept", "x", paste0("Smooth_",  seq_len(p)))
#'  return(basis_matrix)
#'  }
#' set.seed(42)
#' dat <- data.frame(x = runif(100), y = rnorm(100))
#' fit <- gam(y ~ s(x, bs="AMatern", k=12, xt=list(getA=piecewise_linearity, para=c(0.2,0.5))),
#'            data=dat, method="REML")
#' }
#'
#' @export
smooth.construct.AMatern.smooth.spec <- function(object, data, knots) {

nk <- object$p.order[1]
lambda_matern <- object$p.order[2]
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
Omega=as.matrix(bdiag(diag(m)*0,Omega))
# QR decomposition
G <- matrixMultiply(t(C), A)/n
fitqr <- qr(G)
Q <- qr.Q(fitqr, complete = TRUE)[, (m + 1):nrow(G), drop = FALSE]
B = matrixMultiply(C,Q)
fiteigen=matrixEigen(matrixMultiply(t(B),B))
v=fiteigen$vectors[,1:(object$bs.dim-m+1)]
B=matrixMultiply(B,v)
X=cbind(A,B)
Omega <- matrixListProduct(list(t(v),t(Q),Omega,Q,v))
Omega <- Matrix::bdiag(matrix(0, m, m), Omega)
Omega <- as.matrix(Omega)

object$X <- X
object$v = v
if (!object$fixed) {
object$S[[1]] <- (Omega + t(Omega)) / 2
}
object$rank <- object$bs.dim - m + 1
object$null.space.dim <- m
object$null.project <- matrixMultiply(Q,v)
object$df <- ncol(X)
object$smoothfun=matern_basis
class(object) <- c("AMatern.smooth", "mgcv.smooth")
object
}

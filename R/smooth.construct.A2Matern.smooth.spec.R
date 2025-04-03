#' @title Construct a Bivariate Matern-based Smooth Term with Fixed and Random Effects
#'
#' @description
#' This function constructs a smooth term for bivariate data (x1, x2) using Matern splines, incorporating
#' both fixed and random effects to model parameter effects. The fixed effects structure is user-defined
#' via `getA`, stored in `xt$getA`, with corresponding parameters provided in `xt$para` as lists.
#' These functions must take only `x` (data) and `para` (parameters) as inputs, allowing customization of the fixed effect structure.
#'
#' The smooth term dimension (`k`) combines the fixed effect (null space) dimension and the random effect
#' dimension. The `m` parameter specifies two values: the number of knots (`nk`) for generating the Matern
#' spline (based on k-means clustering of `x1` and `x2`), and the scale parameter (`lambda`) for the Matern
#' spline construction. The Matern spline is subjected to PCA, retaining the top principal components such
#' that the number of retained PCs equals `k - fixed effect dimension`.
#'
#' The function returns a smooth term object with the design matrix, penalty matrix, and metadata such as
#' selected Matern knots and scale parameters, suitable for modeling bivariate parameter effects in a GAM.
#'
#' @usage
#' \method{smooth.construct}{A2Matern.smooth.spec}(object, data, knots)
#'
#' @param object A smooth specification object created by `s()`, containing user-defined smoothing parameters.
#' @param data A data frame containing the two covariates (`x1`, `x2`) for the bivariate smooth term.
#' @param knots A list of knots supplied by the user or automatically generated from `x1` and `x2`.
#'
#' @return A smooth term object of class `"A2Matern.smooth"`, `"mgcv.smooth"`, containing:
#' \itemize{
#'   \item `X`: The final design matrix combining fixed and random effect components.
#'   \item `S`: The smoothing penalty matrix.
#'   \item `rank`: The effective rank of the penalty matrix.
#'   \item `null.space.dim`: The dimension of the fixed effect (null space).
#'   \item `null.project`: Projection matrix onto the null space.
#'   \item `smoothfun`: The function used to generate the Matern basis.
#'   \item `kappa_matern`: The selected k-means centroids used as knots for Matern splines.
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
#' dat <- data.frame(x1 = runif(100), x2 = runif(100), y = rnorm(100))
#' fit <- gam(y~s(x1, x2, bs = "A2Matern", k = 12,
#'                xt=list(getA=function(x1,x2,para){cbind(1,x1,x2,x1*x2},
#'                para=NULL)
#'            data = dat, method = "REML")
#' summary(fit)
#' }
#'
#' @export
smooth.construct.A2Matern.smooth.spec <- function(object, data, knots) {

nk <- object$p.order[1]
lambda_matern <- object$p.order[2]
x1 <- data[[object$term[1]]]
x2 <- data[[object$term[2]]]
n <- length(x1)

# Obtain or define default getA
if (is.null(object$xt$getA)) {
getA <- function(x1,x2,para){
A=cbind(1,x1,x2,x1*x2)
return(A)
}  # Default to using x itself
para = 0
A <- getA(x1,x2)
} else {
getA <- object$xt$getA
para =  object$xt$para
A <- getA(x1,x2,para)
}

# Check dimensions
if (!is.matrix(A)) {
  stop("'getA(x)' must return a matrix.")
}
if (nrow(A) != n) {
  stop("The number of rows in matrix 'A' returned by getA(x) must match length of x.")
}

object$getA <- getA # store function for later use
object$para = para
m <- ncol(A)
# Set default bs.dim
if (object$bs.dim < 0) object$bs.dim <- 50
if (object$bs.dim <= m) stop("The smoothing dimension 'k' is too small given the fixed effects.")

object$m <- m
if (is.na(nk)) nk <- min(round(n/3),300)
object$nk <- nk
if(nk<object$bs.dim) stop("The number of knots 'm' is too small.")
if (is.na(lambda_matern)) lambda_matern <- 10*max(sd(x1),sd(x2))
object$lambda_matern <- lambda_matern

kappa_matern=kappa_kmeans_2d(x1=x1,x2=x2,nk)
object$kappa_matern=kappa_matern

fit_design=matern_basis_2d(x1=x1,x2=x2,kappa1=kappa_matern$kappa1,kappa2=kappa_matern$kappa2,lambda=lambda_matern)
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
object$smoothfun=matern_basis_2d
class(object) <- c("A2Matern.smooth", "mgcv.smooth")
object
}

#' @importFrom CppMatrix matrixMultiply matrixEigen
#' @importFrom Matrix bdiag
#' @importFrom mgcv smooth.construct Predict.matrix
#' @export
smooth.construct.Amatern.smooth.spec <- function(object, data, knots) {

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
if (is.na(nk)) nk <- min(round(n/3),100)
object$nk <- nk
if (is.na(lambda_matern)) lambda_matern <- 5*sd(x)
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
v=fiteigen$vectors[,1:(object$bs.dim-m)]
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
object$rank <- object$bs.dim - m
object$null.space.dim <- m
object$null.project <- matrixMultiply(Q,v)
object$df <- ncol(X)
object$smoothfun=matern_basis
class(object) <- c("Amatern.smooth", "mgcv.smooth")
object
}

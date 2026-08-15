#' @usage \method{Predict.matrix}{AMatern.smooth}(object, data)
#' @importFrom mgcv smooth.construct Predict.matrix
#' @importFrom CppMatrix matrixMultiply matrixListProduct
#' @importFrom Matrix bdiag
#' @export
Predict.matrix.AMatern.smooth <- function(object, data) {
x_new <- data[[object$term]]
n_new <- length(x_new)

# Generate A matrix using stored function
A_new <- object$getA(x_new,object$para)
A_new <- taps_A_matrix(A_new, n_new)

m <- ncol(A_new)
Q <- object$null.project
lambda_matern = object$lambda_matern
kappa_matern = object$kappa_matern

fit_design=object$smoothfun(x_new,kappa_matern,lambda_matern)
B=fit_design$DX

# Combine A and B for prediction
C <- cbind(A_new, B)
X <- cbind(A_new, matrixMultiply(C,Q))

X
}

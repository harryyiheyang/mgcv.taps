#' @usage \method{Predict.matrix}{Matern.smooth}(object, data)
#' @importFrom mgcv smooth.construct Predict.matrix
#' @importFrom CppMatrix matrixMultiply matrixListProduct
#' @importFrom Matrix bdiag
#' @export
Predict.matrix.Matern.smooth <- function(object, data) {
  x_new <- data[[object$term]]
  n_new <- length(x_new)
  m <- 2
  Q <- object$null.project
  lambda_matern = object$lambda_matern
  kappa_matern = object$kappa_matern

  fit_design=object$smoothfun(x_new,kappa_matern,lambda_matern)
  B=fit_design$DX

  # Combine A and B for prediction
  C <- cbind(1,x_new,B)
  X <- cbind(1,x_new,matrixMultiply(C,Q))

  X
}

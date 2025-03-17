#' @usage \method{Predict.matrix}{matern.smooth}(object, data)
#' @importFrom mgcv smooth.construct Predict.matrix
#' @importFrom CppMatrix matrixMultiply matrixListProduct
#' @importFrom Matrix bdiag
#' @export
Predict.matrix.matern.smooth <- function(object, data) {
  x_new <- data[[object$term]]
  n_new <- length(x_new)
  A_new <- as.matrix(rep(1,length(x_new)))
  m <- 1
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

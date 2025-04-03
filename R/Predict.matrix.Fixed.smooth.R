#' @title Prediction Matrix for Fixed-Effect-Only Smooth Term
#'
#' @description
#' This function generates the prediction matrix for an `Fixed.smooth` object, using the stored
#' `getA` function and its parameters (`para`) to construct the fixed effect basis for new data.
#'
#' @usage
#' \method{Predict.matrix}{Fixed.smooth}(object, data)
#'
#' @param object A smooth term object of class `"Fixed.smooth"`, containing the stored `getA` function and `para`.
#' @param data A data frame containing the covariate for prediction.
#'
#' @return A matrix representing the fixed effect basis for the new data.
#'
#' @importFrom mgcv Predict.matrix
#'
#' @export
Predict.matrix.Fixed.smooth <- function(object, data) {
# Extract the covariate from the new data
x <- data[[object$term]]

# Use the stored getA function and parameters to construct the prediction matrix
X <- object$getA(x, object$para)

# Return the prediction matrix
return(X)
}

#' @usage \method{Predict.matrix}{cosine.smooth}(object, data)
#' @importFrom mgcv smooth.construct Predict.matrix
#' @export
Predict.matrix.cosine.smooth<-function(object,data) {
  ## prediction method function for the `tr' smooth class
  x <- data[[object$term]]
  funi <- object$transform_fun
  xi = funi(x)
  X<-matrix(0,length(x),object$bs.dim)
  X[,1]=1
  for (i in 2:object$bs.dim){
    X[,i]<- cos((i-1)*pi*xi)
  }
  X # return the prediction matrix
}


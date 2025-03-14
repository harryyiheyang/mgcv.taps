## Adding a penalized truncated power basis class and methods
## as favoured by Ruppert, Wand and Carroll (2003)
## Semiparametric regression CUP. (No advantage to actually
## using this, since mgcv can happily handle non-identity
## penalties.)

#' @export
smooth.construct.cosine.smooth.spec<-function(object,data,knots) {
  ## a truncated power spline constructor method function
  ## object$p.order = null space dimension
  m <- object$p.order[1]
  delta = object$p.order[2]
  if (object$bs.dim<0) object$bs.dim <- 10 ## default
  if (object$bs.dim<=1) stop("k too small for m")
  if (is.na(m)) m <- 1
  object$m=m
  if (is.na(delta)) delta <- 4
  object$delta=delta
  x <- data[[object$term]]  ## the data
  funi = ecdf(x)
  xi = funi(x)
  X<-matrix(0,length(x),object$bs.dim)
  X[,1]=1
  svec=c(1:object$bs.dim)*0
  for (i in 2:object$bs.dim){
  X[,i]<- cos((i-1)*pi*xi)
  svec[i]=(i-1)^delta
  }
  svec[1:m]=0
  object$X<-X # the finished model matrix
  if (!object$fixed) # create the penalty matrix
  {
  object$S[[1]]<-diag(svec)
  }
  object$rank<-object$bs.dim-m  # penalty rank
  object$null.space.dim <- m  # dim. of unpenalized space
  ## store "tr" specific stuff ...

  object$df<-ncol(object$X)     # maximum DoF (if unconstrained)

  class(object)<-"cosine.smooth"  # Give object a class
  object
}

#' @export
Predict.matrix.cosine.smooth<-function(object,data) {
  ## prediction method function for the `tr' smooth class
  x <- data[[object$term]]
  funi = ecdf(x)
  xi = funi(x)
  X<-matrix(0,length(x),object$bs.dim)
  X[,1]=1
  for (i in 2:object$bs.dim){
    X[,i]<- cos((i-1)*pi*xi)
  }
  X # return the prediction matrix
}


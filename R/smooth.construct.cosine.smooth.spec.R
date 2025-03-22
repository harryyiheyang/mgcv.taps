#' @title Construct a Cosine Basis Smooth Term
#'
#' @description
#' This function constructs a smooth term using a cosine basis expansion, where the covariate `x` is first
#' mapped to `[0,1]` using its empirical cumulative distribution function (`ecdf`). The cosine basis functions
#' are then defined as:
#' \deqn{ B_k(x) = \cos((k-1) \pi x) }
#' for `k = 1, ..., k`, where `k` is the basis dimension (`bs.dim`).
#'
#' The smooth term includes both an unpenalized null space of dimension `m` and a penalized component where
#' higher-order basis functions are penalized with an exponent `delta`. The penalty is applied as:
#' \deqn{ s_k = (k-1)^delta, \quad k > m }
#' ensuring that higher-order terms shrink more strongly with increasing `k`.
#'
#' The function constructs a smooth term using a cosine basis expansion, similar to the Fourier series but
#' without sine terms. The covariate `x` is transformed using its empirical cumulative distribution function
#' (`ecdf(x)`) to ensure uniform spacing. The resulting basis functions provide a flexible and smooth representation
#' of periodic or oscillatory effects.
#'
#' The first `m` basis functions form the null space and are unpenalized, while higher-order components
#' are penalized using a diagonal penalty matrix with entries `(k-1)^delta` for `k > m`, ensuring
#' increasing smoothness constraints on higher-frequency components. `m` and `delta` can be specified by users through `s(..,m=c(m,delta),...)`.
#'
#' The final smooth term object contains the design matrix (`X`), penalty matrix (`S`), rank, null space dimension,
#' and degrees of freedom, ensuring compatibility with `mgcv`.
#'
#' @usage
#' \method{smooth.construct}{cosine.smooth.spec}(object, data, knots)
#'
#' @param object A smooth specification object created by `s()`, containing user-defined smoothing parameters.
#' @param data A data frame containing the covariate for the smooth term.
#' @param knots A list of knots supplied by the user or automatically generated from `x` (not used in this method).
#'
#' @return A smooth term object of class `"cosine.smooth"`, `"mgcv.smooth"`, containing:
#' \itemize{
#'   \item `X`: The final design matrix using cosine basis functions.
#'   \item `S`: The smoothing penalty matrix (if applicable).
#'   \item `rank`: The effective rank of the penalty matrix.
#'   \item `null.space.dim`: The dimension of the unpenalized null space.
#'   \item `df`: The effective degrees of freedom of the smooth term.
#'   \item `delta`: The exponent controlling penalty strength.
#' }
#'
#' @importFrom mgcv smooth.construct Predict.matrix
#'
#' @examples
#' \dontrun{
#' library(mgcv)
#' set.seed(42)
#' dat <- data.frame(x = runif(100), y = rnorm(100))
#' fit <- gam(y ~ s(x, bs="cosine", k=10, m=c(2,4), xt=list()),
#'            data=dat, method="REML")
#' }
#'
#' @export
smooth.construct.cosine.smooth.spec<-function(object,data,knots) {
## a truncated power spline constructor method function
## object$p.order = null space dimension
m <- object$p.order[1]
delta = object$p.order[2]
transform_method=object$xt$method
if(is.null(transform_method)){
object$xt$method="empirical quantile"
}
if (object$bs.dim<0) object$bs.dim <- 10 ## default
if (object$bs.dim<=1) stop("k too small for m")
if (is.na(m)) m <- 1
object$m=m
if (is.na(delta)) delta <- 4
object$delta=delta
x <- data[[object$term]]  ## the data
if(object$xt$method=="empirical quantile"){
funi = ecdf(x)
xi = funi(x)
}else{
funi = fit_beta_transform(x,k=object$xt$cluster.size)
xi = funi(x)
}
object$transform_fun=funi
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

class(object)<-c("cosine.smooth", "mgcv.smooth")
object
}

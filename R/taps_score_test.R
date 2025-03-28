#' Perform Score Tests on mgcv Smooth Terms (Null and Penalized Components)
#'
#' This function takes a fitted mgcv `gam` or `bam` object,
#' then performs Score test on the null (unpenalized/fixed) and penalized
#' (smooth) components of each smooth term.
#' It returns a clear summary as a `data.table`.
#'
#' @param fit A fitted `gam` or `bam` model from the `mgcv` package.
#' @param test.component The index of smooth component to be tested.
#' @param null.tol The tolerance of row norm to detect indices of null space. Default to 1e-10.
#' @return A `data.table` summarizing Score statistics and p-values separately
#'         for penalized (smooth) components of each term.
#'
#' @importFrom mgcv vcov.gam
#' @importFrom data.table data.table rbindlist
#' @importFrom ordgam testStat
#' @importFrom stats pchisq
#' @importFrom Matrix bdiag
#' @import CppMatrix
#'
#' @examples
#' \dontrun{
#' fit <- gam(y ~ s(x0) + s(x1), data=dat, method="REML")
#' taps_score_test(fit)
#' }
#'
#' @export
taps_score_test <- function(fit,test.component=1,null.tol=1e-10) {
if (!inherits(fit, "gam")) stop("fit must be a 'gam' or 'bam' object.")
if (is_canonical_link(fit)==F) stop("only canonical link is allowed currently.")
#################### Get basic information #######################
beta <- fit$coefficients
Vb   <- vcov.gam(fit)
X <- predict(fit, fit$model, type = "lpmatrix")
if (nrow(X)>5000) cat("Note: score test can be substantially time-consuming for large sample size.\n")
rdf <- fit$df.residual
smooth_terms <- fit$smooth
p=length(smooth_terms)
phivec=fit$sp

eta <- fit$linear.predictors
mu <- fit$fitted.values
family_info <- fit$family
var_mu <- family_info$variance(mu)
y <- fit$y
pseudo_response <- eta + (y - mu)/var_mu
W_diag <- var_mu
phi0 <- summary(fit)$dispersion
if (is.null(phi0) || !is.numeric(phi0)) {
phi0 <- 1
}
V_phi <- phi0 /W_diag
#################### Get the inverse matrix V_inv #######################

smooth_index_list=list()
random_index_list=list()
S_list=list()

for(i in 1:p){
s <- smooth_terms[[i]]
indices   <- s$first.para:s$last.para
reported_null_dim  <- s$null.space.dim
S_matrix <- s$S[[1]]
if(is.null(s$getA)==0){
detected_null_indices <- indices[1:reported_null_dim]
}else{
col_norms <- apply(S_matrix, 2, function(x) sqrt(sum(x^2)))
detected_null_indices=which(col_norms<null.tol)
}
smooth_indices <- setdiff(indices, detected_null_indices)
if(i!=test.component){
smooth_index_list[[i]]=indices
random_index_list[[i]]=smooth_indices
S_list[[i]]=S_matrix*phivec[i]/phi0
}
if(i==test.component){
Bj=X[,indices]
Thetaj=matrixGeneralizedInverse(S_matrix/norm(S_matrix,"f"))
Gj=matrixListProduct(list(X[,indices],matrixGeneralizedInverse(S_matrix/norm(S_matrix,"f")),t(X[,indices])))
random_index_list[[i]]=smooth_indices
}
}
S_list <- Filter(Negate(is.null), S_list)
smooth_index_list <- Filter(Negate(is.null), smooth_index_list)
random_index_list <- Filter(Negate(is.null), random_index_list)
smooth_index_vec=do.call(c,smooth_index_list)
random_index_vec=do.call(c,random_index_list)
fixed_index_vec=setdiff(1:ncol(X),random_index_vec)
S_All=do.call(bdiag,S_list)
S_All=as.matrix(S_All)
A=X[,fixed_index_vec]
alpha=beta[fixed_index_vec]
B_extend=X[,smooth_index_vec]
## intersect(smooth_index_vec,fixed_index_vec) neq 0
## intersect(random_index_vec,fixed_index_vec) = 0
XtX=matrixMultiply(t(B_extend),B_extend*(1/V_phi))
C=matrixInverse(XtX+S_All)
V_inv=matrixListProduct(list(X[,smooth_index_vec],C,t(X[,smooth_index_vec])))
V_inv=diag(1/V_phi)-t(t(V_inv)*(1/V_phi))*(1/V_phi)
####################### Perform Score Test ########################
V_invA=matrixMultiply(V_inv,A)
P=V_inv-matrixListProduct(list(V_invA,matrixGeneralizedInverse(matrixListProduct(list(t(A),V_invA))),t(V_invA)))
error=pseudo_response-matrixVectorMultiply(A,alpha)
error=matrixVectorMultiply(V_inv,error)
u=max(0,sum(error*(matrixVectorMultiply(Gj,error)))/2)
Cj=matrixListProduct(list(t(Bj),P,Bj))
PGj=matrixMultiply(Cj,Thetaj)
D=matrixMultiply(Thetaj, PGj)
PGj2= matrixMultiply(Cj, D)
e=sum(diag(PGj))/2
h=sum(diag(PGj2))/2
kappa=h/2/e
nu=2*e^2/h
pv=pchisq(u/kappa, nu, lower.tail = FALSE)
out_list=data.table(smooth.term=smooth_terms[[test.component]]$label,smooth.df=nu,smooth.stat=u/kappa,smooth.pvalue=pv)
return(out_list)
}

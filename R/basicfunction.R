matern_basis=function(x,kappa,lambda=10){
DX=outer(x, kappa, FUN = function(mu, a) a - mu)
DK=outer(kappa, kappa, FUN = function(mu, a) a - mu)
DX=abs(DX)/lambda
DK=abs(DK)/lambda
DK=exp(-DK)*(1+DK)
DX=exp(-DX)*(1+DX)
return(list(DX=DX,DK=DK))
}

kappa_quantile=function(x,nk){
n=length(x)
quantile(x[2:(n-1)],seq(0,1,length=nk+2))[2:(nk+1)]
}

mgcv_wald <- function(beta, Vb, indices) {
if (length(indices) == 0) return(c(NA, NA, NA))

est <- beta[indices]
V <- Vb[indices, indices, drop=FALSE]

stat <- as.numeric(t(est) %*% solve(V) %*% est)
df <- length(est)
p_val <- pchisq(stat, df=df, lower.tail=FALSE)

c(statistic=stat, df=df, p.value=p_val)
}

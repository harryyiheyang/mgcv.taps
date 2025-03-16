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

mgcv_wald <- function(est, V, df_fixed = FALSE) {
  # This function handles the Wald test for unpenalized (fixed) components:
  # statistic = b^T V^-1 b, assumed ~ Chi^2(df).
  # If df_fixed = TRUE, we interpret this as a fixed-effect Wald test.
  # Returns (statistic, df, p.value).

  if (length(est) == 0) return(c(NA, NA, NA))

  stat <- as.numeric(t(est) %*% solve(V) %*% est)
  df <- nrow(V)
  p_val <- pchisq(stat, df = df, lower.tail = FALSE)
  c(statistic = stat, df = df, p.value = p_val)
}

testStat_wrapper <- function(p, Xt, V, rank, rdf) {
  # This wrapper calls ordgam::testStat(...) and returns
  # (statistic, df, p.value) in a similar style to mgcv_wald().

  # testStat(...) returns a list with elements:
  #   stat = the test statistic,
  #   pval = the p-value,
  #   rank = the effective rank used.
  # We'll rename them for consistency.

  res <- testStat(p, Xt, V, rank = rank, type = 0, res.df = rdf)
  c(statistic = res$stat, df = res$rank, p.value = res$pval)
}

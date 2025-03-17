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

low_rank_test <- function(beta,X, V, approx.method = "Gamma") {
qrx <- qr(X, tol = 0)
R <- qr.R(qrx)
V <- R %*% V[qrx$pivot, qrx$pivot, drop = FALSE] %*% t(R)
V <- (V + t(V))/2
beta = as.vector(R%*%beta[qrx$pivot])
# Ensure valid method selection
if (!approx.method %in% c("Gamma", "Chisq")) stop("Invalid approximation method. Choose 'Gamma' or 'Chisq'.")

# Compute observed test statistic
Q_obs <- sum(beta^2)

# Compute trace-based approximations
tr_V <- sum(diag(V))           # Trace of V
tr_V2 <- sum(diag(V%*%V))            # Trace of V %*% V

# Expected value of Q
expected_Q <- tr_V

# Variance of Q
var_Q <- 2 * tr_V2

# Chi-square approximation
if (approx.method == "Chisq") {
  p_value <- pchisq(Q_obs / expected_Q, df = 1, lower.tail = FALSE)
  return(list(statistic = Q_obs / expected_Q, df=1,p.value = p_value, method = "Appro. χ²"))
}

# Gamma approximation (default)
k <- expected_Q^2 / var_Q
theta <- var_Q / expected_Q
p_value <- pgamma(Q_obs, shape = k, scale = theta, lower.tail = FALSE)

return(list(statistic = Q_obs, p.value = p_value, df = paste0("shape:", round(k, digits = 2),
                                                              " scale:", round(theta, digits = 2),
                                                              " trQ:", formatC(expected_Q, format = "e", digits = 1)),
            method = "Apprx. Gamma"))
}

matern_basis=function(x,kappa,lambda=10){
DX=outer(x, kappa, FUN = function(mu, a) a - mu)
DK=outer(kappa, kappa, FUN = function(mu, a) a - mu)
DX=abs(DX)/lambda
DK=abs(DK)/lambda
DK=exp(-DK)*(1+DK)
DX=exp(-DX)*(1+DX)
return(list(DX=DX,DK=DK))
}

matern_basis_2d <- function(x1, x2, kappa1, kappa2, lambda = 10) {
n <- length(x1)
m <- length(kappa1)

DX1 <- outer(x1, kappa1, FUN = "-")
DX2 <- outer(x2, kappa2, FUN = "-")
DX <- sqrt(DX1^2 + DX2^2) / lambda

DK1 <- outer(kappa1, kappa1, FUN = "-")
DK2 <- outer(kappa2, kappa2, FUN = "-")
DK <- sqrt(DK1^2 + DK2^2) / lambda

DX <- exp(-DX) * (1 + DX)
DK <- exp(-DK) * (1 + DK)

return(list(DX = DX, DK = DK))
}

kappa_quantile=function(x,nk){
n=length(x)
quantile(x[2:(n-1)],seq(0,1,length=nk+2))[2:(nk+1)]
}

kappa_kmeans_2d <- function(x1, x2, nk) {
  data <- cbind(x1, x2)
  km <- kmeans(data, centers = nk, nstart = 10)
  cluster <- km$cluster
  kappa1 <- tapply(x1, cluster, median)
  kappa2 <- tapply(x2, cluster, median)
  return(list(kappa1 = kappa1, kappa2 = kappa2))
}

mgcv_wald <- function(est, V, df_fixed = FALSE) {
# This function handles the Wald test for unpenalized (fixed) components:
# statistic = b^T V^-1 b, assumed ~ Chi^2(df).
# If df_fixed = TRUE, we interpret this as a fixed-effect Wald test.
# Returns (statistic, df, p.value).

if (length(est) == 0) return(c(NA, NA, NA))

stat <- as.numeric(t(est) %*% matrixGeneralizedInverse(V) %*% est)
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

low_rank_test <- function(beta,X, V) {
qrx <- qr(X, tol = 0)
R <- qr.R(qrx)
V <- R %*% V[qrx$pivot, qrx$pivot, drop = FALSE] %*% t(R)
V <- (V + t(V))/2
beta = as.vector(R%*%beta[qrx$pivot])
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
p_value <- pchisq(Q_obs / expected_Q, df = 1, lower.tail = FALSE)/2
return(list(statistic = Q_obs / expected_Q, df=1,p.value = p_value, method = "Appro. χ²"))
}

fit_beta_transform <- function(x, k = round(length(x)/4)) {
# 1. Data preprocessing: Remove NA values and ensure x is a non-empty numeric vector
x <- na.omit(x)
if (is.null(k)) k = round(length(x)/4)
if (k < 10) stop("Data length is insufficient for fitting. Please use uniform quantile.")

# 2. Discretize data into k intervals and compute the median of each interval
breaks <- quantile(x, probs = seq(0, 1, length.out = k + 1), na.rm = TRUE)
breaks[1] <- breaks[1] - 1e-10  # Avoid boundary issues
breaks[k + 1] <- breaks[k + 1] + 1e-10
intervals <- cut(x, breaks = breaks, include.lowest = TRUE)
medians <- tapply(x, intervals, median, na.rm = TRUE)

# 3. Scale the medians to the [0,1] interval
x_scaled <- (medians - min(medians)) / (max(medians) - min(medians))
x_scaled <- pmax(pmin(x_scaled, 1 - 1e-10), 1e-10)  # Constrain within (0,1)

# 4. Fit a Beta distribution using the method of moments
mu <- mean(x_scaled)
var <- var(x_scaled)
alpha <- mu * (mu * (1 - mu) / var - 1)
beta <- (1 - mu) * (mu * (1 - mu) / var - 1)

if (alpha <= 0 || beta <= 0) {
# If moment estimation fails, fall back to default parameters
warning("Moment estimation failed, using default parameters alpha=1, beta=1")
alpha <- 1
beta <- 1
}

# 5. Return a transformation function: Map empirical quantiles to Beta distribution
function(new_x) {
# Compute empirical quantiles of the input data
ecdf_x <- ecdf(x)(new_x)
# Map to the quantiles of the fitted Beta distribution
qbeta(ecdf_x, shape1 = alpha, shape2 = beta)
}
}

is_canonical_link <- function(fit) {
if (!inherits(fit, "gam")) stop("fit must be a 'gam' object.")

used_link <- fit$family$link

family_name <- fit$family$family
canonical_link <- switch(family_name,
                         "gaussian" = "identity",
                         "poisson" = "log",
                         "binomial" = "logit",
                         "Gamma" = "inverse",
                         "inverse.gaussian" = "1/mu^2",
                         stop("Unknown family"))

return(used_link == canonical_link)
}

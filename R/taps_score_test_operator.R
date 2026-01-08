#' @keywords internal
#' Operator-based implementation for large-n score test
#' Used internally by `taps_score_test()` when sample size is large.
taps_score_test_operator <- function(fit, test.component = 1, null.tol = 1e-10, method = "liu",max_eps=1e-16,max_iter=1e5) {
if (!inherits(fit, "gam")) stop("fit must be a 'gam' or 'bam' object.")

if(fit$family$family!="gaulss"){
beta <- fit$coefficients
X <- predict(fit, newdata = fit$model, type = "lpmatrix")
rdf <- fit$df.residual
smooth_terms <- fit$smooth
p <- length(smooth_terms)
phivec <- fit$sp
eta <- fit$linear.predictors
mu <- fit$fitted.values
mu_eta_func <- fit$family$mu.eta
g_prime_mu <- 1 / mu_eta_func(eta)
family_info <- fit$family
var_mu <- family_info$variance(mu)
y <- fit$y
pseudo_response <- eta + (y - mu) * g_prime_mu
W_diag <- 1 / (var_mu * g_prime_mu^2)
phi0 <- summary(fit)$dispersion
if (is.null(phi0) || !is.numeric(phi0)) phi0 <- 1
V_phi <- phi0 / W_diag
}else{
beta_full <- fit$coefficients
X_full <- predict(fit, newdata = fit$model, type = "lpmatrix")
beta_names <- names(beta_full)
var_start_idx <- which(beta_names == "(Intercept).1")
if (length(var_start_idx) == 0) {
stop("Cannot find '(Intercept).1' in coefficient names. Is this a gaulss model?")
}
mean_indices <- 1:(var_start_idx - 1)
beta <- beta_full[mean_indices]
X <- X_full[, mean_indices, drop = FALSE]
rdf <- fit$df.residual
smooth_terms <- fit$smooth
mean_smooth_indices <- which(sapply(smooth_terms, function(s) {
s$first.para < var_start_idx
}))
p <- length(mean_smooth_indices)
smooth_terms <- smooth_terms[mean_smooth_indices]
phivec <- fit$sp[mean_smooth_indices]
eta <- fit$linear.predictors[, 1]
mu <- fit$fitted.values[, 1]
sigma <- exp(fit$linear.predictors[,2])+fit$coefficients["(Intercept).1"]
var_mu <- sigma^2
y <- fit$y
pseudo_response <- eta + (y - mu)
W_diag <- 1 / (var_mu)
phi0 <- 1
V_phi <- phi0 / W_diag
}

smooth_index_list <- list()
random_index_list <- list()
S_list <- list()

for (i in 1:p) {
s <- smooth_terms[[i]]
indices <- s$first.para:s$last.para
reported_null_dim <- s$null.space.dim
S_matrix <- s$S[[1]]

if (is.null(s$getA) == 0) {
detected_null_indices <- indices[1:reported_null_dim]
} else {
col_norms <- apply(S_matrix, 2, function(x) sqrt(sum(x^2)))
detected_null_indices <- which(col_norms < null.tol)
}

smooth_indices <- setdiff(indices, detected_null_indices)

if (i != test.component) {
smooth_index_list[[i]] <- indices
random_index_list[[i]] <- smooth_indices
S_list[[i]] <- S_matrix * phivec[i] / phi0
}

if (i == test.component) {
Bj <- X[, indices]
Thetaj <- matrixGeneralizedInverse(S_matrix / norm(S_matrix, "f"))
Gj_apply <- function(v) {
matrixVectorMultiply(Bj, matrixVectorMultiply(Thetaj, matrixVectorMultiply(t(Bj), v)))
}
random_index_list[[i]] <- smooth_indices
}
}

S_list <- Filter(Negate(is.null), S_list)
smooth_index_list <- Filter(Negate(is.null), smooth_index_list)
random_index_list <- Filter(Negate(is.null), random_index_list)
smooth_index_vec <- do.call(c, smooth_index_list)
random_index_vec <- do.call(c, random_index_list)
fixed_index_vec <- setdiff(seq_len(ncol(X)), random_index_vec)

A <- X[, fixed_index_vec]
alpha <- beta[fixed_index_vec]
B_extend <- X[, smooth_index_vec]
S_All <- as.matrix(bdiag(S_list))
XtX <- matrixMultiply(t(B_extend), B_extend * (1 / V_phi))
C <- matrixInverse(XtX + S_All)

Vinv_apply <- function(v) {
if (is.matrix(v)) {
return(apply(v, 2, function(col) Vinv_apply(col)))
}
part1 <- v / V_phi
Bt_v <- matrixVectorMultiply(t(B_extend), part1)
C_Bt_v <- matrixVectorMultiply(C, Bt_v)
part2 <- matrixVectorMultiply(B_extend, C_Bt_v)
part1 - part2 / V_phi
}

P_apply <- function(v) {
Vinv_v <- Vinv_apply(v)
Vinv_X <- sapply(1:ncol(A), function(j) Vinv_apply(A[, j]))
AVinv_v <- matrixVectorMultiply(t(A), Vinv_v)
XtVinvX <- matrixMultiply(t(A), Vinv_X)
solve_middle <- matrixVectorMultiply(matrixGeneralizedInverse(XtVinvX), AVinv_v)
correction <- matrixVectorMultiply(Vinv_X, solve_middle)
Vinv_v - correction
}

if (method == "satterthwaite") {
error <- pseudo_response - matrixVectorMultiply(A, alpha)
r <- Vinv_apply(error)
Gj_r <- Gj_apply(r)
u <- sum(r * Gj_r) / 2

q <- ncol(Bj)
Cj <- matrix(0, q, q)
for (i in 1:q) {
Pi <- P_apply(Bj[, i])
for (j in i:q) {
Cj[i, j] <- sum(Bj[, j] * Pi)
if (i != j) Cj[j, i] <- Cj[i, j]
}
}
PGj <- matrixMultiply(Cj, Thetaj)
PGj2 <- matrixMultiply(PGj, PGj)
e <- sum(diag(PGj)) / 2
h <- sum(diag(PGj2)) / 2
kappa <- h / (2 * e)
nu <- 2 * e^2 / h
pv <- pchisq(u / kappa, df = nu, lower.tail = FALSE)
test_stat <- u / kappa
}else if (method == "liu") {
error <- pseudo_response - matrixVectorMultiply(A, alpha)
r <- P_apply(error)
Gj_r <- Gj_apply(r)
u <- sum(r * Gj_r)
q <- ncol(Bj)
eig_theta <- matrixsqrt(Thetaj)
Theta_sqrt <- eig_theta$w
N <- sapply(1:q, function(i) P_apply(Bj[, i]))
BtPB <- crossprod(Bj, N)
Q_small <- matrixListProduct(list(Theta_sqrt, BtPB, Theta_sqrt))
lambda <- eigen(Q_small, symmetric = TRUE, only.values = TRUE)$values
lambda <- lambda[lambda > 1e-16]
pv <- compute_liu_pvalue(q = u, lambda = lambda)
nu <- sum(lambda)^2 / sum(lambda^2)
test_stat <- u / sum(lambda)
}else if (method == "hall") {
error <- pseudo_response - matrixVectorMultiply(A, alpha)
r <- P_apply(error)
Gj_r <- Gj_apply(r)
u <- sum(r * Gj_r)
q <- ncol(Bj)
eig_theta <- matrixsqrt(Thetaj)
Theta_sqrt <- eig_theta$w
N <- sapply(1:q, function(i) P_apply(Bj[, i]))
BtPB <- crossprod(Bj, N)
Q_small <- matrixListProduct(list(Theta_sqrt, BtPB, Theta_sqrt))
lambda <- eigen(Q_small, symmetric = TRUE, only.values = TRUE)$values
lambda <- lambda[lambda > 1e-16]
pv <- compute_hbe_pvalue(q = u, lambda = lambda)
nu <- sum(lambda)^2 / sum(lambda^2)
test_stat <- u / sum(lambda)
}else if (method == "wood") {
error <- pseudo_response - matrixVectorMultiply(A, alpha)
r <- P_apply(error)
Gj_r <- Gj_apply(r)
u <- sum(r * Gj_r)
q <- ncol(Bj)
eig_theta <- matrixsqrt(Thetaj)
Theta_sqrt <- eig_theta$w
N <- sapply(1:q, function(i) P_apply(Bj[, i]))
BtPB <- crossprod(Bj, N)
Q_small <- matrixListProduct(list(Theta_sqrt, BtPB, Theta_sqrt))
lambda <- eigen(Q_small, symmetric = TRUE, only.values = TRUE)$values
lambda <- lambda[lambda > 1e-16]
pv <- compute_wood_pvalue(q = u, lambda = lambda)
nu <- sum(lambda)^2 / sum(lambda^2)
test_stat <- u / sum(lambda)
}else if(method=="davies"){
error <- pseudo_response - matrixVectorMultiply(A, alpha)
r <- P_apply(error)
Gj_r <- Gj_apply(r)
u <- sum(r * Gj_r)
q <- ncol(Bj)
eig_theta <- matrixsqrt(Thetaj)
Theta_sqrt <- eig_theta$w
N <- sapply(1:q, function(i) P_apply(Bj[, i]))
BtPB <- crossprod(Bj, N)
Q_small <- matrixListProduct(list(Theta_sqrt, BtPB, Theta_sqrt))
lambda <- eigen(Q_small, symmetric = TRUE, only.values = TRUE)$values
lambda <- lambda[lambda > 1e-16]
pv=CompQuadForm::davies(q=u,lambda=lambda,lim=max_iter,acc=max_eps)$Qq
nu=sum(lambda)^2/sum(lambda^2)
test_stat=u/sum(lambda)
}else if(method=="imhof"){
error <- pseudo_response - matrixVectorMultiply(A, alpha)
r <- P_apply(error)
Gj_r <- Gj_apply(r)
u <- sum(r * Gj_r)
q <- ncol(Bj)
eig_theta <- matrixsqrt(Thetaj)
Theta_sqrt <- eig_theta$w
N <- sapply(1:q, function(i) P_apply(Bj[, i]))
BtPB <- crossprod(Bj, N)
Q_small <- matrixListProduct(list(Theta_sqrt, BtPB, Theta_sqrt))
lambda <- eigen(Q_small, symmetric = TRUE, only.values = TRUE)$values
lambda <- lambda[lambda > 1e-16]
pv=CompQuadForm::imhof(q=u,lambda=lambda,epsabs=max_eps,epsrel=max_eps,limit=max_iter)$Qq
nu=sum(lambda)^2/sum(lambda^2)
test_stat=u/sum(lambda)
}else {
stop("method must be either 'satterthwaite' or 'davies' or 'liu' or 'imhof' or 'hall' or `wood`")
}

data.table(
smooth.term = smooth_terms[[test.component]]$label,
smooth.df = nu,
smooth.stat = test_stat,
smooth.pvalue = pv,
method = method
)
}

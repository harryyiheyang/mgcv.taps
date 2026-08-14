library(mgcv.taps)

set.seed(20260813)
id <- rep(seq_len(4L), c(2L, 3L, 2L, 3L))
n <- length(id)
x <- seq(-1, 1, length.out = n)
X <- cbind(1, x, x^2)
B <- cbind(1, seq(0.2, 1.1, length.out = n))
G <- matrix(c(0.8, 0.18, 0.18, 0.45), 2L, 2L)
S <- diag(c(0, 0.7, 1.2))
y <- sin(seq_len(n)) + stats::rnorm(n, sd = 0.2)

Z <- matrix(0, n, max(id) * ncol(B))
for (i in seq_len(n)) {
  jj <- (id[i] - 1L) * ncol(B) + seq_len(ncol(B))
  Z[i, jj] <- B[i, ]
}
V <- diag(n) + Z %*% kronecker(diag(max(id)), G) %*% t(Z)
Vinv <- solve(V)
Qinv <- solve(crossprod(X, Vinv %*% X) + S)
expected_beta <- Qinv %*% crossprod(X, Vinv %*% y)
P <- Vinv - Vinv %*% X %*% Qinv %*% t(X) %*% Vinv
P <- (P + t(P)) / 2

expected_moment <- matrix(0, ncol(B), ncol(B))
expected_u <- matrix(0, max(id), ncol(B))
for (g in seq_len(max(id))) {
  ii <- which(id == g)
  ug <- G %*% crossprod(B[ii, , drop = FALSE], P[ii, , drop = FALSE] %*% y)
  posterior <- G - G %*% crossprod(
    B[ii, , drop = FALSE],
    P[ii, ii, drop = FALSE] %*% B[ii, , drop = FALSE]
  ) %*% G
  expected_u[g, ] <- ug
  expected_moment <- expected_moment + tcrossprod(ug) + posterior
}
expected_rss <- sum((P %*% y)^2) + n - sum(diag(P))

observed <- mgcv.taps:::gammfast_projected_moments(
  response = y, X = X, B = B, id = id, G = G, penalty = S,
  return_projection = TRUE, n_threads = 2L
)

stopifnot(
  max(abs(observed$beta - expected_beta)) < 1e-10,
  max(abs(observed$mean_covariance - Qinv)) < 1e-10,
  max(abs(observed$Py - P %*% y)) < 1e-10,
  max(abs(observed$PX - P %*% X)) < 1e-10,
  max(abs(observed$u - expected_u)) < 1e-10,
  max(abs(observed$moment_sum - expected_moment)) < 1e-10,
  abs(observed$trace_P - sum(diag(P))) < 1e-10,
  abs(observed$rss_sum - expected_rss) < 1e-10
)

cat("gammfast projected-moment dense regression passed.\n")

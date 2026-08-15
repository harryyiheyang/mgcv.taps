library(devtools)
load_all(".", recompile = FALSE)

set.seed(20260814)
ng <- 4L
m <- 5L
id <- rep(seq_len(ng), each = m)
X <- cbind(1, rnorm(ng * m), runif(ng * m))
B <- cbind(1, rnorm(ng * m), runif(ng * m))
y <- rnorm(ng * m)
G <- matrix(0, 3L, 3L)
G[1L, 1L] <- 0.7
G[2:3, 2:3] <- matrix(c(0.5, 0.12, 0.12, 0.3), 2L, 2L)
sigma2 <- 1.4
cache <- mgcv.taps:::gammfast_gaussian_cache(cbind(X, y), B, id, 1L)
got <- mgcv.taps:::gammfast_gaussian_projected_cached(
  cache$AtA, cache$BtB, cache$BtA, G, sigma2,
  c(1L, 2L, 2L), fisher = TRUE
)

V <- diag(ng * m)
for (i in seq_len(ng)) {
  ii <- which(id == i)
  V[ii, ii] <- V[ii, ii] + B[ii, ] %*% G %*% t(B[ii, ])
}
Vi <- solve(V)
Xs <- X / sqrt(sigma2)
ys <- y / sqrt(sigma2)
H <- crossprod(Xs, Vi %*% Xs)
P <- Vi - Vi %*% Xs %*% solve(H, crossprod(Xs, Vi))
beta <- solve(H, crossprod(Xs, Vi %*% ys))
moment <- matrix(0, 3L, 3L)
for (i in seq_len(ng)) {
  ii <- which(id == i)
  ti <- crossprod(B[ii, ], P[ii, ] %*% ys)
  Ci <- crossprod(B[ii, ], P[ii, ii] %*% B[ii, ])
  ui <- G %*% ti
  moment <- moment + tcrossprod(ui) + G - G %*% Ci %*% G
}
if (max(abs(got$beta - beta)) > 1e-9 ||
    max(abs(got$moment_sum - moment)) > 1e-8) {
  stop("Cached ordinary-projection moments do not match the dense calculation.")
}

L2 <- t(chol(G[2:3, 2:3]))
E <- vector("list", 4L)
E[[1L]] <- diag(c(2 * G[1L, 1L], 0, 0))
k <- 1L
for (a in seq_len(2L)) {
  for (b in seq_len(a)) {
    k <- k + 1L
    dL <- matrix(0, 2L, 2L)
    dL[a, b] <- if (a == b) L2[a, a] else 1
    Ej <- matrix(0, 3L, 3L)
    Ej[2:3, 2:3] <- dL %*% t(L2) + L2 %*% t(dL)
    E[[k]] <- Ej
  }
}
score <- numeric(4L)
information <- matrix(0, 4L, 4L)
for (j in seq_len(4L)) {
  dVj <- matrix(0, ng * m, ng * m)
  for (i in seq_len(ng)) {
    ii <- which(id == i)
    dVj[ii, ii] <- B[ii, ] %*% E[[j]] %*% t(B[ii, ])
  }
  score[j] <- 0.5 * (
    drop(crossprod(P %*% ys, dVj %*% (P %*% ys))) - sum(diag(P %*% dVj))
  )
  for (k in seq_len(4L)) {
    if (got$information_group[j] != got$information_group[k]) next
    dVk <- matrix(0, ng * m, ng * m)
    for (i in seq_len(ng)) {
      ii <- which(id == i)
      dVk[ii, ii] <- B[ii, ] %*% E[[k]] %*% t(B[ii, ])
    }
    information[j, k] <- 0.5 * sum(diag(P %*% dVj %*% P %*% dVk))
  }
}
if (max(abs(got$score - score)) > 1e-8 ||
    max(abs(got$information - information)) > 1e-8) {
  stop("Cached Fisher score or expected information does not match the dense calculation.")
}

cat("gammfast cached Gaussian projection regression checks passed.\n")

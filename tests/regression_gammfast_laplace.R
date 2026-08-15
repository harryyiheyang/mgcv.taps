library(mgcv)
library(mgcv.taps)
library(gamm4)

set.seed(20260815)
nid <- 80L
m <- 8L
id <- factor(rep(seq_len(nid), each = m))
ii <- as.integer(id)
x <- runif(nid * m, -1, 1)
u0 <- rnorm(nid, sd = 1.2)

for (family_name in c("binary", "poisson")) {
  if (family_name == "binary") {
    eta0 <- -0.3 + 0.65 * sin(pi * x) + u0[ii]
    y <- rbinom(length(eta0), 1, plogis(eta0))
    family <- binomial()
  } else {
    eta0 <- 0.15 + 0.4 * sin(pi * x) + u0[ii]
    y <- rpois(length(eta0), exp(eta0))
    family <- poisson()
  }
  dat <- data.frame(y = y, x = x, id = id)
  fit <- gammfast(
    y ~ s(x, k = 5) + s(id, bs = "re"),
    data = dat, family = family,
    control = list(
      objective.tol = 1e-7, max.outer = 2000L
    )
  )
  fit4 <- gamm4(
    y ~ s(x, k = 5), random = ~(1 | id),
    data = dat, family = family,
    control = lme4::glmerControl(
      optimizer = "bobyqa", tolPwrss = 1e-12,
      optCtrl = list(maxfun = 200000L, rhobeg = 0.1, rhoend = 1e-10)
    )
  )
  vc <- as.data.frame(lme4::VarCorr(fit4$mer))
  G4 <- vc$vcov[vc$grp == "id"]
  eta4 <- stats::predict(fit4$mer)
  if (!fit$converged) {
    stop("A shared-UID Laplace backend did not converge.")
  }
  if (fit$covariance.method !=
      "cached-ordinary-X-Laplace-influence-fixedpoint") {
    stop("The shared-UID Laplace method label is inconsistent.")
  }
  if (abs(log(fit$G[1, 1] / G4)) > 0.08) {
    stop("The shared-UID covariance is not close to gamm4.")
  }
  if (abs(log(fit$sp / fit4$gam$sp)) > 0.25) {
    stop("The shared-UID smoothing parameter is not close to gamm4.")
  }
  if (max(abs(fit$linear.predictors - eta4)) > 0.16) {
    stop("The shared-UID linear predictor is not close to gamm4.")
  }
}

set.seed(20260816)
nid <- 7L
m <- 4L
n <- nid * m
id <- rep(seq_len(nid), each = m)
Xp <- cbind(runif(n, -1, 1), rnorm(n))
B <- cbind(1, runif(n, -0.5, 0.5))
G <- matrix(c(0.8, 0.15, 0.15, 0.45), 2, 2)
phi <- 2.3
w <- runif(n, 0.15, 1.4)
wd <- runif(n, 0.2, 1.2)
dw <- rnorm(n, sd = 0.2)
u <- matrix(rnorm(nid * 2L, sd = 0.4), nid, 2L)
working_cache <- mgcv.taps:::gammfast_gaussian_cache(
  cbind(Xp * sqrt(w), 0), B * sqrt(w), id, n_threads = 2L
)
determinant_cache <- mgcv.taps:::gammfast_gaussian_cache(
  cbind(Xp * sqrt(wd), 0), B * sqrt(wd), id, n_threads = 2L
)
working_mm <- mgcv.taps:::gammfast_gaussian_projected_cached(
  working_cache$AtA, working_cache$BtB, working_cache$BtA,
  G, phi, c(1L, 1L)
)
determinant_crossprod <- mgcv.taps:::gammfast_gaussian_crossprod_cached(
  determinant_cache$AtA, determinant_cache$BtB,
  determinant_cache$BtA, G, n_threads = 2L
)
determinant_mean_covariance <- phi * solve(
  determinant_crossprod$crossprod[seq_len(ncol(Xp)),
                                  seq_len(ncol(Xp)), drop = FALSE]
)
got <- mgcv.taps:::gammfast_laplace_influence_cached(
  X = Xp, B = B, id = id, G = G,
  working_BtB = working_cache$BtB,
  working_BtA = working_cache$BtA,
  working_mean_covariance = working_mm$mean_covariance,
  determinant_BtB = determinant_cache$BtB,
  determinant_BtA = determinant_cache$BtA,
  determinant_mean_covariance = determinant_mean_covariance,
  determinant_derivative = dw, u = u, scale = phi, n_threads = 2L
)

ns <- ncol(Xp)
q <- ncol(B)
Z <- matrix(0, n, ns + nid * q)
Z[, seq_len(ns)] <- Xp
Q <- matrix(0, ncol(Z), ncol(Z))
for (g in seq_len(nid)) {
  jj <- which(id == g)
  kk <- ns + (g - 1L) * q + seq_len(q)
  Z[jj, kk] <- sqrt(phi) * B[jj, , drop = FALSE]
  Q[kk, kk] <- solve(G)
}
working_inverse <- solve(crossprod(Z, (w / phi) * Z) + Q)
determinant_inverse <- solve(crossprod(Z, (wd / phi) * Z) + Q)
leverage <- rowSums((Z %*% determinant_inverse) * Z)
a <- 0.5 * dw * leverage / phi
influence <- drop(working_inverse %*% crossprod(Z, a))
influence_sum <- matrix(0, q, q)
for (g in seq_len(nid)) {
  kk <- ns + (g - 1L) * q + seq_len(q)
  ug <- matrix(u[g, ], q, 1L)
  tg <- sqrt(phi) * matrix(influence[kk], q, 1L)
  influence_g <- tcrossprod(tg, ug) + tcrossprod(ug, tg)
  influence_sum <- influence_sum + influence_g
}
if (max(abs(got$influence_sum - influence_sum)) > 1e-10) {
  stop("The cached Laplace influence does not match the dense Hessians.")
}

cat("gammfast mgcv/Laplace variance regressions passed.\n")

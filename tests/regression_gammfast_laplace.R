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
    data = dat, family = family, inner.max = 5L,
    control = list(
      objective.tol = 1e-7, fixedpoint.tol = 1e-7,
      max.outer = 2000L
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
  if (!fit$converged ||
      fit$covariance.method !=
        "mgcv-fREML-shared-UID-Laplace-fixedpoint" ||
      abs(log(fit$G[1, 1] / G4)) > 0.08 ||
      abs(log(fit$sp / fit4$gam$sp)) > 0.25 ||
      max(abs(fit$linear.predictors - eta4)) > 0.12) {
    stop("The mgcv-driven shared-UID fit is not close to gamm4.")
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
R <- matrix(c(1.2, 0.1, 0.1, 0.9), 2, 2)
w <- runif(n, 0.15, 1.4)
dw <- rnorm(n, sd = 0.2)
u <- matrix(rnorm(nid * 2L, sd = 0.4), nid, 2L)
got <- mgcv.taps:::gammfast_laplace_variance_step(
  X_penalized = Xp, B = B, id = id, G = G,
  smooth_precision = R, working_weight = w,
  weight_derivative = dw, u = u, n_threads = 2L
)

ns <- ncol(Xp)
q <- ncol(B)
Z <- matrix(0, n, ns + nid * q)
Z[, seq_len(ns)] <- Xp
Q <- matrix(0, ncol(Z), ncol(Z))
Q[seq_len(ns), seq_len(ns)] <- R
for (g in seq_len(nid)) {
  jj <- which(id == g)
  kk <- ns + (g - 1L) * q + seq_len(q)
  Z[jj, kk] <- B[jj, , drop = FALSE]
  Q[kk, kk] <- solve(G)
}
Hinv <- solve(crossprod(Z, w * Z) + Q)
leverage <- rowSums((Z %*% Hinv) * Z)
a <- 0.5 * dw * leverage
influence <- drop(Hinv %*% crossprod(Z, a))
moment <- matrix(0, q, q)
conditional <- matrix(0, q, q)
mean <- matrix(0, q, q)
influence_sum <- matrix(0, q, q)
for (g in seq_len(nid)) {
  kk <- ns + (g - 1L) * q + seq_len(q)
  ug <- matrix(u[g, ], q, 1L)
  tg <- matrix(influence[kk], q, 1L)
  Hii <- Hinv[kk, kk, drop = FALSE]
  Cg <- solve(solve(G) + crossprod(B[id == g, , drop = FALSE],
                                   w[id == g] * B[id == g, , drop = FALSE]))
  conditional <- conditional + Cg
  mean <- mean + Hii - Cg
  influence_g <- tcrossprod(tg, ug) + tcrossprod(ug, tg)
  influence_sum <- influence_sum + influence_g
  moment <- moment + tcrossprod(ug) + Hii - influence_g
}
if (max(abs(got$moment_sum - moment)) > 1e-10 ||
    max(abs(got$conditional_sum - conditional)) > 1e-10 ||
    max(abs(got$mean_sum - mean)) > 1e-10 ||
    max(abs(got$influence_sum - influence_sum)) > 1e-10) {
  stop("The shared-UID Laplace variance step does not match the dense Hessian.")
}

cat("gammfast mgcv/Laplace variance regressions passed.\n")

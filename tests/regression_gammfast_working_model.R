library(mgcv)
library(mgcv.taps)

set.seed(20260811)
block_size <- c(1L, 2L, 5L, 3L, 8L)
id_block <- rep(seq_along(block_size), block_size)
n <- length(id_block)
A <- matrix(rnorm(n * 3L), n, 3L)
B <- cbind(1, matrix(rnorm(n * 2L), n, 2L))
L <- matrix(c(0.8, 0, 0, 0.2, 0.6, 0, -0.1, 0.15, 0.5), 3L, 3L)
G <- tcrossprod(L)
phi <- 1.9
W0 <- runif(n, 0.4, 2.2)
R_diag <- phi / W0
got <- mgcv.taps:::gammfast_vinv_apply(
  A, B, id_block, G, R_diag, n_threads = 2L
)
want <- matrix(0, n, ncol(A))
for (i in seq_along(block_size)) {
  ii <- which(id_block == i)
  Vi <- diag(R_diag[ii], length(ii)) + B[ii, , drop = FALSE] %*%
    G %*% t(B[ii, , drop = FALSE])
  want[ii, ] <- solve(Vi, A[ii, , drop = FALSE])
}
if (max(abs(got - want)) > 1e-10) {
  stop("The blockwise Woodbury inverse does not match direct dense solves.")
}

ni <- rep(c(3L, 5L, 7L, 4L), 12L)
id <- factor(rep(seq_along(ni), ni))
time <- unlist(lapply(ni, function(m) seq(0, 1, length.out = m)))
x <- runif(length(id))
off <- 0.3 * sin(2 * pi * x)
u <- rnorm(nlevels(id), sd = 0.2)
y <- 1 + 0.4 * x + 0.5 * sin(2 * pi * x) +
  u[id] + rnorm(length(id), sd = 0.5)
d0 <- data.frame(y = y, x = x, time = time, id = id)
d1 <- data.frame(y = y + off, x = x, time = time, id = id, off = off)
ctl <- list(max.outer = 800L, objective.tol = 1e-6)
f0 <- gammfast(
  y ~ s(x, k = 5L) + s(id, bs = "re") + fs(time, id, k = 3L),
  data = d0, inner.max = 5L, nthreads = 2L, control = ctl
)
f1 <- gammfast(
  y ~ offset(off) + s(x, k = 5L) + s(id, bs = "re") +
    fs(time, id, k = 3L),
  data = d1, inner.max = 5L, nthreads = 2L, control = ctl
)
p0 <- taps_score_test_gamm(f0, method = "liu", n_threads = 2L)$smooth.pvalue
p1 <- taps_score_test_gamm(f1, method = "liu", n_threads = 2L)$smooth.pvalue
if (max(abs(f0$coefficients - f1$coefficients)) > 1e-4 ||
    max(abs(f0$G - f1$G)) > 1e-4 ||
    abs(f0$sig2 - f1$sig2) > 1e-4 || abs(p0 - p1) > 1e-5) {
  stop("Gaussian gammfast is not invariant to an exactly compensated offset.")
}

set.seed(20260812)
ni <- rep(c(3L, 5L, 7L, 4L), 15L)
id <- factor(rep(seq_along(ni), ni))
time <- unlist(lapply(ni, function(m) seq(0, 1, length.out = m)))
x <- runif(length(id))
off <- -0.25 * cos(2 * pi * x)
u <- rnorm(nlevels(id), sd = 0.15)
mu <- exp(
  off + 0.2 + 0.25 * x + 0.5 * sin(2 * pi * x) +
    0.15 * sin(2 * pi * time) + u[id]
)
y <- rTweedie(mu, p = 1.5, phi = 1.8)
d <- data.frame(y = y, x = x, time = time, id = id, off = off)
fit.tw <- gammfast(
  y ~ offset(off) + s(x, k = 5L) + s(id, bs = "re") +
    fs(time, id, k = 3L),
  data = d, family = tw(theta = 1.5, link = "log"),
  inner.max = 5L, nthreads = 2L,
  control = list(max.outer = 1000L, objective.tol = 1e-5)
)
if (!isTRUE(fit.tw$converged)) stop("Tweedie gammfast fixture did not converge.")
variance.tw <- gammfast_variance_test(
  fit.tw, method = "liu", spectrum = "exact", n_threads = 2L
)
if (f0$dispersion.method != "mgcv-fREML" ||
    f0$covariance.method != "mean-Hessian-projected-moment") {
  stop("Gaussian dispersion was not isolated from the covariance fixed point.")
}
if (!is.finite(variance.tw$statistic) ||
    !is.finite(variance.tw$p.value) ||
    variance.tw$working.model != "final PIRLS working Gaussian") {
  stop("Tweedie gammfast variance testing did not use the final PIRLS system.")
}
random.tw <- mgcv.taps:::gammfast_random_info(fit.tw)
E <- mgcv.taps:::extract_pseudo_response(fit.tw)
if (abs(E$phi0 - fit.tw$sig2) > 1e-12 ||
    max(abs(E$V_phi * (fit.tw$weights / fit.tw$sig2) - 1)) > 1e-7) {
  stop("Tweedie gammfast does not expose one consistent actual-scale working model.")
}

X.tw <- predict(
  fit.tw$global, newdata = fit.tw$global$model, type = "lpmatrix"
)
work.tw <- mgcv.taps:::gammfast_working(
  fit.tw$family, fit.tw$y, fit.tw$linear.predictors,
  fit.tw$prior.weights
)
covariance_group <- integer(ncol(random.tw$B))
for (j in seq_along(random.tw$group.index)) {
  covariance_group[random.tw$group.index[[j]]] <- j
}
mapped <- mgcv.taps:::gammfast_non_gaussian_inner_G(
  X = X.tw, B = random.tw$B, id = random.tw$id.index,
  offset = fit.tw$offset, G = fit.tw$G / fit.tw$sig2,
  work = work.tw,
  t_correction = mgcv.taps:::gammfast_t_correction(
    fit.tw$family, fit.tw$y, fit.tw$linear.predictors,
    fit.tw$prior.weights, work.tw$w
  ),
  scale = fit.tw$sig2, covariance_group = covariance_group,
  group_index = random.tw$group.index,
  ng = nrow(fit.tw$random.effects), inner_tol = 0, inner_max = 1L,
  nthreads = 2L
)
G_check <- mapped$G
G_normalized <- fit.tw$G / fit.tw$sig2
G_error <- norm(G_check - G_normalized, "F") /
  (1 + norm(G_normalized, "F"))
if (G_error > 2e-4) {
  stop("Tweedie gammfast covariance is not at the normalized EM fixed point.")
}

V0_inv_apply_actual <- function(v) {
  is_matrix <- is.matrix(v)
  v <- if (is_matrix) v else matrix(v, ncol = 1L)
  ans <- mgcv.taps:::gammfast_vinv_apply(
    v, random.tw$B, random.tw$id.index, fit.tw$G, E$V_phi,
    n_threads = 2L
  )
  if (is_matrix) ans else as.vector(ans)
}
work.score <- work.tw
sw.score <- sqrt(work.score$w / E$phi0)
X.actual <- predict(
  fit.tw$global, newdata = fit.tw$global$model, type = "lpmatrix"
)
X.score <- X.actual * sw.score
B.score <- random.tw$B * sw.score
z.actual <- work.score$z - fit.tw$offset
z.score <- z.actual * sw.score
V0_inv_apply_score <- function(v) {
  is_matrix <- is.matrix(v)
  v <- if (is_matrix) v else matrix(v, ncol = 1L)
  ans <- mgcv.taps:::gammfast_vinv_apply(
    v, B.score, random.tw$id.index, fit.tw$G,
    rep(1, length(fit.tw$y)), n_threads = 2L
  )
  if (is_matrix) ans else as.vector(ans)
}
probe <- cbind(X.score[, seq_len(min(3L, ncol(X.score))), drop = FALSE],
               z.score)
actual_probe <- V0_inv_apply_actual(probe / sw.score) / sw.score
score_probe <- V0_inv_apply_score(probe)
if (max(abs(actual_probe - score_probe)) > 1e-9) {
  stop("Actual-scale and standardized gammfast subject inverses differ.")
}
p_manual <- mgcv.taps:::gammfast_smoothing_score_core(
  g = fit.tw$global, X = X.score, working_response = z.score,
  V0_inv_apply = V0_inv_apply_score, phi0 = E$phi0,
  test.component = 1L,
  null.tol = 1e-10, method = "liu", max_eps = 1e-8,
  max_iter = 1e5
)$smooth.pvalue
p_package <- taps_score_test_gamm(
  fit.tw, method = "liu", n_threads = 2L
)$smooth.pvalue
if (!isTRUE(all.equal(p_package, p_manual, tolerance = 1e-12))) {
  stop("The gammfast score path does not subtract the model offset exactly once.")
}
p_actual <- mgcv.taps:::gammfast_smoothing_score_core(
  g = fit.tw$global, X = X.actual, working_response = z.actual,
  V0_inv_apply = V0_inv_apply_actual, phi0 = E$phi0,
  test.component = 1L, null.tol = 1e-10, method = "liu",
  max_eps = 1e-8, max_iter = 1e5
)$smooth.pvalue
if (!isTRUE(all.equal(p_package, p_actual, tolerance = 1e-10))) {
  stop("Actual-scale and standardized gammfast score tests differ.")
}

p_missing_offset <- mgcv.taps:::gammfast_smoothing_score_core(
  g = fit.tw$global, X = X.score,
  working_response = work.score$z * sw.score,
  V0_inv_apply = V0_inv_apply_score, phi0 = E$phi0,
  test.component = 1L,
  null.tol = 1e-10, method = "liu", max_eps = 1e-8,
  max_iter = 1e5
)$smooth.pvalue
if (abs(p_missing_offset - p_package) < 1e-3) {
  stop("The Tweedie fixture does not detect a missing gammfast score offset.")
}

set.seed(20260813)
mu.nb <- exp(off + 0.1 + 0.2 * x + 0.5 * sin(2 * pi * x) + u[id])
y.nb <- rnbinom(length(id), mu = mu.nb, size = 2.5)
family.nb <- nb(theta = 2.5, link = "log")
eta.nb <- log(mu.nb)
work.nb <- mgcv.taps:::gammfast_working(
  family.nb, y.nb, eta.nb, rep(1, length(y.nb))
)
D.nb <- family.nb$Dd(
  y.nb, mu.nb, family.nb$getTheta(), rep(1, length(y.nb)), level = 0
)
W0.nb <- D.nb$EDmu2 * family.nb$mu.eta(eta.nb)^2 / 2
if (max(abs(work.nb$w - W0.nb)) > 1e-10 ||
    any(!is.finite(work.nb$z))) {
  stop("Negative-binomial gammfast working-scale regression failed.")
}

set.seed(911)
ni.bin <- rep(6L, 40L)
id.bin <- factor(rep(seq_along(ni.bin), ni.bin))
t.bin <- rep(seq(0, 1, length.out = 6L), 40L)
x.bin <- runif(length(id.bin))
u.bin <- rnorm(40L, sd = 0.45)
eta.bin <- -0.2 + sin(2 * pi * x.bin) + u.bin[id.bin] +
  0.2 * cos(pi * t.bin)
y.bin <- rbinom(length(id.bin), 1L, plogis(eta.bin))
d.bin <- data.frame(y = y.bin, x = x.bin, t = t.bin, id = id.bin)
fit.bin <- gammfast(
  y ~ s(x, k = 5L) + s(id, bs = "re") + fs(t, id, k = 2L),
  data = d.bin, family = binomial(), inner.max = 5L, nthreads = 2L,
  control = list(max.outer = 1000L, objective.tol = 1e-5)
)
if (!isTRUE(fit.bin$converged)) {
  stop("Binomial gammfast variance-test fixture did not converge.")
}
variance.bin <- gammfast_variance_test(
  fit.bin, method = "liu", spectrum = "exact", n_threads = 2L
)
if (!is.finite(variance.bin$statistic) ||
    !is.finite(variance.bin$p.value) ||
    variance.bin$family != "binomial" ||
    variance.bin$working.model != "final PIRLS working Gaussian") {
  stop("Binomial variance testing did not report its uncorrected PIRLS calibration.")
}

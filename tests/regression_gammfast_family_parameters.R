library(mgcv)
library(mgcv.taps)

set.seed(20260814)
nid <- 60L
m <- 4L
id <- factor(rep(seq_len(nid), each = m))
x <- stats::runif(nid * m, -1, 1)
u <- stats::rnorm(nid, sd = 0.45)
mu <- exp(0.25 + 0.35 * sin(pi * x) + u[as.integer(id)])
y <- stats::rnbinom(length(mu), mu = mu, size = 1.8)
dat <- data.frame(y = y, x = x, id = id)

fit.qp <- gammfast(
  y ~ s(x, k = 5L) + s(id, bs = "re"),
  data = dat, family = stats::quasipoisson(),
  inner.max = 10L,
  control = list(
    max.outer = 1000L, objective.tol = 1e-5,
    fixedpoint.tol = 1e-5
  )
)
if (!isTRUE(fit.qp$converged) || !is.finite(fit.qp$sig2) ||
    fit.qp$sig2 <= 1 ||
    fit.qp$dispersion.method != "mgcv-fREML" ||
    fit.qp$covariance.method !=
      "mgcv-fREML-shared-UID-Laplace-fixedpoint") {
  stop("Quasi-Poisson scale handling is inconsistent.")
}
if (!is.null(fit.qp$sigma2) || !is.null(fit.qp$G.normalized) ||
    !is.null(fit.qp$fs)) {
  stop("Deprecated duplicate gammfast fields are still present.")
}
if (max(abs(fit.qp$random$covariance[[1L]] - fit.qp$G)) > 1e-12) {
  stop("Quasi-Poisson G is not on the public random-effect scale.")
}

random.qp <- mgcv.taps:::gammfast_random_info(fit.qp)
work.qp <- mgcv.taps:::gammfast_working(
  fit.qp$family, fit.qp$y, fit.qp$linear.predictors,
  fit.qp$prior.weights
)
sw.qp <- sqrt(work.qp$w)
X.qp <- predict(
  fit.qp$global, newdata = fit.qp$global$model, type = "lpmatrix"
)
mm.qp <- mgcv.taps:::gammfast_projected_moments(
  response = sw.qp * (work.qp$z - fit.qp$offset) / sqrt(fit.qp$sig2),
  X = X.qp * sw.qp / sqrt(fit.qp$sig2),
  B = random.qp$B * sw.qp, id = random.qp$id.index,
  G = fit.qp$G / fit.qp$sig2,
  penalty = mgcv.taps:::gammfast_penalty_matrix(
    fit.qp$global, fit.qp$sp, scale = fit.qp$sig2
  )
)
if (max(abs(fit.qp$random.effects - sqrt(fit.qp$sig2) * mm.qp$u)) >
    1e-7) {
  stop("Quasi-Poisson random effects are not on the actual dispersion scale.")
}

family.nb <- mgcv::nb(theta = -1, link = "log")
fit.nb <- gammfast(
  y ~ s(x, k = 5L) + s(id, bs = "re"),
  data = dat, family = family.nb, inner.max = 10L,
  control = list(
    max.outer = 1000L, objective.tol = 1e-5,
    fixedpoint.tol = 1e-5
  )
)
theta.nb <- fit.nb$family$getTheta(TRUE)
if (!isTRUE(fit.nb$converged) || !is.finite(theta.nb) || theta.nb <= 0 ||
    fit.nb$sig2 != 1 ||
    fit.nb$family.parameter.method != "mgcv-estimate.theta") {
  stop("Negative-binomial mgcv family-parameter handling is inconsistent.")
}
if (!is.null(fit.nb$sigma2) || !is.null(fit.nb$G.normalized) ||
    !is.null(fit.nb$fs)) {
  stop("Deprecated duplicate gammfast fields are still present.")
}

cat("gammfast family-parameter regressions passed.\n")

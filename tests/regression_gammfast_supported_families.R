library(mgcv)
library(mgcv.taps)

set.seed(20260814)
nid <- 35L
nr <- 5L
id <- factor(rep(seq_len(nid), each = nr))
x <- stats::runif(nid * nr, -1, 1)
u <- stats::rnorm(nid, sd = 0.35)
eta <- 0.2 + 0.3 * sin(pi * x) + u[as.integer(id)]
mu_log <- exp(eta)
mu_cloglog <- 1 - exp(-exp(eta))
set.seed(20260815)
y_quasibinomial <- stats::rbinom(
  length(eta), 1, stats::plogis(2 * eta)
)
set.seed(20260816)
y_nb_sqrt <- stats::rnbinom(
  length(eta), mu = exp(2 * eta), size = 2
)

lambda <- 3
v <- stats::rnorm(length(mu_log))^2
y0 <- mu_log + mu_log^2 * v / (2 * lambda) -
  mu_log * sqrt(
    4 * mu_log * lambda * v + mu_log^2 * v^2
  ) / (2 * lambda)
z0 <- stats::runif(length(mu_log))
y_inverse_gaussian <- ifelse(
  z0 <= mu_log / (mu_log + y0), y0, mu_log^2 / y0
)

cases <- list(
  list(
    name = "gaussian-identity",
    y = eta + stats::rnorm(length(eta), sd = 0.5),
    family = stats::gaussian()
  ),
  list(
    name = "gaussian-log",
    y = pmax(mu_log + stats::rnorm(length(eta), sd = 0.25), 0.05),
    family = stats::gaussian(link = "log")
  ),
  list(
    name = "binomial-cloglog",
    y = stats::rbinom(length(eta), 1, mu_cloglog),
    family = stats::binomial(link = "cloglog")
  ),
  list(
    name = "quasibinomial-logit",
    y = y_quasibinomial,
    family = stats::quasibinomial(link = "logit")
  ),
  list(
    name = "poisson-sqrt",
    y = stats::rpois(length(eta), mu_log),
    family = stats::poisson(link = "sqrt")
  ),
  list(
    name = "quasipoisson-log",
    y = stats::rnbinom(length(eta), mu = mu_log, size = 2),
    family = stats::quasipoisson(link = "log")
  ),
  list(
    name = "gamma-log",
    y = stats::rgamma(length(eta), shape = 3, scale = mu_log / 3),
    family = stats::Gamma(link = "log")
  ),
  list(
    name = "inverse-gaussian-log",
    y = y_inverse_gaussian,
    family = stats::inverse.gaussian(link = "log")
  ),
  list(
    name = "quasi-mu2",
    y = stats::rgamma(length(eta), shape = 3, scale = mu_log / 3),
    family = stats::quasi(link = "log", variance = "mu^2")
  ),
  list(
    name = "negbin-fixed",
    y = stats::rnbinom(length(eta), mu = mu_log, size = 2),
    family = mgcv::negbin(theta = 2, link = "log")
  ),
  list(
    name = "nb-sqrt",
    y = y_nb_sqrt,
    family = mgcv::nb(theta = 2, link = "sqrt")
  ),
  list(
    name = "Tweedie-fixed",
    y = mgcv::rTweedie(mu_log, p = 1.5, phi = 0.5),
    family = mgcv::Tweedie(p = 1.5, link = "log")
  ),
  list(
    name = "tw-estimated",
    y = mgcv::rTweedie(mu_log, p = 1.5, phi = 0.5),
    family = mgcv::tw(theta = -1.5, link = "log")
  ),
  list(
    name = "beta-probit",
    y = stats::rbeta(
      length(eta), stats::pnorm(eta) * 15,
      (1 - stats::pnorm(eta)) * 15
    ),
    family = mgcv::betar(theta = 15, link = "probit")
  ),
  list(
    name = "scaled-t",
    y = eta + 0.7 * stats::rt(length(eta), df = 5),
    family = mgcv::scat(theta = c(5, 1), link = "identity")
  )
)

results <- data.frame()
for (case in cases) {
  dat <- data.frame(y = case$y, x = x, id = id)
  fit <- gammfast(
    y ~ s(x, k = 5L) + s(id, bs = "re"),
    data = dat, family = case$family,
    control = list(
      max.outer = 200L, objective.tol = 1e-5
    )
  )
  if (!isTRUE(fit$converged) || any(!is.finite(fit$G)) ||
      min(eigen(fit$G, symmetric = TRUE, only.values = TRUE)$values) < -1e-10 ||
      !is.finite(fit$sig2) || fit$sig2 <= 0) {
    stop(case$name, " gammfast smoke fit failed.")
  }
  if (!identical(case$name, "gaussian-identity") &&
      fit$covariance.method !=
        "cached-ordinary-X-Laplace-influence-fixedpoint") {
    stop(case$name, " did not use the Laplace fixed-point path.")
  }
  if (case$name %in% c("gamma-log", "inverse-gaussian-log")) {
    work <- mgcv.taps:::gammfast_working(
      fit$family, fit$y, fit$linear.predictors, fit$prior.weights
    )
    correction <- mgcv.taps:::gammfast_t_correction(
      fit$family, fit$y, fit$linear.predictors,
      fit$prior.weights, work$w
    )
    if (is.null(correction) || any(correction$weight < 0) ||
        max(abs(correction$weight - work$w)) < 1e-8 ||
        any(!is.finite(correction$derivative))) {
      stop(case$name, " determinant correction is invalid.")
    }
  }
  results <- rbind(
    results,
    data.frame(
      family = case$name, G = fit$G[1, 1], sig2 = fit$sig2,
      outer = fit$outer,
      local = if (is.null(fit$local.updates)) 0L else fit$local.updates
    )
  )
}

print(results, row.names = FALSE)
cat("gammfast supported-family regressions passed.\n")

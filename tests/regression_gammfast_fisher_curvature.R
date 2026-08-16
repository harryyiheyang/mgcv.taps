library(mgcv)
library(mgcv.taps)

if ("influence.update" %in% names(formals(gammfast))) {
  stop("The obsolete influence-update API is still exposed.")
}

family <- stats::inverse.gaussian(link = "log")
y_exact <- c(0.2, 0.4, 0.6, 2)
eta_exact <- rep(0, length(y_exact))
prior <- rep(1, length(y_exact))
work <- mgcv.taps:::gammfast_working(family, y_exact, eta_exact, prior)
correction <- mgcv.taps:::gammfast_t_correction(
  family, y_exact, eta_exact, prior, work$w
)

if (!isTRUE(all.equal(correction$weight, work$w, tolerance = 1e-14)) ||
    any(!is.finite(correction$derivative))) {
  stop("The determinant correction did not preserve Fisher curvature.")
}

rinvgauss <- function(mu, phi) {
  lambda <- 1 / phi
  v <- stats::rnorm(length(mu))^2
  candidate <- mu + mu^2 * v / (2 * lambda) -
    mu * sqrt(4 * mu * lambda * v + mu^2 * v^2) / (2 * lambda)
  choose_candidate <- stats::runif(length(mu)) <= mu / (mu + candidate)
  ifelse(choose_candidate, candidate, mu^2 / candidate)
}

set.seed(20260816)
n_id <- 30L
n_per_id <- 6L
id <- factor(rep(seq_len(n_id), each = n_per_id))
x <- stats::runif(length(id), -1, 1)
u <- stats::rnorm(n_id, sd = 0.3)
eta <- 0.1 + 0.25 * sin(pi * x) + u[as.integer(id)]
dat <- data.frame(y = rinvgauss(exp(eta), phi = 0.2), x = x, id = id)

fit <- gammfast(
  y ~ s(x, k = 5L) + s(id, bs = "re"),
  data = dat, family = stats::inverse.gaussian(link = "log"),
  inner.max = 5L,
  control = list(max.outer = 200L, objective.tol = 1e-5)
)
if (!isTRUE(fit$converged) || any(!is.finite(fit$G))) {
  stop("The Fisher-curvature inverse-Gaussian fit failed.")
}

cat("gammfast Fisher-curvature regressions passed.\n")

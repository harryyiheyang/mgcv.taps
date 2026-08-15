library(mgcv)
library(mgcv.taps)

set.seed(20260810)
x <- seq(0, 1, length.out = 240L)
off <- 0.8 * sin(2 * pi * x) + 0.3 * cos(4 * pi * x)
y <- 1 + 0.5 * x + rnorm(length(x), sd = 0.7)
dat0 <- data.frame(y = y, x = x)
dat1 <- data.frame(y = y + off, x = x, off = off)

fit0 <- gam(y ~ s(x, k = 10L), data = dat0, method = "REML")
fit1 <- gam(y ~ offset(off) + s(x, k = 10L), data = dat1, method = "REML")
p0 <- taps_score_test(fit0, method = "davies")$smooth.pvalue
p1 <- taps_score_test(fit1, method = "davies")$smooth.pvalue
if (!isTRUE(all.equal(p0, p1, tolerance = 1e-7))) {
  stop(sprintf(
    paste0("Score p-values are not invariant to an exactly compensated model offset: ",
           "p0=%.16g, p1=%.16g, sig2 diff=%.3g, sp diff=%.3g, coef diff=%.3g."),
    p0, p1, abs(fit0$sig2 - fit1$sig2),
    max(abs(fit0$sp - fit1$sp)),
    max(abs(fit0$coefficients - fit1$coefficients))
  ))
}
fit1.missing.offset <- fit1
fit1.missing.offset$offset <- rep(0, length(x))
p1.missing.offset <- taps_score_test(
  fit1.missing.offset, method = "davies"
)$smooth.pvalue
if (abs(p1.missing.offset - p0) < 0.02) {
  stop("The offset regression fixture does not detect the historical missing-offset bug.")
}

fit0.bam <- bam(
  y ~ s(x, k = 10L), data = dat0,
  method = "fREML", discrete = FALSE
)
fit1.bam <- bam(
  y ~ offset(off) + s(x, k = 10L), data = dat1,
  method = "fREML", discrete = FALSE
)
p0.bam <- taps_score_test(fit0.bam, method = "liu")$smooth.pvalue
p1.bam <- taps_score_test(fit1.bam, method = "liu")$smooth.pvalue
if (!isTRUE(all.equal(p0.bam, p1.bam, tolerance = 1e-7))) {
  stop("Non-discrete bam score p-values are not invariant to a compensated offset.")
}

check_working_identity <- function(fit, label) {
  E <- mgcv.taps:::extract_pseudo_response(fit)
  eta <- as.numeric(fit$linear.predictors)
  mu <- as.numeric(fit$fitted.values)
  if (inherits(fit$family, "extended.family") &&
      is.function(fit$family$Dd)) {
    theta <- fit$family$getTheta()
    Dval <- fit$family$Dd(
      fit$y, mu, theta, fit$prior.weights, level = 0
    )
    W0 <- Dval$EDmu2 * fit$family$mu.eta(eta)^2 / 2
  } else {
    W0 <- fit$prior.weights * fit$family$mu.eta(eta)^2 /
      fit$family$variance(mu)
  }
  phi <- E$phi0
  if (max(abs(fit$weights - W0)) > 1e-7 ||
      max(abs(E$V_phi * fit$weights - phi)) > 1e-7 ||
      max(abs(E$V_phi * (fit$weights / phi) - 1)) > 1e-7) {
    stop(label, " does not satisfy the scale-free working-weight identities.")
  }
}

check_working_identity(fit0, "Gaussian gam")
check_working_identity(fit0.bam, "Gaussian bam")

mu <- exp(-0.7 + 0.4 * x)
y.tw <- rTweedie(mu, p = 1.5, phi = 2.5)
fit.tw <- gam(
  y.tw ~ s(x, k = 10L), family = tw(theta = 1.5, link = "log"),
  method = "REML"
)
E <- mgcv.taps:::extract_pseudo_response(fit.tw)
if (abs(E$phi0 - fit.tw$sig2) > 1e-10 ||
    max(abs(E$V_phi * fit.tw$weights - fit.tw$sig2)) > 1e-7) {
  stop("Extended-family working variance does not include fitted dispersion.")
}
fit.tw.bam <- bam(
  y.tw ~ s(x, k = 10L), family = tw(theta = 1.5, link = "log"),
  method = "fREML", discrete = FALSE
)
check_working_identity(fit.tw, "Tweedie gam")
check_working_identity(fit.tw.bam, "Tweedie bam")

y.nb <- rnbinom(length(x), mu = mu, size = 2.5)
fit.nb <- gam(
  y.nb ~ s(x, k = 10L), family = nb(theta = 2.5, link = "log"),
  method = "REML"
)
fit.nb.bam <- bam(
  y.nb ~ s(x, k = 10L), family = nb(theta = 2.5, link = "log"),
  method = "fREML", discrete = FALSE
)
check_working_identity(fit.nb, "negative-binomial gam")
check_working_identity(fit.nb.bam, "negative-binomial bam")

fit.discrete <- bam(
  y ~ s(x, k = 10L), data = dat0,
  method = "fREML", discrete = TRUE
)
discrete_result <- taps_score_test(fit.discrete)
if (!is.finite(discrete_result$smooth.pvalue)) {
  stop("bam(discrete = TRUE) did not produce a finite post-estimation result.")
}

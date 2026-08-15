library(mgcv)
library(mgcv.taps)

set.seed(20260815)
n <- 220L
dat <- data.frame(
  x = runif(n), z = runif(n), by = runif(n, 0.5, 1.5),
  offset = rnorm(n, sd = 0.1)
)
dat$y <- 1 + dat$x + sin(2 * pi * dat$z) + dat$offset +
  rnorm(n, sd = 0.35)

getA <- function(x, para) cbind(1, x)
form <- y ~ s(x, bs = "AMatern", k = 10,
              xt = list(getA = getA, para = NULL)) +
  s(z, k = 7) + offset(offset)

for (engine in c("gam", "bam", "bam-discrete")) {
  fit <- if (engine == "gam") {
    gam(form, data = dat, method = "REML")
  } else {
    bam(form, data = dat, method = "fREML",
        discrete = engine == "bam-discrete")
  }
  ordinary <- taps_score_test(fit, method = "liu")
  explicit <- taps_score_test(fit, method = "liu", refit = FALSE)
  if (!identical(ordinary, explicit)) {
    stop("refit = FALSE changed the existing score-test path.")
  }

  refitted <- taps_score_test(fit, method = "liu", refit = TRUE)
  fit0 <- mgcv.taps:::taps_score_refit(fit, 1L, FALSE, 1e-10)
  expected_method <- if (engine == "gam") "REML" else "fREML"
  if (!is.finite(refitted$smooth.pvalue) ||
      fit0$call$method != expected_method ||
      (!is.null(fit0$dinfo)) != (engine == "bam-discrete")) {
    stop("The null refit did not preserve the requested mgcv engine.")
  }
  if (engine == "bam-discrete") {
    fit0.cold <- mgcv.taps:::taps_score_refit(fit, 1L, TRUE, 1e-10)
    if (!is.null(fit0.cold$call$coef) || !is.null(fit0.cold$call$in.out)) {
      stop("A discrete bam refit should not use fitted starting values.")
    }
  }
  if (max(abs(fit0$offset - fit$offset)) > 1e-12 ||
      length(fit0$coefficients) != ncol(fit0$.taps_score_X) ||
      max(abs(fit0$sp[-length(fit0$sp)] - fit$sp[2L])) > 1e-8) {
    stop("The null refit did not retain offset, dimensions, or nuisance sp.")
  }
}

getA.vector <- function(x, para) x
fit.vector <- gam(
  y ~ s(x, bs = "AMatern", k = 9,
        xt = list(getA = getA.vector, para = NULL)) + s(z, k = 7),
  data = dat, method = "REML"
)
if (length(predict(fit.vector, newdata = dat[1:5, ])) != 5L) {
  stop("A vector-valued getA was not normalized to a matrix.")
}

getA.scalar <- function(x, para) 1
fit.scalar <- gam(
  y ~ s(x, bs = "AMatern", k = 9,
        xt = list(getA = getA.scalar, para = NULL)),
  data = dat, method = "REML"
)
if (!is.finite(taps_score_test(
  fit.scalar, method = "liu", refit = TRUE
)$smooth.pvalue)) {
  stop("An intercept-only getA refit failed.")
}

fit.by <- gam(
  y ~ s(x, by = by, bs = "AMatern", k = 10,
        xt = list(getA = getA, para = NULL)) + s(z, k = 7),
  data = dat, method = "REML"
)
if (!is.finite(taps_score_test(
  fit.by, method = "liu", refit = TRUE
)$smooth.pvalue)) {
  stop("A numeric-by AMatern refit failed.")
}

getA.complex <- function(x, para) cbind(1, x, x^2)
fit.position <- gam(
  y ~ s(z, k = 7) + s(x, bs = "AMatern", k = 11,
                       xt = list(getA = getA.complex, para = NULL)),
  data = dat, method = "REML"
)
fit.position0 <- mgcv.taps:::taps_score_refit(
  fit.position, 2L, TRUE, 1e-10
)
if (fit.position0$smooth[[2L]]$label != fit.position$smooth[[2L]]$label ||
    !is.finite(taps_score_test(
      fit.position, test.component = 2L, method = "liu", refit = TRUE
    )$smooth.pvalue)) {
  stop("A complex getA was not restored at test.component.")
}

dat2 <- data.frame(x1 = runif(n), x2 = runif(n), z = runif(n))
dat2$y <- dat2$x1 - dat2$x2 + dat2$x1 * dat2$x2 + rnorm(n, sd = 0.4)
getA2 <- function(x1, x2, para) cbind(1, x1, x2, x1 * x2)
fit2 <- gam(
  y ~ s(x1, x2, bs = "A2Matern", k = 13,
        xt = list(getA = getA2, para = NULL)) + s(z, k = 7),
  data = dat2, method = "REML"
)
if (!is.finite(taps_score_test(
  fit2, method = "liu", refit = TRUE
)$smooth.pvalue)) {
  stop("An A2Matern refit failed.")
}

dat$yp <- rpois(n, exp(0.2 + 0.4 * dat$x + 0.3 * sin(2 * pi * dat$z)))
fit.poisson <- bam(
  yp ~ s(x, bs = "AMatern", k = 10,
         xt = list(getA = getA, para = NULL)) + s(z, k = 7),
  data = dat, family = poisson(), method = "fREML", discrete = TRUE
)
if (!is.finite(taps_score_test(
  fit.poisson, method = "liu", refit = TRUE
)$smooth.pvalue)) {
  stop("A non-Gaussian discrete-bam refit failed.")
}

dat$n_trial <- 4L
dat$count <- rbinom(n, dat$n_trial, plogis(dat$x + sin(2 * pi * dat$z)))
fit.count <- gam(
  cbind(count, n_trial - count) ~
    s(x, bs = "AMatern", k = 10,
      xt = list(getA = getA, para = NULL)) + s(z, k = 7),
  data = dat, family = binomial(), method = "REML"
)
fit.count0 <- mgcv.taps:::taps_score_refit(fit.count, 1L, TRUE, 1e-10)
fit.count.explicit <- gam(
  cbind(count, n_trial - count) ~ x + s(z, k = 7),
  data = dat, family = binomial(), method = "REML"
)
count_predictor_difference <- max(abs(
  fit.count0$linear.predictors - fit.count.explicit$linear.predictors
))
count_sp_difference <- max(abs(
  log(fit.count0$sp[1L]) - log(fit.count.explicit$sp)
))
if (count_predictor_difference > 1e-4 || count_sp_difference > 1e-4) {
  stop(
    "The matrix-binomial null refit differs from the explicit fit: predictor ",
    format(count_predictor_difference, digits = 8), ", log-sp ",
    format(count_sp_difference, digits = 8), "."
  )
}

cat("taps score-test null-refit regressions passed.\n")

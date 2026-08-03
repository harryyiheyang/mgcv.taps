library(mgcv)
library(survival)
devtools::load_all(".")
options(warn = 1)

n <- 1000L
ncl <- 50L
mcl <- 20L
sigma_u <- 0.5
nsim <- 100L
lambda <- 0.1
shape <- 2
censor_rate <- 0.704
pv_null_iid <- pv_null_re <- numeric(nsim)
pv_alt_iid <- pv_alt_re <- numeric(nsim)

# NULL seeds: 5910001--5910100
for (i in seq_len(nsim)) {
  set.seed(5910000L + i)
  cl <- factor(rep(seq_len(ncl), each = mcl))
  x1 <- qbeta(pnorm(rnorm(ncl)), 1.5, 1.5)[as.integer(cl)]
  Z <- MASS::mvrnorm(n, rep(0, 3), matrix(0.25, 3, 3) + 0.75 * diag(3))
  x2 <- qbeta(pnorm(Z[, 1]), 1.5, 1.5)
  x3 <- qbeta(pnorm(Z[, 2]), 1.5, 1.5)
  x4 <- qbeta(pnorm(Z[, 3]), 1.5, 1.5)
  u <- rnorm(ncl, 0, sigma_u)
  f1 <- smoothed_linearity(x1, 0)
  t2 <- 2 * pi * x2
  f2 <- 0.4 * sin(t2) + 0.8 * cos(t2) + 1.2 * sin(t2)^2 +
    1.6 * cos(t2)^3 + 2 * sin(t2)^3
  t3 <- 2 * (x3 - 0.5)
  f3 <- 3 * sin(3 * t3) + 6 * exp(-36 * t3^2)
  eta <- f1 + f2 + f3 + u[as.integer(cl)]
  scale <- (lambda * exp(eta))^(-1 / shape)
  event_time <- rweibull(n, shape = shape, scale = scale)
  censor_time <- rexp(n, rate = censor_rate)
  dat <- data.frame(
    time = pmin(event_time, censor_time),
    status = as.integer(event_time <= censor_time), x1, x2, x3, x4, cl
  )
  b0 <- gam(
    Surv(time, status) ~ s(x1, bs = "AMatern", k = 10) +
      s(x2, bs = "cr", k = 10) + s(x3, bs = "cr", k = 15) +
      s(x4, bs = "cr", k = 10),
    data = dat, family = cox.ph(), method = "REML"
  )
  pv_null_iid[i] <- taps_score_test(
    b0, test.component = 1, method = "davies"
  )$smooth.pvalue
  b1 <- gam(
    Surv(time, status) ~ s(x1, bs = "AMatern", k = 10) +
      s(x2, bs = "cr", k = 10) + s(x3, bs = "cr", k = 15) +
      s(x4, bs = "cr", k = 10) + s(cl, bs = "re"),
    data = dat, family = cox.ph(), method = "REML"
  )
  pv_null_re[i] <- taps_score_test(
    b1, test.component = 1, method = "davies"
  )$smooth.pvalue
}

# Alternative seeds: 5920001--5920100
for (i in seq_len(nsim)) {
  set.seed(5920000L + i)
  cl <- factor(rep(seq_len(ncl), each = mcl))
  x1 <- qbeta(pnorm(rnorm(ncl)), 1.5, 1.5)[as.integer(cl)]
  Z <- MASS::mvrnorm(n, rep(0, 3), matrix(0.25, 3, 3) + 0.75 * diag(3))
  x2 <- qbeta(pnorm(Z[, 1]), 1.5, 1.5)
  x3 <- qbeta(pnorm(Z[, 2]), 1.5, 1.5)
  x4 <- qbeta(pnorm(Z[, 3]), 1.5, 1.5)
  u <- rnorm(ncl, 0, sigma_u)
  f1 <- smoothed_linearity(x1, 0.5)
  t2 <- 2 * pi * x2
  f2 <- 0.4 * sin(t2) + 0.8 * cos(t2) + 1.2 * sin(t2)^2 +
    1.6 * cos(t2)^3 + 2 * sin(t2)^3
  t3 <- 2 * (x3 - 0.5)
  f3 <- 3 * sin(3 * t3) + 6 * exp(-36 * t3^2)
  eta <- f1 + f2 + f3 + u[as.integer(cl)]
  scale <- (lambda * exp(eta))^(-1 / shape)
  event_time <- rweibull(n, shape = shape, scale = scale)
  censor_time <- rexp(n, rate = censor_rate)
  dat <- data.frame(
    time = pmin(event_time, censor_time),
    status = as.integer(event_time <= censor_time), x1, x2, x3, x4, cl
  )
  b0 <- gam(
    Surv(time, status) ~ s(x1, bs = "AMatern", k = 10) +
      s(x2, bs = "cr", k = 10) + s(x3, bs = "cr", k = 15) +
      s(x4, bs = "cr", k = 10),
    data = dat, family = cox.ph(), method = "REML"
  )
  pv_alt_iid[i] <- taps_score_test(
    b0, test.component = 1, method = "davies"
  )$smooth.pvalue
  b1 <- gam(
    Surv(time, status) ~ s(x1, bs = "AMatern", k = 10) +
      s(x2, bs = "cr", k = 10) + s(x3, bs = "cr", k = 15) +
      s(x4, bs = "cr", k = 10) + s(cl, bs = "re"),
    data = dat, family = cox.ph(), method = "REML"
  )
  pv_alt_re[i] <- taps_score_test(
    b1, test.component = 1, method = "davies"
  )$smooth.pvalue
}

out <- list(null_iid = pv_null_iid, null_re = pv_null_re,
            alternative_iid = pv_alt_iid, alternative_re = pv_alt_re)
if (any(!is.finite(unlist(out, use.names = FALSE)))) {
  stop("Cox RE simulation produced non-finite p-values.")
}
saveRDS(out, file.path("example", "re", "cox_re_pvalues.rds"))

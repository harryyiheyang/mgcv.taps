library(mgcv)
library(survival)
devtools::load_all(".")
options(warn = 1)

n <- 1000L
nsim <- 100L
lambda <- 0.1
shape <- 2
censor_rate <- 0.704
pv_null <- numeric(nsim)
pv_alt <- numeric(nsim)
censor_null <- numeric(nsim)
censor_alt <- numeric(nsim)

# NULL seeds: 910001--910100
for (i in seq_len(nsim)) {
  set.seed(910000L + i)
  X <- MASS::mvrnorm(n, rep(0, 4), matrix(0.25, 4, 4) + 0.75 * diag(4))
  x1 <- qbeta(pnorm(X[, 1]), 1.5, 1.5)
  x2 <- qbeta(pnorm(X[, 2]), 1.5, 1.5)
  x3 <- qbeta(pnorm(X[, 3]), 1.5, 1.5)
  x4 <- qbeta(pnorm(X[, 4]), 1.5, 1.5)
  f1 <- smoothed_linearity(x1, 0)
  t2 <- 2 * pi * x2
  f2 <- 0.4 * sin(t2) + 0.8 * cos(t2) + 1.2 * sin(t2)^2 +
    1.6 * cos(t2)^3 + 2 * sin(t2)^3
  t3 <- 2 * (x3 - 0.5)
  f3 <- 3 * sin(3 * t3) + 6 * exp(-36 * t3^2)
  eta <- f1 + f2 + f3
  scale <- (lambda * exp(eta))^(-1 / shape)
  event_time <- rweibull(n, shape = shape, scale = scale)
  censor_time <- rexp(n, rate = censor_rate)
  censor_null[i] <- mean(event_time > censor_time)
  dat <- data.frame(
    time = pmin(event_time, censor_time),
    status = as.integer(event_time <= censor_time), x1, x2, x3, x4
  )
  b <- gam(
    Surv(time, status) ~ s(x1, bs = "AMatern", k = 10) +
      s(x2, bs = "cr", k = 10) + s(x3, bs = "cr", k = 15) +
      s(x4, bs = "cr", k = 10),
    data = dat, family = cox.ph(), method = "REML"
  )
  pv_null[i] <- taps_score_test(
    b, test.component = 1, method = "davies"
  )$smooth.pvalue
}

# Alternative seeds: 920001--920100
for (i in seq_len(nsim)) {
  set.seed(920000L + i)
  X <- MASS::mvrnorm(n, rep(0, 4), matrix(0.25, 4, 4) + 0.75 * diag(4))
  x1 <- qbeta(pnorm(X[, 1]), 1.5, 1.5)
  x2 <- qbeta(pnorm(X[, 2]), 1.5, 1.5)
  x3 <- qbeta(pnorm(X[, 3]), 1.5, 1.5)
  x4 <- qbeta(pnorm(X[, 4]), 1.5, 1.5)
  f1 <- smoothed_linearity(x1, 0.5)
  t2 <- 2 * pi * x2
  f2 <- 0.4 * sin(t2) + 0.8 * cos(t2) + 1.2 * sin(t2)^2 +
    1.6 * cos(t2)^3 + 2 * sin(t2)^3
  t3 <- 2 * (x3 - 0.5)
  f3 <- 3 * sin(3 * t3) + 6 * exp(-36 * t3^2)
  eta <- f1 + f2 + f3
  scale <- (lambda * exp(eta))^(-1 / shape)
  event_time <- rweibull(n, shape = shape, scale = scale)
  censor_time <- rexp(n, rate = censor_rate)
  censor_alt[i] <- mean(event_time > censor_time)
  dat <- data.frame(
    time = pmin(event_time, censor_time),
    status = as.integer(event_time <= censor_time), x1, x2, x3, x4
  )
  b <- gam(
    Surv(time, status) ~ s(x1, bs = "AMatern", k = 10) +
      s(x2, bs = "cr", k = 10) + s(x3, bs = "cr", k = 15) +
      s(x4, bs = "cr", k = 10),
    data = dat, family = cox.ph(), method = "REML"
  )
  pv_alt[i] <- taps_score_test(
    b, test.component = 1, method = "davies"
  )$smooth.pvalue
}

out <- list(null = pv_null, alternative = pv_alt)
if (!identical(lengths(out), c(null = 100L, alternative = 100L)) ||
    any(!is.finite(unlist(out, use.names = FALSE)))) {
  stop("Cox simulation must produce 100 finite p-values per scenario.")
}
saveRDS(out, file.path("example", "cox_pvalues.rds"))

censoring <- data.frame(
  scenario = c("null", "alternative"),
  n = n,
  simulations = nsim,
  censoring_distribution = "exponential",
  censoring_rate = censor_rate,
  mean_censoring_fraction = c(mean(censor_null), mean(censor_alt)),
  sd_censoring_fraction = c(sd(censor_null), sd(censor_alt))
)
write.csv(censoring, file.path("example", "cox_censoring_summary.csv"),
          row.names = FALSE)

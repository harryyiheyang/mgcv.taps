library(mgcv)
devtools::load_all(".")

n <- 1000L
nsim <- 100L
theta <- c(-2, 0.3)
pv_null <- numeric(nsim)
pv_alt <- numeric(nsim)

rzip <- function(gamma, theta) {
  lambda <- exp(gamma)
  eta <- theta[1] + exp(theta[2]) * gamma
  presence <- 1 - exp(-exp(eta))
  y <- numeric(length(gamma))
  ii <- presence > runif(length(gamma))
  y[ii] <- qpois(runif(sum(ii), dpois(0, lambda[ii]), 1), lambda[ii])
  y
}

# NULL seeds: 1010001--1010100
for (i in seq_len(nsim)) {
  set.seed(1010000L + i)
  X <- MASS::mvrnorm(n, rep(0, 4), matrix(0.25, 4, 4) + 0.75 * diag(4))
  x0 <- qbeta(pnorm(X[, 1]), 1.5, 1.5)
  x1 <- qbeta(pnorm(X[, 2]), 1.5, 1.5)
  x2 <- qbeta(pnorm(X[, 3]), 1.5, 1.5)
  x3 <- qbeta(pnorm(X[, 4]), 1.5, 1.5)
  f1 <- 2 * x0
  t2 <- 2 * pi * x1
  f2 <- 0.4 * sin(t2) + 0.8 * cos(t2) + 1.2 * sin(t2)^2 +
    1.6 * cos(t2)^3 + 2 * sin(t2)^3
  t3 <- 2 * (x2 - 0.5)
  f3 <- 3 * sin(3 * t3) + 6 * exp(-36 * t3^2)
  y <- rzip((f1 + f2 + f3) / 3 - 1, theta)
  dat <- data.frame(x0, x1, x2, x3, y)
  b <- gam(
    y ~ s(x0, bs = "AMatern") + s(x1, bs = "cr") +
      s(x2, bs = "cr") + s(x3, bs = "cr"),
    data = dat, family = ziP(), method = "REML"
  )
  pv_null[i] <- taps_score_test(
    b, test.component = 1, method = "davies"
  )$smooth.pvalue
}

# Alternative seeds: 1020001--1020100
for (i in seq_len(nsim)) {
  set.seed(1020000L + i)
  X <- MASS::mvrnorm(n, rep(0, 4), matrix(0.25, 4, 4) + 0.75 * diag(4))
  x0 <- qbeta(pnorm(X[, 1]), 1.5, 1.5)
  x1 <- qbeta(pnorm(X[, 2]), 1.5, 1.5)
  x2 <- qbeta(pnorm(X[, 3]), 1.5, 1.5)
  x3 <- qbeta(pnorm(X[, 4]), 1.5, 1.5)
  f1 <- smoothed_linearity(x0, 0.5)
  t2 <- 2 * pi * x1
  f2 <- 0.4 * sin(t2) + 0.8 * cos(t2) + 1.2 * sin(t2)^2 +
    1.6 * cos(t2)^3 + 2 * sin(t2)^3
  t3 <- 2 * (x2 - 0.5)
  f3 <- 3 * sin(3 * t3) + 6 * exp(-36 * t3^2)
  y <- rzip((f1 + f2 + f3) / 3 - 1, theta)
  dat <- data.frame(x0, x1, x2, x3, y)
  b <- gam(
    y ~ s(x0, bs = "AMatern") + s(x1, bs = "cr") +
      s(x2, bs = "cr") + s(x3, bs = "cr"),
    data = dat, family = ziP(), method = "REML"
  )
  pv_alt[i] <- taps_score_test(
    b, test.component = 1, method = "davies"
  )$smooth.pvalue
}

out <- list(null = pv_null, alternative = pv_alt)
if (!identical(lengths(out), c(null = 100L, alternative = 100L)) ||
    any(!is.finite(unlist(out, use.names = FALSE)))) {
  stop("ZIP simulation must produce 100 finite p-values per scenario.")
}
saveRDS(out, file.path("example", "zip_pvalues.rds"))

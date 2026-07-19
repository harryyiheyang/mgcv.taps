library(mgcv)
devtools::load_all(".")

n <- 1000L
nsim <- 100L
J <- 5L
cutpoints <- c(-2, -0.5, 0.5, 2)
pv_null <- numeric(nsim)
pv_alt <- numeric(nsim)

# NULL seeds: 810001--810100
for (i in seq_len(nsim)) {
  set.seed(810000L + i)
  X <- MASS::mvrnorm(n, rep(0, 4), matrix(0.25, 4, 4) + 0.75 * diag(4))
  x1 <- qbeta(pnorm(X[, 1]), 1.5, 1.5)
  x2 <- qbeta(pnorm(X[, 2]), 1.5, 1.5)
  x3 <- qbeta(pnorm(X[, 3]), 1.5, 1.5)
  x4 <- qbeta(pnorm(X[, 4]), 1.5, 1.5)
  f1 <- 2 * x1
  t2 <- 2 * pi * x2
  f2 <- 0.4 * sin(t2) + 0.8 * cos(t2) + 1.2 * sin(t2)^2
  t3 <- 2 * (x3 - 0.5)
  f3 <- 3 * sin(3 * t3) + 6 * exp(-36 * t3^2)
  eta <- f1 + f2 + f3
  F <- cbind(0, sapply(cutpoints, function(a) plogis(a - eta)), 1)
  prob <- F[, -1, drop = FALSE] - F[, -ncol(F), drop = FALSE]
  y <- apply(prob, 1, function(p) sample.int(J, 1L, prob = p))
  dat <- data.frame(y = as.integer(y), x1, x2, x3, x4)
  b <- gam(
    y ~ s(x1, bs = "AMatern", k = 10) + s(x2, bs = "cr", k = 10) +
      s(x3, bs = "cr", k = 15) + s(x4, bs = "cr", k = 10),
    data = dat, family = ocat(R = J), method = "REML"
  )
  pv_null[i] <- taps_score_test(
    b, test.component = 1, method = "davies"
  )$smooth.pvalue
}

# Alternative seeds: 820001--820100
for (i in seq_len(nsim)) {
  set.seed(820000L + i)
  X <- MASS::mvrnorm(n, rep(0, 4), matrix(0.25, 4, 4) + 0.75 * diag(4))
  x1 <- qbeta(pnorm(X[, 1]), 1.5, 1.5)
  x2 <- qbeta(pnorm(X[, 2]), 1.5, 1.5)
  x3 <- qbeta(pnorm(X[, 3]), 1.5, 1.5)
  x4 <- qbeta(pnorm(X[, 4]), 1.5, 1.5)
  f1 <- smoothed_linearity(x1, 0.5)
  t2 <- 2 * pi * x2
  f2 <- 0.4 * sin(t2) + 0.8 * cos(t2) + 1.2 * sin(t2)^2
  t3 <- 2 * (x3 - 0.5)
  f3 <- 3 * sin(3 * t3) + 6 * exp(-36 * t3^2)
  eta <- f1 + f2 + f3
  F <- cbind(0, sapply(cutpoints, function(a) plogis(a - eta)), 1)
  prob <- F[, -1, drop = FALSE] - F[, -ncol(F), drop = FALSE]
  y <- apply(prob, 1, function(p) sample.int(J, 1L, prob = p))
  dat <- data.frame(y = as.integer(y), x1, x2, x3, x4)
  b <- gam(
    y ~ s(x1, bs = "AMatern", k = 10) + s(x2, bs = "cr", k = 10) +
      s(x3, bs = "cr", k = 15) + s(x4, bs = "cr", k = 10),
    data = dat, family = ocat(R = J), method = "REML"
  )
  pv_alt[i] <- taps_score_test(
    b, test.component = 1, method = "davies"
  )$smooth.pvalue
}

out <- list(null = pv_null, alternative = pv_alt)
if (!identical(lengths(out), c(null = 100L, alternative = 100L)) ||
    any(!is.finite(unlist(out, use.names = FALSE)))) {
  stop("ocat simulation must produce 100 finite p-values per scenario.")
}
saveRDS(out, file.path("example", "ocat_pvalues.rds"))

library(mgcv)
devtools::load_all(".")

n <- 1000L
ncl <- 50L
mcl <- 20L
sigma_u <- 0.5
nsim <- 100L
n_trial <- 4L
pv_null_iid <- pv_null_re <- numeric(nsim)
pv_alt_iid <- pv_alt_re <- numeric(nsim)

# NULL seeds: 5310001--5310100
for (i in seq_len(nsim)) {
  set.seed(5310000L + i)
  cl <- factor(rep(seq_len(ncl), each = mcl))
  x0 <- qbeta(pnorm(rnorm(ncl)), 1.5, 1.5)[as.integer(cl)]
  Z <- MASS::mvrnorm(n, rep(0, 3), matrix(0.25, 3, 3) + 0.75 * diag(3))
  x1 <- qbeta(pnorm(Z[, 1]), 1.5, 1.5)
  x2 <- qbeta(pnorm(Z[, 2]), 1.5, 1.5)
  x3 <- qbeta(pnorm(Z[, 3]), 1.5, 1.5)
  u <- rnorm(ncl, 0, sigma_u)
  t2 <- 2 * pi * x1
  f2 <- 0.4 * sin(t2) + 0.8 * cos(t2) + 1.2 * sin(t2)^2 +
    1.6 * cos(t2)^3 + 2 * sin(t2)^3
  t3 <- 2 * (x2 - 0.5)
  f3 <- 3 * sin(3 * t3) + 6 * exp(-36 * t3^2)
  prob <- binomial()$linkinv(f2 + f3 - 1 + u[as.integer(cl)])
  y <- rbinom(n, size = n_trial, prob = prob)
  dat <- data.frame(x0, x1, x2, x3, y, n_trial, cl)
  b0 <- gam(
    cbind(y, n_trial - y) ~ s(x0, bs = "AMatern") +
      s(x1, bs = "cr") + s(x2, bs = "cr") + s(x3, bs = "cr"),
    data = dat, family = binomial(), method = "REML"
  )
  pv_null_iid[i] <- taps_score_test(
    b0, test.component = 1, method = "davies"
  )$smooth.pvalue
  b1 <- gam(
    cbind(y, n_trial - y) ~ s(x0, bs = "AMatern") +
      s(x1, bs = "cr") + s(x2, bs = "cr") + s(x3, bs = "cr") +
      s(cl, bs = "re"),
    data = dat, family = binomial(), method = "REML"
  )
  pv_null_re[i] <- taps_score_test(
    b1, test.component = 1, method = "davies"
  )$smooth.pvalue
}

# Alternative seeds: 5320001--5320100
for (i in seq_len(nsim)) {
  set.seed(5320000L + i)
  cl <- factor(rep(seq_len(ncl), each = mcl))
  x0 <- qbeta(pnorm(rnorm(ncl)), 1.5, 1.5)[as.integer(cl)]
  Z <- MASS::mvrnorm(n, rep(0, 3), matrix(0.25, 3, 3) + 0.75 * diag(3))
  x1 <- qbeta(pnorm(Z[, 1]), 1.5, 1.5)
  x2 <- qbeta(pnorm(Z[, 2]), 1.5, 1.5)
  x3 <- qbeta(pnorm(Z[, 3]), 1.5, 1.5)
  u <- rnorm(ncl, 0, sigma_u)
  f1 <- smoothed_linearity(x0, 0.5)
  t2 <- 2 * pi * x1
  f2 <- 0.4 * sin(t2) + 0.8 * cos(t2) + 1.2 * sin(t2)^2 +
    1.6 * cos(t2)^3 + 2 * sin(t2)^3
  t3 <- 2 * (x2 - 0.5)
  f3 <- 3 * sin(3 * t3) + 6 * exp(-36 * t3^2)
  prob <- binomial()$linkinv(f1 + f2 + f3 - 1 + u[as.integer(cl)])
  y <- rbinom(n, size = n_trial, prob = prob)
  dat <- data.frame(x0, x1, x2, x3, y, n_trial, cl)
  b0 <- gam(
    cbind(y, n_trial - y) ~ s(x0, bs = "AMatern") +
      s(x1, bs = "cr") + s(x2, bs = "cr") + s(x3, bs = "cr"),
    data = dat, family = binomial(), method = "REML"
  )
  pv_alt_iid[i] <- taps_score_test(
    b0, test.component = 1, method = "davies"
  )$smooth.pvalue
  b1 <- gam(
    cbind(y, n_trial - y) ~ s(x0, bs = "AMatern") +
      s(x1, bs = "cr") + s(x2, bs = "cr") + s(x3, bs = "cr") +
      s(cl, bs = "re"),
    data = dat, family = binomial(), method = "REML"
  )
  pv_alt_re[i] <- taps_score_test(
    b1, test.component = 1, method = "davies"
  )$smooth.pvalue
}

out <- list(null_iid = pv_null_iid, null_re = pv_null_re,
            alternative_iid = pv_alt_iid, alternative_re = pv_alt_re)
if (any(!is.finite(unlist(out, use.names = FALSE)))) {
  stop("Binomial RE simulation produced non-finite p-values.")
}
saveRDS(out, file.path("example", "re", "binomial_re_pvalues.rds"))

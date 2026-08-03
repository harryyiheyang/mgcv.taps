library(mgcv)
devtools::load_all(".")

n <- 1000L
ncl <- 50L
mcl <- 20L
sigma_u <- 0.5
nsim <- 100L
J <- 5L
cutpoints <- c(-2, -0.5, 0.5, 2)
pv_null_iid <- pv_null_re <- numeric(nsim)
pv_alt_iid <- pv_alt_re <- numeric(nsim)

# NULL seeds: 5810001--5810100
for (i in seq_len(nsim)) {
  set.seed(5810000L + i)
  cl <- factor(rep(seq_len(ncl), each = mcl))
  x1 <- qbeta(pnorm(rnorm(ncl)), 1.5, 1.5)[as.integer(cl)]
  Z <- MASS::mvrnorm(n, rep(0, 3), matrix(0.25, 3, 3) + 0.75 * diag(3))
  x2 <- qbeta(pnorm(Z[, 1]), 1.5, 1.5)
  x3 <- qbeta(pnorm(Z[, 2]), 1.5, 1.5)
  x4 <- qbeta(pnorm(Z[, 3]), 1.5, 1.5)
  u <- rnorm(ncl, 0, sigma_u)
  f1 <- 2 * x1
  t2 <- 2 * pi * x2
  f2 <- 0.4 * sin(t2) + 0.8 * cos(t2) + 1.2 * sin(t2)^2
  t3 <- 2 * (x3 - 0.5)
  f3 <- 3 * sin(3 * t3) + 6 * exp(-36 * t3^2)
  eta <- f1 + f2 + f3 + u[as.integer(cl)]
  F <- cbind(0, sapply(cutpoints, function(a) plogis(a - eta)), 1)
  prob <- F[, -1, drop = FALSE] - F[, -ncol(F), drop = FALSE]
  y <- apply(prob, 1, function(p) sample.int(J, 1L, prob = p))
  dat <- data.frame(y = as.integer(y), x1, x2, x3, x4, cl)
  b0 <- gam(
    y ~ s(x1, bs = "AMatern", k = 10) + s(x2, bs = "cr", k = 10) +
      s(x3, bs = "cr", k = 15) + s(x4, bs = "cr", k = 10),
    data = dat, family = ocat(R = J), method = "REML"
  )
  pv_null_iid[i] <- taps_score_test(
    b0, test.component = 1, method = "davies"
  )$smooth.pvalue
  b1 <- gam(
    y ~ s(x1, bs = "AMatern", k = 10) + s(x2, bs = "cr", k = 10) +
      s(x3, bs = "cr", k = 15) + s(x4, bs = "cr", k = 10) +
      s(cl, bs = "re"),
    data = dat, family = ocat(R = J), method = "REML"
  )
  pv_null_re[i] <- taps_score_test(
    b1, test.component = 1, method = "davies"
  )$smooth.pvalue
}

# Alternative seeds: 5820001--5820100
for (i in seq_len(nsim)) {
  set.seed(5820000L + i)
  cl <- factor(rep(seq_len(ncl), each = mcl))
  x1 <- qbeta(pnorm(rnorm(ncl)), 1.5, 1.5)[as.integer(cl)]
  Z <- MASS::mvrnorm(n, rep(0, 3), matrix(0.25, 3, 3) + 0.75 * diag(3))
  x2 <- qbeta(pnorm(Z[, 1]), 1.5, 1.5)
  x3 <- qbeta(pnorm(Z[, 2]), 1.5, 1.5)
  x4 <- qbeta(pnorm(Z[, 3]), 1.5, 1.5)
  u <- rnorm(ncl, 0, sigma_u)
  f1 <- smoothed_linearity(x1, 0.5)
  t2 <- 2 * pi * x2
  f2 <- 0.4 * sin(t2) + 0.8 * cos(t2) + 1.2 * sin(t2)^2
  t3 <- 2 * (x3 - 0.5)
  f3 <- 3 * sin(3 * t3) + 6 * exp(-36 * t3^2)
  eta <- f1 + f2 + f3 + u[as.integer(cl)]
  F <- cbind(0, sapply(cutpoints, function(a) plogis(a - eta)), 1)
  prob <- F[, -1, drop = FALSE] - F[, -ncol(F), drop = FALSE]
  y <- apply(prob, 1, function(p) sample.int(J, 1L, prob = p))
  dat <- data.frame(y = as.integer(y), x1, x2, x3, x4, cl)
  b0 <- gam(
    y ~ s(x1, bs = "AMatern", k = 10) + s(x2, bs = "cr", k = 10) +
      s(x3, bs = "cr", k = 15) + s(x4, bs = "cr", k = 10),
    data = dat, family = ocat(R = J), method = "REML"
  )
  pv_alt_iid[i] <- taps_score_test(
    b0, test.component = 1, method = "davies"
  )$smooth.pvalue
  b1 <- gam(
    y ~ s(x1, bs = "AMatern", k = 10) + s(x2, bs = "cr", k = 10) +
      s(x3, bs = "cr", k = 15) + s(x4, bs = "cr", k = 10) +
      s(cl, bs = "re"),
    data = dat, family = ocat(R = J), method = "REML"
  )
  pv_alt_re[i] <- taps_score_test(
    b1, test.component = 1, method = "davies"
  )$smooth.pvalue
}

out <- list(null_iid = pv_null_iid, null_re = pv_null_re,
            alternative_iid = pv_alt_iid, alternative_re = pv_alt_re)
if (any(!is.finite(unlist(out, use.names = FALSE)))) {
  stop("ocat RE simulation produced non-finite p-values.")
}
saveRDS(out, file.path("example", "re", "ocat_re_pvalues.rds"))

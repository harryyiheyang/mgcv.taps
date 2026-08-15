library(mgcv)
library(mgcv.taps)

set.seed(20260819)
nid <- 80L
m <- 8L
n <- nid * m
id <- factor(rep(seq_len(nid), each = m))
ii <- as.integer(id)
x <- runif(n, -1, 1)
time <- runif(n)
u0 <- rnorm(nid, sd = 0.55)
eta <- 0.1 + 0.45 * sin(pi * x) + u0[ii]
y <- rpois(n, exp(eta))
dat <- data.frame(y = y, x = x, time = time, id = id)

fit <- gammfast(
  y ~ s(x, k = 6L) + s(id, bs = "re") + fs(time, id, k = 3L),
  data = dat, family = poisson(), nthreads = 2L,
  control = list(objective.tol = 1e-7, max.outer = 200L)
)
if (!fit$converged) stop("The zero 3D-fs Poisson fit did not converge.")

intercept_re <- fit$random$group.index[[1L]]
fs_re <- fit$random$group.index[[2L]]
Gfs <- fit$G[fs_re, fs_re, drop = FALSE]
efs <- eigen(Gfs, symmetric = TRUE, only.values = TRUE)$values
if (any(!is.finite(Gfs)) || min(efs) <= 0) {
  stop("The zero 3D-fs covariance estimate is not positive definite.")
}
if (any(fit$trace$inner_limit)) {
  stop("The zero 3D-fs fit reached the covariance-update cap.")
}
cat(
  "zero 3D-fs Poisson: elapsed =", format(fit$elapsed, digits = 5),
  ", outer =", fit$outer,
  ", intercept variance =", format(fit$G[intercept_re, intercept_re], digits = 5),
  ", fs Frobenius norm =", format(norm(Gfs, "F"), digits = 5),
  ", fs eigenvalues =", paste(format(efs, digits = 5), collapse = ","),
  "\n"
)
print(Gfs)
print(
  fit$trace[, c(
    "outer", "dobjective", "pirls", "pirls_converged",
    "fixedpoint_G", "variance_evaluations", "inner_last", "inner_limit"
  )],
  row.names = FALSE
)

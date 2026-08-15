library(mgcv)
library(mgcv.taps)

set.seed(20260817)
nid <- 80L
m <- 8L
n <- nid * m
id <- factor(rep(seq_len(nid), each = m))
ii <- as.integer(id)
x <- runif(n, -1, 1)
z <- rnorm(n)
time <- runif(n)
u0 <- rnorm(nid, sd = 0.55)
uf1 <- rnorm(nid, sd = 0.30)
uf2 <- rnorm(nid, sd = 0.18)
eta <- 0.1 + 0.45 * sin(pi * x) + u0[ii] +
  sqrt(2) * cos(pi * time) * uf1[ii] +
  sqrt(2) * cos(2 * pi * time) * uf2[ii]
y <- rpois(n, exp(eta))
dat <- data.frame(y = y, x = x, z = z, time = time, id = id)

fit <- gammfast(
  y ~ s(x, k = 6L) + s(id, bs = "re") +
    s(z, id, bs = "re") + fs(time, id, k = 2L),
  data = dat, family = poisson(), nthreads = 2L,
  control = list(objective.tol = 1e-7, max.outer = 200L)
)
if (!fit$converged) stop("The zero scalar-re Poisson fit did not converge.")

zero_re <- fit$random$group.index[[2L]]
intercept_re <- fit$random$group.index[[1L]]
fs_re <- fit$random$group.index[[3L]]
if (fit$G[zero_re, zero_re] >= 0.01) {
  stop("The true-zero scalar re variance did not approach the boundary.")
}
if (any(fit$trace$inner_limit)) {
  stop("The true-zero scalar re fit reached the covariance-update cap.")
}

B <- fit$random$B
Lfs <- t(chol(fit$G[fs_re, fs_re]))
Bfs <- B[, fs_re, drop = FALSE] %*% Lfs
Z0 <- matrix(0, n, nid)
Zs <- matrix(0, n, nid)
Zfs <- matrix(0, n, nid * length(fs_re))
for (g in seq_len(nid)) {
  rows <- which(ii == g)
  Z0[rows, g] <- B[rows, intercept_re]
  Zs[rows, g] <- B[rows, zero_re]
  cols <- (g - 1L) * length(fs_re) + seq_along(fs_re)
  Zfs[rows, cols] <- Bfs[rows, , drop = FALSE]
}
para <- list(
  Z0 = list(diag(nid), rank = nid),
  Zs = list(diag(nid), rank = nid),
  Zfs = list(diag(ncol(Zfs)), sp = 1, rank = ncol(Zfs))
)
t0 <- proc.time()[[3]]
native <- gam(
  y ~ s(x, k = 6L) + Z0 + Zs + Zfs,
  data = dat, family = poisson(), method = "REML",
  sp = c(-1, -1, 1, -1), paraPen = para
)
native_elapsed <- proc.time()[[3]] - t0
native_intercept <- 1 / native$full.sp[1L]
native_zero <- 1 / native$full.sp[2L]
native_global_sp <- native$full.sp[4L]
eta_error <- max(abs(predict(native, type = "link") - fit$linear.predictors))
if (!native$converged || native_zero >= 0.01) {
  stop("The dense mgcv zero scalar-re reference did not converge to the boundary.")
}
if (abs(log(native_intercept / fit$G[intercept_re, intercept_re])) > 0.05 ||
    abs(log(native_global_sp / fit$sp)) > 0.05 || eta_error > 0.02) {
  stop("The dense mgcv scalar-re reference does not match gammfast.")
}
cat(
  "zero scalar-re Poisson: elapsed =", format(fit$elapsed, digits = 5),
  ", outer =", fit$outer,
  ", intercept variance =", format(fit$G[intercept_re, intercept_re], digits = 5),
  ", zero-re variance =", format(fit$G[zero_re, zero_re], digits = 5),
  ", fs eigenvalues =",
  paste(format(eigen(fit$G[fs_re, fs_re], symmetric = TRUE,
                     only.values = TRUE)$values, digits = 5), collapse = ","),
  "\n"
)
cat(
  "dense mgcv conditional reference: elapsed =",
  format(native_elapsed, digits = 5),
  ", intercept variance =", format(native_intercept, digits = 5),
  ", zero-re variance =", format(native_zero, digits = 5),
  ", max eta error =", format(eta_error, digits = 5), "\n"
)
print(
  fit$trace[, c(
    "outer", "dobjective", "pirls", "pirls_converged",
    "fixedpoint_G", "variance_evaluations", "inner_last", "inner_limit"
  )],
  row.names = FALSE
)

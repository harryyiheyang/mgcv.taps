library(mgcv)
library(mgcv.taps)

legacy_formals <- c(
  "formula", "data", "family", "weights", "inner.max", "nthreads",
  "discrete", "control", "verbose"
)
if (!identical(head(names(formals(gammfast)), length(legacy_formals)),
               legacy_formals)) {
  stop("gammfast no longer preserves its existing positional API.")
}

set.seed(20260815)
nid <- 40L
m <- 7L
n <- nid * m
id <- factor(rep(seq_len(nid), each = m))
ii <- as.integer(id)
x <- runif(n, -1, 1)
time <- runif(n)
u0 <- rnorm(nid, sd = 0.45)
u1 <- rnorm(nid, sd = 0.25)
eta <- 0.2 + 0.45 * sin(pi * x) + u0[ii] +
  sqrt(2) * cos(pi * time) * u1[ii]
y <- rpois(n, exp(eta))
dat <- data.frame(y = y, x = x, time = time, id = id)

fit <- gammfast(
  y ~ s(x, k = 5L) + s(id, bs = "re") + fs(time, id, k = 2L),
  data = dat, family = poisson(), nthreads = 2L,
  control = list(objective.tol = 1e-7, max.outer = 200L)
)
if (!fit$converged) stop("The reordered Poisson profile did not converge.")
if (fit$inner.tol != 1e-5 || fit$inner.max != 300L ||
    fit$pirls.tol != 1e-6 || fit$pirls.max != 100L) {
  stop("The non-Gaussian inner defaults are inconsistent.")
}
if (fit$optimization.method != "profiled-PIRLS-G-blockwise-fREML") {
  stop("The reordered non-Gaussian optimization label is missing.")
}
B <- fit$random$B
L <- t(chol(fit$G))
Bw <- B %*% L
Z <- matrix(0, n, nid * ncol(B))
for (g in seq_len(nid)) {
  rows <- which(ii == g)
  cols <- (g - 1L) * ncol(B) + seq_len(ncol(B))
  Z[rows, cols] <- Bw[rows, , drop = FALSE]
}
colnames(Z) <- paste0(
  "id", rep(seq_len(nid), each = ncol(B)),
  "_z", rep(seq_len(ncol(B)), nid)
)
form_native <- update(fit$global.formula, y ~ . + Z)
para <- list(Z = list(diag(ncol(Z)), sp = 1, rank = ncol(Z)))
native <- gam(
  form_native, data = dat, family = poisson(), method = "REML",
  sp = c(1, fit$sp), paraPen = para
)
eta_native <- as.numeric(predict(native, type = "link"))
eta_error <- max(abs(eta_native - fit$linear.predictors))
if (eta_error > 2e-5) {
  stop("The fixed-G, fixed-smoothing mgcv mode does not match gammfast.")
}
native_outer <- gam(
  form_native, data = dat, family = poisson(), method = "REML",
  sp = c(1, -1), paraPen = para
)
lambda_error <- abs(log(native_outer$sp / fit$sp))
if (lambda_error > 0.03) {
  stop("The blockwise smoothing update is not close to nested mgcv at fixed G.")
}

cat("Poisson reordered profile passed; outer =", fit$outer,
    ", elapsed =", format(fit$elapsed, digits = 4),
    ", max eta error =", format(eta_error, digits = 4),
    ", lambda log error =", format(lambda_error, digits = 4), "\n")

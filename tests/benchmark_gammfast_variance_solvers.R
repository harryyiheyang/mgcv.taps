library(mgcv)
library(mgcv.taps)

set.seed(20260821)
nid <- as.integer(Sys.getenv("GAMMFAST_BENCHMARK_NID", "300"))
family_name <- tolower(Sys.getenv("GAMMFAST_BENCHMARK_FAMILY", "binary"))
if (!family_name %in% c("binary", "poisson")) {
  stop("GAMMFAST_BENCHMARK_FAMILY must be binary or poisson.")
}
m <- 5L
id <- factor(rep(seq_len(nid), each = m))
ii <- as.integer(id)
x <- runif(nid * m, -1, 1)
t <- runif(nid * m)
Sigma <- matrix(
  c(0.8, 0.12, -0.05, 0.08,
    0.12, 0.40, 0.08, -0.04,
    -0.05, 0.08, 0.24, 0.05,
    0.08, -0.04, 0.05, 0.16),
  4L, 4L
)
u <- MASS::mvrnorm(nid, rep(0, 4L), Sigma)
B <- cbind(
  1,
  sqrt(2) * cos(pi * t),
  sqrt(2) * cos(2 * pi * t),
  sqrt(2) * cos(3 * pi * t)
)
eta <- -0.2 + 0.7 * sin(pi * x) + rowSums(B * u[ii, ])
if (family_name == "binary") {
  y <- rbinom(length(eta), 1, plogis(eta))
  family <- binomial()
} else {
  y <- rpois(length(eta), exp(pmin(eta, 2)))
  family <- poisson()
}
dat <- data.frame(
  y = y, x = x, t = t, id = id,
  fs1 = B[, 2L], fs2 = B[, 3L], fs3 = B[, 4L]
)

solvers <- "fixedpoint"
max_outer <- as.integer(Sys.getenv("GAMMFAST_BENCHMARK_MAX_OUTER", "5000"))
inner_max <- as.integer(Sys.getenv("GAMMFAST_BENCHMARK_INNER", "30"))
fit_tolerance <- as.numeric(Sys.getenv("GAMMFAST_BENCHMARK_TOL", "1e-7"))
out <- vector("list", length(solvers))
fits <- vector("list", length(solvers))
for (j in seq_along(solvers)) {
  tm <- system.time(
    fit <- gammfast(
      y ~ s(x, k = 10L) + s(id, bs = "re") + fs(t, id, k = 3L),
      data = dat, family = family, inner.max = inner_max,
      nthreads = 1L,
      control = list(
        objective.tol = fit_tolerance, fixedpoint.tol = fit_tolerance,
        max.outer = max_outer
      )
    )
  )
  fits[[j]] <- fit
  out[[j]] <- data.frame(
    solver = solvers[j], converged = fit$converged,
    outer = fit$outer,
    variance_evaluations = fit$variance.evaluations,
    local_updates = fit$local.updates,
    elapsed = unname(tm["elapsed"]),
    fixedpoint_G = tail(fit$trace$fixedpoint_G, 1L),
    eta_change = tail(fit$trace$eta_change, 1L),
    min_eigen_G = min(eigen(fit$G, symmetric = TRUE)$values),
    G11 = fit$G[1L, 1L]
  )
  print(out[[j]])
}
out <- do.call(rbind, out)
reference <- fits[[1L]]
out$G_difference <- vapply(
  fits,
  function(fit) norm(fit$G - reference$G, "F") /
    (1 + norm(reference$G, "F")),
  numeric(1L)
)
out$eta_difference <- vapply(
  fits,
  function(fit) max(abs(fit$linear.predictors - reference$linear.predictors)),
  numeric(1L)
)
print(out)
print(reference$G)

if (identical(
    tolower(Sys.getenv("GAMMFAST_BENCHMARK_GAMM4", "false")), "true"
)) {
  library(gamm4)
  tm4 <- system.time(
    fit4 <- gamm4(
      y ~ s(x, k = 10L),
      random = ~(1 | id) + (0 + fs1 + fs2 + fs3 | id),
      data = dat, family = family,
      control = lme4::glmerControl(
        optimizer = "bobyqa", tolPwrss = 1e-10,
        optCtrl = list(maxfun = 200000L)
      )
    )
  )
  vc4 <- as.data.frame(lme4::VarCorr(fit4$mer))
  print(vc4)
  fs4 <- vc4$vcov[
    !is.na(vc4$var1) & vc4$var1 %in% c("fs1", "fs2", "fs3") &
      is.na(vc4$var2)
  ]
  fs_fast <- diag(fits[[1L]]$G)[2:4]
  comparison <- data.frame(
    component = c("fs1", "fs2", "fs3"),
    truth = diag(Sigma)[2:4], gammfast = fs_fast, lme4 = fs4,
    gammfast_shrinkage_from_lme4 = 1 - fs_fast / fs4,
    gammfast_shrinkage_from_truth = 1 - fs_fast / diag(Sigma)[2:4],
    lme4_shrinkage_from_truth = 1 - fs4 / diag(Sigma)[2:4]
  )
  print(data.frame(method = "gamm4/lme4", elapsed = unname(tm4["elapsed"])))
  print(comparison)
}

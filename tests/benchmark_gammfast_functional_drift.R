library(mgcv)
library(mgcv.taps)

n_rep <- as.integer(Sys.getenv("GAMMFAST_FUNCTIONAL_REPS", "5"))
n_id <- as.integer(Sys.getenv("GAMMFAST_FUNCTIONAL_NID", "70"))
min_per_id <- as.integer(Sys.getenv("GAMMFAST_FUNCTIONAL_MIN_PER_ID", "5"))
max_per_id <- as.integer(Sys.getenv("GAMMFAST_FUNCTIONAL_MAX_PER_ID", "10"))
pirls_max <- as.integer(Sys.getenv("GAMMFAST_FUNCTIONAL_PIRLS_MAX", "1000"))
if (n_rep < 1L || n_id < 20L || min_per_id < 4L ||
    max_per_id < min_per_id || pirls_max < 1L) {
  stop("The functional-drift benchmark dimensions are invalid.")
}

available_families <- c("gaussian", "binary", "nb", "sbeta")
family_names <- trimws(strsplit(Sys.getenv(
  "GAMMFAST_FUNCTIONAL_FAMILIES", paste(available_families, collapse = ",")
), ",", fixed = TRUE)[[1L]])
if (!length(family_names) || any(!family_names %in% available_families)) {
  stop("GAMMFAST_FUNCTIONAL_FAMILIES contains an unsupported family.")
}

available_drift_scales <- c(low = 0.02, medium = 0.12, high = 0.40)
drift_names <- trimws(strsplit(Sys.getenv(
  "GAMMFAST_FUNCTIONAL_DRIFTS",
  paste(names(available_drift_scales), collapse = ",")
), ",", fixed = TRUE)[[1L]])
if (!length(drift_names) ||
    any(!drift_names %in% names(available_drift_scales))) {
  stop("GAMMFAST_FUNCTIONAL_DRIFTS contains an unsupported level.")
}
drift_scales <- available_drift_scales[drift_names]

intercept_variance <- 0.16
D_scale <- diag(c(1, sqrt(0.5), sqrt(0.1)))
R_ar1 <- outer(1:3, 1:3, function(i, j) 0.5^abs(i - j))
rownames(D_scale) <- colnames(D_scale) <- paste0("cos", 1:3)
rownames(R_ar1) <- colnames(R_ar1) <- paste0("cos", 1:3)

truth_covariance <- function(drift_scale) {
  drift_scale * D_scale %*% R_ar1 %*% D_scale
}

family_factory <- function(name) {
  switch(name,
    gaussian = stats::gaussian(), binary = stats::binomial(),
    nb = mgcv::nb(theta = 3, link = "log"),
    sbeta = mgcv::betar(theta = 25, link = "logit"),
    stop("Unknown benchmark family.")
  )
}

gamm4_family_factory <- function(name) {
  switch(name,
    gaussian = stats::gaussian(), binary = stats::binomial(),
    nb = MASS::negative.binomial(theta = 3, link = "log"),
    stop("gamm4 is not used for this benchmark family.")
  )
}

example_mean_function <- function(time) {
  angle <- 2 * pi * time
  0.4 * sin(angle) + 0.8 * cos(angle) +
    1.2 * sin(angle)^2 + 1.6 * cos(angle)^3 +
    2 * sin(angle)^3
}

mean_grid <- seq(0, 1, length.out = 10001L)
mean_center <- mean(example_mean_function(mean_grid))
family_mean_scale <- c(gaussian = 1, binary = 0.45, nb = 0.35, sbeta = 0.45)
family_intercept <- c(gaussian = 0.2, binary = -0.2, nb = -0.5, sbeta = -0.2)

cosine_basis <- function(time) {
  sqrt(2) * cbind(
    cos1 = cos(pi * time), cos2 = cos(2 * pi * time),
    cos3 = cos(3 * pi * time)
  )
}

simulate_data <- function(family_name, drift_scale, seed) {
  set.seed(seed)
  n_per_id <- sample(
    min_per_id:max_per_id, n_id, replace = TRUE,
    prob = seq(max_per_id - min_per_id + 1L, 1L)
  )
  id <- factor(rep(seq_len(n_id), n_per_id))
  n <- length(id)
  mixture <- stats::runif(n) < 0.7
  raw_time <- numeric(n)
  raw_time[mixture] <- stats::rbeta(sum(mixture), 1.4, 2.6)
  raw_time[!mixture] <- stats::rbeta(sum(!mixture), 3.8, 1.5)
  time <- (raw_time - min(raw_time)) / diff(range(raw_time))
  if (length(unique(time)) != n) {
    stop("The sparse functional-time design unexpectedly contains ties.")
  }

  basis <- cosine_basis(time)
  subject_index <- as.integer(id)
  random_intercept <- stats::rnorm(n_id, sd = sqrt(intercept_variance))
  G_fs <- truth_covariance(drift_scale)
  random_cosine <- matrix(stats::rnorm(3L * n_id), n_id, 3L) %*%
    chol(G_fs)
  eta_mean <- family_intercept[family_name] +
    family_mean_scale[family_name] *
      (example_mean_function(time) - mean_center)
  eta_random <- random_intercept[subject_index] +
    rowSums(basis * random_cosine[subject_index, , drop = FALSE])
  eta <- eta_mean + eta_random

  y <- switch(family_name,
    gaussian = eta + stats::rnorm(n, sd = 0.55),
    binary = stats::rbinom(n, 1L, stats::plogis(eta)),
    nb = stats::rnbinom(n, mu = exp(eta), size = 3),
    sbeta = {
      mu <- pmin(pmax(stats::plogis(eta), 1e-8), 1 - 1e-8)
      stats::rbeta(n, mu * 25, (1 - mu) * 25)
    }
  )
  data.frame(
    y = y, time = time, id = id,
    cos1 = basis[, 1L], cos2 = basis[, 2L], cos3 = basis[, 3L],
    eta_mean = unname(eta_mean)
  )
}

gammfast_formula <- y ~ s(time, k = 12L, bs = "cr") +
  s(id, bs = "re") + fs(time, id, k = 3L)
gamm4_formula <- y ~ s(time, k = 12L, bs = "cr")
gamm4_random <- ~(1 | id) + (0 + cos1 + cos2 + cos3 | id)

fit_gammfast <- function(dat, family_name, influence_update) {
  fit <- gammfast(
    gammfast_formula, data = dat, family = family_factory(family_name),
    inner.max = 5L, inner.tol = 1e-5, pirls.max = pirls_max,
    nthreads = 1L, influence.update = influence_update,
    control = list(max.outer = 1000L, objective.tol = 1e-5)
  )
  group_types <- vapply(
    fit$random$groups, function(group) group$type, character(1L)
  )
  intercept_index <- fit$random$group.index[[match(
    "random_intercept", group_types
  )]]
  fs_index <- fit$random$group.index[[match("fs_cosine", group_types)]]
  list(
    G = fit$G[fs_index, fs_index, drop = FALSE],
    intercept = fit$G[intercept_index, intercept_index],
    mean_rmse = sqrt(mean((fit$global.fitted - dat$eta_mean)^2)),
    converged = isTRUE(fit$converged), singular = FALSE, outer = fit$outer
  )
}

fit_gamm4 <- function(dat, family_name) {
  if (!requireNamespace("gamm4", quietly = TRUE) ||
      !requireNamespace("lme4", quietly = TRUE)) {
    stop("The corrected covariance reference requires gamm4 and lme4.")
  }
  fit <- gamm4::gamm4(
    gamm4_formula, random = gamm4_random, data = dat,
    family = gamm4_family_factory(family_name)
  )
  vc <- lme4::VarCorr(fit$mer)
  fs_index <- which(vapply(vc, function(x) {
    identical(colnames(x), paste0("cos", 1:3))
  }, logical(1L)))
  intercept_index <- which(vapply(vc, function(x) {
    identical(colnames(x), "(Intercept)")
  }, logical(1L)))
  if (length(fs_index) != 1L || length(intercept_index) != 1L) {
    stop("gamm4 did not return the expected separated covariance blocks.")
  }
  messages <- fit$mer@optinfo$conv$lme4$messages
  list(
    G = unclass(vc[[fs_index]]),
    intercept = unname(vc[[intercept_index]][1L, 1L]),
    mean_rmse = sqrt(mean((as.numeric(stats::predict(
      fit$gam, newdata = dat, type = "link"
    )) - dat$eta_mean)^2)),
    converged = is.null(messages), singular = lme4::isSingular(fit$mer),
    outer = NA_integer_
  )
}

inla_hyper_mode <- function(fit, pattern) {
  index <- grep(pattern, rownames(fit$summary.hyperpar))
  if (length(index) != 1L) {
    stop("INLA did not return one hyperparameter matching: ", pattern)
  }
  value <- fit$summary.hyperpar[index, "mode"]
  if (!is.finite(value)) value <- fit$summary.hyperpar[index, "0.5quant"]
  unname(value)
}

fit_inla <- function(dat, family_name) {
  if (!requireNamespace("INLA", quietly = TRUE)) {
    stop("The beta covariance reference requires INLA.")
  }
  subject_index <- as.integer(dat$id)
  n_subject <- nlevels(dat$id)
  smooth <- mgcv::smoothCon(
    mgcv::s(time, k = 12L, bs = "cr"), data = dat,
    absorb.cons = TRUE
  )[[1L]]
  penalty_eigen <- eigen(smooth$S[[1L]], symmetric = TRUE)
  positive <- penalty_eigen$values >
    max(penalty_eigen$values) * 1e-10
  penalized_design <- smooth$X %*%
    penalty_eigen$vectors[, positive, drop = FALSE] %*%
    diag(1 / sqrt(penalty_eigen$values[positive]))
  null_design <- smooth$X %*%
    penalty_eigen$vectors[, !positive, drop = FALSE]
  colnames(null_design) <- paste0("mean_null", seq_len(ncol(null_design)))
  inla_data <- transform(
    dat, intercept = 1, random_intercept = subject_index,
    random_fs1 = subject_index,
    random_fs2 = subject_index + n_subject,
    random_fs3 = subject_index + 2L * n_subject
  )
  inla_data <- cbind(inla_data, null_design)
  stack <- INLA::inla.stack(
    data = list(y = dat$y), A = list(1, penalized_design),
    effects = list(
      inla_data[setdiff(names(inla_data), "y")],
      mean_penalized = seq_len(ncol(penalized_design))
    )
  )
  likelihood <- switch(family_name,
    gaussian = "gaussian", binary = "binomial",
    nb = "nbinomial", sbeta = "beta"
  )
  control_family <- switch(family_name,
    nb = list(hyper = list(theta = list(
      initial = log(3), fixed = TRUE
    ))),
    sbeta = list(hyper = list(theta = list(
      initial = log(25), fixed = TRUE
    ))),
    list()
  )
  fixed_terms <- paste(c("intercept", colnames(null_design)), collapse = " + ")
  formula <- stats::as.formula(paste0(
    "y ~ 0 + ", fixed_terms,
    " + f(mean_penalized, model = 'iid')",
    " + f(random_intercept, model = 'iid')",
    " + f(random_fs1, cos1, model = 'iid3d', n = ",
    3L * n_subject,
    ", hyper = list(theta1 = list(prior = 'wishart3d', ",
    "param = c(7, 0.6, 0.3, 0.06, 0, 0, 0))))",
    " + f(random_fs2, cos2, copy = 'random_fs1')",
    " + f(random_fs3, cos3, copy = 'random_fs1')"
  ))
  environment(formula) <- environment()
  fit <- INLA::inla(
    formula, data = INLA::inla.stack.data(stack), family = likelihood,
    control.family = control_family,
    control.predictor = list(A = INLA::inla.stack.A(stack)),
    control.compute = list(config = TRUE), verbose = FALSE
  )
  precision <- vapply(1:3, function(j) inla_hyper_mode(
    fit, paste0("^Precision for random_fs1 \\(component ", j, "\\)$")
  ), numeric(1L))
  correlation <- diag(3L)
  correlation[1L, 2L] <- correlation[2L, 1L] <-
    inla_hyper_mode(fit, "^Rho1:2 for random_fs1$")
  correlation[1L, 3L] <- correlation[3L, 1L] <-
    inla_hyper_mode(fit, "^Rho1:3 for random_fs1$")
  correlation[2L, 3L] <- correlation[3L, 2L] <-
    inla_hyper_mode(fit, "^Rho2:3 for random_fs1$")
  standard_deviation <- 1 / sqrt(precision)
  G <- tcrossprod(standard_deviation) * correlation
  penalized_mean <- fit$summary.random$mean_penalized$mean[
    match(seq_len(ncol(penalized_design)),
      fit$summary.random$mean_penalized$ID)
  ]
  fixed_names <- c("intercept", colnames(null_design))
  fixed_mean <- fit$summary.fixed[fixed_names, "mean"]
  eta_mean <- drop(cbind(1, null_design) %*% fixed_mean) +
    drop(penalized_design %*% penalized_mean)
  list(
    G = G,
    intercept = 1 / inla_hyper_mode(
      fit, "^Precision for random_intercept$"
    ),
    mean_rmse = sqrt(mean((eta_mean - dat$eta_mean)^2)),
    converged = all(is.finite(G)), singular = FALSE, outer = NA_integer_
  )
}

reference_method <- function(family_name) {
  if (identical(family_name, "sbeta")) "inla" else "gamm4"
}

run_timed <- function(fun) {
  gc(verbose = FALSE)
  value <- NULL
  error <- NULL
  elapsed <- system.time({
    value <- tryCatch(fun(), error = function(e) {
      error <<- conditionMessage(e)
      NULL
    })
  })[["elapsed"]]
  list(value = value, error = error, elapsed = unname(elapsed))
}

pack_covariance <- function(G) {
  c(
    g11 = G[1L, 1L], g22 = G[2L, 2L], g33 = G[3L, 3L],
    g12 = G[1L, 2L], g13 = G[1L, 3L], g23 = G[2L, 3L]
  )
}

unpack_covariance <- function(x) {
  matrix(c(
    x["g11"], x["g12"], x["g13"],
    x["g12"], x["g22"], x["g23"],
    x["g13"], x["g23"], x["g33"]
  ), 3L, 3L, byrow = TRUE)
}

records <- list()
record_index <- 0L
method_orders <- list(
  c("reference", "frozen", "refreshed"),
  c("refreshed", "reference", "frozen"),
  c("frozen", "refreshed", "reference")
)
for (family_name in family_names) {
  family_seed_index <- match(family_name, available_families)
  for (drift_name in names(drift_scales)) {
    drift_seed_index <- match(drift_name, names(available_drift_scales))
    drift_scale <- unname(drift_scales[drift_name])
    G_truth <- truth_covariance(drift_scale)
    for (replicate_index in seq_len(n_rep)) {
      seed <- 20261000L + family_seed_index * 1000L +
        drift_seed_index * 100L + replicate_index
      dat <- simulate_data(family_name, drift_scale, seed)
      order_index <- 1L + (replicate_index + drift_seed_index - 2L) %% 3L
      for (method_key in method_orders[[order_index]]) {
        method <- if (identical(method_key, "reference")) {
          reference_method(family_name)
        } else {
          method_key
        }
        result <- switch(method,
          gamm4 = run_timed(function() fit_gamm4(dat, family_name)),
          inla = run_timed(function() fit_inla(dat, family_name)),
          frozen = run_timed(function() fit_gammfast(
            dat, family_name, "frozen"
          )),
          refreshed = run_timed(function() fit_gammfast(
            dat, family_name, "refreshed"
          ))
        )
        value <- result$value
        covariance <- if (is.null(value)) {
          setNames(rep(NA_real_, 6L), names(pack_covariance(G_truth)))
        } else {
          pack_covariance(value$G)
        }
        record_index <- record_index + 1L
        records[[record_index]] <- data.frame(
          family = family_name, drift = drift_name,
          drift_scale = drift_scale, replicate = replicate_index,
          method = method, reference = reference_method(family_name),
          intercept = if (is.null(value)) NA_real_ else value$intercept,
          t(covariance), elapsed = result$elapsed,
          mean_rmse = if (is.null(value)) NA_real_ else value$mean_rmse,
          converged = !is.null(value) && isTRUE(value$converged),
          singular = if (is.null(value)) NA else value$singular,
          outer = if (is.null(value)) NA_integer_ else value$outer,
          error = if (is.null(result$error)) NA_character_ else result$error,
          check.names = FALSE
        )
      }
    }
  }
}

results <- do.call(rbind, records)
results$relative_G_error_truth <- NA_real_
results$diagonal_log_rmse_truth <- NA_real_
results$correlation_rmse_truth <- NA_real_
results$relative_G_error_reference <- NA_real_
results$relative_difference_vs_frozen <- NA_real_

covariance_names <- names(pack_covariance(diag(3L)))
key <- paste(results$family, results$drift, results$replicate)
for (i in seq_len(nrow(results))) {
  G_estimate <- unpack_covariance(unlist(results[i, covariance_names]))
  G_truth <- truth_covariance(results$drift_scale[i])
  if (any(!is.finite(G_estimate))) next
  results$relative_G_error_truth[i] <- norm(G_estimate - G_truth, "F") /
    norm(G_truth, "F")
  results$diagonal_log_rmse_truth[i] <- sqrt(mean(log(
    pmax(diag(G_estimate), 1e-10) / diag(G_truth)
  )^2))
  if (all(diag(G_estimate) > 0)) {
    results$correlation_rmse_truth[i] <- sqrt(mean(
      (cov2cor(G_estimate)[upper.tri(G_estimate)] -
        R_ar1[upper.tri(R_ar1)])^2
    ))
  }

  reference_index <- key == key[i] & results$method == results$reference[i]
  frozen_index <- key == key[i] & results$method == "frozen"
  if (sum(reference_index) == 1L) {
    G_reference <- unpack_covariance(unlist(
      results[reference_index, covariance_names]
    ))
    results$relative_G_error_reference[i] <-
      norm(G_estimate - G_reference, "F") /
        pmax(norm(G_reference, "F"), 1e-10)
  }
  if (identical(results$method[i], "refreshed") &&
      sum(frozen_index) == 1L) {
    G_frozen <- unpack_covariance(unlist(
      results[frozen_index, covariance_names]
    ))
    results$relative_difference_vs_frozen[i] <-
      norm(G_estimate - G_frozen, "F") /
        pmax(norm(G_frozen, "F"), 1e-10)
  } else {
    results$relative_difference_vs_frozen[i] <- 0
  }
}

median_finite <- function(x) {
  x <- x[is.finite(x)]
  if (length(x)) stats::median(x) else NA_real_
}
summary_accuracy <- aggregate(
  cbind(
    intercept, relative_G_error_truth, diagonal_log_rmse_truth,
    correlation_rmse_truth, relative_G_error_reference, mean_rmse,
    relative_difference_vs_frozen
  ) ~ family + drift + drift_scale + method,
  data = results, FUN = median_finite, na.action = stats::na.pass
)
summary_speed <- aggregate(
  elapsed ~ family + drift + method, data = results,
  FUN = function(x) c(median = median(x), mean = mean(x))
)
summary_convergence <- aggregate(
  converged ~ family + drift + method, data = results,
  FUN = function(x) paste(sum(x), "/", length(x))
)
summary_singularity <- aggregate(
  singular ~ family + drift + method, data = results,
  FUN = function(x) paste(sum(x, na.rm = TRUE), "/", sum(!is.na(x)))
)

cat("\nUnstructured functional-covariance accuracy:\n")
print(summary_accuracy, row.names = FALSE)
cat("\nElapsed seconds by family, drift, and method:\n")
print(summary_speed, row.names = FALSE)
cat("\nConvergence summary:\n")
print(summary_convergence, row.names = FALSE)
cat("\nSingular-fit summary:\n")
print(summary_singularity, row.names = FALSE)
if (any(!results$converged)) {
  cat("\nFit errors/non-convergence:\n")
  print(results[!results$converged, c(
    "family", "drift", "replicate", "method", "outer", "error"
  )], row.names = FALSE)
}

output_file <- Sys.getenv("GAMMFAST_FUNCTIONAL_BENCHMARK_OUTPUT", "")
if (nzchar(output_file)) {
  utils::write.csv(results, output_file, row.names = FALSE)
}

invisible(list(
  results = results, accuracy = summary_accuracy,
  speed = summary_speed, convergence = summary_convergence,
  singularity = summary_singularity
))

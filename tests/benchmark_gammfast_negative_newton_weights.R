library(mgcv)
library(mgcv.taps)

n_rep <- as.integer(Sys.getenv("GAMMFAST_NEGATIVE_W_REPS", "5"))
n_id <- as.integer(Sys.getenv("GAMMFAST_NEGATIVE_W_NID", "60"))
min_per_id <- as.integer(Sys.getenv("GAMMFAST_NEGATIVE_W_MIN_PER_ID", "5"))
max_per_id <- as.integer(Sys.getenv("GAMMFAST_NEGATIVE_W_MAX_PER_ID", "10"))
phi <- as.numeric(Sys.getenv("GAMMFAST_NEGATIVE_W_PHI", "1"))
if (n_rep < 1L || n_id < 20L || min_per_id < 4L ||
    max_per_id < min_per_id || !is.finite(phi) || phi <= 0) {
  stop("The negative-Newton-weight benchmark dimensions are invalid.")
}

available_drift_scales <- c(low = 0.02, medium = 0.12, high = 0.40)
drift_names <- trimws(strsplit(Sys.getenv(
  "GAMMFAST_NEGATIVE_W_DRIFTS",
  paste(names(available_drift_scales), collapse = ",")
), ",", fixed = TRUE)[[1L]])
if (!length(drift_names) ||
    any(!drift_names %in% names(available_drift_scales))) {
  stop("GAMMFAST_NEGATIVE_W_DRIFTS contains an unsupported level.")
}
drift_scales <- available_drift_scales[drift_names]

available_methods <- c("gammfast", "gamm4", "mgcv_observed_newton")
method_names <- trimws(strsplit(Sys.getenv(
  "GAMMFAST_NEGATIVE_W_METHODS", paste(available_methods, collapse = ",")
), ",", fixed = TRUE)[[1L]])
if (!length(method_names) || any(!method_names %in% available_methods)) {
  stop("GAMMFAST_NEGATIVE_W_METHODS contains an unsupported method.")
}
strict <- identical(Sys.getenv("GAMMFAST_NEGATIVE_W_STRICT", "false"), "true")

intercept_variance <- 0.16
D_scale <- diag(c(1, sqrt(0.5), sqrt(0.1)))
R_ar1 <- outer(1:3, 1:3, function(i, j) 0.5^abs(i - j))
truth_covariance <- function(drift_scale) {
  drift_scale * D_scale %*% R_ar1 %*% D_scale
}

example_mean_function <- function(time) {
  angle <- 2 * pi * time
  0.4 * sin(angle) + 0.8 * cos(angle) +
    1.2 * sin(angle)^2 + 1.6 * cos(angle)^3 +
    2 * sin(angle)^3
}
mean_grid <- seq(0, 1, length.out = 10001L)
mean_center <- mean(example_mean_function(mean_grid))

cosine_basis <- function(time) {
  sqrt(2) * cbind(
    cos1 = cos(pi * time), cos2 = cos(2 * pi * time),
    cos3 = cos(3 * pi * time)
  )
}

# For inverse-Gaussian responses with log link,
# V(mu) = mu^3 and g(mu) = log(mu). Hence the exact observed-Newton
# multiplier is alpha = 2 * y / mu - 1. The corresponding observed
# weight is negative exactly when y < mu / 2, while the Fisher weight
# remains positive.
inverse_gaussian_weights <- function(y, mu, prior_weights = 1) {
  list(
    observed = prior_weights * (2 * y - mu) / mu^2,
    fisher = prior_weights / mu,
    alpha = 2 * y / mu - 1
  )
}

exact_mu <- rep(1, 4L)
exact_y <- c(0.2, 0.4, 0.6, 2)
exact_weights <- inverse_gaussian_weights(exact_y, exact_mu)
if (!identical(exact_weights$observed < 0, exact_y < exact_mu / 2) ||
    any(exact_weights$fisher <= 0)) {
  stop("The exact inverse-Gaussian negative-weight construction failed.")
}
cat("Exact inverse-Gaussian/log construction:\n")
print(data.frame(
  y = exact_y, mu = exact_mu, alpha = exact_weights$alpha,
  observed_W = exact_weights$observed,
  fisher_W = exact_weights$fisher
), row.names = FALSE)

rinvgauss <- function(mu, phi) {
  # Michael-Schucany-Haas generator in the GLM parameterization
  # Var(Y) = phi * mu^3, so the conventional shape is lambda = 1 / phi.
  lambda <- 1 / phi
  v <- stats::rnorm(length(mu))^2
  candidate <- mu + mu^2 * v / (2 * lambda) -
    mu * sqrt(4 * mu * lambda * v + mu^2 * v^2) / (2 * lambda)
  choose_candidate <- stats::runif(length(mu)) <=
    mu / (mu + candidate)
  ifelse(choose_candidate, candidate, mu^2 / candidate)
}

simulate_data <- function(drift_scale, seed) {
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

  basis <- cosine_basis(time)
  subject_index <- as.integer(id)
  random_intercept <- stats::rnorm(
    n_id, sd = sqrt(intercept_variance)
  )
  G_fs <- truth_covariance(drift_scale)
  random_cosine <- matrix(stats::rnorm(3L * n_id), n_id, 3L) %*%
    chol(G_fs)
  eta_mean <- -0.15 + 0.32 *
    (example_mean_function(time) - mean_center)
  eta_random <- random_intercept[subject_index] +
    rowSums(basis * random_cosine[subject_index, , drop = FALSE])
  eta <- eta_mean + eta_random
  mu <- exp(eta)
  y <- rinvgauss(mu, phi)

  data.frame(
    y = y, time = time, id = id,
    cos1 = basis[, 1L], cos2 = basis[, 2L], cos3 = basis[, 3L],
    eta_mean = eta_mean, eta = eta, mu = mu
  )
}

gammfast_formula <- y ~ s(time, k = 12L, bs = "cr") +
  s(id, bs = "re") + fs(time, id, k = 3L)
gamm4_formula <- y ~ s(time, k = 12L, bs = "cr")
gamm4_random <- ~(1 | id) + (0 + cos1 + cos2 + cos3 | id)

relative_G_error <- function(estimate, truth) {
  sqrt(sum((estimate - truth)^2)) / sqrt(sum(truth^2))
}

diagnose_weights <- function(y, eta, prior_weights = 1) {
  mu <- exp(eta)
  weight <- inverse_gaussian_weights(y, mu, prior_weights)
  c(
    negative_fraction = mean(weight$observed < 0),
    minimum_observed = min(weight$observed),
    minimum_fisher = min(weight$fisher),
    alpha_q01 = unname(stats::quantile(weight$alpha, 0.01)),
    alpha_median = stats::median(weight$alpha)
  )
}

run_timed <- function(code) {
  start <- proc.time()[["elapsed"]]
  value <- code()
  list(value = value, elapsed = proc.time()[["elapsed"]] - start)
}

fit_gammfast <- function(dat, truth_G) {
  fit <- gammfast(
    gammfast_formula, data = dat,
    family = stats::inverse.gaussian(link = "log"),
    inner.max = 5L, inner.tol = 1e-5, pirls.max = 1000L,
    nthreads = 1L,
    control = list(max.outer = 1000L, objective.tol = 1e-5)
  )
  group_types <- vapply(
    fit$random$groups, function(group) group$type, character(1L)
  )
  intercept_index <- fit$random$group.index[[match(
    "random_intercept", group_types
  )]]
  fs_index <- fit$random$group.index[[match("fs_cosine", group_types)]]
  G_estimate <- fit$G[fs_index, fs_index, drop = FALSE]
  final_diagnostic <- diagnose_weights(
    fit$y, fit$linear.predictors, fit$prior.weights
  )
  fisher_work <- mgcv.taps:::gammfast_working(
    fit$family, fit$y, fit$linear.predictors, fit$prior.weights
  )
  determinant <- mgcv.taps:::gammfast_t_correction(
    fit$family, fit$y, fit$linear.predictors,
    fit$prior.weights, fisher_work$w
  )
  if (!isTRUE(all.equal(
    determinant$weight, fisher_work$w, tolerance = 1e-14
  ))) {
    stop("gammfast did not retain Fisher determinant curvature.")
  }

  c(
    converged = isTRUE(fit$converged),
    G_relative_error = relative_G_error(G_estimate, truth_G),
    intercept_variance = fit$G[intercept_index, intercept_index],
    mean_rmse = sqrt(mean((fit$global.fitted - dat$eta_mean)^2)),
    eta_rmse = sqrt(mean((fit$linear.predictors - dat$eta)^2)),
    outer = fit$outer,
    final_diagnostic
  )
}

fit_gamm4 <- function(dat, truth_G) {
  if (!requireNamespace("gamm4", quietly = TRUE) ||
      !requireNamespace("lme4", quietly = TRUE)) {
    stop("gamm4 and lme4 are required for the full-covariance reference.")
  }
  fit <- gamm4::gamm4(
    gamm4_formula, random = gamm4_random, data = dat,
    family = stats::inverse.gaussian(link = "log")
  )
  vc <- lme4::VarCorr(fit$mer)
  fs_index <- which(vapply(vc, function(x) {
    identical(colnames(x), paste0("cos", 1:3))
  }, logical(1L)))
  intercept_index <- which(vapply(vc, function(x) {
    identical(colnames(x), "(Intercept)")
  }, logical(1L)))
  if (length(fs_index) != 1L || length(intercept_index) != 1L) {
    stop("gamm4 did not return the expected covariance blocks.")
  }
  fitted_eta <- as.numeric(stats::predict(fit$mer, type = "link"))
  mean_eta <- as.numeric(stats::predict(fit$gam, type = "link"))
  messages <- fit$mer@optinfo$conv$lme4$messages
  c(
    converged = is.null(messages),
    G_relative_error = relative_G_error(unclass(vc[[fs_index]]), truth_G),
    intercept_variance = unname(vc[[intercept_index]][1L, 1L]),
    mean_rmse = sqrt(mean((mean_eta - dat$eta_mean)^2)),
    eta_rmse = sqrt(mean((fitted_eta - dat$eta)^2)),
    outer = NA_real_,
    diagnose_weights(dat$y, fitted_eta)
  )
}

fit_mgcv_newton <- function(dat) {
  fit <- mgcv::gam(
    y ~ s(time, k = 12L, bs = "cr") + s(id, bs = "re") +
      s(id, by = cos1, bs = "re") +
      s(id, by = cos2, bs = "re") +
      s(id, by = cos3, bs = "re"),
    data = dat, family = stats::inverse.gaussian(link = "log"),
    method = "REML"
  )
  random_labels <- vapply(fit$smooth, function(smooth) {
    if (grepl("id", smooth$label, fixed = TRUE)) smooth$label else ""
  }, character(1L))
  random_labels <- random_labels[nzchar(random_labels)]
  mean_eta <- as.numeric(stats::predict(
    fit, type = "link", exclude = random_labels
  ))
  fitted_eta <- as.numeric(stats::predict(fit, type = "link"))
  diagnostic <- diagnose_weights(dat$y, fitted_eta)
  if (length(fit$working.weights) == nrow(dat)) {
    diagnostic["mgcv_negative_fraction"] <-
      mean(fit$working.weights < 0)
    diagnostic["mgcv_minimum_working"] <- min(fit$working.weights)
  }
  c(
    converged = isTRUE(fit$converged), G_relative_error = NA_real_,
    intercept_variance = NA_real_,
    mean_rmse = sqrt(mean((mean_eta - dat$eta_mean)^2)),
    eta_rmse = sqrt(mean((fitted_eta - dat$eta)^2)),
    outer = fit$outer.info$iter, diagnostic
  )
}

results <- data.frame()
data_diagnostics <- data.frame()
for (drift_name in names(drift_scales)) {
  drift_scale <- unname(drift_scales[drift_name])
  truth_G <- truth_covariance(drift_scale)
  for (replicate_index in seq_len(n_rep)) {
    seed <- 20260816L + match(drift_name, names(drift_scales)) * 1000L +
      replicate_index
    dat <- simulate_data(drift_scale, seed)
    truth_diagnostic <- diagnose_weights(dat$y, dat$eta)
    data_diagnostics <- rbind(data_diagnostics, data.frame(
      drift = drift_name, replicate = replicate_index,
      negative_fraction = truth_diagnostic["negative_fraction"],
      minimum_observed = truth_diagnostic["minimum_observed"],
      alpha_q01 = truth_diagnostic["alpha_q01"], row.names = NULL
    ))

    methods <- list(
      gammfast = function() fit_gammfast(dat, truth_G),
      gamm4 = function() fit_gamm4(dat, truth_G),
      mgcv_observed_newton = function() fit_mgcv_newton(dat)
    )[method_names]
    for (method_name in names(methods)) {
      timed <- try(run_timed(methods[[method_name]]), silent = TRUE)
      if (inherits(timed, "try-error")) {
        warning(method_name, " failed: ", as.character(timed))
        results <- rbind(results, data.frame(
          method = method_name, drift = drift_name,
          replicate = replicate_index, elapsed = NA_real_,
          converged = FALSE, G_relative_error = NA_real_,
          intercept_variance = NA_real_, mean_rmse = NA_real_,
          eta_rmse = NA_real_, outer = NA_real_,
          negative_fraction = NA_real_,
          minimum_observed = NA_real_, minimum_fisher = NA_real_,
          alpha_q01 = NA_real_, alpha_median = NA_real_,
          mgcv_negative_fraction = NA_real_,
          mgcv_minimum_working = NA_real_
        ))
        next
      }
      value <- timed$value
      results <- rbind(results, data.frame(
        method = method_name, drift = drift_name,
        replicate = replicate_index, elapsed = timed$elapsed,
        converged = as.logical(value["converged"]),
        G_relative_error = value["G_relative_error"],
        intercept_variance = value["intercept_variance"],
        mean_rmse = value["mean_rmse"], eta_rmse = value["eta_rmse"],
        outer = value["outer"],
        negative_fraction = value["negative_fraction"],
        minimum_observed = value["minimum_observed"],
        minimum_fisher = value["minimum_fisher"],
        alpha_q01 = value["alpha_q01"],
        alpha_median = value["alpha_median"],
        mgcv_negative_fraction = if (
          "mgcv_negative_fraction" %in% names(value)
        ) value["mgcv_negative_fraction"] else NA_real_,
        mgcv_minimum_working = if (
          "mgcv_minimum_working" %in% names(value)
        ) value["mgcv_minimum_working"] else NA_real_,
        row.names = NULL
      ))
    }
  }
}

cat("\nObserved negative weights at the true simulated means:\n")
print(data_diagnostics, row.names = FALSE)
cat("\nFit-level diagnostics:\n")
print(results, row.names = FALSE)

aggregate_mean <- function(variable) {
  stats::aggregate(
    results[[variable]],
    list(method = results$method, drift = results$drift),
    function(x) mean(x, na.rm = TRUE)
  )
}
for (variable in c(
  "elapsed", "G_relative_error", "mean_rmse", "eta_rmse",
  "negative_fraction"
)) {
  cat("\nMean ", variable, ":\n", sep = "")
  print(aggregate_mean(variable), row.names = FALSE)
}

if (strict &&
    any(results$method == "gammfast" & !results$converged)) {
  stop("At least one gammfast Fisher fit failed to converge.")
}
cat("\nNegative-Newton-weight GAMM benchmark completed.\n")

library(mgcv)
library(mgcv.taps)

n_rep <- as.integer(Sys.getenv("GAMMFAST_INFLUENCE_REPS", "5"))
n_id <- as.integer(Sys.getenv("GAMMFAST_INFLUENCE_NID", "70"))
n_per_id <- as.integer(Sys.getenv("GAMMFAST_INFLUENCE_NPER", "8"))
if (n_rep < 1L || n_id < 10L || n_per_id < 4L) {
  stop("The influence benchmark dimensions are invalid.")
}

family_names <- c("gaussian", "binary", "nb", "sbeta")
component_names <- c("intercept", "slope1", "slope2")
mgcv_labels <- c("s(id)", "s(r1,id)", "s(r2,id)")
truth_by_family <- rbind(
  gaussian = c(intercept = 0.64, slope1 = 0.36, slope2 = 0.25),
  binary = c(intercept = 0.64, slope1 = 0.36, slope2 = 0.25),
  nb = c(intercept = 0.64, slope1 = 0.36, slope2 = 0.25),
  sbeta = c(intercept = 0.36, slope1 = 0.16, slope2 = 0.09)
)

family_factory <- function(name) {
  switch(name,
    gaussian = stats::gaussian(),
    binary = stats::binomial(),
    nb = mgcv::nb(theta = 3, link = "log"),
    sbeta = mgcv::betar(theta = 20, link = "logit"),
    stop("Unknown benchmark family.")
  )
}

simulate_data <- function(name, seed) {
  set.seed(seed)
  truth <- truth_by_family[name, ]
  id <- factor(rep(seq_len(n_id), each = n_per_id))
  id_index <- as.integer(id)
  grid <- seq(-1, 1, length.out = n_per_id)
  r1 <- rep(grid, n_id)
  r1 <- r1 / sqrt(mean(r1^2))
  r2 <- rep(grid^2 - mean(grid^2), n_id)
  r2 <- r2 / sqrt(mean(r2^2))
  x <- stats::runif(n_id * n_per_id, -1, 1)
  random <- cbind(
    stats::rnorm(n_id, sd = sqrt(truth[1L])),
    stats::rnorm(n_id, sd = sqrt(truth[2L])),
    stats::rnorm(n_id, sd = sqrt(truth[3L]))
  )
  eta <- 0.15 + 0.45 * sin(pi * x) +
    random[id_index, 1L] + r1 * random[id_index, 2L] +
    r2 * random[id_index, 3L]
  y <- switch(name,
    gaussian = eta + stats::rnorm(length(eta), sd = 0.55),
    binary = stats::rbinom(length(eta), 1L, stats::plogis(eta)),
    nb = stats::rnbinom(length(eta), mu = exp(eta), size = 3),
    sbeta = {
      mu <- stats::plogis(eta)
      stats::rbeta(length(eta), mu * 20, (1 - mu) * 20)
    }
  )
  data.frame(y = y, x = x, r1 = r1, r2 = r2, id = id)
}

model_formula <- y ~ s(x, k = 6L, bs = "cr") +
  s(id, bs = "re") + s(r1, id, bs = "re") +
  s(r2, id, bs = "re")

extract_mgcv_variance <- function(fit) {
  labels <- vapply(fit$smooth, function(s) s$label, character(1L))
  if (!all(mgcv_labels %in% labels)) {
    stop("mgcv::gam did not return the expected random-effect smooths.")
  }
  out <- vapply(mgcv_labels, function(label) {
    smooth <- fit$smooth[[match(label, labels)]]
    if (is.null(smooth$first.sp) || length(smooth$S.scale) != 1L) {
      stop("An mgcv random-effect smooth has invalid variance metadata.")
    }
    fit$sig2 * smooth$S.scale / fit$sp[smooth$first.sp]
  }, numeric(1L))
  setNames(out, component_names)
}

fit_mgcv <- function(dat, family_name) {
  fit <- mgcv::gam(
    model_formula, data = dat, family = family_factory(family_name),
    method = "REML"
  )
  list(
    variance = extract_mgcv_variance(fit), converged = isTRUE(fit$converged),
    outer = if (is.null(fit$outer.info$iter)) NA_integer_ else fit$outer.info$iter
  )
}

fit_gammfast <- function(dat, family_name, influence_update) {
  fit <- gammfast(
    model_formula, data = dat, family = family_factory(family_name),
    inner.max = 5L, inner.tol = 1e-5, nthreads = 1L,
    influence.update = influence_update,
    control = list(max.outer = 500L, objective.tol = 1e-5)
  )
  list(
    variance = setNames(diag(fit$G), component_names),
    converged = isTRUE(fit$converged), outer = fit$outer,
    variance_evaluations = if (is.null(fit$variance.evaluations)) {
      NA_integer_
    } else {
      fit$variance.evaluations
    }
  )
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

rows <- list()
row_index <- 0L
for (family_index in seq_along(family_names)) {
  family_name <- family_names[family_index]
  for (replicate_index in seq_len(n_rep)) {
    seed <- 20260900L + family_index * 100L + replicate_index
    dat <- simulate_data(family_name, seed)
    truth <- truth_by_family[family_name, ]
    gammfast_order <- if (replicate_index %% 2L) {
      c("frozen", "refreshed")
    } else {
      c("refreshed", "frozen")
    }
    method_order <- c("mgcv", gammfast_order)
    for (method in method_order) {
      result <- if (identical(method, "mgcv")) {
        run_timed(function() fit_mgcv(dat, family_name))
      } else {
        run_timed(function() fit_gammfast(dat, family_name, method))
      }
      if (is.null(result$value)) {
        row_index <- row_index + 1L
        rows[[row_index]] <- data.frame(
          family = family_name, replicate = replicate_index,
          method = method, component = component_names,
          estimate = NA_real_, truth = unname(truth),
          elapsed = result$elapsed, converged = FALSE,
          outer = NA_integer_, variance_evaluations = NA_integer_,
          error = result$error
        )
        next
      }
      for (component_index in seq_along(component_names)) {
        row_index <- row_index + 1L
        rows[[row_index]] <- data.frame(
          family = family_name, replicate = replicate_index,
          method = method, component = component_names[component_index],
          estimate = unname(result$value$variance[component_index]),
          truth = unname(truth[component_index]),
          elapsed = result$elapsed,
          converged = result$value$converged,
          outer = result$value$outer,
          variance_evaluations = if (is.null(
            result$value$variance_evaluations
          )) {
            NA_integer_
          } else {
            result$value$variance_evaluations
          },
          error = NA_character_
        )
      }
    }
  }
}

results <- do.call(rbind, rows)
key <- paste(results$family, results$replicate, results$component)
mgcv_rows <- results$method == "mgcv"
mgcv_lookup <- setNames(results$estimate[mgcv_rows], key[mgcv_rows])
results$mgcv <- unname(mgcv_lookup[key])
results$log_error_vs_mgcv <- ifelse(
  results$method == "mgcv", 0,
  abs(log(pmax(results$estimate, 1e-10) / pmax(results$mgcv, 1e-10)))
)
results$absolute_error_vs_mgcv <- abs(results$estimate - results$mgcv)
results$relative_error_vs_truth <-
  abs(results$estimate - results$truth) / results$truth
frozen_rows <- results$method == "frozen"
frozen_lookup <- setNames(results$estimate[frozen_rows], key[frozen_rows])
results$frozen <- unname(frozen_lookup[key])
results$relative_difference_vs_frozen <- ifelse(
  results$method == "refreshed",
  abs(results$estimate - results$frozen) / pmax(abs(results$frozen), 1e-8),
  0
)

summary_accuracy <- aggregate(
  cbind(
    estimate, absolute_error_vs_mgcv, log_error_vs_mgcv,
    relative_error_vs_truth, relative_difference_vs_frozen
  ) ~
    family + method + component,
  data = results, FUN = median, na.rm = TRUE
)
timing_rows <- results[!duplicated(results[c(
  "family", "replicate", "method"
)]), ]
summary_speed <- aggregate(
  elapsed ~ family + method, data = timing_rows,
  FUN = function(x) c(median = median(x), mean = mean(x))
)

cat("\nMedian variance-component accuracy:\n")
print(summary_accuracy, row.names = FALSE)
cat("\nElapsed seconds by family and method:\n")
print(summary_speed, row.names = FALSE)
cat("\nConvergence/error summary:\n")
print(aggregate(
  converged ~ family + method, data = timing_rows,
  FUN = function(x) paste(sum(x), "/", length(x))
), row.names = FALSE)

output_file <- Sys.getenv("GAMMFAST_INFLUENCE_BENCHMARK_OUTPUT", "")
if (nzchar(output_file)) {
  utils::write.csv(results, output_file, row.names = FALSE)
}

invisible(list(
  results = results, accuracy = summary_accuracy, speed = summary_speed
))

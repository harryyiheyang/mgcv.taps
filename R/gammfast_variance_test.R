#' Post-estimation variance-component test for gammfast
#'
#' Evaluates the fitted subject-level random structure using one completed
#' `gammfast()` fit. The fitted covariance `G`, global smoothing
#' parameters, dispersion, coefficients, and BLUPs are frozen. No null model
#' is fitted and no PIRLS, REML, BLUP, or smoothing-parameter iteration is run.
#'
#' The statistic is the Wood-type random-effect quadratic form obtained after
#' profiling the global mean coefficients while retaining their fitted
#' penalties. Its conditional reference distribution is a weighted sum of
#' independent chi-squared variables. For Gaussian fits this is the fitted
#' Gaussian working model. For every supported non-Gaussian family it is the
#' final PIRLS working-Gaussian model, using the family's final working weights
#' and fitted scale. The exact spectrum or randomized spectral moments are
#' evaluated in C++.
#'
#' This is a Wood-type conditional post-estimation evaluation of the joint
#' fitted random structure. It freezes all fitted parameters, performs no null
#' refit or Schur projection, and is not a Rao or Gaussian `U / sqrt(I)` score
#' test. Binary fits currently use this same conditional working-model
#' calibration; no CPQL small-sample correction is applied.
#'
#' @param fit A converged `gammfast` fit from any family supported by
#'   [gammfast()].
#' @param method Quadratic-form p-value method. `"auto"` uses Davies for an
#'   exact spectrum and randomized Liu for a stochastic spectrum.
#' @param spectrum Spectrum evaluation, `"auto"`, `"exact"`, or
#'   `"randomized"`. The automatic path switches above `q_threshold`.
#' @param q_threshold Maximum random-effect coefficient dimension for the
#'   automatic exact path.
#' @param n_probe Number of Rademacher probes for randomized spectral moments.
#' @param seed Reproducible random-projection seed.
#' @param max_eps Absolute error tolerance for Davies' method.
#' @param max_iter Maximum integration steps for Davies' method.
#' @param n_threads Number of subject-level OpenMP threads.
#'
#' @return A one-row `data.table` containing the statistic, p-value, reference
#'   rank, and computation details.
#' @export
gammfast_variance_test <- function(fit,
                                   method = c("auto", "davies", "liu"),
                                   spectrum = c("auto", "exact", "randomized"),
                                   q_threshold = 1000L,
                                   n_probe = 100L, seed = 20260810L,
                                   max_eps = 1e-8, max_iter = 1e5,
                                   n_threads = 1L) {
  method <- match.arg(method)
  spectrum <- match.arg(spectrum)
  requested_method <- method
  requested_spectrum <- spectrum
  if (!inherits(fit, "gammfast")) {
    stop("fit must be a 'gammfast' object.")
  }
  if (!isTRUE(fit$converged)) {
    stop("The gammfast fit must converge before variance-component testing.")
  }
  phi0 <- fit$sigma2
  if (length(phi0) != 1L || !is.finite(phi0) || phi0 <= 0) {
    stop("The gammfast fit has an invalid dispersion estimate.")
  }
  if (is.null(fit$G) || !is.matrix(fit$G) ||
      any(!is.finite(fit$G))) {
    stop("The gammfast fit has an invalid random-effect covariance.")
  }
  if (length(n_threads) != 1L || !is.finite(n_threads) || n_threads < 1) {
    stop("n_threads must be a positive integer.")
  }
  n_threads <- as.integer(n_threads)
  if (length(q_threshold) != 1L || !is.finite(q_threshold) ||
      q_threshold < 1 || q_threshold != as.integer(q_threshold)) {
    stop("q_threshold must be a positive integer.")
  }
  if (length(n_probe) != 1L || !is.finite(n_probe) || n_probe < 2 ||
      n_probe != as.integer(n_probe)) {
    stop("n_probe must be an integer of at least two.")
  }
  if (length(seed) != 1L || !is.finite(seed) || seed != as.integer(seed)) {
    stop("seed must be one finite integer.")
  }
  q_threshold <- as.integer(q_threshold)
  n_probe <- as.integer(n_probe)
  seed <- as.integer(seed)

  g <- fit$global
  X <- stats::predict(g, newdata = g$model, type = "lpmatrix") / sqrt(phi0)
  random <- gammfast_random_info(fit)
  B <- random$B / sqrt(phi0)
  gaussian_family <- identical(tolower(fit$family$family[1]), "gaussian")
  if (gaussian_family) {
    R_diag <- rep(1, nrow(X))
    working_model <- "Gaussian"
  } else {
    work <- gammfast_working(
      fit$family, fit$y, fit$linear.predictors, fit$prior.weights,
      nthreads = n_threads
    )
    if (any(work$w <= 0)) {
      stop("The final PIRLS working weights must be positive for variance testing.")
    }
    R_diag <- 1 / work$w
    working_model <- "final PIRLS working Gaussian"
  }
  p <- ncol(X)
  penalty <- matrix(0, p, p)
  if (length(g$S)) {
    if (length(g$S) != length(g$off) || length(g$S) != length(g$sp)) {
      stop("The fitted global penalties and smoothing parameters are inconsistent.")
    }
    for (j in seq_along(g$S)) {
      ii <- g$off[j] + seq_len(nrow(g$S[[j]])) - 1L
      penalty[ii, ii] <- penalty[ii, ii] + g$sp[j] * g$S[[j]] / phi0
    }
  }

  q <- nrow(fit$random.effects) * ncol(B)
  spectrum_used <- spectrum
  if (spectrum_used == "auto") {
    spectrum_used <- if (q <= q_threshold) "exact" else "randomized"
  }
  method_used <- method
  if (method_used == "auto") {
    method_used <- if (spectrum_used == "exact") "davies" else "liu"
  }
  if (spectrum_used == "randomized" && method_used == "davies") {
    stop("Davies calibration requires spectrum = 'exact'.")
  }

  if (spectrum_used == "exact") {
    spectrum_result <- gammfast_re_quadratic_spectrum(
      X = X, B = B, id = random$id.index, G = fit$G,
      R_diag = R_diag, penalty = penalty,
      random_effects = as.matrix(fit$random.effects),
      n_threads = n_threads
    )
    lambda <- spectrum_result$lambda
    lambda_scale <- max(lambda)
    if (!is.finite(lambda_scale) || lambda_scale <= 0) {
      stop("The variance-component reference distribution has zero rank.")
    }
    lambda <- lambda[lambda > lambda_scale * .Machine$double.eps^0.8]
    statistic <- spectrum_result$statistic
    moment_relative_mcse <- NA_real_
    reference_rank <- length(lambda)
    information_rank <- spectrum_result$information_rank
    probe_count <- NA_integer_
  } else {
    had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv)
    on.exit({
      if (had_seed) {
        assign(".Random.seed", old_seed, envir = .GlobalEnv)
      } else if (exists(".Random.seed", envir = .GlobalEnv,
                        inherits = FALSE)) {
        rm(".Random.seed", envir = .GlobalEnv)
      }
    }, add = TRUE)
    set.seed(seed)
    probes <- matrix(
      sample(c(-1, 1), q * n_probe, replace = TRUE),
      nrow = q, ncol = n_probe
    )
    spectrum_result <- gammfast_re_quadratic_random_moments(
      X = X, B = B, id = random$id.index, G = fit$G,
      R_diag = R_diag, penalty = penalty,
      random_effects = as.matrix(fit$random.effects), probes = probes,
      n_threads = n_threads
    )
    moments <- spectrum_result$moments
    moment_mcse <- apply(spectrum_result$probe_moments, 2L, stats::sd) /
      sqrt(n_probe)
    moment_relative_mcse <- max(moment_mcse / abs(moments))
    statistic <- spectrum_result$statistic
    reference_rank <- NA_integer_
    information_rank <- NA_integer_
    probe_count <- n_probe
  }

  p_method <- method_used
  fallback <- FALSE
  fallback_from <- NA_character_
  davies_ifault <- NA_integer_
  if (method_used == "davies") {
    calibration <- gammfast_quadratic_pvalue(
      statistic = statistic, lambda = lambda, method = "davies",
      max_eps = max_eps, max_iter = max_iter
    )
    p_value <- calibration$p.value
    p_method <- calibration$method
    fallback <- calibration$fallback
    fallback_from <- calibration$fallback.from
    davies_ifault <- calibration$davies.ifault
  } else if (spectrum_used == "exact") {
    p_value <- compute_liu_pvalue(statistic, lambda)
  } else {
    p_value <- gammfast_liu_moment_pvalue(statistic, moments)
    p_method <- "randomized-liu"
  }

  data.table::data.table(
    component = "joint fitted random structure",
    statistic = statistic,
    p.value = p_value,
    requested.method = requested_method,
    method = p_method,
    fallback = fallback,
    fallback.from = fallback_from,
    davies.ifault = davies_ifault,
    requested.spectrum = requested_spectrum,
    spectrum = spectrum_used,
    reference.rank = reference_rank,
    information.rank = information_rank,
    n.probe = probe_count,
    moment.relative.mcse = moment_relative_mcse,
    n.subject = spectrum_result$n_subject,
    basis.dimension = spectrum_result$basis_dimension,
    random.groups = length(random$groups),
    family = fit$family$family[1],
    working.model = working_model,
    cpql.corrected = FALSE,
    conditional = TRUE,
    null.refit = FALSE,
    post.estimation = TRUE,
    fitted.parameters.frozen = TRUE,
    gaussian.score = FALSE,
    schur = FALSE,
    full.random.design = FALSE
  )
}

gammfast_quadratic_pvalue <- function(statistic, lambda,
                                      method = c("davies", "liu"),
                                      max_eps = 1e-8, max_iter = 1e5) {
  method <- match.arg(method)
  if (length(statistic) != 1L || !is.finite(statistic) || statistic < 0) {
    stop("statistic must be one finite non-negative value.")
  }
  if (!is.numeric(lambda) || !length(lambda) ||
      any(!is.finite(lambda)) || any(lambda <= 0)) {
    stop("lambda must contain positive finite weights.")
  }
  if (method == "liu") {
    return(list(
      p.value = compute_liu_pvalue(statistic, lambda),
      method = "liu", fallback = FALSE,
      fallback.from = NA_character_, davies.ifault = NA_integer_
    ))
  }

  davies_result <- CompQuadForm::davies(
    q = statistic, lambda = lambda, lim = max_iter, acc = max_eps
  )
  p_value <- davies_result$Qq
  fallback <- !is.finite(p_value) || p_value <= 0 || p_value > 1 ||
    davies_result$ifault != 0
  if (fallback) {
    p_value <- compute_liu_pvalue(statistic, lambda)
  }
  list(
    p.value = p_value,
    method = if (fallback) "liu-fallback" else "davies",
    fallback = fallback,
    fallback.from = if (fallback) "davies" else NA_character_,
    davies.ifault = as.integer(davies_result$ifault)
  )
}

gammfast_liu_moment_pvalue <- function(q, moments) {
  c1 <- moments[1L]
  c2 <- moments[2L]
  c3 <- moments[3L]
  c4 <- moments[4L]
  if (length(moments) != 4L || any(!is.finite(moments)) ||
      any(moments <= 0)) {
    stop("The randomized spectral moments must be positive and finite.")
  }
  s1 <- c3 / c2^(3 / 2)
  s2 <- c4 / c2^2
  if (s1^2 > s2) {
    a <- 1 / (s1 - sqrt(s1^2 - s2))
    delta <- s1 * a^3 - a^2
    l <- a^2 - 2 * delta
  } else {
    delta <- 0
    l <- 1 / s2
    a <- sqrt(l)
  }
  muQ <- c1
  sigmaQ <- sqrt(2 * c2)
  tstar <- (q - muQ) / sigmaQ
  muX <- l + delta
  sigmaX <- sqrt(2) * a
  stats::pchisq(
    tstar * sigmaX + muX, df = l, ncp = delta, lower.tail = FALSE
  )
}

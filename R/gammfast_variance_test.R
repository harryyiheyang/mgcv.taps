#' Post-estimation subject-covariance quadratic evaluation for gammfast
#'
#' Evaluates the entire fitted subject-level covariance after projecting the
#' unpenalized fixed-effect space. Global penalized smooths enter the fitted
#' null working covariance as low-rank covariance components. The fitted
#' coefficients, smoothing parameters, dispersion, family parameters, and
#' subject covariance are frozen; no null model or fitting iteration is run.
#'
#' The statistic is `r' P0 K P0 r`, where `K` is the complete fitted
#' subject-level covariance and `P0` is the fixed-effect projection under the
#' target-free covariance `V0 = R + K_global`. Thus the tested subject
#' covariance supplies the quadratic kernel but is excluded from the null
#' projection covariance. Its reference distribution is a
#' weighted sum of independent chi-squared variables, evaluated by Davies or
#' randomized Liu. This is not a Rao score test and does not test `vec(G)`.
#'
#' @param fit A converged `gammfast` fit.
#' @param method Quadratic-form p-value method.
#' @param spectrum Exact or randomized spectrum evaluation.
#' @param q_threshold Maximum fitted subject-covariance factor dimension for
#'   automatic exact evaluation.
#' @param n_probe Number of Rademacher probes for randomized moments.
#' @param seed Reproducible randomized-spectrum seed.
#' @param max_eps Davies absolute error tolerance.
#' @param max_iter Davies maximum integration steps.
#' @param n_threads Number of OpenMP threads.
#'
#' @return A one-row `data.table` with the quadratic statistic, p-value, and
#'   computation diagnostics.
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
  if (!inherits(fit, "gammfast")) stop("fit must be a 'gammfast' object.")
  if (!isTRUE(fit$converged)) {
    stop("The gammfast fit must converge before variance-component evaluation.")
  }
  phi0 <- fit$sig2
  if (length(phi0) != 1L || !is.finite(phi0) || phi0 <= 0) {
    stop("The gammfast fit has an invalid dispersion estimate.")
  }
  if (is.null(fit$G) || !is.matrix(fit$G) || any(!is.finite(fit$G))) {
    stop("The gammfast fit has an invalid subject covariance.")
  }
  if (length(n_threads) != 1L || !is.finite(n_threads) || n_threads < 1) {
    stop("n_threads must be a positive integer.")
  }
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
  n_threads <- as.integer(n_threads)
  q_threshold <- as.integer(q_threshold)
  n_probe <- as.integer(n_probe)
  seed <- as.integer(seed)

  g <- fit$global
  X <- stats::predict(g, newdata = g$model, type = "lpmatrix") / sqrt(phi0)
  random <- gammfast_random_info(fit)
  B <- random$B / sqrt(phi0)
  gaussian_family <- identical(tolower(fit$family$family[1]), "gaussian")
  if (gaussian_family) {
    response <- (fit$y - fit$offset) / sqrt(phi0)
    R_diag <- rep(1, length(response))
    working_model <- "Gaussian"
  } else {
    work <- gammfast_working(
      fit$family, fit$y, fit$linear.predictors, fit$prior.weights,
      nthreads = n_threads
    )
    if (any(work$w <= 0)) {
      stop("The final PIRLS working weights must be positive.")
    }
    response <- (work$z - fit$offset) / sqrt(phi0)
    R_diag <- 1 / work$w
    working_model <- "final PIRLS working Gaussian"
  }
  penalty <- gammfast_penalty_matrix(g, g$sp, scale = phi0)

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

  probes <- matrix(numeric(), 0L, 0L)
  if (spectrum_used == "randomized") {
    had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv)
    on.exit({
      if (had_seed) {
        assign(".Random.seed", old_seed, envir = .GlobalEnv)
      } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
        rm(".Random.seed", envir = .GlobalEnv)
      }
    }, add = TRUE)
    set.seed(seed)
    probes <- matrix(sample(c(-1, 1), q * n_probe, replace = TRUE),
                     nrow = q, ncol = n_probe)
  }
  result <- gammfast_variance_quadratic(
    response = response, X = X, beta = fit$coefficients,
    B = B, id = random$id.index, G = fit$G, R_diag = R_diag,
    penalty = penalty, probes = probes,
    exact = identical(spectrum_used, "exact"),
    eigen_tol = 1e-8, n_threads = n_threads
  )
  if (result$g_negative > 0L) {
    warning(result$g_negative, " negative fitted-covariance eigenvalue(s) were removed; minimum = ",
            format(result$g_min_eigen, digits = 4), ".")
  }
  if (result$penalty_negative > 0L) {
    warning(result$penalty_negative, " negative global-penalty eigenvalue(s) were removed; minimum = ",
            format(result$penalty_min_eigen, digits = 4), ".")
  }
  if (result$reference_negative > 0L) {
    warning(result$reference_negative, " negative reference eigenvalue(s) were removed; minimum = ",
            format(result$reference_min_eigen, digits = 4), ".")
  }

  statistic <- result$statistic
  moment_relative_mcse <- NA_real_
  reference_rank <- NA_integer_
  probe_count <- NA_integer_
  reference_mean <- NA_real_
  reference_sd <- NA_real_
  if (spectrum_used == "exact") {
    lambda <- result$lambda
    lambda <- lambda[is.finite(lambda) & lambda > 0]
    if (!length(lambda)) {
      warning("The fitted subject covariance has zero positive reference rank.")
      p_value <- NA_real_
      p_method <- NA_character_
      fallback <- FALSE
      fallback_from <- NA_character_
      davies_ifault <- NA_integer_
    } else if (method_used == "davies") {
      calibration <- gammfast_quadratic_pvalue(
        statistic, lambda, "davies", max_eps, max_iter
      )
      p_value <- calibration$p.value
      p_method <- calibration$method
      fallback <- calibration$fallback
      fallback_from <- calibration$fallback.from
      davies_ifault <- calibration$davies.ifault
    } else {
      p_value <- compute_liu_pvalue(statistic, lambda)
      p_method <- "liu"
      fallback <- FALSE
      fallback_from <- NA_character_
      davies_ifault <- NA_integer_
    }
    reference_rank <- length(lambda)
    if (length(lambda)) {
      reference_mean <- sum(lambda)
      reference_sd <- sqrt(2 * sum(lambda^2))
    }
  } else {
    moments <- result$moments
    probe_moments <- result$probe_moments
    moment_mcse <- apply(probe_moments, 2L, stats::sd) / sqrt(n_probe)
    moment_relative_mcse <- max(moment_mcse / pmax(abs(moments), .Machine$double.eps))
    p_value <- gammfast_liu_moment_pvalue(statistic, moments)
    p_method <- "randomized-liu"
    fallback <- FALSE
    fallback_from <- NA_character_
    davies_ifault <- NA_integer_
    probe_count <- n_probe
    reference_mean <- moments[1L]
    reference_sd <- sqrt(2 * moments[2L])
  }

  data.table::data.table(
    component = "joint fitted subject covariance",
    statistic = statistic, p.value = p_value,
    requested.method = requested_method, method = p_method,
    fallback = fallback, fallback.from = fallback_from,
    davies.ifault = davies_ifault,
    requested.spectrum = requested_spectrum, spectrum = spectrum_used,
    reference.rank = reference_rank, n.probe = probe_count,
    reference.mean = reference_mean, reference.sd = reference_sd,
    statistic.to.reference.mean = statistic / reference_mean,
    moment.relative.mcse = moment_relative_mcse,
    n.subject = result$n_subject,
    basis.dimension = result$basis_dimension,
    fixed.effect.rank = result$fixed_rank,
    global.smooth.rank = result$smooth_rank,
    random.groups = length(random$groups), family = fit$family$family[1],
    working.model = working_model,
    quadratic.form = "r' P0 K_subject P0 r",
    projection.covariance = "R + K_global_smooth",
    tested.covariance.in.projection = FALSE,
    fixed.effect.projected = TRUE,
    global.smooth.as.covariance = TRUE,
    vec.G.test = FALSE, score.test = FALSE,
    conditional = TRUE, null.refit = FALSE, post.estimation = TRUE,
    fitted.parameters.frozen = TRUE, gaussian.score = FALSE,
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
    return(list(p.value = compute_liu_pvalue(statistic, lambda),
                method = "liu", fallback = FALSE,
                fallback.from = NA_character_, davies.ifault = NA_integer_))
  }
  davies_result <- CompQuadForm::davies(
    q = statistic, lambda = lambda, lim = max_iter, acc = max_eps
  )
  p_value <- davies_result$Qq
  fallback <- !is.finite(p_value) || p_value <= 0 || p_value > 1 ||
    davies_result$ifault != 0
  if (fallback) p_value <- compute_liu_pvalue(statistic, lambda)
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
  if (length(moments) != 4L || any(!is.finite(moments)) || any(moments <= 0)) {
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
  stats::pchisq(tstar * sigmaX + muX, df = l, ncp = delta,
                lower.tail = FALSE)
}

#' Unified Conditional Mean-Structure Evaluation for GAMMs
#'
#' Evaluates a penalized mean-structure smooth in a fitted model that contains an
#' `re` or `fs` random-effect structure. The function dispatches to the fitted
#' covariance representation supplied by `mgcv::gam()` or `mgcv::bam()`,
#' `gamm4::gamm4()`, or [gammfast()]. It never refits a null model.
#'
#' For native `gam` and supported `bam` fits, fitted `re`/`fs` smooths are
#' retained as penalized nuisance covariance components. For `gamm4`, the
#' sparse `lme4` covariance factor is combined with the fitted mgcv smooth
#' representation. For `gammfast`, the subject-block covariance
#' `B_i G B_i^T + R_i` is applied by the blockwise C++ inverse kernel; it is not
#' replaced by a diagonal working covariance.
#'
#' This is a conditional post-estimation quadratic evaluation. Fitted
#' coefficients, smoothing parameters, dispersion, and random-effect
#' covariance parameters are frozen; no null model is refitted. The historical
#' score-test name is retained for API compatibility, but the result is not a
#' Rao score test and is not calibrated as `U / sqrt(I)` by a Gaussian limit.
#'
#' @param fit A fitted `gam`, `bam`, `gamm4`, or `gammfast` object containing
#'   an `re` or `fs` structure.
#' @param test.component Integer index of the mean-structure smooth to test.
#'   An `re` or `fs` smooth cannot itself be selected.
#' @param null.tol Penalty null-space detection tolerance.
#' @param method Quadratic-form p-value method.
#' @param max_eps Absolute error tolerance for Davies' method.
#' @param max_iter Maximum integration steps for Davies' method.
#' @param eps_mu Pseudo-response tolerance for native mgcv fits.
#' @param n_threads Number of threads for applicable working-model and
#'   subject-block operations.
#'
#' @return A `data.table` containing the tested mean smooth, p-value, method,
#'   and explicit conditional/no-refit provenance. Attribute `backend` records
#'   the selected covariance backend.
#' @export
taps_score_test_gamm <- function(fit, test.component = 1,
                                 null.tol = 1e-10,
                                 method = "davies",
                                 max_eps = 1e-8, max_iter = 1e5,
                                 eps_mu = 1e-12, n_threads = 1L) {
  backend <- gamm_score_backend(fit)
  g <- if (backend == "gamm4") fit$gam else if (backend == "gammfast") {
    fit$global
  } else {
    fit
  }
  gamm_validate_mean_component(g, test.component)

  if (backend == "gammfast") {
    ans <- taps_score_test_gammfast_impl(
      fit = fit, test.component = test.component, null.tol = null.tol,
      method = method, max_eps = max_eps, max_iter = max_iter,
      n_threads = n_threads
    )
  } else if (backend == "gamm4") {
    ans <- taps_score_test_gamm4_re(
      fit = fit, test.component = test.component, null.tol = null.tol,
      method = method, max_eps = max_eps, max_iter = max_iter
    )
  } else {
    ans <- taps_score_test_re(
      fit = fit, test.component = test.component, null.tol = null.tol,
      method = method, max_eps = max_eps, max_iter = max_iter,
      eps_mu = eps_mu, n_threads = n_threads
    )
  }
  ans[["conditional"]] <- TRUE
  ans[["null.refit"]] <- FALSE
  ans[["post.estimation"]] <- TRUE
  ans[["calibration"]] <- "weighted quadratic"
  ans[["gaussian.score"]] <- FALSE
  ans[["schur"]] <- FALSE
  attr(ans, "backend") <- backend
  attr(ans, "conditional") <- TRUE
  attr(ans, "null.refit") <- FALSE
  attr(ans, "post.estimation") <- TRUE
  if (backend == "gammfast") {
    attr(ans, "full.random.design") <- FALSE
  }
  ans
}

gamm_score_backend <- function(fit) {
  if (inherits(fit, "gammfast")) return("gammfast")
  if (inherits(fit, "gamm4")) return("gamm4")
  if (!inherits(fit, "gam")) {
    stop("fit must be a fitted 'gam', 'bam', 'gamm4', or 'gammfast' object.")
  }
  has_random <- any(vapply(fit$smooth, gamm_is_random_smooth, logical(1)))
  if (!has_random) {
    stop("A native gam or bam fit must contain an re or fs smooth.")
  }
  "mgcv"
}

gamm_is_random_smooth <- function(s) {
  any(grepl("random.effect|fs.interaction", class(s)))
}

gamm_validate_mean_component <- function(g, test.component) {
  p <- length(g$smooth)
  if (p < 1L || length(test.component) != 1L || is.na(test.component)) {
    stop("test.component must be a valid mean-smooth index.")
  }
  test_index <- suppressWarnings(as.integer(test.component))
  if (is.na(test_index) || test.component != test_index ||
      test_index < 1L || test_index > p) {
    stop("test.component must be a valid mean-smooth index.")
  }
  if (gamm_is_random_smooth(g$smooth[[test_index]])) {
    stop("test.component must select a mean-structure smooth, not an re or fs smooth.")
  }
  invisible(test_index)
}

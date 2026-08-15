#' Conditional Post-Estimation Evaluation for a Specified Smooth Term
#'
#' Computes a fitted-model quadratic evaluation for the penalized component of
#' a selected smooth term in a fitted \code{gam} or \code{bam} object. The
#' model is handled exactly through its mgcv smooth representation. Hence
#' mgcv \code{re}/\code{fs} terms, when present, are ordinary fitted smoothing
#' components in this interface; no separate random-effect covariance backend
#' is activated. The remaining smooth terms are treated as
#' nuisance components and are profiled out via a structured covariance
#' approximation. P-values can be computed by several methods for quadratic
#' forms in normal variables.
#'
#' The historical score-test name is retained for compatibility. With
#' \code{refit = FALSE}, the fitted coefficients, smoothing parameters, and
#' dispersion are frozen, so the result is a conditional post-estimation
#' evaluation rather than a Rao score test. With \code{refit = TRUE}, the tested
#' AMatern or A2Matern penalized component is removed and the null model is
#' refitted from the current fit. Calibration uses the weighted quadratic-form
#' spectrum; it is not a Gaussian `U / sqrt(I)` test.
#'
#' @param fit A fitted \code{gam} or \code{bam} model object from \pkg{mgcv}.
#'   Use [taps_score_test_gamm()] for \code{gamm4} or \code{gammfast} fits, or
#'   when an explicit fitted random-effect covariance backend is required.
#' @param test.component Integer. Index of the smooth term to be tested. Default is \code{1}.
#' @param null.tol Numeric. Row-norm threshold used to detect null-space basis
#'   columns of the penalty matrix when \code{getA} is unavailable. Default is \code{1e-10}.
#' @param method Character. Method for computing the p-value of the score
#'   statistic. One of \code{"davies"} (default), \code{"liu"}, \code{"hbe"},
#'   \code{"wood"}, \code{"satterthwaite"}, or \code{"saddlepoint"}.
#' @param max_eps Numeric. Absolute error tolerance passed to
#'   \code{CompQuadForm::davies}. Default is \code{1e-8}.
#' @param max_iter Integer. Maximum number of integration steps passed to
#'   \code{CompQuadForm::davies}. Default is \code{1e5}.
#' @param eps_mu Numeric. Tolerance passed to \code{extract_pseudo_response}.
#'   Default is \code{1e-12}.
#' @param n_threads Integer. Number of threads used by supported pseudo-response
#'   and structured random-effect calculations. Default is \code{1}.
#' @param refit Logical. If \code{TRUE}, refit the null mean model by replacing
#'   the tested AMatern or A2Matern smooth with its \code{getA} covariates.
#'   Default is \code{FALSE}.
#' @param sp.refit Logical. If \code{TRUE}, re-estimate nuisance smoothing
#'   parameters by REML for \code{gam} or fREML for \code{bam} during a null
#'   refit. If \code{FALSE}, retain their fitted values. Default is \code{TRUE}.
#'   Discrete \code{bam} refits use mgcv's default starting values because
#'   fitted coefficient or smoothing-parameter starts can change fREML
#'   convergence.
#'
#' @return A \code{data.table} with three columns:
#'   \describe{
#'     \item{smooth.term}{Label of the tested smooth term.}
#'     \item{smooth.pvalue}{Score test p-value.}
#'     \item{method}{The p-value method used.}
#'   }
#'
#' @importFrom mgcv vcov.gam
#' @importFrom data.table data.table rbindlist
#' @importFrom ordgam testStat
#' @importFrom Matrix bdiag
#' @import CppMatrix
#' @importFrom survey pchisqsum
#' @import CompQuadForm
#'
#' @examples
#' \dontrun{
#' library(mgcv)
#' fit <- gam(y ~ s(x0) + s(x1), data = dat, method = "REML")
#' taps_score_test(fit, test.component = 1)
#' }
#'
#' @export
taps_score_test <- function(fit, test.component = 1, null.tol = 1e-10,
                            method = "davies", max_eps = 1e-8, max_iter = 1e5,
                            eps_mu = 1e-12, n_threads = 1,
                            refit = FALSE, sp.refit = TRUE) {
  if (inherits(fit, "gammfast")) {
    stop("taps_score_test does not accept gammfast fits; use taps_score_test_gamm().")
  }
  if (inherits(fit, "gamm4")) {
    stop("taps_score_test does not accept gamm4 fits; use taps_score_test_gamm().")
  }
  if (!inherits(fit, "gam")) stop("fit must be a 'gam' or 'bam' object.")
  if (length(refit) != 1L || is.na(refit) || !is.logical(refit)) {
    stop("refit must be a single TRUE or FALSE value.")
  }
  if (length(sp.refit) != 1L || is.na(sp.refit) || !is.logical(sp.refit)) {
    stop("sp.refit must be a single TRUE or FALSE value.")
  }

  if (refit) fit <- taps_score_refit(fit, test.component, sp.refit, null.tol)

  if (identical(fit$family$family, "Cox PH")) {
    return(taps_score_test_cox(
      fit = fit, test.component = test.component, null.tol = null.tol,
      method = method, max_eps = max_eps, max_iter = max_iter
    ))
  }

  if (grepl("^zero inflated poisson", tolower(fit$family$family))) {
    return(taps_score_test_zip(
      fit = fit, test.component = test.component, null.tol = null.tol,
      method = method, max_eps = max_eps, max_iter = max_iter
    ))
  }

  res <- extract_pseudo_response(
    fit, eps_mu = eps_mu, n_threads = n_threads
  )

  pseudo_response <- res$pseudo_response
  V_phi           <- res$V_phi
  phi0            <- res$phi0

  model_offset <- fit$offset
  if (is.null(model_offset)) model_offset <- rep(0, length(pseudo_response))
  model_offset <- as.numeric(model_offset)
  if (length(model_offset) != length(pseudo_response) ||
      any(!is.finite(model_offset))) {
    stop("fit has an invalid model offset.")
  }
  pseudo_response <- pseudo_response - model_offset

  beta         <- fit$coefficients
  X            <- fit$.taps_score_X
  if (is.null(X)) X <- predict(fit, newdata = fit$model, type = "lpmatrix")
  smooth_terms <- fit$smooth
  p            <- length(smooth_terms)

  if (p < 1L) stop("fit must contain at least one smooth term.")
  if (length(test.component) != 1L || is.na(test.component)) {
    stop("test.component must be a valid smooth-term index.")
  }
  test_component_int <- suppressWarnings(as.integer(test.component))
  if (is.na(test_component_int) || test.component != test_component_int ||
      test_component_int < 1L || test_component_int > p) {
    stop("test.component must be a valid smooth-term index.")
  }
  test.component <- test_component_int

  if (!is.null(res$valid_idx)) {
    idx             <- res$valid_idx
    if (sum(idx) == 0) stop("No valid observations for testing")
    pseudo_response <- pseudo_response[idx]
    V_phi           <- V_phi[idx]
    X               <- X[idx, , drop = FALSE]
  }

  smooth_index_list <- list()
  random_index_list <- list()
  S_list            <- list()

  get_scaled_penalty <- function(s) {
    if (is.null(s$first.sp) || is.null(s$last.sp)) {
      stop("Penalized smooth has no smoothing-parameter index.")
    }
    sp_idx <- seq.int(s$first.sp, s$last.sp)
    if (length(sp_idx) != length(s$S)) {
      stop("Penalized smooth has inconsistent penalty and smoothing-parameter indexes.")
    }
    sp <- fit$sp[sp_idx]
    if (anyNA(sp)) {
      stop("Penalized smooth has missing smoothing-parameter values.")
    }
    Reduce(`+`, Map(function(S_matrix, sp_value) S_matrix * sp_value / phi0,
                    s$S, sp))
  }

  for (i in seq_len(p)) {
    s       <- smooth_terms[[i]]
    indices <- s$first.para:s$last.para
    reported_null_dim <- s$null.space.dim
    is_zero_rank <- !is.null(s$rank) && isTRUE(all(s$rank == 0))
    is_fixed_smooth <- isTRUE(s$fixed) || is.null(s$S) || length(s$S) == 0 || is_zero_rank

    if (i == test.component && is_fixed_smooth) {
      stop("test.component must be a penalized smooth with fx = FALSE.")
    }

    if (is_fixed_smooth) next

    S_matrix <- get_scaled_penalty(s)
    S_norm <- norm(S_matrix, "f")

    if (!is.finite(S_norm) || S_norm <= 0) {
      if (i == test.component) {
        stop("test.component must have a non-zero penalty matrix.")
      }
      next
    }

    if (is.null(s$getA) == 0) {
      detected_null_indices <- if (reported_null_dim > 0) indices[seq_len(reported_null_dim)] else integer(0)
    } else {
      col_norms <- sqrt(colSums(S_matrix^2))
      detected_null_indices <- indices[which(col_norms < null.tol)]
    }

    smooth_indices <- setdiff(indices, detected_null_indices)

    if (i != test.component) {
      smooth_index_list[[i]] <- indices
      random_index_list[[i]] <- smooth_indices
      S_list[[i]]            <- S_matrix
    }

    if (i == test.component) {
      Bj       <- X[, indices]
      Thetaj   <- matrixGeneralizedInverse(S_matrix / S_norm)
      Gj_apply <- function(v) {
        v_is_matrix <- is.matrix(v)
        v_mat       <- if (v_is_matrix) v else matrix(v, ncol = 1)
        Bt_v        <- matrixMultiply(Bj, v_mat, transA = TRUE)
        out         <- matrixMultiply(Bj, matrixMultiply(Thetaj, Bt_v))
        if (v_is_matrix) out else as.vector(out)
      }
      random_index_list[[i]] <- smooth_indices
    }
  }

  S_list            <- Filter(Negate(is.null), S_list)
  smooth_index_list <- Filter(Negate(is.null), smooth_index_list)
  random_index_list <- Filter(Negate(is.null), random_index_list)
  smooth_index_vec  <- do.call(c, smooth_index_list)
  random_index_vec  <- do.call(c, random_index_list)
  fixed_index_vec   <- setdiff(seq_len(ncol(X)), random_index_vec)

  A        <- as.matrix(X[, fixed_index_vec])
  alpha    <- beta[fixed_index_vec]
  B_extend <- X[, smooth_index_vec]
  S_All    <- as.matrix(bdiag(S_list))
  XtX      <- matrixMultiply(B_extend, B_extend * (1 / V_phi), transA = TRUE)
  C        <- matrixInverse(XtX + S_All)

  Vinv_apply <- function(v) {
    v_is_matrix <- is.matrix(v)
    v_mat       <- if (v_is_matrix) v else matrix(v, ncol = 1)
    part1       <- v_mat / V_phi
    Bt_v        <- matrixMultiply(B_extend, part1, transA = TRUE)
    C_Bt_v      <- matrixMultiply(C, Bt_v)
    part2       <- matrixMultiply(B_extend, C_Bt_v)
    out         <- part1 - part2 / V_phi
    if (v_is_matrix) out else as.vector(out)
  }

  Vinv_A      <- Vinv_apply(A)
  XtVinvX     <- matrixMultiply(A, Vinv_A, transA = TRUE)
  XtVinvX_inv <- matrixGeneralizedInverse(XtVinvX)

  P_apply <- function(v) {
    v_is_matrix  <- is.matrix(v)
    v_mat        <- if (v_is_matrix) v else matrix(v, ncol = 1)
    Vinv_v       <- Vinv_apply(v_mat)
    AVinv_v      <- matrixMultiply(A, Vinv_v, transA = TRUE)
    solve_middle <- matrixMultiply(XtVinvX_inv, AVinv_v)
    out          <- Vinv_v - matrixMultiply(Vinv_A, solve_middle)
    if (v_is_matrix) out else as.vector(out)
  }

  if (method == "satterthwaite") {
    error  <- pseudo_response - matrixVectorMultiply(A, alpha)
    r      <- Vinv_apply(error)
    Gj_r   <- Gj_apply(r)
    u      <- max(0,sum(r * Gj_r) / 2)

    q  <- ncol(Bj)
    Cj <- matrix(0, q, q)
    for (i in 1:q) {
      Pi <- P_apply(Bj[, i])
      for (j in i:q) {
        Cj[i, j] <- sum(Bj[, j] * Pi)
        if (i != j) Cj[j, i] <- Cj[i, j]
      }
    }
    PGj  <- matrixMultiply(Cj, Thetaj)
    PGj2 <- matrixMultiply(PGj, PGj)
    e    <- sum(diag(PGj)) / 2
    h    <- sum(diag(PGj2)) / 2
    kappa <- h / (2 * e)
    nu    <- 2 * e^2 / h
    pv    <- pchisq(u / kappa, df = nu, lower.tail = FALSE)
  } else {
    error  <- pseudo_response - matrixVectorMultiply(A, alpha)
    r      <- P_apply(error)
    Gj_r   <- Gj_apply(r)
    u      <- max(0,sum(r * Gj_r))
    q      <- ncol(Bj)
    eig_theta <- matrixsqrt(Thetaj)
    Theta_sqrt <- eig_theta$w
    N      <- P_apply(Bj)
    BtPB   <- matrixMultiply(Bj, N, transA = TRUE)
    Q_small <- matrixListProduct(list(Theta_sqrt, BtPB, Theta_sqrt))
    lambda  <- eigen(Q_small, symmetric = TRUE, only.values = TRUE)$values
    lambda  <- lambda[lambda > 1e-15]

    pv <- switch(method,
                 liu        = compute_liu_pvalue(q = u, lambda = lambda),
                 hbe        = compute_hbe_pvalue(q = u, lambda = lambda),
                 wood       = compute_wood_pvalue(q = u, lambda = lambda),
                 saddlepoint = survey::pchisqsum(
                   x = u, df = rep(1, length(lambda)), a = lambda,
                   lower.tail = FALSE, method = "saddlepoint"
                 ),
                 davies = {
                   davies_res <- CompQuadForm::davies(q = u, lambda = lambda,
                                                      lim = max_iter, acc = max_eps)
                   if (davies_res$Qq <= 0 || davies_res$Qq > 1.0) {
                     compute_liu_pvalue(u, lambda)
                   } else {
                     davies_res$Qq
                   }
                 },
                 stop("method must be one of 'satterthwaite', 'davies', 'liu', 'hbe', 'wood', 'saddlepoint'")
    )
  }

  data.table(
    smooth.term    = smooth_terms[[test.component]]$label,
    smooth.pvalue  = pv,
    method         = method
  )
}

#' Experimental Schur-Projected Quadratic Diagnostic
#'
#' Computes an experimental fitted-model quadratic diagnostic for one penalized smooth after projecting its
#' variance-component score off the score span of selected nuisance smooths.
#' Only ordinary mgcv \code{gam}/\code{bam} fits are supported. Fitted mgcv
#' \code{re}/\code{fs} terms are permitted and remain ordinary smoothing
#' components; this function does not activate a separate random-effect
#' covariance backend.
#' Cox and zero-inflated Poisson models are intentionally not handled by this
#' function.
#'
#' All fitted quantities are frozen and no null model is refitted. This
#' diagnostic is separate from the primary conditional definition used by
#' [taps_score_test_gamm()] and is not called by that interface. Its reference
#' law is evaluated from the signed eigenvalues of the projected quadratic
#' kernel; it is never standardized as `U / sqrt(I)` or calibrated by a
#' Gaussian score limit.
#'
#' Each entry of \code{schur.component} identifies a whole smooth term. If a
#' selected smooth has several penalties, its fitted, smoothing-parameter
#' weighted penalty is treated as one local covariance direction.
#'
#' @param fit A fitted mgcv \code{gam} or \code{bam} object. \code{gamm4} and
#'   \code{gammfast} objects are not supported.
#' @param test.component Integer index of the smooth term to be tested.
#' @param schur.component Integer vector indexing nuisance smooth terms. Use
#'   \code{integer(0)} for no Schur projection.
#' @param null.tol Numeric tolerance used to detect unpenalized columns.
#' @param method Null calibration. Currently only \code{"davies"}, using the
#'   signed eigenvalues of the Schur-projected quadratic form, is accepted.
#'   A labeled saddlepoint fallback is used if Davies returns an invalid tail
#'   probability for a signed kernel.
#' @param max_eps Absolute error tolerance passed to
#'   \code{CompQuadForm::davies}.
#' @param max_iter Maximum integration steps passed to
#'   \code{CompQuadForm::davies}.
#' @param eps_mu Tolerance passed to \code{extract_pseudo_response} for
#'   ordinary \code{gam}/\code{bam} fits.
#' @param n_threads Number of threads passed to \code{extract_pseudo_response}.
#'
#' @return A one-row \code{data.table} containing the p-value, raw and
#'   efficient scores, raw and efficient information, and the fraction of
#'   target information removed by the Schur projection. Per-component
#'   projection coefficients, standardized Fisher couplings, and signed
#'   information-loss contributions are returned as labeled strings separated
#'   by \code{", "}, so multiple nuisance terms remain distinguishable in a
#'   data-frame or CSV cell.
#'
#' @examples
#' \dontrun{
#' library(mgcv)
#' fit <- gam(y ~ s(x0) + s(x1), data = dat, method = "REML")
#' taps_score_test_schur(fit, test.component = 1, schur.component = 2)
#' }
#'
#' @export
taps_score_test_schur <- function(fit, test.component = 1,
                                  schur.component = integer(0),
                                  null.tol = 1e-10,
                                  method = "davies",
                                  max_eps = 1e-8, max_iter = 1e5,
                                  eps_mu = 1e-12, n_threads = 1) {
  if (length(method) != 1L || !identical(method, "davies")) {
    stop("method must be 'davies'.")
  }

  if (inherits(fit, "gammfast") || inherits(fit, "gamm4")) {
    stop("taps_score_test_schur accepts only mgcv gam or bam fits.")
  }

  if (!inherits(fit, "gam")) {
    stop("fit must be an mgcv 'gam' or 'bam' object.")
  }
  if (identical(fit$family$family, "Cox PH") ||
      grepl("^zero inflated poisson", tolower(fit$family$family))) {
    stop("taps_score_test_schur currently supports only ordinary mgcv gam/bam GLM-style fits.")
  }

  res <- extract_pseudo_response(
    fit, eps_mu = eps_mu, n_threads = n_threads
  )
  pseudo_response <- res$pseudo_response
  V_phi <- res$V_phi
  X_keep <- NULL
  if (!is.null(res$valid_idx)) {
    X_keep <- res$valid_idx
    if (sum(X_keep) == 0L) stop("No valid observations for testing.")
    pseudo_response <- pseudo_response[X_keep]
    V_phi <- V_phi[X_keep]
  }

  V0_inv_apply <- function(v) {
    v_is_matrix <- is.matrix(v)
    v_mat <- if (v_is_matrix) v else matrix(v, ncol = 1L)
    out <- v_mat / V_phi
    if (v_is_matrix) out else as.vector(out)
  }

  offset <- stats::model.offset(fit$model)
  if (is.null(offset)) offset <- numeric(length(res$pseudo_response))
  if (!is.null(X_keep)) offset <- offset[X_keep]

  taps_score_test_schur_core(
    g = fit,
    pseudo_response = pseudo_response,
    V_phi = V_phi,
    phi0 = res$phi0,
    offset = offset,
    V0_inv_apply = V0_inv_apply,
    test.component = test.component,
    schur.component = schur.component,
    null.tol = null.tol,
    method = method,
    max_eps = max_eps,
    max_iter = max_iter,
    row_keep = X_keep
  )
}

validate_schur_components <- function(smooth_terms, test.component,
                                      schur.component) {
  p <- length(smooth_terms)
  if (p < 1L) stop("fit must contain at least one smooth term.")

  if (length(test.component) != 1L || is.na(test.component)) {
    stop("test.component must be a valid smooth-term index.")
  }
  test_int <- suppressWarnings(as.integer(test.component))
  if (is.na(test_int) || test.component != test_int ||
      test_int < 1L || test_int > p) {
    stop("test.component must be a valid smooth-term index.")
  }

  if (is.null(schur.component)) schur.component <- integer(0)
  if (!is.atomic(schur.component) || anyNA(schur.component)) {
    stop("schur.component must be an integer vector of smooth-term indexes.")
  }
  schur_int <- suppressWarnings(as.integer(schur.component))
  if (length(schur_int) > 0L &&
      (any(schur.component != schur_int) || any(schur_int < 1L) ||
       any(schur_int > p))) {
    stop("schur.component must contain valid smooth-term indexes.")
  }
  if (anyDuplicated(schur_int)) {
    stop("schur.component must not contain duplicated indexes.")
  }
  if (test_int %in% schur_int) {
    stop("test.component cannot also be a schur.component.")
  }

  list(test = test_int, schur = schur_int)
}

expand_penalty_to_coefficients <- function(S_matrix, n_coef, smooth_label) {
  S_matrix <- as.matrix(S_matrix)
  if (nrow(S_matrix) != ncol(S_matrix)) {
    stop("Penalty matrix for '", smooth_label, "' must be square.")
  }
  if (nrow(S_matrix) == n_coef) return(S_matrix)
  if (nrow(S_matrix) == 0L || n_coef %% nrow(S_matrix) != 0L) {
    stop(
      "Penalty dimension for '", smooth_label,
      "' is incompatible with its model-matrix columns."
    )
  }

  n_block <- n_coef %/% nrow(S_matrix)
  kronecker(diag(n_block), S_matrix)
}

schur_scaled_penalty <- function(s, sp, phi0, n_coef) {
  if (is.null(s$first.sp) || is.null(s$last.sp)) {
    stop("Penalized smooth has no smoothing-parameter index.")
  }
  sp_idx <- seq.int(s$first.sp, s$last.sp)
  if (length(sp_idx) != length(s$S)) {
    stop("Penalized smooth has inconsistent penalty and smoothing-parameter indexes.")
  }
  sp_value <- sp[sp_idx]
  if (anyNA(sp_value)) {
    stop("Penalized smooth has missing smoothing-parameter values.")
  }
  Reduce(`+`, Map(
    function(S_matrix, value) {
      expand_penalty_to_coefficients(S_matrix, n_coef, s$label) *
        value / phi0
    },
    s$S, sp_value
  ))
}

format_schur_report <- function(labels, values) {
  if (!length(values)) return("")
  value_text <- formatC(values, digits = 8L, format = "fg", flag = "#")
  paste(paste0(labels, "=", value_text), collapse = ", ")
}

taps_score_test_schur_core <- function(g, pseudo_response, V_phi, phi0,
                                       offset, V0_inv_apply,
                                       test.component, schur.component,
                                       null.tol, method, max_eps, max_iter,
                                       row_keep = NULL) {
  smooth_terms <- g$smooth
  component_index <- validate_schur_components(
    smooth_terms, test.component, schur.component
  )
  test.component <- component_index$test
  schur.component <- component_index$schur

  X <- stats::predict(g, newdata = g$model, type = "lpmatrix")
  if (!is.null(row_keep)) X <- X[row_keep, , drop = FALSE]
  if (nrow(X) != length(pseudo_response) ||
      length(offset) != length(pseudo_response) ||
      any(!is.finite(offset))) {
    stop("Model matrix, offset, and pseudo-response have incompatible lengths.")
  }

  beta <- g$coefficients
  p <- length(smooth_terms)
  smooth_index_list <- list()
  random_index_list <- list()
  S_list <- list()
  score_components <- vector("list", 1L + length(schur.component))
  score_term_index <- c(test.component, schur.component)

  for (i in seq_len(p)) {
    s <- smooth_terms[[i]]
    indices <- s$first.para:s$last.para
    reported_null_dim <- s$null.space.dim
    is_zero_rank <- !is.null(s$rank) && isTRUE(all(s$rank == 0))
    is_fixed_smooth <- isTRUE(s$fixed) || is.null(s$S) ||
      length(s$S) == 0L || is_zero_rank

    if (i %in% score_term_index && is_fixed_smooth) {
      stop("test.component and schur.component terms must be penalized smooths with fx = FALSE.")
    }
    if (is_fixed_smooth) next

    S_matrix <- schur_scaled_penalty(
      s, g$sp, phi0, n_coef = length(indices)
    )
    S_norm <- norm(S_matrix, "f")
    if (!is.finite(S_norm) || S_norm <= 0) {
      if (i %in% score_term_index) {
        stop("test.component and schur.component terms must have non-zero penalty matrices.")
      }
      next
    }

    if (!is.null(s$getA)) {
      detected_null_indices <- if (reported_null_dim > 0L) {
        indices[seq_len(reported_null_dim)]
      } else {
        integer(0)
      }
    } else {
      col_norms <- sqrt(colSums(S_matrix^2))
      detected_null_indices <- indices[col_norms < null.tol]
    }
    smooth_indices <- setdiff(indices, detected_null_indices)

    if (i != test.component) {
      smooth_index_list[[i]] <- indices
      random_index_list[[i]] <- smooth_indices
      S_list[[i]] <- S_matrix
    } else {
      random_index_list[[i]] <- smooth_indices
    }

    component_position <- match(i, score_term_index)
    if (!is.na(component_position)) {
      score_components[[component_position]] <- list(
        B = as.matrix(X[, indices, drop = FALSE]),
        theta = matrixGeneralizedInverse(S_matrix / S_norm),
        label = s$label,
        index = i
      )
    }
  }

  if (any(vapply(score_components, is.null, logical(1)))) {
    stop("Failed to construct one or more requested score components.")
  }

  S_list <- Filter(Negate(is.null), S_list)
  smooth_index_list <- Filter(Negate(is.null), smooth_index_list)
  random_index_list <- Filter(Negate(is.null), random_index_list)
  smooth_index_vec <- if (length(smooth_index_list)) {
    unlist(smooth_index_list, use.names = FALSE)
  } else {
    integer(0)
  }
  random_index_vec <- if (length(random_index_list)) {
    unlist(random_index_list, use.names = FALSE)
  } else {
    integer(0)
  }
  fixed_index_vec <- setdiff(seq_len(ncol(X)), random_index_vec)

  A <- as.matrix(X[, fixed_index_vec, drop = FALSE])
  alpha <- beta[fixed_index_vec]
  B_extend <- as.matrix(X[, smooth_index_vec, drop = FALSE])

  if (ncol(B_extend) > 0L) {
    S_All <- as.matrix(Matrix::bdiag(S_list))
    V0_B <- V0_inv_apply(B_extend)
    XtV0X <- matrixMultiply(B_extend, V0_B, transA = TRUE)
    C <- matrixInverse(XtV0X + S_All)

    Vinv_apply <- function(v) {
      v_is_matrix <- is.matrix(v)
      v_mat <- if (v_is_matrix) v else matrix(v, ncol = 1L)
      V0_v <- V0_inv_apply(v_mat)
      Bt_v <- matrixMultiply(B_extend, V0_v, transA = TRUE)
      out <- V0_v - matrixMultiply(V0_B, matrixMultiply(C, Bt_v))
      if (v_is_matrix) out else as.vector(out)
    }
  } else {
    Vinv_apply <- V0_inv_apply
  }

  Vinv_A <- Vinv_apply(A)
  XtVinvX <- matrixMultiply(A, Vinv_A, transA = TRUE)
  XtVinvX_inv <- matrixGeneralizedInverse(XtVinvX)

  P_apply <- function(v) {
    v_is_matrix <- is.matrix(v)
    v_mat <- if (v_is_matrix) v else matrix(v, ncol = 1L)
    Vinv_v <- Vinv_apply(v_mat)
    AVinv_v <- matrixMultiply(A, Vinv_v, transA = TRUE)
    solve_middle <- matrixMultiply(XtVinvX_inv, AVinv_v)
    out <- Vinv_v - matrixMultiply(Vinv_A, solve_middle)
    if (v_is_matrix) out else as.vector(out)
  }

  error <- pseudo_response - offset - matrixVectorMultiply(A, alpha)
  r <- P_apply(error)
  B_score <- do.call(cbind, lapply(score_components, `[[`, "B"))
  PB_score <- P_apply(B_score)
  M <- matrixMultiply(B_score, PB_score, transA = TRUE)
  M <- (M + t(M)) / 2

  component_width <- vapply(score_components, function(x) ncol(x$B), integer(1))
  component_end <- cumsum(component_width)
  component_start <- component_end - component_width + 1L
  component_columns <- Map(seq.int, component_start, component_end)
  n_component <- length(score_components)

  score <- numeric(n_component)
  trace_score <- numeric(n_component)
  information <- matrix(0, n_component, n_component)
  Bt_r <- matrixMultiply(B_score, r, transA = TRUE)

  for (i in seq_len(n_component)) {
    ii <- component_columns[[i]]
    theta_i <- score_components[[i]]$theta
    z_i <- Bt_r[ii, , drop = FALSE]
    q_i <- as.numeric(matrixMultiply(
      z_i, matrixMultiply(theta_i, z_i), transA = TRUE
    ))
    trace_i <- sum(diag(matrixMultiply(M[ii, ii, drop = FALSE], theta_i)))
    score[i] <- (q_i - trace_i) / 2
    trace_score[i] <- trace_i

    for (j in i:n_component) {
      jj <- component_columns[[j]]
      theta_j <- score_components[[j]]$theta
      info_ij <- matrixListProduct(list(
        theta_i,
        M[ii, jj, drop = FALSE],
        theta_j,
        M[jj, ii, drop = FALSE]
      ))
      information[i, j] <- sum(diag(info_ij)) / 2
      information[j, i] <- information[i, j]
    }
  }

  I_FF <- information[1L, 1L]
  if (!is.finite(I_FF) || I_FF <= 0) {
    stop("The target score has non-positive Fisher information.")
  }

  if (n_component == 1L) {
    schur_coef <- numeric(0)
    schur_coupling <- numeric(0)
    schur_contribution <- numeric(0)
    schur_rank <- 0L
    efficient_score <- score[1L]
    efficient_information <- I_FF
  } else {
    I_FS <- information[1L, -1L, drop = FALSE]
    I_SS <- information[-1L, -1L, drop = FALSE]
    schur_rank <- qr(I_SS, tol = 1e-10)$rank
    I_SS_inv <- matrixGeneralizedInverse(I_SS)
    schur_coef <- as.numeric(matrixMultiply(I_FS, I_SS_inv))
    schur_scale <- sqrt(I_FF * pmax(0, diag(I_SS)))
    schur_coupling <- rep(NA_real_, length(schur_scale))
    positive_scale <- is.finite(schur_scale) & schur_scale > 0
    schur_coupling[positive_scale] <-
      as.numeric(I_FS)[positive_scale] / schur_scale[positive_scale]
    schur_contribution <- schur_coef * as.numeric(I_FS) / I_FF
    efficient_score <- score[1L] - sum(schur_coef * score[-1L])
    efficient_information <- as.numeric(
      I_FF - matrixListProduct(list(I_FS, I_SS_inv, t(I_FS)))
    )
  }

  information_tol <- 1e-10 * max(1, I_FF)
  if (!is.finite(efficient_information) ||
      efficient_information <= information_tol) {
    stop("Schur projection leaves no positive target information.")
  }
  schur_fraction <- max(0, min(1, 1 - efficient_information / I_FF))
  schur_labels <- vapply(
    score_components[-1L], `[[`, character(1), "label"
  )

  kernel_weight <- c(1, -schur_coef)
  q_component <- 2 * score + trace_score
  efficient_q <- sum(kernel_weight * q_component)
  C_eff <- as.matrix(Matrix::bdiag(Map(
    function(component, weight) component$theta * weight,
    score_components, kernel_weight
  )))
  M_sqrt <- matrixsqrt(M)$w
  Q_small <- matrixListProduct(list(M_sqrt, C_eff, M_sqrt))
  Q_small <- (Q_small + t(Q_small)) / 2
  lambda <- eigen(Q_small, symmetric = TRUE, only.values = TRUE)$values
  lambda_tol <- 1e-12 * max(1, max(abs(lambda)))
  lambda <- lambda[abs(lambda) > lambda_tol]
  if (length(lambda) == 0L) {
    stop("Schur-projected score kernel has no non-zero eigenvalues.")
  }

  davies_ifault <- NA_integer_
  davies_res <- CompQuadForm::davies(
    q = efficient_q,
    lambda = lambda,
    lim = max_iter,
    acc = max_eps
  )
  davies_ifault <- as.integer(davies_res$ifault)
  pv <- davies_res$Qq
  if (!is.finite(pv) || pv <= 0 || pv > 1) {
    if (all(lambda > 0)) {
      pv <- compute_liu_pvalue(efficient_q, lambda)
      method <- "liu-fallback"
    } else {
      pv <- survey::pchisqsum(
        x = efficient_q, df = rep(1, length(lambda)), a = lambda,
        lower.tail = FALSE, method = "saddlepoint"
      )
      method <- "saddlepoint-fallback"
      if (!is.finite(pv) || pv <= 0 || pv > 1) {
        stop(
          "Davies and saddlepoint calibration failed for the signed ",
          "Schur-projected kernel."
        )
      }
    }
  }

  data.table::data.table(
    smooth.term = score_components[[1L]]$label,
    schur.terms = paste(
      schur_labels, collapse = ", "
    ),
    schur.coefficients = format_schur_report(schur_labels, schur_coef),
    schur.couplings = format_schur_report(schur_labels, schur_coupling),
    schur.contributions = format_schur_report(
      schur_labels, schur_contribution
    ),
    smooth.pvalue = pv,
    score = score[1L],
    efficient.score = efficient_score,
    information = I_FF,
    efficient.information = efficient_information,
    schur.fraction = schur_fraction,
    schur.rank = schur_rank,
    min.eigenvalue = min(lambda),
    max.eigenvalue = max(lambda),
    davies.ifault = davies_ifault,
    method = method,
    conditional = TRUE,
    null.refit = FALSE,
    post.estimation = TRUE,
    gaussian.score = FALSE,
    experimental = TRUE
  )
}

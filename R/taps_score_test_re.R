taps_score_test_re <- function(fit, test.component = 1, null.tol = 1e-10,
                               method = "davies", max_eps = 1e-8,
                               max_iter = 1e5, eps_mu = 1e-12,
                               n_threads = 1) {
  if (length(n_threads) != 1L || is.na(n_threads) ||
      !is.numeric(n_threads) || !is.finite(n_threads) ||
      n_threads < 1 || n_threads > .Machine$integer.max ||
      n_threads != as.integer(n_threads)) {
    stop("n_threads must be a positive integer for include_re = TRUE.")
  }
  n_threads <- as.integer(n_threads)
  if (is.null(fit$sig2) || length(fit$sig2) != 1L ||
      !is.numeric(fit$sig2) || !is.finite(fit$sig2) || fit$sig2 <= 0) {
    stop("include_re = TRUE requires a finite positive fitted dispersion.")
  }
  res <- extract_pseudo_response(
    fit, eps_mu = eps_mu, n_threads = n_threads
  )

  pseudo_response <- res$pseudo_response
  V_phi           <- res$V_phi
  phi0            <- res$phi0

  beta         <- fit$coefficients
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

  if (!any(vapply(smooth_terms, native_re_is_smooth, logical(1)))) {
    stop("include_re = TRUE requires at least one mgcv re or fs smooth.")
  }
  if (native_re_is_smooth(smooth_terms[[test.component]])) {
    stop("include_re = TRUE supports re and fs smooths as nuisance terms only.")
  }

  data <- native_re_design_data(fit)

  if (!is.null(res$valid_idx)) {
    idx             <- res$valid_idx
    if (sum(idx) == 0) stop("No valid observations for testing")
    pseudo_response <- pseudo_response[idx]
    V_phi           <- V_phi[idx]
    data            <- data[idx, , drop = FALSE]
  }

  fixed_blocks      <- list()
  dense_blocks      <- list()
  local_blocks      <- list()
  random_index_list <- list()
  dense_S_list      <- list()

  if (fit$nsdf > 0L) {
    B_param <- stats::model.matrix(fit$pterms, data)
    if (ncol(B_param) != fit$nsdf) {
      stop("Parametric design has inconsistent coefficient dimensions.")
    }
    fixed_blocks[["parametric"]] <- list(
      indices = seq_len(fit$nsdf), X = B_param
    )
  }

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
    is_fixed_smooth <- isTRUE(s$fixed) || is.null(s$S) || length(s$S) == 0 ||
      is_zero_rank

    if (i == test.component && is_fixed_smooth) {
      stop("test.component must be a penalized smooth with fx = FALSE.")
    }

    B_dense <- if (native_re_is_smooth(s)) NULL else {
      native_re_dense_smooth_design(s, data)
    }

    if (is_fixed_smooth) {
      if (native_re_is_smooth(s)) {
        stop("include_re = TRUE requires re and fs nuisance smooths to be penalized.")
      }
      fixed_blocks[[as.character(i)]] <- list(indices = indices, X = B_dense)
      next
    }

    S_matrix <- get_scaled_penalty(s)
    S_norm <- norm(S_matrix, "f")

    if (!is.finite(S_norm) || S_norm <= 0) {
      if (i == test.component) {
        stop("test.component must have a non-zero penalty matrix.")
      }
      if (native_re_is_smooth(s)) {
        stop("Native re and fs nuisance smooths must have a non-zero penalty.")
      }
      fixed_blocks[[as.character(i)]] <- list(indices = indices, X = B_dense)
      next
    }

    if (is.null(s$getA) == 0) {
      detected_null_indices <- if (reported_null_dim > 0) {
        indices[seq_len(reported_null_dim)]
      } else integer(0)
    } else {
      col_norms <- sqrt(colSums(S_matrix^2))
      detected_null_indices <- indices[which(col_norms < null.tol)]
    }

    smooth_indices <- setdiff(indices, detected_null_indices)

    if (i != test.component) {
      if (native_re_is_smooth(s)) {
        if (length(detected_null_indices)) {
          stop("Native re and fs nuisance penalties must be full rank within uid.")
        }
        local_design <- native_re_local_smooth_design(s, data)
        local_penalty <- native_re_local_penalty(
          S_matrix, n_group = length(local_design$levels),
          block_size = ncol(local_design$X)
        )
        local_blocks[[as.character(i)]] <- list(
          indices = indices, X = local_design$X,
          group = local_design$group, levels = local_design$levels,
          S = local_penalty
        )
      } else {
        dense_blocks[[as.character(i)]] <- list(indices = indices, X = B_dense)
        dense_S_list[[i]] <- S_matrix
      }
      random_index_list[[i]] <- smooth_indices
    }

    if (i == test.component) {
      Bj       <- B_dense
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

    if (length(detected_null_indices)) {
      fixed_blocks[[as.character(i)]] <- list(
        indices = detected_null_indices,
        X = B_dense[, match(detected_null_indices, indices), drop = FALSE]
      )
    }
  }

  fixed_blocks      <- Filter(Negate(is.null), fixed_blocks)
  dense_blocks      <- Filter(Negate(is.null), dense_blocks)
  local_blocks      <- Filter(Negate(is.null), local_blocks)
  random_index_list <- Filter(Negate(is.null), random_index_list)
  dense_S_list      <- Filter(Negate(is.null), dense_S_list)

  fixed_design <- native_re_bind_dense(fixed_blocks, length(pseudo_response))
  dense_design <- native_re_bind_dense(dense_blocks, length(pseudo_response))
  local_design <- native_re_bind_local(local_blocks, length(pseudo_response))
  random_index_vec <- if (length(random_index_list)) {
    do.call(c, random_index_list)
  } else integer(0)
  fixed_index_vec <- setdiff(seq_along(beta), random_index_vec)
  if (!identical(fixed_design$indices, fixed_index_vec)) {
    stop("Direct design reconstruction did not recover the fixed coefficient columns.")
  }

  A        <- fixed_design$X
  alpha    <- beta[fixed_design$indices]
  dense_S <- native_re_block_diag(dense_S_list)
  profile <- native_re_diagonal_profile(
    B_dense = dense_design$X, B_local = local_design$X,
    group = local_design$group, V_phi = V_phi,
    dense_S = dense_S, local_S = local_design$S,
    n_threads = n_threads
  )
  Vinv_apply <- profile$apply

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
    u      <- max(0, sum(r * Gj_r) / 2)

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
    u      <- max(0, sum(r * Gj_r))
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
                 stop("method must be one of 'satterthwaite', 'davies', 'liu', 'hbe', 'wood', 'saddlepoint'"))
  }

  data.table(
    smooth.term    = smooth_terms[[test.component]]$label,
    smooth.pvalue  = pv,
    method         = method
  )
}

taps_score_test_gammfast_impl <- function(fit, test.component = 1,
                                           null.tol = 1e-10,
                                           method = "davies",
                                           max_eps = 1e-8, max_iter = 1e5,
                                           n_threads = 1L) {
  if (!isTRUE(fit$converged)) {
    stop("The gammfast fit must converge before a TAPS score test.")
  }
  g <- fit$global
  work <- gammfast_working(
    fit$family, fit$y, fit$linear.predictors, fit$prior.weights,
    nthreads = n_threads
  )
  model_offset <- fit$offset
  if (is.null(model_offset)) model_offset <- g$offset
  if (is.null(model_offset)) model_offset <- numeric(length(work$z))
  model_offset <- as.numeric(model_offset)
  if (length(model_offset) != length(work$z) ||
      any(!is.finite(model_offset))) {
    stop("The gammfast fit has an invalid model offset.")
  }
  phi0 <- fit$sigma2
  if (length(phi0) != 1L || !is.finite(phi0) || phi0 <= 0) {
    stop("The gammfast fit has an invalid dispersion estimate.")
  }
  if (is.null(fit$G) || !is.matrix(fit$G) || any(!is.finite(fit$G))) {
    stop("The gammfast fit has an invalid random-effect covariance.")
  }

  X <- stats::predict(g, newdata = g$model, type = "lpmatrix")
  sw <- sqrt(work$w / phi0)
  X_work <- X * sw
  z_work <- (work$z - model_offset) * sw
  random <- gammfast_random_info(fit)
  B_work <- random$B * sw
  R_work <- rep(1, length(z_work))

  V0_inv_apply <- function(v) {
    v_is_matrix <- is.matrix(v)
    v_mat <- if (v_is_matrix) v else matrix(v, ncol = 1L)
    out <- gammfast_vinv_apply(
      v_mat, B_work, random$id.index, fit$G, R_work,
      n_threads = n_threads
    )
    if (v_is_matrix) out else as.vector(out)
  }

  gammfast_smoothing_score_core(
    g = g, X = X_work, working_response = z_work,
    V0_inv_apply = V0_inv_apply, phi0 = phi0,
    test.component = test.component, null.tol = null.tol,
    method = method, max_eps = max_eps, max_iter = max_iter
  )
}

gammfast_smoothing_score_core <- function(g, X, working_response,
                                          V0_inv_apply, phi0, test.component,
                                          null.tol, method, max_eps,
                                          max_iter) {
  beta <- g$coefficients
  smooth_terms <- g$smooth
  p <- length(smooth_terms)
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

  get_scaled_penalty <- function(s) {
    if (is.null(s$first.sp) || is.null(s$last.sp)) {
      stop("Penalized smooth has no smoothing-parameter index.")
    }
    sp_idx <- seq.int(s$first.sp, s$last.sp)
    if (length(sp_idx) != length(s$S)) {
      stop("Penalized smooth has inconsistent penalty and smoothing-parameter indexes.")
    }
    sp <- g$sp[sp_idx]
    if (anyNA(sp)) {
      stop("Penalized smooth has missing smoothing-parameter values.")
    }
    Reduce(`+`, Map(function(S_matrix, sp_value) {
      S_matrix * sp_value / phi0
    }, s$S, sp))
  }

  smooth_index_list <- list()
  random_index_list <- list()
  S_list <- list()
  for (i in seq_len(p)) {
    s <- smooth_terms[[i]]
    indices <- s$first.para:s$last.para
    is_zero_rank <- !is.null(s$rank) && isTRUE(all(s$rank == 0))
    is_fixed_smooth <- isTRUE(s$fixed) || is.null(s$S) ||
      length(s$S) == 0L || is_zero_rank
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
      detected_null_indices <- if (s$null.space.dim > 0) {
        indices[seq_len(s$null.space.dim)]
      } else {
        integer()
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
      Bj <- X[, indices, drop = FALSE]
      Thetaj <- matrixGeneralizedInverse(S_matrix / S_norm)
      Gj_apply <- function(v) {
        v_is_matrix <- is.matrix(v)
        v_mat <- if (v_is_matrix) v else matrix(v, ncol = 1L)
        Bt_v <- matrixMultiply(Bj, v_mat, transA = TRUE)
        out <- matrixMultiply(Bj, matrixMultiply(Thetaj, Bt_v))
        if (v_is_matrix) out else as.vector(out)
      }
      random_index_list[[i]] <- smooth_indices
    }
  }

  S_list <- Filter(Negate(is.null), S_list)
  smooth_index_list <- Filter(Negate(is.null), smooth_index_list)
  random_index_list <- Filter(Negate(is.null), random_index_list)
  smooth_index_vec <- if (length(smooth_index_list)) {
    unlist(smooth_index_list, use.names = FALSE)
  } else {
    integer()
  }
  random_index_vec <- if (length(random_index_list)) {
    unlist(random_index_list, use.names = FALSE)
  } else {
    integer()
  }
  fixed_index_vec <- setdiff(seq_len(ncol(X)), random_index_vec)
  A <- as.matrix(X[, fixed_index_vec, drop = FALSE])
  alpha <- beta[fixed_index_vec]
  B_extend <- X[, smooth_index_vec, drop = FALSE]

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
    out <- Vinv_v - matrixMultiply(
      Vinv_A, matrixMultiply(XtVinvX_inv, AVinv_v)
    )
    if (v_is_matrix) out else as.vector(out)
  }

  error <- working_response - matrixVectorMultiply(A, alpha)
  if (method == "satterthwaite") {
    r <- P_apply(error)
    u <- max(0, sum(r * Gj_apply(r)) / 2)
    q <- ncol(Bj)
    Cj <- matrix(0, q, q)
    for (i in seq_len(q)) {
      Pi <- P_apply(Bj[, i])
      for (j in i:q) {
        Cj[i, j] <- sum(Bj[, j] * Pi)
        if (i != j) Cj[j, i] <- Cj[i, j]
      }
    }
    PGj <- matrixMultiply(Cj, Thetaj)
    PGj2 <- matrixMultiply(PGj, PGj)
    e <- sum(diag(PGj)) / 2
    h <- sum(diag(PGj2)) / 2
    kappa <- h / (2 * e)
    nu <- 2 * e^2 / h
    pv <- stats::pchisq(u / kappa, df = nu, lower.tail = FALSE)
  } else {
    r <- P_apply(error)
    u <- max(0, sum(r * Gj_apply(r)))
    Theta_sqrt <- matrixsqrt(Thetaj)$w
    N <- P_apply(Bj)
    BtPB <- matrixMultiply(Bj, N, transA = TRUE)
    Q_small <- matrixListProduct(list(Theta_sqrt, BtPB, Theta_sqrt))
    lambda <- eigen(Q_small, symmetric = TRUE, only.values = TRUE)$values
    lambda <- lambda[lambda > 1e-15]
    pv <- switch(method,
      liu = compute_liu_pvalue(q = u, lambda = lambda),
      hbe = compute_hbe_pvalue(q = u, lambda = lambda),
      wood = compute_wood_pvalue(q = u, lambda = lambda),
      saddlepoint = survey::pchisqsum(
        x = u, df = rep(1, length(lambda)), a = lambda,
        lower.tail = FALSE, method = "saddlepoint"
      ),
      davies = {
        davies_res <- CompQuadForm::davies(
          q = u, lambda = lambda, lim = max_iter, acc = max_eps
        )
        if (davies_res$Qq <= 0 || davies_res$Qq > 1) {
          compute_liu_pvalue(u, lambda)
        } else {
          davies_res$Qq
        }
      },
      stop("method must be one of 'satterthwaite', 'davies', 'liu', 'hbe', 'wood', 'saddlepoint'")
    )
  }

  data.table::data.table(
    smooth.term = smooth_terms[[test.component]]$label,
    smooth.pvalue = pv,
    method = method
  )
}

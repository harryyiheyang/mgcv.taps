gamm4_v0_inv_apply_re <- function(V_phi, Zt, Lambdat, scale) {
  if (nrow(Zt) == 0L) {
    return(function(v) {
      v_is_matrix <- is.matrix(v)
      v_mat <- if (v_is_matrix) v else matrix(v, ncol = 1)
      out <- v_mat / V_phi
      if (v_is_matrix) out else as.vector(out)
    })
  }

  F <- sqrt(scale) * Matrix::t(Lambdat %*% Zt)
  R_inv <- Matrix::Diagonal(x = 1 / V_phi)
  R_inv_F <- R_inv %*% F
  K <- Matrix::Diagonal(n = ncol(F)) + Matrix::crossprod(F, R_inv_F)
  K <- Matrix::forceSymmetric(K, uplo = "L")
  K_chol <- Matrix::Cholesky(K, LDL = FALSE, perm = TRUE)

  function(v) {
    v_is_matrix <- is.matrix(v)
    v_mat <- if (v_is_matrix) v else matrix(v, ncol = 1)
    R_inv_v <- v_mat / V_phi
    F_R_inv_v <- Matrix::crossprod(F, R_inv_v)
    K_inv_F_R_inv_v <- Matrix::solve(K_chol, F_R_inv_v)
    out <- R_inv_v - R_inv_F %*% K_inv_F_R_inv_v
    out <- as.matrix(out)
    if (v_is_matrix) out else as.vector(out)
  }
}

extract_random_block_gamm4_re <- function(fit) {
  if (!inherits(fit$gam, "gam") || !inherits(fit$mer, "merMod")) {
    stop("gamm4 fit must contain compatible $gam and $mer components.")
  }
  re_block <- extract_random_block_gamm4(fit)
  if (!inherits(re_block$Zt, "sparseMatrix") ||
      !inherits(re_block$Lambdat, "sparseMatrix") ||
      ncol(re_block$Zt) != nrow(fit$gam$model)) {
    stop("gamm4 lme4 random-effect matrices are inconsistent or not sparse.")
  }
  re_block
}

taps_score_test_gamm4_re <- function(fit, test.component = 1,
                                      null.tol = 1e-10,
                                      method = "davies", max_eps = 1e-8,
                                      max_iter = 1e5) {
  if (!inherits(fit, "gamm4")) {
    stop("fit must be a 'gamm4' object.")
  }
  if (!requireNamespace("lme4", quietly = TRUE)) {
    stop("Package 'lme4' is required for gamm4 TAPS tests.")
  }
  if (utils::packageVersion("lme4") < numeric_version("2.0.6")) {
    stop("gamm4 TAPS requires lme4 >= 2.0.6.")
  }

  g <- fit$gam
  res <- extract_pseudo_response_gamm4(fit)
  pseudo_response <- res$pseudo_response
  V_phi <- res$V_phi
  phi0 <- res$phi0
  offset <- stats::model.offset(g$model)
  if (is.null(offset)) offset <- numeric(length(pseudo_response))
  if (length(offset) != length(pseudo_response) ||
      any(!is.finite(offset))) {
    stop("gamm4 model offset is missing, non-finite, or has the wrong length.")
  }
  beta <- g$coefficients
  X <- predict(g, newdata = g$model, type = "lpmatrix")
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

  re_block <- extract_random_block_gamm4_re(fit)
  if (nrow(re_block$Zt) == 0L) {
    stop("include_re = TRUE requires a gamm4 lme4 random-effect structure.")
  }
  V0_inv_apply <- gamm4_v0_inv_apply_re(
    V_phi = V_phi, Zt = re_block$Zt,
    Lambdat = re_block$Lambdat, scale = re_block$scale
  )

  smooth_index_list <- list()
  random_index_list <- list()
  S_list <- list()

  for (i in seq_len(p)) {
    s <- smooth_terms[[i]]
    indices <- s$first.para:s$last.para
    reported_null_dim <- s$null.space.dim
    is_zero_rank <- !is.null(s$rank) && isTRUE(all(s$rank == 0))
    is_fixed_smooth <- isTRUE(s$fixed) || is.null(s$S) ||
      length(s$S) == 0 || is_zero_rank

    if (i == test.component && is_fixed_smooth) {
      stop("test.component must be a penalized smooth with fx = FALSE.")
    }
    if (is_fixed_smooth) next

    S_matrix <- gamm4_scaled_penalty(s, g$sp, phi0)
    S_norm <- norm(S_matrix, "f")
    if (!is.finite(S_norm) || S_norm <= 0) {
      if (i == test.component) {
        stop("test.component must have a non-zero penalty matrix.")
      }
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
      smooth_index_list[[i]] <- indices
      random_index_list[[i]] <- smooth_indices
      S_list[[i]] <- S_matrix
    }
    if (i == test.component) {
      Bj <- X[, indices, drop = FALSE]
      Thetaj <- matrixGeneralizedInverse(S_matrix / S_norm)
      Gj_apply <- function(v) {
        v_is_matrix <- is.matrix(v)
        v_mat <- if (v_is_matrix) v else matrix(v, ncol = 1)
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
    do.call(c, smooth_index_list)
  } else integer(0)
  random_index_vec <- if (length(random_index_list)) {
    do.call(c, random_index_list)
  } else integer(0)
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
      v_mat <- if (v_is_matrix) v else matrix(v, ncol = 1)
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
    v_mat <- if (v_is_matrix) v else matrix(v, ncol = 1)
    Vinv_v <- Vinv_apply(v_mat)
    AVinv_v <- matrixMultiply(A, Vinv_v, transA = TRUE)
    solve_middle <- matrixMultiply(XtVinvX_inv, AVinv_v)
    out <- Vinv_v - matrixMultiply(Vinv_A, solve_middle)
    if (v_is_matrix) out else as.vector(out)
  }

  if (method == "satterthwaite") {
    error <- pseudo_response - offset - matrixVectorMultiply(A, alpha)
    r <- Vinv_apply(error)
    Gj_r <- Gj_apply(r)
    u <- max(0, sum(r * Gj_r) / 2)
    q <- ncol(Bj)
    Cj <- matrix(0, q, q)
    for (i in 1:q) {
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
    pv <- pchisq(u / kappa, df = nu, lower.tail = FALSE)
  } else {
    error <- pseudo_response - offset - matrixVectorMultiply(A, alpha)
    r <- P_apply(error)
    Gj_r <- Gj_apply(r)
    u <- max(0, sum(r * Gj_r))
    q <- ncol(Bj)
    eig_theta <- matrixsqrt(Thetaj)
    Theta_sqrt <- eig_theta$w
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
        if (davies_res$Qq <= 0 || davies_res$Qq > 1.0) {
          compute_liu_pvalue(u, lambda)
        } else {
          davies_res$Qq
        }
      },
      stop("method must be one of 'satterthwaite', 'davies', 'liu', 'hbe', 'wood', 'saddlepoint'", call. = FALSE)
    )
  }

  data.table(
    smooth.term = smooth_terms[[test.component]]$label,
    smooth.pvalue = pv,
    method = method
  )
}

native_re_is_smooth <- function(s) {
  any(c("random.effect", "fs.interaction") %in% class(s))
}

native_re_design_data <- function(fit) {
  data <- fit$model
  if (!inherits(fit, "bam") || is.null(fit$dinfo)) return(data)

  discrete_mf <- utils::getFromNamespace("discrete.mf", "mgcv")
  dk <- discrete_mf(
    fit$dinfo$gp, mf = data, names.pmf = fit$dinfo$pmf.names,
    full = TRUE
  )
  for (name in names(dk$mf)) {
    lookup_column <- dk$ks[name, 1]
    values <- dk$mf[[name]]
    if (is.matrix(values)) {
      data[[name]] <- values[dk$k[, lookup_column], , drop = FALSE]
    } else {
      data[[name]] <- values[dk$k[, lookup_column]]
    }
  }
  data
}

native_re_dense_smooth_design <- function(s, data) {
  B <- mgcv::PredictMat(s, data)
  n_coef <- s$last.para - s$first.para + 1L
  if (ncol(B) != n_coef) {
    stop("Smooth design has inconsistent coefficient dimensions.")
  }
  as.matrix(B)
}

native_re_local_smooth_design <- function(s, data) {
  n_coef <- s$last.para - s$first.para + 1L

  if (inherits(s, "random.effect")) {
    mf <- stats::model.frame(s$form, data, na.action = stats::na.pass)
    factor_index <- which(vapply(mf, is.factor, logical(1)))
    if (length(factor_index) != 1L) {
      stop("Native re smooths must have exactly one factor grouping variable.")
    }
    group_factor <- mf[[factor_index]]
    group_levels <- levels(group_factor)
    if (!length(group_levels) || n_coef != length(group_levels)) {
      stop("Native re smooths must have one coefficient per grouping level.")
    }
    X <- matrix(1, nrow(mf), 1L)
    value_index <- setdiff(seq_along(mf), factor_index)
    for (j in value_index) {
      value <- mf[[j]]
      if ((!is.numeric(value) && !is.logical(value)) || is.matrix(value)) {
        stop("Native re covariates must be scalar numeric variables.")
      }
      X[, 1L] <- X[, 1L] * as.numeric(value)
    }
    group <- as.character(group_factor)
  } else if (inherits(s, "fs.interaction")) {
    group_factor <- data[[s$fterm]]
    group_levels <- as.character(s$flev)
    group <- as.character(group_factor)
    if (anyNA(match(group, group_levels))) {
      stop("Native fs smooth has an unrecognised factor level.")
    }
    data0 <- data
    data0[[s$fterm]] <- NULL
    base_smooth <- s
    class(base_smooth) <- base_smooth$base$bs
    base_smooth$rank <- base_smooth$base$rank
    base_smooth$null.space.dim <- base_smooth$base$null.space.dim
    base_smooth$bs.dim <- base_smooth$base$bs.dim
    base_smooth$term <- base_smooth$base$term
    X <- mgcv::Predict.matrix(base_smooth, data0) %*% s$P
    if (ncol(X) * length(group_levels) != n_coef) {
      stop("Native fs smooth has inconsistent local coefficient dimensions.")
    }
  } else {
    stop("Uid-local design is only available for mgcv re and fs smooths.")
  }

  if (anyNA(group) || anyNA(X) || any(!is.finite(X))) {
    stop("Native random-effect design contains missing or non-finite values.")
  }
  list(X = as.matrix(X), group = group, levels = group_levels)
}

native_re_local_penalty <- function(S, n_group, block_size,
                                    tolerance = 1e-10) {
  expected <- n_group * block_size
  if (nrow(S) != expected || ncol(S) != expected) {
    stop("Random-effect penalty has inconsistent uid-block dimensions.")
  }
  reference <- S[seq_len(block_size), seq_len(block_size), drop = FALSE]
  scale <- max(1, max(abs(reference)))
  for (g in seq_len(n_group)) {
    index <- (g - 1L) * block_size + seq_len(block_size)
    if (max(abs(S[index, index, drop = FALSE] - reference)) >
        tolerance * scale) {
      stop("Random-effect penalty blocks are not identical across uid.")
    }
    if (n_group > 1L &&
        max(abs(S[index, -index, drop = FALSE])) > tolerance * scale) {
      stop("Random-effect penalty is not block diagonal by uid.")
    }
  }
  reference
}

native_re_bind_dense <- function(blocks, n) {
  if (!length(blocks)) {
    return(list(indices = integer(0), X = matrix(0, n, 0L)))
  }
  indices <- unlist(lapply(blocks, `[[`, "indices"), use.names = FALSE)
  X <- do.call(cbind, lapply(blocks, `[[`, "X"))
  order_index <- order(indices)
  list(indices = indices[order_index], X = X[, order_index, drop = FALSE])
}

native_re_bind_local <- function(blocks, n) {
  if (!length(blocks)) {
    stop("The native GAMM backend requires a penalized re or fs nuisance smooth.")
  }
  block_order <- order(vapply(blocks, function(block) {
    min(block$indices)
  }, numeric(1)))
  blocks <- blocks[block_order]
  group_levels <- blocks[[1]]$levels
  group <- blocks[[1]]$group
  for (block in blocks) {
    if (!setequal(block$levels, group_levels) ||
        !identical(block$group, group)) {
      stop("All re and fs nuisance smooths must use the same uid grouping.")
    }
  }
  group_index <- match(group, group_levels)
  if (length(group_index) != n || anyNA(group_index)) {
    stop("Uid grouping is inconsistent with the model rows.")
  }
  list(
    indices = unlist(lapply(blocks, `[[`, "indices"), use.names = FALSE),
    X = do.call(cbind, lapply(blocks, `[[`, "X")),
    group = as.integer(group_index),
    levels = group_levels,
    S = as.matrix(Matrix::bdiag(lapply(blocks, `[[`, "S")))
  )
}

native_re_block_diag <- function(mats) {
  if (!length(mats)) return(matrix(0, 0L, 0L))
  as.matrix(Matrix::bdiag(mats))
}

native_re_diagonal_profile <- function(B_dense, B_local, group, V_phi,
                                       dense_S, local_S, n_threads) {
  W <- 1 / as.numeric(V_phi)
  cpp_profile <- native_re_uid_profile_create_cpp(
    dense = B_dense, local = B_local, group = group, weight = W,
    dense_penalty = dense_S, local_penalty = local_S,
    n_threads = n_threads
  )
  C_schur <- as.matrix(cpp_profile$schur)
  C_schur_inv <- if (nrow(C_schur)) {
    matrixGeneralizedInverse((C_schur + t(C_schur)) / 2)
  } else {
    matrix(0, 0L, 0L)
  }
  pointer <- cpp_profile$pointer

  Vinv_apply <- function(v) {
    v_is_matrix <- is.matrix(v)
    v_mat <- if (v_is_matrix) v else matrix(v, ncol = 1L)
    out <- native_re_uid_profile_apply_cpp(pointer, v_mat, C_schur_inv)
    if (v_is_matrix) out else as.vector(out)
  }

  list(
    apply = Vinv_apply,
    n_group = cpp_profile$n_group,
    block_size = cpp_profile$block_size
  )
}

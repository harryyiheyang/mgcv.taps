gammfast_weight_cache <- function(X, B, id, response, weight, nthreads) {
  sw <- sqrt(weight)
  gammfast_gaussian_cache(
    cbind(X * sw, response * sw),
    B * sw, id, n_threads = nthreads
  )
}

gammfast_working_cache <- function(X, B, id, work, offset, nthreads) {
  gammfast_weight_cache(
    X, B, id, work$z - offset, work$w, nthreads
  )
}

gammfast_fixed_mode_cached <- function(G0, X, B, id, cache, G, scale,
                                       sp, nthreads) {
  cp <- gammfast_gaussian_crossprod_cached(
    cache$AtA, cache$BtB, cache$BtA, G, n_threads = nthreads
  )
  H <- (cp$crossprod + t(cp$crossprod)) / 2
  p <- ncol(X)
  XtX <- H[seq_len(p), seq_len(p), drop = FALSE]
  Xtz <- H[seq_len(p), p + 1L]
  penalty <- gammfast_penalty_matrix(G0, sp)
  precision <- (XtX + penalty + t(XtX + penalty)) / 2
  R <- chol(precision)
  beta <- backsolve(R, forwardsolve(t(R), Xtz))
  precision_inverse <- chol2inv(R)
  u <- sqrt(scale) * gammfast_gaussian_blup_cached(
    cache$BtB, cache$BtA, G, beta, scale
  )
  list(
    beta = beta, u = u,
    Vp = scale * precision_inverse,
    Ve = scale * precision_inverse %*% XtX %*% precision_inverse,
    logdet = cp$logdet
  )
}

gammfast_cached_moment <- function(cache, G, scale, covariance_group,
                                   group_index, ng, correction = NULL,
                                   mm = NULL) {
  if (is.null(mm)) {
    mm <- gammfast_gaussian_projected_cached(
      cache$AtA, cache$BtB, cache$BtA, G, scale,
      covariance_group, fisher = FALSE
    )
  }
  moment <- mm$moment_sum
  if (!is.null(correction)) moment <- moment - correction / scale
  G_new <- gammfast_project_covariance(moment / ng, group_index)
  list(G = G_new, mm = mm)
}

gammfast_prepare_influence_caches <- function(
    X, B, id, work, offset, t_correction, nthreads) {
  working <- gammfast_working_cache(X, B, id, work, offset, nthreads)
  same_curvature <- is.null(t_correction) ||
    !any(t_correction$weight != work$w)
  determinant <- if (same_curvature) {
    working
  } else {
    gammfast_weight_cache(
      X, B, id, work$z - offset, t_correction$weight, nthreads
    )
  }
  list(
    working = working, determinant = determinant,
    same.curvature = same_curvature
  )
}

gammfast_cached_influence <- function(
    X, B, id, G, scale, covariance_group, t_correction, nthreads,
    caches) {
  cache <- caches$working
  mm <- gammfast_gaussian_projected_cached(
    cache$AtA, cache$BtB, cache$BtA, G, scale,
    covariance_group, fisher = FALSE
  )
  if (is.null(t_correction)) {
    return(list(mm = mm, correction = NULL))
  }
  determinant_cache <- caches$determinant
  determinant_mean_covariance <- mm$mean_covariance
  if (!caches$same.curvature) {
    determinant_crossprod <- gammfast_gaussian_crossprod_cached(
      determinant_cache$AtA, determinant_cache$BtB,
      determinant_cache$BtA, G, n_threads = nthreads
    )
    p <- ncol(X)
    determinant_precision <- determinant_crossprod$crossprod[
      seq_len(p), seq_len(p), drop = FALSE
    ] / scale
    determinant_mean_covariance <- chol2inv(chol(
      (determinant_precision + t(determinant_precision)) / 2
    ))
  }
  influence <- gammfast_laplace_influence_cached(
    X = X, B = B, id = id, G = G,
    working_BtB = cache$BtB, working_BtA = cache$BtA,
    working_mean_covariance = mm$mean_covariance,
    determinant_BtB = determinant_cache$BtB,
    determinant_BtA = determinant_cache$BtA,
    determinant_mean_covariance = determinant_mean_covariance,
    determinant_derivative = t_correction$derivative,
    u = sqrt(scale) * mm$u, scale = scale,
    same_curvature = caches$same.curvature,
    n_threads = nthreads
  )
  list(mm = mm, correction = influence$influence_sum)
}

gammfast_inner_limit <- function(inner, inner_max) {
  is.finite(inner_max) && inner >= inner_max
}

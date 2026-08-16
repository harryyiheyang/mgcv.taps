gammfast_non_gaussian_inner_G <- function(
    X, B, id, offset, G, work, t_correction, mean_penalty,
    scale, group_index, ng, inner_tol, inner_max,
    nthreads) {
  caches <- gammfast_prepare_influence_caches(
    X, B, id, work, offset, t_correction, nthreads
  )
  cache <- caches$working
  state <- gammfast_cached_influence(
    X, B, id, G, scale, mean_penalty,
    t_correction, nthreads,
    caches
  )
  inner <- 0L
  local_change <- Inf
  repeat {
    mapped <- gammfast_cached_moment(
      cache, G, scale, group_index, ng,
      mean_penalty = mean_penalty,
      correction = state$correction, mm = state$mm
    )
    state$mm <- NULL
    local_change <- norm(mapped$G - G, "F") / (1 + norm(G, "F"))
    G <- mapped$G
    inner <- inner + 1L
    if (local_change < inner_tol || gammfast_inner_limit(inner, inner_max)) break
  }
  list(
    G = G, cache = cache, inner = inner,
    local_change = local_change,
    inner_limit = gammfast_inner_limit(inner, inner_max) &&
      local_change >= inner_tol
  )
}

gammfast_non_gaussian_profile <- function(
    family, y, prior_weights, offset, G0, X, B, id, G, eta, sp, scale,
    group_index, ng, inner_tol, inner_max, nthreads,
    pirls_control) {
  mu <- family$linkinv(eta)
  deviance_old <- sum(family$dev.resids(y, mu, prior_weights))
  if (!is.finite(deviance_old)) {
    stop("The initial profiled PIRLS deviance is not finite.")
  }
  deviance_change <- Inf
  eta_change <- Inf
  theta_change <- Inf
  inner_total <- 0L
  inner_limit <- FALSE
  mode <- NULL
  pirls_converged <- FALSE
  damp_eta <- identical(tolower(family$family[1]), "inverse.gaussian")

  for (pirls in seq_len(pirls_control$maxit)) {
    work <- gammfast_working(family, y, eta, prior_weights)
    t_correction <- gammfast_t_correction(
      family, y, eta, prior_weights, work$w
    )
    mean_penalty <- gammfast_penalty_matrix(G0, sp, scale = scale)
    local <- gammfast_non_gaussian_inner_G(
      X = X, B = B, id = id, offset = offset, G = G,
      work = work, t_correction = t_correction,
      mean_penalty = mean_penalty, scale = scale,
      group_index = group_index,
      ng = ng, inner_tol = inner_tol, inner_max = inner_max,
      nthreads = nthreads
    )
    G <- local$G
    mode <- gammfast_conditional_mode_cached(
      G0, local$cache, G, scale, sp, nthreads
    )
    eta_new <- offset + drop(X %*% mode$beta) +
      rowSums(B * mode$u[id, , drop = FALSE])
    if (damp_eta) eta_new <- (eta + eta_new) / 2
    mu_new <- family$linkinv(eta_new)
    deviance <- sum(family$dev.resids(y, mu_new, prior_weights))
    if (!is.finite(deviance)) {
      stop("The profiled PIRLS deviance is not finite.")
    }
    deviance_change <- abs(deviance - deviance_old) /
      (0.1 + abs(deviance))
    eta_change <- max(abs(eta_new - eta)) / (1 + max(abs(eta)))
    theta_update <- gammfast_update_theta(
      family, y, mu_new, prior_weights, scale
    )
    updated_scale <- theta_update$scale
    scale_change <- abs(log(updated_scale / scale))
    if (!isTRUE(all.equal(updated_scale, scale))) {
      G <- G * scale / updated_scale
      scale <- updated_scale
    }
    theta_change <- max(theta_update$change, scale_change)
    eta <- eta_new
    deviance_old <- deviance
    inner_total <- inner_total + local$inner
    inner_limit <- inner_limit || local$inner_limit
    pirls_converged <- pirls > 1L &&
      deviance_change < pirls_control$epsilon &&
      eta_change < pirls_control$epsilon &&
      theta_change < pirls_control$epsilon
    if (pirls_converged) break
  }
  list(
    G = G, scale = scale, beta = mode$beta, u = mode$u, eta = eta,
    mode = mode, pirls = pirls, pirls_converged = pirls_converged,
    deviance = deviance_old, deviance_change = deviance_change,
    eta_change = eta_change, theta_change = theta_change,
    inner_total = inner_total, inner_last = local$inner,
    fixedpoint_G = local$local_change, inner_limit = inner_limit
  )
}

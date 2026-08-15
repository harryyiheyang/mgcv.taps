gammfast_non_gaussian <- function(formula, global_formula, shell, G0, Sl,
                                  X, y, prior_weights, offset, B, id,
                                  id_factor, random_structure, inner.max,
                                  inner.tol, pirls.max, pirls.tol,
                                  covariance_group, nthreads, control,
                                  verbose, call) {
  family <- shell$family
  ng <- nlevels(id_factor)
  eta <- as.numeric(shell$linear.predictors)
  initial_variance <- stats::var(eta)
  if (!is.finite(initial_variance) || initial_variance <= 0) {
    initial_variance <- 1
  }
  G <- gammfast_initial_covariance(
    random_structure$group.index, initial_variance * 0.15
  )
  family_scale <- if (is.numeric(shell$sig2) && length(shell$sig2) == 1L) {
    shell$sig2
  } else {
    1
  }
  estimate_working_phi <- gammfast_estimate_working_scale(family)
  pirls_control <- list(epsilon = pirls.tol, maxit = pirls.max)
  work <- gammfast_working(family, y, eta, prior_weights)
  initial <- gammfast_general_step(
    G0, Sl, X, work$z - offset, B, id, G,
    sqrt(work$w), NULL, nthreads,
    phi = family_scale, estimate_phi = estimate_working_phi
  )
  rho_start <- pmin(initial$rho, 25)
  sp_current <- initial$sp
  family_scale <- initial$phi
  sp_boundary <- isTRUE(initial$sp_boundary) || any(log(sp_current) >= 25)
  initial_cache <- gammfast_working_cache(
    X, B, id, work, offset, nthreads
  )
  initial_mode <- gammfast_fixed_mode_cached(
    G0, X, B, id, initial_cache, G, family_scale,
    sp_current, nthreads
  )
  beta_state <- initial_mode$beta
  u_state <- initial_mode$u
  eta <- offset + drop(X %*% beta_state) +
    rowSums(B * u_state[id, , drop = FALSE])
  objective_old <- NULL
  trace <- data.frame()
  converged <- FALSE
  variance_evaluations <- 0L
  start_time <- proc.time()[[3]]

  for (outer in seq_len(control$max.outer)) {
    profile <- gammfast_non_gaussian_profile(
      family = family, y = y, prior_weights = prior_weights,
      offset = offset, G0 = G0, X = X, B = B, id = id, G = G,
      eta = eta, sp = sp_current, scale = family_scale,
      covariance_group = covariance_group,
      group_index = random_structure$group.index, ng = ng,
      inner_tol = inner.tol, inner_max = inner.max,
      nthreads = nthreads, pirls_control = pirls_control
    )
    G <- profile$G
    family_scale <- profile$scale
    eta <- profile$eta
    variance_evaluations <- variance_evaluations + profile$inner_total

    work <- gammfast_working(family, y, eta, prior_weights)
    step <- gammfast_general_step(
      G0, Sl, X, work$z - offset, B, id, G,
      sqrt(work$w), rho_start, nthreads,
      optimize_sp = !sp_boundary,
      sp_fixed = if (sp_boundary) sp_current else NULL,
      phi = family_scale,
      estimate_phi = estimate_working_phi && !sp_boundary
    )
    sp_new <- step$sp
    rho_start <- pmin(step$rho, 25)
    if (step$phi != family_scale) {
      G <- G * family_scale / step$phi
      family_scale <- step$phi
    }
    sp_boundary <- sp_boundary || isTRUE(step$sp_boundary) ||
      any(log(sp_new) >= 25)
    if (!is.null(step$fr$outer.info$conv) &&
        grepl("divergence", step$fr$outer.info$conv, ignore.case = TRUE)) {
      sp_boundary <- TRUE
    }
    dsp <- max(abs(log(sp_new / sp_current)))
    objective <- step$fr$reml + 0.5 * step$logdet
    dobjective <- if (is.null(objective_old)) {
      Inf
    } else {
      abs(objective - objective_old) / (1 + abs(objective_old))
    }
    trace <- rbind(trace, data.frame(
      outer = outer, objective = objective, dobjective = dobjective,
      dsp = dsp, sig2 = family_scale,
      pirls = profile$pirls,
      pirls_converged = profile$pirls_converged,
      pirls_deviance = profile$deviance,
      pirls_deviance_change = profile$deviance_change,
      eta_change = profile$eta_change,
      theta_change = profile$theta_change,
      fixedpoint_G = profile$fixedpoint_G,
      variance_evaluations = profile$inner_total,
      variance_evaluations_total = variance_evaluations,
      inner_last = profile$inner_last,
      inner_limit = profile$inner_limit
    ))
    if (verbose) {
      cat(
        "outer", outer, "working fREML/Laplace criterion",
        format(objective, digits = 8), "relative change",
        format(dobjective, digits = 3), "PIRLS", profile$pirls,
        "G updates", profile$inner_total, "\n"
      )
    }
    converged <- outer > 2L && profile$pirls_converged &&
      dobjective < control$objective.tol
    sp_current <- sp_new
    objective_old <- objective
    if (converged) break
  }

  final_profile <- gammfast_non_gaussian_profile(
    family = family, y = y, prior_weights = prior_weights,
    offset = offset, G0 = G0, X = X, B = B, id = id, G = G,
    eta = eta, sp = sp_current, scale = family_scale,
    covariance_group = covariance_group,
    group_index = random_structure$group.index, ng = ng,
    inner_tol = inner.tol, inner_max = inner.max,
    nthreads = nthreads, pirls_control = pirls_control
  )
  G <- final_profile$G
  family_scale <- final_profile$scale
  beta <- final_profile$beta
  names(beta) <- names(shell$coefficients)
  u <- final_profile$u
  rownames(u) <- levels(id_factor)
  colnames(u) <- random_structure$column.names
  eta <- final_profile$eta
  eta_global <- offset + drop(X %*% beta)
  eta_random <- eta - eta_global
  mu <- family$linkinv(eta)
  work_final <- gammfast_working(family, y, eta, prior_weights)

  shell$coefficients <- beta
  shell$sp <- sp_current
  shell$sig2 <- family_scale
  shell$scale <- family_scale
  shell$weights <- work_final$w
  shell$linear.predictors <- eta_global
  shell$fitted.values <- family$linkinv(eta_global)
  shell$family <- family
  shell$Vp <- final_profile$mode$Vp
  shell$Ve <- final_profile$mode$Ve

  G_actual <- family_scale * G
  random_structure$covariance <- lapply(
    random_structure$group.index,
    function(jj) G_actual[jj, jj, drop = FALSE]
  )

  fit <- list(
    coefficients = beta, Vp = final_profile$mode$Vp,
    Ve = final_profile$mode$Ve, random.effects = u, G = G_actual,
    dispersion.method = if (estimate_working_phi) {
      "mgcv-fREML"
    } else if (inherits(family, "extended.family") &&
               is.numeric(family$scale) && length(family$scale) == 1L &&
               family$scale < 0) {
      "mgcv-estimate.theta"
    } else {
      "family-fixed"
    },
    family.parameter.method = if (inherits(family, "extended.family") &&
                                  !is.null(family$n.theta) &&
                                  family$n.theta > 0L) {
      "mgcv-estimate.theta"
    } else {
      "family-fixed"
    },
    covariance.method = "cached-ordinary-X-Laplace-influence-fixedpoint",
    optimization.method = "profiled-PIRLS-G-blockwise-fREML",
    sp = sp_current, fitted.values = mu, linear.predictors = eta,
    global.fitted = eta_global, random.fitted = eta_random,
    residuals = y - mu, y = y, weights = work_final$w,
    prior.weights = prior_weights, offset = offset, sig2 = family_scale,
    objective = objective_old,
    objective.type = "working-fREML-Laplace-block-criterion",
    converged = converged && final_profile$pirls_converged,
    outer = nrow(trace), inner.max = inner.max, inner.tol = inner.tol,
    pirls.tol = pirls_control$epsilon, pirls.max = pirls_control$maxit,
    variance.evaluations = variance_evaluations + final_profile$inner_total,
    elapsed = proc.time()[[3]] - start_time, trace = trace,
    formula = formula, global.formula = global_formula, global = shell,
    family = family, random = random_structure,
    call = call, control = control
  )
  class(fit) <- "gammfast"
  if (!fit$converged) {
    warning("gammfast reached control$max.outer before blockwise convergence.")
  }
  fit
}

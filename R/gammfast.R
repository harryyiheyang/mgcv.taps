utils::globalVariables(".gammfast_prior_weights")

#' Fast additive model with subject-specific trajectories
#'
#' Fits an additive model with independent subject-level random groups. A
#' mandatory `s(id, bs = "re")` supplies the random intercept; additional
#' `s(x, id, bs = "re")` terms supply independent random slopes and custom
#' `fs(x, id, k)` markers supply independent cosine random trajectories. All
#' groups share one subject ID and are marginalized by subject, while global
#' smooth terms use the `mgcv` fREML machinery.
#'
#' @param formula An mgcv-style formula containing `s(id, bs = "re")`, zero or
#'   more independent random-slope `s(..., id, bs = "re")` terms, and zero or
#'   more custom `fs(x, id, k)` markers. Standard mgcv
#'   `s(..., bs = "fs")` terms are not supported by `gammfast()`.
#' @param data A data frame.
#' @param family An mgcv family object. Supported families are Gaussian,
#'   binomial and quasi-binomial, Poisson and quasi-Poisson, Gamma, inverse
#'   Gaussian, quasi-likelihood, negative binomial, scaled t, beta regression,
#'   and Tweedie. Standard links supported by mgcv are accepted.
#' @param weights Optional prior weights.
#' @param inner.max Optional safety limit for shared-UID covariance updates per
#'   fixed PIRLS working model. The default is 300 iterations.
#' @param inner.tol Relative covariance fixed-point tolerance.
#' @param pirls.max Maximum PIRLS iterations for each fixed global
#'   smoothing-parameter vector.
#' @param pirls.tol Relative PIRLS deviance tolerance.
#' @param influence.update Non-Gaussian Laplace-influence update strategy.
#'   `"frozen"` preserves the original implementation by evaluating the
#'   influence once per PIRLS working model. `"refreshed"` reevaluates it at
#'   every subject-covariance update while reusing the working caches.
#' @param nthreads Number of OpenMP threads used for subject-level operations.
#' @param discrete Must be `FALSE`. The current solver does not implement the
#'   mgcv/bam discrete model-matrix representation.
#' @param control A list with `objective.tol` and `max.outer`. Outer iterations
#'   stop on relative convergence of the branch-specific reported objective.
#'   For non-Gaussian fits this is explicitly a blockwise working
#'   fREML/Laplace criterion. Each covariance iteration stops at `inner.tol`
#'   or reaches `inner.max`. The former `fixedpoint.tol` entry is accepted as
#'   an alias for `inner.tol`.
#' @param verbose Whether to print outer-iteration diagnostics.
#'
#' @return An object of class `"gammfast"` containing the global fit, the full
#'   random-trajectory covariance `G` on the linear-predictor scale, subject
#'   BLUPs, fitted values, and convergence diagnostics.
#'
#' @details `gammfast()` deliberately does not use mgcv's discrete fitting
#'   path. mgcv constructs the global design and penalties (`fit = FALSE`) and
#'   its fREML machinery updates global smoothing parameters and dispersion.
#'   `gammfast()` replaces only the working cross-products by their current
#'   shared-UID marginal forms and estimates the subject covariance without
#'   expanding an `n` by `n_subject * q` random-effect design. For models for
#'   which the expected working-curvature derivative is available, the
#'   covariance update combines the conditional BLUP moment, the ordinary
#'   full-`X` mean-estimation Schur correction, and the Laplace determinant
#'   influence correction. It does not eigendecompose or reparameterize the
#'   global mean penalty for covariance updates. Gamma and inverse-Gaussian
#'   models retain mgcv Fisher weights for P-IRLS and use Fisher curvature
#'   throughout the determinant leverage. Extended-family nuisance
#'   parameters are updated by mgcv's family-parameter optimizer. Formula
#'   offsets are retained in fitting and prediction. For non-Gaussian models,
#'   each current global smoothing-parameter vector is held fixed while PIRLS
#'   converges; within every PIRLS working model, `W` and `z` are held fixed
#'   while the subject covariance converges. Only after that profile step does
#'   the blockwise mgcv fREML smoothing update run. This is a blockwise profile
#'   approximation, not Wood's exact nested-Newton LAML derivative with the
#'   covariance response Schur-complemented from the outer Hessian.
#'
#' Each custom `fs(x, id, k)` covariate is affinely scaled over its fitted
#' range to `[0, 1]`, reflecting the working assumption of approximately
#' uniform functional-domain coverage. Its basis is
#' `sqrt(2) * cos(h * pi * x)`, `h = 1, ..., k`, with no intercept column.
#' Every fs block has an unrestricted positive-definite covariance. Random
#' groups are mutually independent, and the explicit random-intercept block is
#' independent of every slope and fs block.
#'
#' @examples
#' \dontrun{
#' library(jmcm)
#' data(aids)
#' aids$id <- factor(aids$id)
#' fit <- gammfast(
#'   sqrt(cd4) ~ s(time, k = 8) + age + s(id, bs = "re") +
#'     fs(time, id, k = 3),
#'   data = aids
#' )
#' }
#'
#' @export
gammfast <- function(formula, data, family = stats::gaussian(), weights = NULL,
                     inner.max = 300L, nthreads = 1L,
                     discrete = FALSE, control = list(), verbose = FALSE,
                     inner.tol = 1e-5, pirls.max = 100L,
                     pirls.tol = 1e-6,
                     influence.update = c("frozen", "refreshed")) {
  if (!inherits(formula, "formula")) stop("formula must be a model formula.")
  if (!is.data.frame(data)) stop("data must be a data frame.")
  influence.update <- match.arg(influence.update)
  family_info <- gammfast_validate_family(family)
  family <- family_info$family
  if (is.null(weights)) weights <- rep(1, nrow(data))
  if (!is.numeric(weights) || length(weights) != nrow(data) ||
      any(!is.finite(weights)) || any(weights < 0)) {
    stop("weights must be a finite non-negative numeric vector with one value per row.")
  }
  if (!any(weights > 0)) stop("At least one weight must be positive.")
  if (!is.list(control)) stop("control must be a list.")
  if (!is.null(control$fixedpoint.tol)) {
    if (!missing(inner.tol)) {
      stop("Specify inner.tol or control$fixedpoint.tol, not both.")
    }
    inner.tol <- control$fixedpoint.tol
    control$fixedpoint.tol <- NULL
  }
  weight_name <- ".gammfast_prior_weights"
  if (weight_name %in% names(data)) {
    stop("data contains the reserved column '.gammfast_prior_weights'.")
  }
  data[[weight_name]] <- weights
  if (!is.numeric(inner.max) || length(inner.max) != 1L ||
      !is.finite(inner.max) || inner.max < 1 ||
      inner.max != floor(inner.max)) {
    stop("inner.max must be a positive finite integer.")
  }
  if (!is.numeric(inner.tol) || length(inner.tol) != 1L ||
      !is.finite(inner.tol) || inner.tol <= 0) {
    stop("inner.tol must be one positive finite value.")
  }
  if (!is.numeric(pirls.max) || length(pirls.max) != 1L ||
      !is.finite(pirls.max) || pirls.max < 1 ||
      pirls.max != floor(pirls.max)) {
    stop("pirls.max must be a positive integer.")
  }
  if (!is.numeric(pirls.tol) || length(pirls.tol) != 1L ||
      !is.finite(pirls.tol) || pirls.tol <= 0) {
    stop("pirls.tol must be one positive finite value.")
  }
  if (!is.numeric(nthreads) || length(nthreads) != 1L ||
      !is.finite(nthreads) || nthreads < 1 || nthreads != floor(nthreads)) {
    stop("nthreads must be a positive integer.")
  }
  inner.max <- as.integer(inner.max)
  pirls.max <- as.integer(pirls.max)
  nthreads <- as.integer(nthreads)
  if (length(discrete) != 1L || is.na(discrete) || !is.logical(discrete)) {
    stop("discrete must be a single TRUE or FALSE value.")
  }
  if (isTRUE(discrete)) {
    stop("gammfast does not support discrete = TRUE; use discrete = FALSE.")
  }

  defaults <- list(
    objective.tol = 1e-7,
    max.outer = 5000L
  )
  unknown <- setdiff(names(control), names(defaults))
  if (length(unknown)) {
    stop("Unknown control entries: ", paste(unknown, collapse = ", "), ".")
  }
  defaults[names(control)] <- control
  control <- defaults
  if (!is.finite(control$objective.tol) || control$objective.tol <= 0) {
    stop("control$objective.tol must be positive.")
  }
  if (!is.finite(control$max.outer) || control$max.outer < 1) {
    stop("control$max.outer must be a positive integer.")
  }
  control$max.outer <- as.integer(control$max.outer)

  parsed <- gammfast_parse_formula(formula, data)
  id_name <- parsed$id
  global_formula <- parsed$formula
  if (!id_name %in% names(data)) stop("The fs id variable is absent from data.")
  if (anyNA(data[[id_name]])) {
    stop("gammfast does not allow missing subject IDs.")
  }

  G0 <- mgcv::gam(
    global_formula, data = data, family = family,
    weights = .gammfast_prior_weights,
    method = "REML", fit = FALSE, na.action = stats::na.fail
  )
  family <- G0$family
  shell <- gammfast_prefit_shell(G0, family)
  X <- as.matrix(G0$X)
  if (!length(G0$S)) {
    stop("gammfast currently requires at least one global smooth term.")
  }
  gammfast_validate_smooths(G0)
  y <- as.numeric(G0$y)
  if (length(y) != nrow(data)) stop("The response length does not match data.")
  prior_weights <- as.numeric(shell$prior.weights)
  offset <- stats::model.offset(shell$model)
  if (is.null(offset)) offset <- numeric(length(y))
  if (length(offset) != length(y) || any(!is.finite(offset))) {
    stop("The model offset must be finite with one value per observation.")
  }
  id_factor <- factor(data[[id_name]])
  id <- as.integer(id_factor)
  ng <- nlevels(id_factor)
  random_structure <- gammfast_random_structure(
    parsed$groups, data, id_name, levels(id_factor)
  )
  B <- random_structure$B
  covariance_group <- integer(ncol(B))
  for (j in seq_along(random_structure$group.index)) {
    covariance_group[random_structure$group.index[[j]]] <- j
  }
  random_structure$id <- id_name
  random_structure$id.levels <- levels(id_factor)
  random_structure$id.index <- id

  Sl.setup <- utils::getFromNamespace("Sl.setup", "mgcv")
  Sl.Xprep <- utils::getFromNamespace("Sl.Xprep", "mgcv")
  Sl.postproc <- utils::getFromNamespace("Sl.postproc", "mgcv")
  initial.sp <- utils::getFromNamespace("initial.sp", "mgcv")
  fast.REML.fit <- utils::getFromNamespace("fast.REML.fit", "mgcv")
  Sl <- Sl.setup(G0)

  if (!family_info$gaussian) {
    return(gammfast_non_gaussian(
      formula = formula, global_formula = global_formula,
      shell = shell, G0 = G0, Sl = Sl, X = X, y = y,
      prior_weights = prior_weights, offset = offset, B = B, id = id,
      id_factor = id_factor, random_structure = random_structure,
      inner.max = inner.max, inner.tol = inner.tol,
      pirls.max = pirls.max, pirls.tol = pirls.tol,
      covariance_group = covariance_group,
      nthreads = nthreads, control = control, verbose = verbose,
      call = match.call(), influence.update = influence.update
    ))
  }
  if (any(prior_weights != 1)) {
    stop("Gaussian gammfast currently requires unit prior weights.")
  }

  y_work <- y - offset
  sigma2 <- stats::var(y_work) * 0.5
  G <- gammfast_initial_covariance(random_structure$group.index, 0.3)
  A <- cbind(X, y_work)
  crossprod_cache <- gammfast_gaussian_cache(
    A, B, id, n_threads = nthreads
  )
  rho_start <- NULL
  sp_old <- NULL
  sigma2_old <- NULL
  objective_old <- NULL
  trace <- data.frame()
  converged <- FALSE
  start_time <- proc.time()[[3]]

  for (outer in seq_len(control$max.outer)) {
    cp <- gammfast_gaussian_crossprod_cached(
      crossprod_cache$AtA, crossprod_cache$BtB,
      crossprod_cache$BtA, G, n_threads = nthreads
    )
    H <- (cp$crossprod + t(cp$crossprod)) / 2
    p <- ncol(X)
    XtX <- H[seq_len(p), seq_len(p), drop = FALSE]
    Xty <- H[seq_len(p), p + 1L]
    yty <- H[p + 1L, p + 1L]
    R1 <- chol(XtX)
    f <- forwardsolve(t(R1), Xty)
    rss_extra <- max(0, yty - sum(f^2))
    um <- Sl.Xprep(Sl, R1, nt = nthreads)
    if (is.null(rho_start)) {
      rho_start <- log(initial.sp(R1, G0$S, G0$off))
    }
    fr <- fast.REML.fit(
      um$Sl, um$X, f, rho = rho_start,
      L = G0$L, rho.0 = G0$lsp0,
      log.phi = log(sigma2), phi.fixed = FALSE,
      rss.extra = rss_extra, nobs = length(y), Mp = um$Mp,
      nt = nthreads, gamma = 1
    )
    pars <- gammfast_reml_parameters(fr, G0, estimate_phi = TRUE)
    sigma2 <- pars$phi
    pp <- Sl.postproc(
      Sl, fr, um$undrop, R1, cov = TRUE,
      scale = sigma2, L = G0$L, nt = nthreads
    )
    beta <- pp$beta
    rho_start <- pars$rho
    eta_global_work <- drop(X %*% beta)

    inner <- 0L
    fixedpoint_G <- Inf
    repeat {
      inner <- inner + 1L
      mm <- gammfast_gaussian_projected_cached(
        crossprod_cache$AtA, crossprod_cache$BtB,
        crossprod_cache$BtA, G, sigma2, covariance_group,
        fisher = FALSE
      )
      G_new <- gammfast_project_covariance(
        mm$moment_sum / ng, random_structure$group.index
      )
      fixedpoint_G <- norm(G_new - G, "F") / (1 + norm(G, "F"))
      G <- G_new
      if (fixedpoint_G < inner.tol ||
          gammfast_inner_limit(inner, inner.max)) break
    }
    dphi <- if (is.null(sigma2_old)) Inf else abs(log(sigma2 / sigma2_old))
    sp_now <- pars$sp
    dsp <- if (is.null(sp_old)) Inf else max(abs(log(sp_now / sp_old)))
    objective <- fr$reml + 0.5 * cp$logdet
    dobjective <- if (is.null(objective_old)) {
      Inf
    } else {
      abs(objective - objective_old) / (1 + abs(objective_old))
    }
    trace <- rbind(trace, data.frame(
      outer = outer, objective = objective,
      dobjective = dobjective, fixedpoint_G = fixedpoint_G,
      dphi = dphi, dsp = dsp,
      sig2 = sigma2, fixed_updates = inner
    ))
    if (verbose) {
      cat(
        "outer", outer, "total REML", format(objective, digits = 8),
        "relative change", format(dobjective, digits = 3),
        "G updates", inner, "\n"
      )
    }
    converged <- outer > 2L && dobjective < control$objective.tol
    sp_old <- sp_now
    sigma2_old <- sigma2
    objective_old <- objective
    if (converged) break
  }

  final <- gammfast_global_step(
    G0, Sl, X, y_work, B, id, G, sigma2, rho_start, nthreads,
    crossprod_cache = crossprod_cache, estimate_phi = TRUE
  )
  sigma2 <- final$phi
  beta <- final$beta
  names(beta) <- names(shell$coefficients)
  eta_global_work <- drop(X %*% beta)
  eta_global <- offset + eta_global_work
  u <- sqrt(sigma2) * gammfast_gaussian_blup_cached(
    crossprod_cache$BtB, crossprod_cache$BtA, G, beta, sigma2
  )
  rownames(u) <- levels(id_factor)
  colnames(u) <- random_structure$column.names
  eta_random <- rowSums(B * u[id, , drop = FALSE])
  eta <- eta_global + eta_random
  random_structure$covariance <- lapply(
    random_structure$group.index,
    function(jj) sigma2 * G[jj, jj, drop = FALSE]
  )

  shell$coefficients <- beta
  shell$sp <- final$sp
  shell$sig2 <- sigma2
  shell$scale <- sigma2
  shell$weights <- prior_weights
  shell$Vp <- final$pp$Vp
  shell$Ve <- final$pp$Ve
  shell$linear.predictors <- eta_global
  shell$fitted.values <- eta_global

  fit <- list(
    coefficients = beta,
    Vp = final$pp$Vp,
    Ve = final$pp$Ve,
    random.effects = u,
    G = sigma2 * G,
    dispersion.method = "mgcv-fREML",
    family.parameter.method = "family-fixed",
    covariance.method = "mean-Hessian-projected-moment",
    influence.update = "not-applicable",
    sp = final$sp,
    fitted.values = eta,
    linear.predictors = eta,
    global.fitted = eta_global,
    random.fitted = eta_random,
    residuals = y - eta,
    y = y, weights = prior_weights, offset = offset,
    sig2 = sigma2,
    prior.weights = prior_weights,
    objective = final$fr$reml + 0.5 * final$logdet,
    converged = converged,
    outer = nrow(trace),
    inner.max = inner.max, inner.tol = inner.tol,
    elapsed = proc.time()[[3]] - start_time,
    trace = trace,
    formula = formula,
    global.formula = global_formula,
    global = shell,
    family = family,
    random = random_structure,
    call = match.call(),
    control = control
  )
  class(fit) <- "gammfast"
  if (!converged) {
    warning("gammfast reached control$max.outer before convergence.")
  }
  fit
}

gammfast_validate_family <- function(family) {
  if (is.character(family)) {
    if (length(family) != 1L) stop("family must have length one.")
    family <- switch(tolower(family),
      gaussian = stats::gaussian(),
      binomial = stats::binomial(),
      quasibinomial = stats::quasibinomial(),
      poisson = stats::poisson(),
      quasipoisson = stats::quasipoisson(),
      gamma = stats::Gamma(),
      inverse.gaussian = stats::inverse.gaussian(),
      quasi = stats::quasi(),
      stop("Unsupported character family; use an mgcv family object for custom settings.")
    )
  }
  if (!inherits(family, "family")) {
    stop("family must be an mgcv family object.")
  }
  family_name <- tolower(family$family[1])
  if (identical(family$family[1], "Cox PH")) {
    stop("gammfast does not support Cox PH because its working Hessian is not observation-diagonal.")
  }
  if (grepl("^zero inflated poisson", family_name)) {
    stop("gammfast does not support zero-inflated Poisson because its working Hessian is not observation-diagonal.")
  }
  gaussian_family <- identical(family_name, "gaussian")
  gaussian <- gaussian_family && identical(family$link, "identity")
  supported <- gaussian_family ||
    identical(family_name, "binomial") ||
    identical(family_name, "quasibinomial") ||
    identical(family_name, "poisson") ||
    identical(family_name, "quasipoisson") ||
    identical(family_name, "gamma") ||
    identical(family_name, "inverse.gaussian") ||
    identical(family_name, "quasi") ||
    grepl("^negative binomial", family_name) ||
    identical(family_name, "scaled t") ||
    identical(family_name, "beta regression") ||
    grepl("^tweedie", family_name)
  if (!supported) {
    stop(
      "Unsupported family. See ?gammfast for the supported mgcv families."
    )
  }
  list(family = family, gaussian = gaussian)
}

gammfast_prefit_shell <- function(G0, family) {
  shell <- G0
  shell$model <- G0$mf
  shell$y <- as.numeric(G0$y)
  shell$prior.weights <- as.numeric(G0$w)
  shell$weights <- as.numeric(G0$w)
  shell$offset <- as.numeric(G0$offset)
  shell$coefficients <- rep(0, ncol(G0$X))
  names(shell$coefficients) <- G0$term.names
  family_name <- tolower(family$family[1])
  if (identical(family_name, "gaussian") &&
      identical(family$link, "identity")) {
    eta <- shell$y
  } else if (family_name %in% c("binomial", "quasibinomial")) {
    eta <- rep(family$linkfun(0.5), length(shell$y))
  } else if (family_name %in% c("poisson", "quasipoisson")) {
    eta <- family$linkfun(shell$y + 0.1)
  } else {
    mu <- pmax(shell$y, sqrt(.Machine$double.eps))
    valid <- is.finite(mu) & family$validmu(mu)
    if (!all(valid)) {
      center <- mean(mu[valid])
      if (!is.finite(center)) center <- 1
      mu[!valid] <- center
    }
    eta <- family$linkfun(mu)
    eta[!is.finite(eta)] <- 0
  }
  shell$linear.predictors <- as.numeric(eta)
  shell$fitted.values <- family$linkinv(shell$linear.predictors)
  shell$residuals <- shell$y - shell$fitted.values
  shell$family <- family
  shell$sig2 <- 1
  shell$scale <- 1
  shell$sp <- numeric(0)
  shell$method <- "fREML"
  shell$xlevels <- stats::.getXlevels(G0$terms, G0$mf)
  shell$na.action <- attr(G0$mf, "na.action")
  class(shell) <- c("gam", "glm", "lm")
  shell
}

gammfast_validate_smooths <- function(G0) {
  for (s in G0$smooth) {
    cls <- class(s)
    if (any(grepl("random.effect|fs.interaction", cls))) {
      stop(
        "Subject-level random effects must be supplied as s(..., bs = 're') ",
        "or custom fs(x, id, k) terms in the gammfast formula."
      )
    }
    if (length(s$S) > 1L) {
      stop("Overlapping or multi-penalty global smooths are not currently supported.")
    }
  }
  invisible(NULL)
}

gammfast_working <- function(family, y, eta, prior_weights) {
  if (inherits(family, "extended.family") && is.function(family$Dd)) {
    mu <- family$linkinv(eta)
    theta <- family$getTheta()
    Dval <- family$Dd(y, mu, theta, prior_weights, level = 0)
    mu_eta <- family$mu.eta(eta)
    D_eta <- Dval$Dmu * mu_eta
    w <- Dval$EDmu2 * mu_eta^2 / 2
    z <- eta - D_eta / (2 * w)
  } else {
    mu <- family$linkinv(eta)
    mu_eta <- family$mu.eta(eta)
    variance <- family$variance(mu)
    z <- eta + (y - mu) / mu_eta
    w <- prior_weights * mu_eta^2 / variance
  }
  if (length(z) != length(y) || length(w) != length(y) ||
      any(!is.finite(z)) || any(!is.finite(w)) || any(w < 0) ||
      !any(w > 0)) {
    stop("The family produced an invalid diagonal PIRLS working system.")
  }
  list(z = as.numeric(z), w = as.numeric(w))
}

gammfast_t_correction <- function(family, y, eta, prior_weights,
                                  working_weight) {
  family_name <- tolower(family$family[1])
  mu <- family$linkinv(eta)
  determinant_weight <- working_weight
  determinant_derivative <- NULL

  if (inherits(family, "extended.family")) {
    theta <- family$getTheta(TRUE)
    family <- utils::getFromNamespace("fix.family.link", "mgcv")(family)
    link_term <- -2 * family$g2g(mu)
    mu_eta <- family$mu.eta(eta)
    if (grepl("^negative binomial", family_name)) {
      determinant_derivative <- working_weight * (
        link_term - mu_eta * (1 / mu + 1 / (theta + mu))
      )
    } else if (identical(family_name, "tweedie")) {
      determinant_derivative <- working_weight * (
        link_term - theta * mu_eta / mu
      )
    } else if (identical(family_name, "scaled t")) {
      determinant_derivative <- working_weight * link_term
    } else if (identical(family_name, "beta regression")) {
      information <- base::psigamma(mu * theta, deriv = 1) +
        base::psigamma((1 - mu) * theta, deriv = 1)
      information_derivative <- base::psigamma(mu * theta, deriv = 2) -
        base::psigamma((1 - mu) * theta, deriv = 2)
      determinant_derivative <- working_weight * (
        link_term +
          theta * mu_eta * information_derivative / information
      )
    }
  } else {
    fix_link <- utils::getFromNamespace("fix.family.link", "mgcv")
    fix_variance <- utils::getFromNamespace("fix.family.var", "mgcv")
    family <- fix_variance(fix_link(family))
    mu_eta <- family$mu.eta(eta)
    mu_eta2 <- -family$d2link(mu) * mu_eta^3
    variance <- family$variance(mu)
    variance_derivative <- family$dvar(mu)
    fisher_derivative <- working_weight * (
      2 * mu_eta2 / mu_eta -
        mu_eta * variance_derivative / variance
    )
    determinant_derivative <- fisher_derivative
  }

  if (is.null(determinant_derivative)) return(NULL)
  if (length(determinant_weight) != length(y) ||
      length(determinant_derivative) != length(y) ||
      any(!is.finite(determinant_weight)) ||
      any(!is.finite(determinant_derivative)) ||
      any(determinant_weight < 0)) {
    stop("The family produced an invalid Laplace determinant correction.")
  }
  list(
    weight = as.numeric(determinant_weight),
    derivative = as.numeric(determinant_derivative)
  )
}

gammfast_update_theta <- function(family, y, mu, prior_weights, scale) {
  n_theta <- if (is.null(family$n.theta)) 0L else family$n.theta
  scale_control <- if (is.null(family$scale)) 1 else family$scale
  estimate_scale <- is.numeric(scale_control) && length(scale_control) == 1L &&
    scale_control < 0
  if (n_theta < 1L && !estimate_scale) {
    return(list(scale = scale, change = 0))
  }
  theta_old <- family$getTheta()
  state_old <- c(theta_old, if (estimate_scale) log(scale) else numeric())
  estimate.theta <- utils::getFromNamespace("estimate.theta", "mgcv")
  theta_new <- estimate.theta(
    theta_old, family, y, mu, scale = scale_control,
    wt = prior_weights, tol = 1e-8
  )
  if (estimate_scale) {
    scale <- exp(theta_new[length(theta_new)])
    theta_new <- theta_new[-length(theta_new)]
  }
  if (n_theta > 0L) family$putTheta(theta_new)
  state_new <- c(theta_new, if (estimate_scale) log(scale) else numeric())
  change <- if (length(state_new)) {
    max(abs(state_new - state_old)) / (1 + max(abs(state_old)))
  } else {
    0
  }
  list(scale = scale, change = change)
}

gammfast_estimate_working_scale <- function(family) {
  family_name <- tolower(family$family[1])
  family_name %in% c(
    "gaussian", "gamma", "inverse.gaussian", "quasi", "quasibinomial",
    "quasipoisson"
  ) || grepl("^tweedie\\(", family_name)
}

gammfast_reml_parameters <- function(fr, G0, estimate_phi, phi = 1) {
  nsp <- length(G0$S)
  expected <- nsp + as.integer(estimate_phi)
  if (length(fr$rho.full) != expected ||
      any(!is.finite(fr$rho.full))) {
    stop("mgcv returned inconsistent fREML parameters.")
  }
  sp <- exp(fr$rho.full[seq_len(nsp)])
  if (estimate_phi) {
    if (length(fr$rho) < 2L) stop("mgcv did not return the fREML scale.")
    phi <- exp(fr$rho.full[expected])
    rho <- fr$rho[-length(fr$rho)]
  } else {
    rho <- fr$rho
  }
  if (!is.finite(phi) || phi <= 0) stop("mgcv returned an invalid scale.")
  list(rho = rho, sp = sp, phi = phi)
}

gammfast_general_step <- function(G0, Sl, X, z, B, id, G, sw,
                                  rho_start, nthreads, optimize_sp = TRUE,
                                  sp_fixed = NULL, phi = 1,
                                  estimate_phi = FALSE) {
  Sl.Xprep <- utils::getFromNamespace("Sl.Xprep", "mgcv")
  Sl.postproc <- utils::getFromNamespace("Sl.postproc", "mgcv")
  initial.sp <- utils::getFromNamespace("initial.sp", "mgcv")
  fast.REML.fit <- utils::getFromNamespace("fast.REML.fit", "mgcv")
  A <- cbind(X, z) * sw
  Bw <- B * sw
  cp <- gammfast_gaussian_crossprod(
    A, Bw, id, G, n_threads = nthreads
  )
  H <- (cp$crossprod + t(cp$crossprod)) / 2
  p <- ncol(X)
  XtX <- H[seq_len(p), seq_len(p), drop = FALSE]
  Xtz <- H[seq_len(p), p + 1L]
  ztz <- H[p + 1L, p + 1L]
  if (!optimize_sp) {
    if (estimate_phi) {
      stop("A fixed-smoothing step cannot estimate the working scale.")
    }
    if (is.null(sp_fixed) || length(sp_fixed) != length(G0$S) ||
        any(!is.finite(sp_fixed)) || any(sp_fixed <= 0)) {
      stop("Fixed global smoothing parameters are invalid.")
    }
    penalty <- gammfast_penalty_matrix(G0, sp_fixed)
    Q <- XtX
    Q <- (Q + penalty + t(Q + penalty)) / 2
    Rq <- chol(Q)
    beta <- backsolve(Rq, forwardsolve(t(Rq), Xtz))
    Qinv <- chol2inv(Rq)
    Vp <- phi * Qinv
    Ve <- phi * Qinv %*% XtX %*% Qinv
    rss <- max(0, ztz - 2 * sum(beta * Xtz) +
      sum(beta * drop(XtX %*% beta)))
    criterion <- 0.5 *
      (rss + sum(beta * drop(penalty %*% beta))) / phi
    return(list(
      beta = beta,
      fr = list(
        rho = log(sp_fixed), rho.full = log(sp_fixed), reml = criterion
      ),
      pp = list(beta = beta, Vp = Vp, Ve = Ve),
      logdet = cp$logdet, sp_boundary = FALSE,
      rho = log(sp_fixed), sp = sp_fixed, phi = phi
    ))
  }
  R1 <- chol(XtX)
  f <- forwardsolve(t(R1), Xtz)
  rss_extra <- max(0, ztz - sum(f^2))
  um <- Sl.Xprep(Sl, R1, nt = nthreads)
  if (is.null(rho_start)) {
    rho_start <- log(initial.sp(R1, G0$S, G0$off))
  }
  boundary_warning <- FALSE
  fr <- withCallingHandlers(
    fast.REML.fit(
      um$Sl, um$X, f, rho = rho_start,
      L = G0$L, rho.0 = G0$lsp0,
      log.phi = log(phi), phi.fixed = !estimate_phi,
      rss.extra = rss_extra, nobs = length(z), Mp = um$Mp,
      nt = nthreads, gamma = 1
    ),
    warning = function(w) {
      if (grepl("Possible divergence detected in fast.REML.fit",
                conditionMessage(w), fixed = TRUE)) {
        boundary_warning <<- TRUE
        invokeRestart("muffleWarning")
      }
    }
  )
  pars <- gammfast_reml_parameters(
    fr, G0, estimate_phi = estimate_phi, phi = phi
  )
  pp <- Sl.postproc(
    Sl, fr, um$undrop, R1, cov = TRUE,
    scale = pars$phi, L = G0$L, nt = nthreads
  )
  list(
    beta = pp$beta, fr = fr, pp = pp, logdet = cp$logdet,
    sp_boundary = boundary_warning,
    rho = pars$rho, sp = pars$sp, phi = pars$phi
  )
}

gammfast_parse_formula <- function(formula, data) {
  tt <- stats::terms(formula, specials = "fs", keep.order = TRUE)
  labels <- attr(tt, "term.labels")
  drop_index <- integer()
  groups <- list()
  id_name <- NULL
  has_intercept <- FALSE
  formula_environment <- environment(formula)
  if (is.null(formula_environment)) formula_environment <- parent.frame()

  for (i in seq_along(labels)) {
    term_call <- str2lang(labels[i])
    if (!is.call(term_call)) next
    call_name <- as.character(term_call[[1L]])
    args <- as.list(term_call)[-1L]
    arg_names <- names(args)
    if (is.null(arg_names)) arg_names <- rep("", length(args))
    positional <- which(!nzchar(arg_names))

    if (identical(call_name, "s") && "bs" %in% arg_names) {
      bs <- eval(args[[which(arg_names == "bs")[1L]]],
                 envir = formula_environment)
      if (length(bs) != 1L || !is.character(bs)) {
        stop("Each s(..., bs = ...) basis name must evaluate to one character value.")
      }
      if (identical(bs, "fs")) {
        stop(
          "gammfast does not support mgcv s(..., bs = 'fs'); use the custom ",
          "fs(x, id, k) marker."
        )
      }
      if (!identical(bs, "re")) next
      variables <- vapply(args[positional], function(x) {
        if (!is.symbol(x)) {
          stop("gammfast s(..., bs = 're') variables must be bare column names.")
        }
        as.character(x)
      }, character(1))
      missing_variables <- setdiff(variables, names(data))
      if (length(missing_variables)) {
        stop("Random-effect variables are absent from data: ",
             paste(missing_variables, collapse = ", "), ".")
      }
      factor_variables <- variables[vapply(
        variables, function(variable) is.factor(data[[variable]]), logical(1)
      )]
      if (length(factor_variables) != 1L) {
        stop("Each s(..., bs = 're') must contain exactly one factor ID variable.")
      }
      numeric_variables <- setdiff(variables, factor_variables)
      if (length(numeric_variables) &&
          any(!vapply(numeric_variables, function(variable) {
            x <- data[[variable]]
            is.numeric(x) || is.logical(x)
          }, logical(1)))) {
        stop("Random-slope covariates must be scalar numeric variables.")
      }
      this_id <- unname(factor_variables[1L])
      if (is.null(id_name)) id_name <- this_id
      if (!identical(id_name, this_id)) {
        stop("All gammfast random groups must share the same factor ID variable.")
      }
      is_intercept <- length(numeric_variables) == 0L
      has_intercept <- has_intercept || is_intercept
      groups[[length(groups) + 1L]] <- list(
        type = if (is_intercept) "random_intercept" else "random_slope",
        label = labels[i], id = this_id, variables = numeric_variables,
        k = 1L
      )
      drop_index <- c(drop_index, i)
      next
    }

    if (identical(call_name, "fs")) {
      if (length(positional) < 2L ||
          !is.symbol(args[[positional[1L]]]) ||
          !is.symbol(args[[positional[2L]]])) {
        stop("fs must be written as fs(x, id, k = ...).")
      }
      variable <- as.character(args[[positional[1L]]])
      this_id <- as.character(args[[positional[2L]]])
      missing_variables <- setdiff(c(variable, this_id), names(data))
      if (length(missing_variables)) {
        stop("fs variables are absent from data: ",
             paste(missing_variables, collapse = ", "), ".")
      }
      if (!is.numeric(data[[variable]])) {
        stop("The custom fs covariate must be numeric: ", variable, ".")
      }
      if (!is.factor(data[[this_id]])) {
        stop("The custom fs ID variable must be a factor: ", this_id, ".")
      }
      if (anyNA(data[[variable]]) || anyNA(data[[this_id]])) {
        stop("gammfast does not allow missing custom fs covariates or IDs.")
      }
      if (is.null(id_name)) id_name <- this_id
      if (!identical(id_name, this_id)) {
        stop("All gammfast random groups must share the same factor ID variable.")
      }
      k_position <- which(arg_names == "k")
      k <- if (length(k_position)) {
        eval(args[[k_position[1L]]], envir = formula_environment)
      } else {
        5L
      }
      if (length(k) != 1L || !is.finite(k) || k < 1L ||
          k != as.integer(k)) {
        stop("The fs basis dimension k must be a positive integer.")
      }
      groups[[length(groups) + 1L]] <- list(
        type = "fs_cosine", label = labels[i], id = this_id,
        variables = variable, k = as.integer(k)
      )
      drop_index <- c(drop_index, i)
    }
  }

  if (is.null(id_name) || !length(groups)) {
    stop("formula must contain subject-level random groups.")
  }
  if (!has_intercept) {
    stop("formula must contain an explicit random intercept s(id, bs = 're').")
  }
  global_formula <- formula(stats::drop.terms(
    tt, sort(unique(drop_index)), keep.response = TRUE
  ))
  environment(global_formula) <- environment(formula)
  list(formula = global_formula, id = id_name, groups = groups)
}

gammfast_cosine_basis <- function(time01, k) {
  B <- matrix(0, length(time01), k)
  for (j in seq_len(k)) {
    B[, j] <- sqrt(2) * cos(j * pi * time01)
  }
  B
}

gammfast_random_structure <- function(groups, data, id_name, id_levels) {
  B_list <- vector("list", length(groups))
  group_index <- vector("list", length(groups))
  column_names <- character()
  column_end <- 0L
  for (j in seq_along(groups)) {
    group <- groups[[j]]
    if (identical(group$type, "fs_cosine")) {
      x <- data[[group$variables]]
      xr <- range(x)
      if (any(!is.finite(xr)) || diff(xr) <= 0) {
        stop("The fs covariate must have a non-zero finite range: ",
             group$variables, ".")
      }
      x01 <- (x - xr[1L]) / diff(xr)
      Bj <- gammfast_cosine_basis(x01, group$k)
      group$range <- xr
      names_j <- paste0("fs(", group$variables, "):cos", seq_len(group$k))
    } else {
      Bj <- rep(1, nrow(data))
      if (length(group$variables)) {
        for (variable in group$variables) {
          value <- data[[variable]]
          if (anyNA(value) || any(!is.finite(value))) {
            stop("Random-slope covariates must be finite: ", variable, ".")
          }
          Bj <- Bj * as.numeric(value)
        }
      }
      Bj <- matrix(Bj, ncol = 1L)
      names_j <- if (identical(group$type, "random_intercept")) {
        "re:(Intercept)"
      } else {
        paste0("re:", paste(group$variables, collapse = ":"))
      }
    }
    jj <- column_end + seq_len(ncol(Bj))
    column_end <- max(jj)
    group$columns <- jj
    groups[[j]] <- group
    B_list[[j]] <- Bj
    group_index[[j]] <- jj
    column_names <- c(column_names, names_j)
  }
  B <- do.call(cbind, B_list)
  colnames(B) <- column_names
  list(
    B = B, groups = groups, group.index = group_index,
    column.names = column_names, id = id_name, id.levels = id_levels
  )
}

gammfast_initial_covariance <- function(group_index, value) {
  K <- max(unlist(group_index, use.names = FALSE))
  G <- matrix(0, K, K)
  for (jj in group_index) G[jj, jj] <- diag(value, length(jj))
  G
}

gammfast_penalty_matrix <- function(G0, sp, scale = 1) {
  p <- ncol(G0$X)
  if (is.null(p) || p < 1L) stop("The global model matrix is invalid.")
  if (length(G0$S) != length(G0$off) || length(G0$S) != length(sp) ||
      any(!is.finite(sp)) || any(sp <= 0)) {
    stop("The global penalties and smoothing parameters are inconsistent.")
  }
  if (length(scale) != 1L || !is.finite(scale) || scale <= 0) {
    stop("The penalty scale must be one positive finite value.")
  }
  penalty <- matrix(0, p, p)
  for (j in seq_along(G0$S)) {
    ii <- G0$off[j] + seq_len(nrow(G0$S[[j]])) - 1L
    penalty[ii, ii] <- penalty[ii, ii] + sp[j] * G0$S[[j]]
  }
  (penalty + t(penalty)) / (2 * scale)
}

gammfast_project_covariance <- function(G, group_index) {
  out <- matrix(0, nrow(G), ncol(G))
  for (jj in group_index) {
    A <- (G[jj, jj, drop = FALSE] + t(G[jj, jj, drop = FALSE])) / 2
    E <- eigen(A, symmetric = TRUE)
    scale <- max(1, max(abs(E$values)))
    E$values <- pmax(E$values, scale * 1e-10)
    out[jj, jj] <- E$vectors %*% (E$values * t(E$vectors))
  }
  (out + t(out)) / 2
}

gammfast_random_info <- function(object) {
  if (is.null(object$random)) {
    stop("The gammfast object has no random structure.")
  }
  object$random
}

gammfast_random_design <- function(random_structure, newdata) {
  B_list <- vector("list", length(random_structure$groups))
  clamped <- rep(FALSE, nrow(newdata))
  for (j in seq_along(random_structure$groups)) {
    group <- random_structure$groups[[j]]
    missing_variables <- setdiff(group$variables, names(newdata))
    if (length(missing_variables)) {
      stop("newdata is missing random-effect variables: ",
           paste(missing_variables, collapse = ", "), ".")
    }
    if (identical(group$type, "fs_cosine")) {
      x <- newdata[[group$variables]]
      if (!is.numeric(x) || any(!is.finite(x))) {
        stop("Newdata fs covariates must be finite and numeric.")
      }
      x01 <- (x - group$range[1L]) / diff(group$range)
      clamped <- clamped | x01 < 0 | x01 > 1
      B_list[[j]] <- gammfast_cosine_basis(pmax(0, pmin(1, x01)), group$k)
    } else {
      Bj <- rep(1, nrow(newdata))
      for (variable in group$variables) {
        value <- newdata[[variable]]
        if ((!is.numeric(value) && !is.logical(value)) ||
            anyNA(value) || any(!is.finite(value))) {
          stop("Newdata random-slope covariates must be finite and numeric.")
        }
        Bj <- Bj * as.numeric(value)
      }
      B_list[[j]] <- matrix(Bj, ncol = 1L)
    }
  }
  B <- do.call(cbind, B_list)
  colnames(B) <- random_structure$column.names
  list(B = B, clamped = clamped)
}

gammfast_global_step <- function(G0, Sl, X, y, B, id, G, sigma2,
                                 rho_start, nthreads,
                                 crossprod_cache = NULL,
                                 estimate_phi = FALSE) {
  Sl.Xprep <- utils::getFromNamespace("Sl.Xprep", "mgcv")
  Sl.postproc <- utils::getFromNamespace("Sl.postproc", "mgcv")
  fast.REML.fit <- utils::getFromNamespace("fast.REML.fit", "mgcv")
  if (is.null(crossprod_cache)) {
    A <- cbind(X, y)
    cp <- gammfast_gaussian_crossprod(
      A, B, id, G, n_threads = nthreads
    )
  } else {
    cp <- gammfast_gaussian_crossprod_cached(
      crossprod_cache$AtA, crossprod_cache$BtB,
      crossprod_cache$BtA, G, n_threads = nthreads
    )
  }
  H <- (cp$crossprod + t(cp$crossprod)) / 2
  p <- ncol(X)
  XtX <- H[seq_len(p), seq_len(p), drop = FALSE]
  Xty <- H[seq_len(p), p + 1L]
  yty <- H[p + 1L, p + 1L]
  R1 <- chol(XtX)
  f <- forwardsolve(t(R1), Xty)
  rss_extra <- max(0, yty - sum(f^2))
  um <- Sl.Xprep(Sl, R1, nt = nthreads)
  fr <- fast.REML.fit(
    um$Sl, um$X, f, rho = rho_start,
    L = G0$L, rho.0 = G0$lsp0,
    log.phi = log(sigma2), phi.fixed = !estimate_phi,
    rss.extra = rss_extra, nobs = length(y), Mp = um$Mp,
    nt = nthreads, gamma = 1
  )
  pars <- gammfast_reml_parameters(
    fr, G0, estimate_phi = estimate_phi, phi = sigma2
  )
  pp <- Sl.postproc(
    Sl, fr, um$undrop, R1, cov = TRUE,
    scale = pars$phi, L = G0$L, nt = nthreads
  )
  list(
    beta = pp$beta, fr = fr, pp = pp, logdet = cp$logdet,
    rho = pars$rho, sp = pars$sp, phi = pars$phi
  )
}

#' @export
predict.gammfast <- function(object, newdata = NULL,
                             type = c("link", "response"),
                             include.random = TRUE, ...) {
  type <- match.arg(type)
  if (is.null(newdata)) {
    eta <- if (include.random) object$linear.predictors else object$global.fitted
  } else {
    eta <- as.numeric(stats::predict(object$global, newdata = newdata, type = "link"))
    if (include.random) {
      random <- gammfast_random_info(object)
      if (!random$id %in% names(newdata)) {
        stop("newdata must contain the gammfast subject ID variable.")
      }
      design <- gammfast_random_design(random, newdata)
      ii <- match(as.character(newdata[[random$id]]), random$id.levels)
      u <- matrix(0, length(ii), ncol(design$B))
      known <- !is.na(ii)
      u[known, ] <- object$random.effects[ii[known], , drop = FALSE]
      eta <- eta + rowSums(design$B * u)
    }
  }
  if (type == "response") object$family$linkinv(eta) else eta
}

#' @export
print.gammfast <- function(x, ...) {
  cat(x$family$family, "gammfast fit\n")
  cat("Formula: ")
  print(x$formula)
  cat("Subjects:", nrow(x$random.effects), "\n")
  random <- gammfast_random_info(x)
  cat("Independent random groups:", length(random$groups), "\n")
  cat("Subject random dimension:", ncol(x$random.effects), "\n")
  cat("Outer iterations:", x$outer, "\n")
  cat("Converged:", x$converged, "\n")
  invisible(x)
}

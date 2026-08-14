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
#' @param family An mgcv family object. Supported families are Gaussian with
#'   identity link, binomial with logit or probit link, Poisson and
#'   quasi-Poisson with log link, negative binomial, scaled t, beta regression,
#'   Tweedie, and ordered categorical.
#' @param weights Optional prior weights.
#' @param inner.max Maximum shared-UID covariance updates per mgcv P-IRLS outer
#'   iteration, between 1 and 30.
#' @param nthreads Number of OpenMP threads used for subject-level operations.
#' @param discrete Must be `FALSE`. The current solver does not implement the
#'   mgcv/bam discrete model-matrix representation.
#' @param control A list with `objective.tol`, `fixedpoint.tol`, and
#'   `max.outer`. Local covariance updates stop early at `fixedpoint.tol`.
#' @param verbose Whether to print outer-iteration diagnostics.
#'
#' @return An object of class `"gammfast"` containing the global fit, the full
#'   random-trajectory covariance `G` on the linear-predictor scale, its
#'   unit-dispersion representation `G.normalized`, subject BLUPs, fitted values,
#'   and convergence diagnostics.
#'
#' @details `gammfast()` deliberately does not use mgcv's discrete fitting
#'   path. mgcv constructs the global design and penalties (`fit = FALSE`) and
#'   its fREML machinery estimates global smoothing parameters and dispersion.
#'   `gammfast()` replaces only the working cross-products by their current
#'   shared-UID marginal forms and estimates the subject covariance without
#'   expanding an `n` by `n_subject * q` random-effect design. For canonical
#'   binomial-logit and Poisson-log models, the covariance score includes the
#'   conditional, penalized-mean, and Laplace conditional-mode corrections.
#'   Global smooths remain mgcv mean coefficients with a `paraPen` precision;
#'   their covariance coupling is evaluated by a low-dimensional Schur
#'   complement. Extended-family nuisance parameters are updated by mgcv's
#'   family-parameter optimizer. Formula offsets are retained in fitting and
#'   prediction.
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
                     inner.max = 30L, nthreads = 1L,
                     discrete = FALSE, control = list(), verbose = FALSE) {
  if (!inherits(formula, "formula")) stop("formula must be a model formula.")
  if (!is.data.frame(data)) stop("data must be a data frame.")
  family_info <- gammfast_validate_family(family)
  family <- family_info$family
  if (is.null(weights)) weights <- rep(1, nrow(data))
  if (!is.numeric(weights) || length(weights) != nrow(data) ||
      any(!is.finite(weights)) || any(weights < 0)) {
    stop("weights must be a finite non-negative numeric vector with one value per row.")
  }
  if (!any(weights > 0)) stop("At least one weight must be positive.")
  weight_name <- ".gammfast_prior_weights"
  if (weight_name %in% names(data)) {
    stop("data contains the reserved column '.gammfast_prior_weights'.")
  }
  data[[weight_name]] <- weights
  if (length(inner.max) != 1L || !is.finite(inner.max) ||
      inner.max < 1 || inner.max > 30) {
    stop("inner.max must be an integer from 1 to 30.")
  }
  if (length(nthreads) != 1L || !is.finite(nthreads) || nthreads < 1) {
    stop("nthreads must be a positive integer.")
  }
  inner.max <- as.integer(inner.max)
  nthreads <- as.integer(nthreads)
  if (length(discrete) != 1L || is.na(discrete) || !is.logical(discrete)) {
    stop("discrete must be a single TRUE or FALSE value.")
  }
  if (isTRUE(discrete)) {
    stop("gammfast does not support discrete = TRUE; use discrete = FALSE.")
  }

  defaults <- list(
    objective.tol = 1e-7,
    fixedpoint.tol = 1e-7,
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
  if (!is.finite(control$fixedpoint.tol) || control$fixedpoint.tol <= 0) {
    stop("control$fixedpoint.tol must be positive.")
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
      inner.max = inner.max,
      nthreads = nthreads, control = control, verbose = verbose,
      call = match.call()
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
    sp_inner <- pars$sp
    root_sigma <- sqrt(sigma2)

    for (inner in seq_len(inner.max)) {
      mm <- gammfast_projected_moments(
        response = y_work / root_sigma,
        X = X / root_sigma, B = B, id = id, G = G,
        penalty = gammfast_penalty_matrix(G0, sp_inner, scale = sigma2),
        n_threads = nthreads
      )
      G <- gammfast_project_covariance(
        mm$moment_sum / ng, random_structure$group.index
      )
    }
    mm_check <- gammfast_projected_moments(
      response = y_work / root_sigma,
      X = X / root_sigma, B = B, id = id, G = G,
      penalty = gammfast_penalty_matrix(G0, sp_inner, scale = sigma2),
      n_threads = nthreads
    )
    G_check <- gammfast_project_covariance(
      mm_check$moment_sum / ng, random_structure$group.index
    )
    beta <- mm_check$beta
    eta_global_work <- drop(X %*% beta)
    fixedpoint_G <- norm(G_check - G, "F") / (1 + norm(G, "F"))
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
      fixedpoint_sigma = 0, dphi = dphi, dsp = dsp,
      sigma2 = sigma2
    ))
    if (verbose) {
      cat(
        "outer", outer, "objective", format(objective, digits = 8),
        "fixedpoint", format(fixedpoint_G, digits = 3),
        "\n"
      )
    }
    converged <- outer > 2L &&
      dobjective < control$objective.tol &&
      max(fixedpoint_G, dphi, dsp) < control$fixedpoint.tol
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
  root_sigma <- sqrt(sigma2)
  mm <- gammfast_projected_moments(
    response = y_work / root_sigma,
    X = X / root_sigma, B = B, id = id, G = G,
    penalty = gammfast_penalty_matrix(
      G0, final$sp, scale = sigma2
    ),
    n_threads = nthreads
  )
  u <- sqrt(sigma2) * mm$u
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
    G.normalized = G,
    sigma2 = sigma2,
    dispersion.method = "mgcv-fREML",
    family.parameter.method = "family-fixed",
    covariance.method = "mean-Hessian-projected-moment",
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
    inner.max = inner.max,
    elapsed = proc.time()[[3]] - start_time,
    trace = trace,
    formula = formula,
    global.formula = global_formula,
    global = shell,
    family = family,
    random = random_structure,
    fs = gammfast_legacy_fs(random_structure, sigma2 * G),
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
      poisson = stats::poisson(),
      quasipoisson = stats::quasipoisson(),
      stop("Character family must be 'gaussian', 'binomial', 'poisson', or 'quasipoisson'; use an mgcv family object otherwise.")
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
  gaussian <- identical(family_name, "gaussian")
  supported <- gaussian ||
    identical(family_name, "binomial") ||
    identical(family_name, "poisson") ||
    identical(family_name, "quasipoisson") ||
    grepl("^negative binomial", family_name) ||
    identical(family_name, "scaled t") ||
    identical(family_name, "beta regression") ||
    identical(family_name, "tweedie") ||
    identical(family_name, "ordered categorical")
  if (!supported) {
    stop(
      "Unsupported family. gammfast currently supports Gaussian, binomial, ",
      "Poisson, quasi-Poisson, negative binomial, scaled t, beta regression, ",
      "Tweedie, and ordered categorical."
    )
  }
  if (gaussian && !identical(family$link, "identity")) {
    stop("Gaussian gammfast currently requires the identity link.")
  }
  if (identical(family_name, "binomial") &&
      !family$link %in% c("logit", "probit")) {
    stop("Binomial gammfast currently supports only logit and probit links.")
  }
  if (family_name %in% c("poisson", "quasipoisson") &&
      !identical(family$link, "log")) {
    stop("Poisson and quasi-Poisson gammfast currently require the log link.")
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
  if (identical(family_name, "gaussian")) {
    eta <- shell$y
  } else if (identical(family_name, "binomial")) {
    eta <- rep(family$linkfun(0.5), length(shell$y))
  } else if (identical(family_name, "ordered categorical")) {
    eta <- rep(0, length(shell$y))
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

gammfast_non_gaussian <- function(formula, global_formula, shell, G0, Sl,
                                  X, y, prior_weights, offset, B, id,
                                  id_factor, random_structure, inner.max,
                                  nthreads, control, verbose, call) {
  family <- shell$family
  family_name <- tolower(family$family[1])
  laplace_variance <-
    (identical(family_name, "binomial") &&
       identical(family$link, "logit")) ||
    (identical(family_name, "poisson") &&
       identical(family$link, "log"))
  laplace_penalty <- if (laplace_variance) {
    gammfast_laplace_penalty_setup(G0)
  } else {
    NULL
  }
  X_penalized <- if (laplace_variance) {
    X %*% laplace_penalty$penalized_vectors
  } else {
    NULL
  }
  ng <- nlevels(id_factor)
  eta <- as.numeric(shell$linear.predictors)
  eta_global <- eta
  eta_random <- numeric(length(y))
  initial_variance <- stats::var(eta)
  if (!is.finite(initial_variance) || initial_variance <= 0) {
    initial_variance <- 1
  }
  G <- gammfast_initial_covariance(
    random_structure$group.index, initial_variance * 0.15
  )
  rho_start <- if (length(shell$sp) == length(G0$S)) {
    pmin(log(shell$sp), 25)
  } else {
    NULL
  }
  sp_old <- NULL
  sp_current <- NULL
  objective_old <- NULL
  trace <- data.frame()
  converged <- FALSE
  family_scale <- if (is.numeric(shell$sig2) && length(shell$sig2) == 1L) {
    shell$sig2
  } else {
    1
  }
  scale_old <- NULL
  estimate_working_phi <- gammfast_estimate_working_scale(family)
  sp_boundary <- !is.null(rho_start) && any(rho_start >= 20)
  variance_evaluations <- 0L
  local_updates <- 0L
  start_time <- proc.time()[[3]]

  for (outer in seq_len(control$max.outer)) {
    evaluations_start <- variance_evaluations
    updates_start <- local_updates
    work <- gammfast_working(
      family, y, eta, prior_weights, nthreads = nthreads
    )
    step <- gammfast_general_step(
      G0, Sl, X, work$z - offset, B, id, G,
      sqrt(work$w), rho_start, nthreads,
      optimize_sp = !sp_boundary,
      sp_fixed = sp_current, phi = family_scale,
      estimate_phi = estimate_working_phi
    )
    beta <- step$beta
    rho_start <- pmin(step$rho, 25)
    sp_current <- step$sp
    family_scale <- step$phi
    sp_boundary <- sp_boundary || isTRUE(step$sp_boundary) ||
      any(log(step$sp) >= 25)
    if (!is.null(step$fr$outer.info$conv) &&
        grepl("divergence", step$fr$outer.info$conv, ignore.case = TRUE)) {
      sp_boundary <- TRUE
    }
    sw <- sqrt(work$w)
    root_scale <- sqrt(family_scale)
    working_response <- sw * (work$z - offset) / root_scale
    Xw <- X * sw / root_scale
    Bw <- B * sw
    penalty <- gammfast_penalty_matrix(
      G0, step$sp, scale = family_scale
    )

    if (laplace_variance) {
      smooth_precision <- crossprod(
        laplace_penalty$penalized_vectors,
        penalty %*% laplace_penalty$penalized_vectors
      )
      covariance_map <- function(G_now) {
        variance_evaluations <<- variance_evaluations + 1L
        mm_now <- gammfast_projected_moments(
          response = working_response, X = Xw, B = Bw, id = id,
          G = G_now, penalty = penalty, n_threads = nthreads
        )
        vcm_now <- gammfast_laplace_variance_step(
          X_penalized = X_penalized, B = B, id = id, G = G_now,
          smooth_precision = smooth_precision,
          working_weight = work$w / family_scale,
          weight_derivative = work$dw / family_scale,
          u = root_scale * mm_now$u, n_threads = nthreads
        )
        moment_now <- vcm_now$moment_sum / (family_scale * ng)
        list(
          G = gammfast_project_covariance(
            moment_now, random_structure$group.index
          ),
          moment = moment_now, mm = mm_now, vcm = vcm_now
        )
      }
      for (inner in seq_len(inner.max)) {
        mapped <- covariance_map(G)
        local_change <- norm(mapped$G - G, "F") /
          (1 + norm(G, "F"))
        G <- mapped$G
        beta_inner <- mapped$mm$beta
        u_inner <- root_scale * mapped$mm$u
        eta_inner <- offset + drop(X %*% beta_inner) +
          rowSums(B * u_inner[id, , drop = FALSE])
        local_updates <- local_updates + 1L
        work <- gammfast_working(
          family, y, eta_inner, prior_weights, nthreads = nthreads
        )
        sw <- sqrt(work$w)
        working_response <- sw * (work$z - offset) / root_scale
        Xw <- X * sw / root_scale
        Bw <- B * sw
        if (local_change <= control$fixedpoint.tol) break
      }
    } else {
      for (inner in seq_len(inner.max)) {
        variance_evaluations <- variance_evaluations + 1L
        mm <- gammfast_projected_moments(
          response = working_response, X = Xw, B = Bw, id = id, G = G,
          penalty = penalty, n_threads = nthreads
        )
        G <- gammfast_project_covariance(
          mm$moment_sum / ng, random_structure$group.index
        )
      }
    }

    moment_check <- if (laplace_variance) {
      check <- covariance_map(G)
      mm <- check$mm
      check$moment * ng
    } else {
      variance_evaluations <- variance_evaluations + 1L
      mm <- gammfast_projected_moments(
        response = working_response, X = Xw, B = Bw, id = id, G = G,
        penalty = penalty, n_threads = nthreads
      )
      mm$moment_sum
    }
    G_check <- gammfast_project_covariance(
      moment_check / ng, random_structure$group.index
    )
    fixedpoint_G <- norm(G_check - G, "F") / (1 + norm(G, "F"))
    beta <- mm$beta
    eta_global <- offset + drop(X %*% beta)
    u <- root_scale * mm$u
    eta_random <- rowSums(B * u[id, , drop = FALSE])
    eta_new <- eta_global + eta_random
    mu_new <- family$linkinv(eta_new)
    theta_update <- gammfast_update_theta(
      family, y, mu_new, prior_weights, family_scale
    )
    updated_scale <- theta_update$scale
    if (!isTRUE(all.equal(updated_scale, family_scale))) {
      G <- G * family_scale / updated_scale
      family_scale <- updated_scale
    }

    sp_now <- step$sp
    dsp <- if (is.null(sp_old)) Inf else max(abs(log(sp_now / sp_old)))
    if (sp_boundary) dsp <- 0
    dphi <- if (is.null(scale_old)) Inf else
      abs(log(family_scale / scale_old))
    eta_change <- max(abs(eta_new - eta)) / (1 + max(abs(eta)))
    objective <- step$fr$reml + 0.5 * step$logdet
    dobjective <- if (is.null(objective_old)) {
      Inf
    } else {
      abs(objective - objective_old) / (1 + abs(objective_old))
    }
    trace <- rbind(trace, data.frame(
      outer = outer, objective = objective, dobjective = dobjective,
      fixedpoint_G = fixedpoint_G, fixedpoint_sigma = 0,
      dphi = dphi, dsp = dsp, sigma2 = family_scale,
      eta_change = eta_change, theta_change = theta_update$change,
      variance_evaluations = variance_evaluations - evaluations_start,
      variance_evaluations_total = variance_evaluations,
      local_updates = local_updates - updates_start,
      local_updates_total = local_updates
    ))
    if (verbose) {
      cat(
        "outer", outer, "objective", format(objective, digits = 8),
        "eta", format(eta_change, digits = 3),
        "fixedpoint", format(fixedpoint_G, digits = 3), "\n"
      )
    }
    converged <- outer > 2L &&
      dobjective < control$objective.tol &&
      max(fixedpoint_G, dphi, dsp, theta_update$change) <
        control$fixedpoint.tol
    eta <- eta_new
    sp_old <- sp_now
    scale_old <- family_scale
    objective_old <- objective
    if (converged) break
  }

  work <- gammfast_working(
    family, y, eta, prior_weights, nthreads = nthreads
  )
  final <- gammfast_general_step(
    G0, Sl, X, work$z - offset, B, id, G,
    sqrt(work$w), rho_start, nthreads,
    optimize_sp = FALSE, sp_fixed = sp_current,
    phi = family_scale, estimate_phi = FALSE
  )
  sw <- sqrt(work$w)
  root_scale <- sqrt(family_scale)
  mm <- gammfast_projected_moments(
    response = sw * (work$z - offset) / root_scale,
    X = X * sw / root_scale, B = B * sw,
    id = id, G = G,
    penalty = gammfast_penalty_matrix(
      G0, final$sp, scale = family_scale
    ),
    n_threads = nthreads
  )
  beta <- mm$beta
  names(beta) <- names(shell$coefficients)
  eta_global <- offset + drop(X %*% beta)
  u <- root_scale * mm$u
  rownames(u) <- levels(id_factor)
  colnames(u) <- random_structure$column.names
  eta_random <- rowSums(B * u[id, , drop = FALSE])
  eta <- eta_global + eta_random
  mu <- family$linkinv(eta)
  work_final <- gammfast_working(
    family, y, eta, prior_weights, nthreads = nthreads
  )

  shell$coefficients <- beta
  shell$sp <- final$sp
  shell$sig2 <- family_scale
  shell$scale <- family_scale
  shell$weights <- work_final$w
  shell$linear.predictors <- eta_global
  shell$fitted.values <- family$linkinv(eta_global)
  shell$family <- family
  shell$Vp <- final$pp$Vp
  shell$Ve <- final$pp$Ve

  G_normalized <- G
  G_actual <- family_scale * G_normalized
  random_structure$covariance <- lapply(
    random_structure$group.index,
    function(jj) G_actual[jj, jj, drop = FALSE]
  )

  fit <- list(
    coefficients = beta,
    Vp = final$pp$Vp,
    Ve = final$pp$Ve,
    random.effects = u, G = G_actual, G.normalized = G_normalized,
    sigma2 = family_scale,
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
    covariance.method = if (laplace_variance) {
      "mgcv-fREML-shared-UID-Laplace-fixedpoint"
    } else {
      "mean-Hessian-projected-moment"
    },
    sp = final$sp, fitted.values = mu,
    linear.predictors = eta, global.fitted = eta_global,
    random.fitted = eta_random, residuals = y - mu,
    y = y, weights = work_final$w, prior.weights = prior_weights,
    offset = offset, sig2 = family_scale,
    objective = final$fr$reml + 0.5 * final$logdet,
    converged = converged, outer = nrow(trace), inner.max = inner.max,
    variance.evaluations = variance_evaluations,
    local.updates = local_updates,
    elapsed = proc.time()[[3]] - start_time, trace = trace,
    formula = formula, global.formula = global_formula, global = shell,
    family = family,
    random = random_structure,
    fs = gammfast_legacy_fs(random_structure, G_actual),
    call = call, control = control
  )
  class(fit) <- "gammfast"
  if (!converged) {
    warning("gammfast reached control$max.outer before convergence.")
  }
  fit
}

gammfast_working <- function(family, y, eta, prior_weights, nthreads) {
  family_name <- tolower(family$family[1])
  if (identical(family_name, "ordered categorical")) {
    out <- ocat_folded(
      eta = eta, y_int = as.integer(y), alpha = family$getTheta(TRUE),
      eps_mu = 1e-12, n_threads = nthreads
    )
    z <- out$z_star
    w <- out$w_star
  } else if (inherits(family, "extended.family") &&
             is.function(family$Dd)) {
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
  dw <- if (identical(family_name, "binomial") &&
            identical(family$link, "logit")) {
    w * (1 - 2 * family$linkinv(eta))
  } else if (identical(family_name, "poisson") &&
             identical(family$link, "log")) {
    w
  } else {
    NULL
  }
  list(
    z = as.numeric(z), w = as.numeric(w),
    dw = if (is.null(dw)) NULL else as.numeric(dw)
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
  tolower(family$family[1]) %in% c("quasipoisson", "quasibinomial")
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
        data[variables], is.factor, logical(1)
      )]
      if (length(factor_variables) != 1L) {
        stop("Each s(..., bs = 're') must contain exactly one factor ID variable.")
      }
      numeric_variables <- setdiff(variables, factor_variables)
      if (length(numeric_variables) &&
          any(!vapply(data[numeric_variables], function(x) {
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

gammfast_legacy_fs <- function(random_structure, G = NULL) {
  fs_index <- which(vapply(random_structure$groups, function(x) {
    identical(x$type, "fs_cosine")
  }, logical(1)))
  if (length(fs_index) != 1L) return(NULL)
  group <- random_structure$groups[[fs_index]]
  out <- list(
    time = group$variables, id = random_structure$id, k = group$k,
    time.range = group$range, id.levels = random_structure$id.levels,
    B = random_structure$B[, group$columns, drop = FALSE],
    id.index = random_structure$id.index
  )
  if (!is.null(G)) out$G <- G[group$columns, group$columns, drop = FALSE]
  out
}

gammfast_random_info <- function(object) {
  if (!is.null(object$random)) return(object$random)
  if (is.null(object$fs)) stop("The gammfast object has no random structure.")
  group <- list(
    type = "fs_cosine", label = paste0("fs(", object$fs$time, ")"),
    id = object$fs$id, variables = object$fs$time, k = object$fs$k,
    range = object$fs$time.range, columns = seq_len(object$fs$k)
  )
  list(
    B = object$fs$B, groups = list(group),
    group.index = list(seq_len(object$fs$k)),
    column.names = colnames(object$random.effects), id = object$fs$id,
    id.levels = object$fs$id.levels, id.index = object$fs$id.index
  )
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

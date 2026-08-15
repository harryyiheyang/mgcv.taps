taps_A_matrix <- function(A, n) {
  if (length(A) == 0L) stop("getA must not return an empty object.")
  if (is.atomic(A) && is.null(dim(A))) {
    if (length(A) == 1L) A <- rep(A, n)
    if (length(A) == n) A <- matrix(A, ncol = 1L)
  }
  if (!is.matrix(A) || !is.numeric(A)) {
    stop("getA must return a numeric vector or matrix.")
  }
  if (nrow(A) != n) stop("The number of rows returned by getA must match the data.")
  if (any(!is.finite(A))) stop("getA must return finite values.")
  A
}

taps_score_refit <- function(fit, test.component, sp.refit, null.tol) {
  p <- length(fit$smooth)
  if (length(test.component) != 1L || !is.numeric(test.component) ||
      !is.finite(test.component) || test.component != floor(test.component) ||
      test.component < 1L || test.component > p) {
    stop("test.component must be a valid smooth-term index.")
  }
  test.component <- as.integer(test.component)

  target <- fit$smooth[[test.component]]
  is_taps <- inherits(target, "AMatern.smooth") ||
    inherits(target, "A2Matern.smooth")
  if (!is_taps || is.null(target$getA)) {
    stop("refit = TRUE requires an AMatern or A2Matern smooth with getA.")
  }
  if (length(target$S) != 1L || target$first.sp != target$last.sp) {
    stop("refit = TRUE currently requires a one-penalty AMatern smooth.")
  }
  if (!is.null(target$id)) {
    stop("refit = TRUE does not support a tested smooth with a shared id.")
  }
  if (!inherits(fit$formula, "formula")) {
    stop("refit = TRUE does not support multiple-formula models.")
  }

  refit_data <- as.data.frame(fit$model)
  n <- nrow(refit_data)
  refit_data$.taps_response <- stats::model.response(fit$model)
  x_args <- lapply(target$term, function(term) refit_data[[term]])
  if (any(vapply(x_args, is.null, logical(1)))) {
    stop("The fitted model frame does not contain the tested smooth covariates.")
  }
  A <- do.call(target$getA, c(x_args, list(target$para)))
  A <- taps_A_matrix(A, n)

  if (!identical(target$by, "NA")) {
    by <- refit_data[[target$by]]
    if (is.null(by) || !is.numeric(by) || is.matrix(by) ||
        !is.null(target$by.level)) {
      stop("refit = TRUE supports only a numeric by variable.")
    }
    A <- A * as.numeric(by)
  }

  X_full <- predict(fit, newdata = fit$model, type = "lpmatrix")
  target_idx <- target$first.para:target$last.para
  m <- target$null.space.dim
  if (length(m) != 1L || is.na(m) || m < 0L || m >= length(target_idx)) {
    stop("The tested smooth has no penalized component.")
  }
  pen_local <- seq.int(m + 1L, length(target_idx))
  B_test <- X_full[, target_idx[pen_local], drop = FALSE]
  S_test <- as.matrix(target$S[[1]][pen_local, pen_local, drop = FALSE])
  S_norm <- norm(S_test, "F")
  if (!is.finite(S_norm) || S_norm <= null.tol) {
    stop("The tested smooth has an invalid penalized block.")
  }

  Xp <- X_full[, seq_len(fit$nsdf), drop = FALSE]
  X_keep <- Xp
  keep <- logical(ncol(A))
  rank0 <- qr(X_keep, tol = null.tol)$rank
  for (j in seq_len(ncol(A))) {
    rank1 <- qr(cbind(X_keep, A[, j]), tol = null.tol)$rank
    if (rank1 > rank0) {
      keep[j] <- TRUE
      X_keep <- cbind(X_keep, A[, j])
      rank0 <- rank1
    }
  }
  A <- A[, keep, drop = FALSE]
  A_names <- if (ncol(A)) paste0(".taps_A_", seq_len(ncol(A))) else character()
  if (length(A_names)) refit_data[A_names] <- as.data.frame(A)

  tt <- stats::terms(fit$formula, specials = c("s", "te", "ti", "t2"))
  special <- sort(unique(unlist(attr(tt, "specials"), use.names = FALSE)))
  variables <- attr(tt, "variables")
  smooth_calls <- lapply(special, function(i) variables[[i + 1L]])
  split <- mgcv::interpret.gam(fit$formula)
  if (length(smooth_calls) != p || length(split$smooth.spec) != p) {
    stop("refit = TRUE does not support factor-by smooth expansion.")
  }

  nuisance_calls <- smooth_calls[-test.component]
  param_terms <- attr(stats::terms(split$pf), "term.labels")
  smooth_terms <- vapply(nuisance_calls, function(x) {
    paste(deparse(x), collapse = "")
  }, character(1))
  null_formula <- stats::reformulate(
    c(param_terms, smooth_terms, A_names), response = ".taps_response",
    intercept = attr(stats::terms(split$pf), "intercept")
  )
  environment(null_formula) <- environment(fit$formula)

  offset <- fit$offset
  if (is.null(offset)) offset <- rep(0, n)
  weights <- stats::model.weights(fit$model)
  if (is.null(weights)) weights <- rep(1, n)
  refit_data$.taps_offset <- as.numeric(offset)
  refit_data$.taps_weights <- as.numeric(weights)

  full_sp <- fit$sp
  if (!is.null(fit$full.sp) && length(fit$full.sp)) full_sp <- fit$full.sp
  target_sp <- target$first.sp
  if (length(full_sp) < target_sp) {
    stop("The fitted smoothing parameters do not match the tested smooth.")
  }
  nuisance_sp <- full_sp[-target_sp]
  sp <- if (sp.refit) rep(-1, length(nuisance_sp)) else nuisance_sp
  if (!length(sp)) sp <- NULL

  is_bam <- inherits(fit, "bam")
  is_discrete <- is_bam && !is.null(fit$dinfo)
  refit_method <- if (is_bam) "fREML" else "REML"
  refit_env <- new.env(parent = environment(fit$formula))
  refit_env$refit_data <- refit_data
  refit_env$refit_family <- fit$family

  refit_call <- fit$call
  refit_call[[1L]] <- if (is_bam) quote(mgcv::bam) else quote(mgcv::gam)
  refit_call$formula <- null_formula
  refit_call$data <- quote(refit_data)
  refit_call$family <- quote(refit_family)
  refit_call$weights <- quote(.taps_weights)
  refit_call$offset <- quote(.taps_offset)
  refit_call$subset <- NULL
  refit_call$na.action <- quote(stats::na.fail)
  refit_call$method <- refit_method
  refit_call$sp <- sp
  refit_call$G <- NULL
  refit_call$start <- NULL
  refit_call$coef <- NULL
  refit_call$in.out <- NULL
  if (is_bam) refit_call$discrete <- is_discrete
  refit_call$fit <- FALSE
  G0 <- eval(refit_call, envir = refit_env)

  eta0 <- fit$linear.predictors -
    as.vector(B_test %*% fit$coefficients[target_idx[pen_local]]) - offset
  working_weights <- fit$weights
  if (is.null(working_weights)) working_weights <- weights
  w <- sqrt(pmax(working_weights, 0))
  qrx <- qr(G0$X * w, tol = null.tol)
  if (qrx$rank != ncol(G0$X)) {
    stop("The null refit design is rank deficient after adding getA.")
  }
  start <- qr.coef(qrx, eta0 * w)
  if (length(start) != ncol(G0$X) || any(!is.finite(start))) {
    stop("The null refit produced invalid starting coefficients.")
  }

  refit_call$fit <- TRUE
  if (is_bam) {
    if (!is_discrete) refit_call$coef <- start
  } else {
    refit_call$start <- start
  }
  if (sp.refit && length(nuisance_sp) && !is_discrete) {
    refit_call$in.out <- list(sp = nuisance_sp, scale = fit$sig2)
  }
  fit0 <- eval(refit_call, envir = refit_env)

  X0 <- predict(fit0, newdata = fit0$model, type = "lpmatrix")
  sp0 <- fit0$sp
  if (!is.null(fit0$full.sp) && length(fit0$full.sp)) sp0 <- fit0$full.sp
  fit0$sp <- c(sp0, 1)
  fit0$.taps_score_X <- cbind(X0, B_test)
  fit0$coefficients <- c(fit0$coefficients, rep(0, ncol(B_test)))

  target0 <- target
  target0$S <- list(S_test)
  target0$rank <- ncol(B_test)
  target0$null.space.dim <- 0L
  target0$first.para <- ncol(X0) + 1L
  target0$last.para <- ncol(X0) + ncol(B_test)
  target0$first.sp <- length(sp0) + 1L
  target0$last.sp <- target0$first.sp
  target0$by <- "NA"
  target0$by.level <- NULL
  fit0$smooth <- append(fit0$smooth, list(target0), after = test.component - 1L)
  fit0
}

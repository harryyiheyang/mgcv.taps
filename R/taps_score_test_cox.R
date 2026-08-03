cox_peto_information <- function(X, eta, time, status, strata) {
  X <- as.matrix(X)
  eta <- as.numeric(eta)
  time <- as.numeric(time)
  status <- as.integer(status)
  strata <- as.integer(factor(strata))

  if (nrow(X) != length(eta) || length(time) != length(eta) ||
      length(status) != length(eta) || length(strata) != length(eta)) {
    stop("Cox design, eta, time, status, and strata must have matching rows.")
  }

  H <- matrix(0, ncol(X), ncol(X))
  U <- numeric(ncol(X))
  events <- 0L
  for (s in sort(unique(strata))) {
    ii <- strata == s
    if (!any(status[ii] == 1L)) next
    ss <- cox_peto_suffstat(
      X = X[ii, , drop = FALSE], eta = eta[ii], time = time[ii],
      status = status[ii]
    )
    H <- H + as.matrix(ss$information)
    U <- U + as.numeric(ss$score)
    events <- events + as.integer(ss$events)
  }

  if (events == 0L) stop("Cox score test requires at least one event.")
  list(information = (H + t(H)) / 2, score = U, events = events)
}

taps_score_test_cox <- function(fit, test.component = 1, null.tol = 1e-10,
                                method = "davies", max_eps = 1e-8,
                                max_iter = 1e5) {
  X <- predict(fit, newdata = fit$model, type = "lpmatrix")
  beta <- as.numeric(fit$coefficients)
  smooth_terms <- fit$smooth
  p <- length(smooth_terms)

  if (p < 1L) stop("fit must contain at least one smooth term.")
  if (length(test.component) != 1L || is.na(test.component)) {
    stop("test.component must be a valid smooth-term index.")
  }
  test_component_int <- as.integer(test.component)
  if (is.na(test_component_int) || test.component != test_component_int ||
      test_component_int < 1L || test_component_int > p) {
    stop("test.component must be a valid smooth-term index.")
  }
  test.component <- test_component_int

  all_random <- integer(0)
  test_random <- integer(0)
  theta_test <- NULL
  S_full <- matrix(0, ncol(X), ncol(X))

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
    S_matrix <- as.matrix(Reduce(`+`, Map(
      function(S, sp_value) S * sp_value, s$S, sp
    )))
    S_norm <- norm(S_matrix, "f")
    if (!is.finite(S_norm) || S_norm <= 0) {
      if (i == test.component) {
        stop("test.component must have a non-zero penalty matrix.")
      }
      next
    }

    reported_null_dim <- s$null.space.dim
    if (!is.null(s$getA)) {
      null_indices <- if (reported_null_dim > 0L) {
        indices[seq_len(reported_null_dim)]
      } else {
        integer(0)
      }
    } else {
      col_norms <- sqrt(colSums(S_matrix^2))
      null_indices <- indices[col_norms < null.tol]
    }
    random_indices <- setdiff(indices, null_indices)
    all_random <- c(all_random, random_indices)

    if (i == test.component) {
      local_random <- match(random_indices, indices)
      theta <- matrixGeneralizedInverse(S_matrix / S_norm)
      theta_test <- as.matrix(theta[local_random, local_random, drop = FALSE])
      test_random <- random_indices
    } else {
      S_full[indices, indices] <- S_full[indices, indices] + S_matrix
    }
  }

  all_random <- sort(unique(all_random))
  test_random <- sort(unique(test_random))
  if (!length(test_random) || is.null(theta_test)) {
    stop("Could not identify the penalized coefficients of test.component.")
  }

  fixed_indices <- setdiff(seq_len(ncol(X)), all_random)
  other_random <- setdiff(all_random, test_random)
  null_indices <- sort(c(fixed_indices, other_random))
  X_null <- X[, null_indices, drop = FALSE]
  S_null <- S_full[null_indices, null_indices, drop = FALSE]
  b_null <- beta[null_indices]
  offset <- as.numeric(fit$linear.predictors) -
    matrixVectorMultiply(X, beta)

  if (is.matrix(fit$y)) {
    time <- fit$y[, 1]
    strata <- fit$y[, 2]
  } else {
    time <- fit$y
    strata <- rep(1L, length(time))
  }
  status <- fit$prior.weights

  converged <- FALSE
  for (iter in seq_len(50L)) {
    eta_null <- offset + matrixVectorMultiply(X_null, b_null)
    stat_null <- cox_peto_information(
      X_null, eta_null, time, status, strata
    )
    g <- stat_null$score - matrixVectorMultiply(S_null, b_null)
    H_pen <- stat_null$information + S_null
    step <- matrixVectorMultiply(matrixGeneralizedInverse(H_pen), g)
    if (any(!is.finite(step))) {
      stop("Cox null fit produced a non-finite Newton step.")
    }
    b_null <- b_null + step
    d_eta <- matrixVectorMultiply(X_null, step)
    if (sqrt(mean(d_eta^2)) < 1e-8) {
      converged <- TRUE
      break
    }
  }
  if (!converged) stop("Cox null fit did not converge in 50 Newton iterations.")

  eta_null <- offset + matrixVectorMultiply(X_null, b_null)
  # The nuisance block holds every coefficient fitted by the null model, so the
  # penalized coefficients of the other smooths must be profiled out alongside
  # the unpenalized ones; their penalty enters the nuisance score and
  # information.
  Z_indices <- null_indices
  if (length(Z_indices)) {
    Z <- X[, Z_indices, drop = FALSE]
    S_Z <- S_full[Z_indices, Z_indices, drop = FALSE]
    keep_Z <- colSums((Z - matrix(colMeans(Z), nrow(Z), ncol(Z),
                                  byrow = TRUE))^2) > null.tol |
      diag(S_Z) > 0
    Z <- Z[, keep_Z, drop = FALSE]
    S_Z <- S_Z[keep_Z, keep_Z, drop = FALSE]
    b_Z <- b_null[keep_Z]
  } else {
    Z <- matrix(nrow = nrow(X), ncol = 0L)
    S_Z <- matrix(0, 0L, 0L)
    b_Z <- numeric(0)
  }
  B1 <- X[, test_random, drop = FALSE]
  D <- cbind(B1, Z)
  stat <- cox_peto_information(D, eta_null, time, status, strata)
  q <- ncol(B1)
  U1 <- stat$score[seq_len(q)]
  P1 <- stat$information[seq_len(q), seq_len(q), drop = FALSE]

  if (ncol(Z)) {
    iz <- q + seq_len(ncol(Z))
    H1Z <- stat$information[seq_len(q), iz, drop = FALSE]
    HZZ <- stat$information[iz, iz, drop = FALSE] + S_Z
    UZ <- stat$score[iz] - matrixVectorMultiply(S_Z, b_Z)
    HZZ_inv <- matrixGeneralizedInverse(HZZ)
    U1 <- U1 - matrixVectorMultiply(
      H1Z, matrixVectorMultiply(HZZ_inv, UZ)
    )
    P1 <- P1 - matrixMultiply(
      H1Z, matrixMultiply(HZZ_inv, H1Z, transB = TRUE)
    )
  }
  P1 <- (P1 + t(P1)) / 2

  u <- max(0, as.numeric(crossprod(
    U1, matrixVectorMultiply(theta_test, U1)
  )))
  theta_sqrt <- matrixsqrt(theta_test)$w
  Q_small <- matrixListProduct(list(theta_sqrt, P1, theta_sqrt))
  lambda <- eigen((Q_small + t(Q_small)) / 2,
                  symmetric = TRUE, only.values = TRUE)$values
  lambda <- lambda[lambda > 1e-15]
  if (!length(lambda)) stop("Cox score test has no positive mixture weights.")

  if (method == "satterthwaite") {
    e <- sum(lambda)
    h <- sum(lambda^2)
    kappa <- h / e
    nu <- e^2 / h
    pv <- pchisq(u / kappa, df = nu, lower.tail = FALSE)
  } else {
    pv <- switch(
      method,
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
      stop(paste0(
        "method must be one of 'satterthwaite', 'davies', 'liu', ",
        "'hbe', 'wood', 'saddlepoint'"
      ))
    )
  }

  data.table(
    smooth.term = smooth_terms[[test.component]]$label,
    smooth.pvalue = pv,
    method = method
  )
}

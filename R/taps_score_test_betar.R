betar_joint_information <- function(X, eta, y, theta, family, weights) {
  X <- as.matrix(X)
  eta <- as.numeric(eta)
  y <- as.numeric(y)
  theta <- as.numeric(theta)
  weights <- as.numeric(weights)

  if (length(theta) != 1L) {
    stop("Beta regression must have one precision parameter.")
  }
  if (nrow(X) != length(eta) || length(y) != length(eta) ||
      length(weights) != length(eta)) {
    stop("Beta design, linear predictor, response, and weights must have matching rows.")
  }
  if (!grepl("^beta regression", tolower(family$family))) {
    stop("The Beta score test supports only mgcv::betar().")
  }

  mu <- family$linkinv(eta)
  d <- mgcv::dDeta(y, mu, weights, theta, family, deriv = 2)
  if (!all(d$good)) {
    stop("Beta score calculation produced non-finite derivatives.")
  }

  s <- -as.numeric(d$Deta) / 2
  h <- as.numeric(d$Deta2) / 2
  cross_theta <- crossprod(X, as.numeric(d$Detath) / 2)
  H <- crossprod(X, X * h)

  list(
    score = as.numeric(crossprod(X, s)),
    information = (H + t(H)) / 2,
    cross_theta = cross_theta,
    theta_score = -sum(d$Dth) / 2,
    theta_information = sum(d$Dth2) / 2
  )
}

taps_score_test_betar <- function(fit, test.component = 1,
                                  null.tol = 1e-10,
                                  method = "davies", max_eps = 1e-8,
                                  max_iter = 1e5) {
  if (!grepl("^beta regression", tolower(fit$family$family))) {
    stop("The Beta score test supports only mgcv::betar().")
  }
  if (!isTRUE(fit$family$n.theta > 0L)) {
    stop("The Beta precision parameter must have been estimated by mgcv.")
  }

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
  Theta1 <- NULL
  S_full <- matrix(0, ncol(X), ncol(X))

  for (i in seq_len(p)) {
    sm <- smooth_terms[[i]]
    indices <- sm$first.para:sm$last.para
    is_zero_rank <- !is.null(sm$rank) && isTRUE(all(sm$rank == 0))
    is_fixed_smooth <- isTRUE(sm$fixed) || is.null(sm$S) ||
      length(sm$S) == 0L || is_zero_rank

    if (i == test.component && is_fixed_smooth) {
      stop("test.component must be a penalized smooth with fx = FALSE.")
    }
    if (is_fixed_smooth) next
    if (is.null(sm$first.sp) || is.null(sm$last.sp)) {
      stop("Penalized smooth has no smoothing-parameter index.")
    }

    sp_idx <- seq.int(sm$first.sp, sm$last.sp)
    if (length(sp_idx) != length(sm$S)) {
      stop("Penalized smooth has inconsistent penalty and smoothing-parameter indexes.")
    }
    sp <- fit$sp[sp_idx]
    if (anyNA(sp)) {
      stop("Penalized smooth has missing smoothing-parameter values.")
    }
    S_matrix <- as.matrix(Reduce(`+`, Map(
      function(S, sp_value) S * sp_value, sm$S, sp
    )))
    S_norm <- norm(S_matrix, "f")
    if (!is.finite(S_norm) || S_norm <= 0) {
      if (i == test.component) {
        stop("test.component must have a non-zero penalty matrix.")
      }
      next
    }

    reported_null_dim <- sm$null.space.dim
    if (!is.null(sm$getA)) {
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
      Theta <- matrixGeneralizedInverse(S_matrix / S_norm)
      Theta1 <- as.matrix(Theta[local_random, local_random, drop = FALSE])
      test_random <- random_indices
    } else {
      S_full[indices, indices] <- S_full[indices, indices] + S_matrix
    }
  }

  all_random <- sort(unique(all_random))
  test_random <- sort(unique(test_random))
  if (!length(test_random) || is.null(Theta1)) {
    stop("Could not identify the penalized coefficients of test.component.")
  }

  fixed_indices <- setdiff(seq_len(ncol(X)), all_random)
  other_random <- setdiff(all_random, test_random)
  null_indices <- sort(c(fixed_indices, other_random))
  X_null <- X[, null_indices, drop = FALSE]
  S_null <- S_full[null_indices, null_indices, drop = FALSE]
  b_null <- beta[null_indices]
  theta_null <- as.numeric(fit$family$getTheta(FALSE))
  y <- as.numeric(fit$y)
  weights <- if (is.null(fit$prior.weights)) {
    rep(1, length(y))
  } else {
    as.numeric(fit$prior.weights)
  }
  offset <- as.numeric(fit$linear.predictors) -
    matrixVectorMultiply(X, beta)

  converged <- FALSE
  for (iter in seq_len(50L)) {
    eta_null <- offset + matrixVectorMultiply(X_null, b_null)
    stat_null <- betar_joint_information(
      X_null, eta_null, y, theta_null, fit$family, weights
    )
    g <- c(
      stat_null$score - matrixVectorMultiply(S_null, b_null),
      stat_null$theta_score
    )
    H_pen <- rbind(
      cbind(stat_null$information + S_null, stat_null$cross_theta),
      c(t(stat_null$cross_theta), stat_null$theta_information)
    )
    step <- matrixVectorMultiply(matrixGeneralizedInverse(H_pen), g)
    if (any(!is.finite(step))) {
      stop("Beta null fit produced a non-finite Newton step.")
    }

    nb <- length(b_null)
    step_b <- step[seq_len(nb)]
    step_theta <- step[nb + 1L]
    b_null <- b_null + step_b
    theta_null <- theta_null + step_theta
    d_eta <- matrixVectorMultiply(X_null, step_b)
    if (max(sqrt(mean(d_eta^2)), abs(step_theta)) < 1e-8) {
      converged <- TRUE
      break
    }
  }
  if (!converged) {
    stop("Beta null fit did not converge in 50 Newton iterations.")
  }

  eta_null <- offset + matrixVectorMultiply(X_null, b_null)
  Z <- X[, fixed_indices, drop = FALSE]
  B1 <- X[, test_random, drop = FALSE]
  D <- cbind(B1, Z)
  stat <- betar_joint_information(
    D, eta_null, y, theta_null, fit$family, weights
  )
  q <- ncol(B1)
  U1 <- stat$score[seq_len(q)]
  H11 <- stat$information[seq_len(q), seq_len(q), drop = FALSE]
  H1theta <- stat$cross_theta[seq_len(q), , drop = FALSE]

  iz <- q + seq_len(ncol(Z))
  UZ <- stat$score[iz]
  H1Z <- stat$information[seq_len(q), iz, drop = FALSE]
  HZZ <- stat$information[iz, iz, drop = FALSE]
  HZtheta <- stat$cross_theta[iz, , drop = FALSE]
  UN <- c(UZ, stat$theta_score)
  H1N <- cbind(H1Z, H1theta)
  HNN <- rbind(
    cbind(HZZ, HZtheta),
    c(t(HZtheta), stat$theta_information)
  )
  HNN_inv <- matrixGeneralizedInverse(HNN)
  U1 <- U1 - matrixVectorMultiply(
    H1N, matrixVectorMultiply(HNN_inv, UN)
  )
  P1 <- H11 - matrixMultiply(
    H1N, matrixMultiply(HNN_inv, H1N, transB = TRUE)
  )
  P1 <- (P1 + t(P1)) / 2

  u <- max(0, as.numeric(crossprod(
    U1, matrixVectorMultiply(Theta1, U1)
  )))
  Theta_sqrt <- matrixsqrt(Theta1)$w
  Q_small <- matrixListProduct(list(Theta_sqrt, P1, Theta_sqrt))
  lambda <- eigen((Q_small + t(Q_small)) / 2,
                  symmetric = TRUE, only.values = TRUE)$values
  lambda <- lambda[lambda > 1e-15]
  if (!length(lambda)) {
    stop("Beta score test has no positive mixture weights.")
  }

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

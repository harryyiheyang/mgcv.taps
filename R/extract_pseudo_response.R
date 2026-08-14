extract_pseudo_response <- function(fit, ...) {
  args <- list(...)

  if (inherits(fit, "qgam")) {
    eta <- fit$linear.predictors
    y   <- fit$y
    theta <- tryCatch(log(fit$family$getTheta(TRUE)), error = function(e) fit$lsig)
    if (is.null(theta)) stop("Could not extract theta (lsig) from qgam object.")
    wt    <- if (is.null(fit$prior.weights)) rep(1, length(y)) else fit$prior.weights
    Dval  <- fit$family$Dd(y, eta, theta, wt)
    pseudo_response <- eta - Dval$Dmu / Dval$Dmu2
    W_diag <- Dval$Dmu2 / 2
    phi0   <- 1
    valid_idx <- NULL
  } else if (grepl("Ordered Categorical", fit$family$family)) {
    eta       <- as.numeric(fit$linear.predictors)
    y_int     <- as.integer(fit$y)
    alpha_cuts <- fit$family$getTheta(TRUE)
    cpp_result <- ocat_folded(
      eta       = eta,
      y_int     = y_int,
      alpha     = alpha_cuts,
      eps_mu    = args$eps_mu    %||% 1e-12,
      n_threads = args$n_threads %||% 1
    )
    pseudo_response <- cpp_result$z_star
    W_diag  <- cpp_result$w_star
    phi0    <- 1.0
    valid_idx <- !is.na(pseudo_response) & W_diag > 0
  } else if (grepl("^cnorm", fit$family$family)) {
    eta <- fit$linear.predictors
    censor_attr <- attr(fit$y, "censor")
    y_lower <- ifelse(censor_attr == -Inf, -Inf, fit$y)
    y_upper <- ifelse(censor_attr == -Inf, fit$y, censor_attr)

    sig <- tryCatch(fit$family$getTheta(TRUE), error = function(e) sqrt(fit$sig2))

    u_L <- (y_lower - eta) / sig
    u_U <- (y_upper - eta) / sig

    l_prime <- numeric(length(eta))
    W_diag  <- numeric(length(eta))

    idx_exact <- which(y_lower == y_upper)
    if (length(idx_exact) > 0) {
      l_prime[idx_exact] <- u_L[idx_exact] / sig
      W_diag[idx_exact]  <- 1 / (sig^2)
    }

    idx_rc <- which(y_upper == Inf)
    if (length(idx_rc) > 0) {
      u <- u_L[idx_rc]
      h_u <- exp(dnorm(u, log = TRUE) - pnorm(u, lower.tail = FALSE, log.p = TRUE))
      l_prime[idx_rc] <- h_u / sig
      W_diag[idx_rc]  <- h_u * (h_u - u) / (sig^2)
    }

    idx_lc <- which(y_lower == -Inf)
    if (length(idx_lc) > 0) {
      u <- u_U[idx_lc]
      g_u <- exp(dnorm(u, log = TRUE) - pnorm(u, log.p = TRUE))
      l_prime[idx_lc] <- -g_u / sig
      W_diag[idx_lc]  <- g_u * (u + g_u) / (sig^2)
    }

    idx_int <- which(y_lower != y_upper & y_upper != Inf & y_lower != -Inf)
    if (length(idx_int) > 0) {
      u1 <- u_L[idx_int]
      u2 <- u_U[idx_int]
      p_diff <- pmax(pnorm(u2) - pnorm(u1), 1e-15)
      d1 <- dnorm(u1)
      d2 <- dnorm(u2)
      term1 <- (d1 - d2) / p_diff
      term2 <- (u1 * d1 - u2 * d2) / p_diff
      l_prime[idx_int] <- term1 / sig
      W_diag[idx_int]  <- (term1^2 - term2) / (sig^2)
    }

    pseudo_response <- eta + l_prime / W_diag
    phi0 <- 1.0
    valid_idx <- !is.na(pseudo_response) & W_diag > 1e-12
  } else if (inherits(fit$family, "extended.family") &&
             is.function(fit$family$Dd)) {
    eta <- as.numeric(fit$linear.predictors)
    y <- fit$y
    wt <- if (is.null(fit$prior.weights)) {
      rep(1, length(y))
    } else {
      as.numeric(fit$prior.weights)
    }
    work <- gammfast_working(
      fit$family, y, eta, wt
    )
    pseudo_response <- work$z
    W_diag <- work$w
    phi0 <- fit$sig2
    if (is.null(phi0) || length(phi0) != 1L ||
        !is.numeric(phi0) || !is.finite(phi0) || phi0 <= 0) {
      phi0 <- 1.0
    }
    phi0 <- as.numeric(phi0)
    valid_idx <- is.finite(pseudo_response) & is.finite(W_diag) &
      W_diag > 1e-12
  } else {
    eta <- fit$linear.predictors
    y <- fit$y
    wt <- as.numeric(fit$prior.weights)
    work <- gammfast_working(
      fit$family, y, eta, wt
    )
    pseudo_response <- work$z
    W_diag <- work$w
    phi0 <- fit$sig2
    if (is.null(phi0) || length(phi0) != 1L ||
        !is.numeric(phi0) || !is.finite(phi0) || phi0 <= 0) {
      phi0 <- 1.0
    }
    phi0 <- as.numeric(phi0)
    valid_idx <- NULL
  }

  list(
    pseudo_response = pseudo_response,
    V_phi           = phi0 / W_diag,
    phi0            = phi0,
    valid_idx       = valid_idx
  )
}

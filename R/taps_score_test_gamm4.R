gamm4_allowed_families <- c(
  "gaussian", "binomial", "poisson", "Gamma", "inverse.gaussian"
)

gamm4_is_negative_binomial <- function(family_name) {
  grepl("^Negative Binomial\\(", family_name)
}

extract_pseudo_response_gamm4 <- function(fit) {
  g <- fit$gam
  mer <- fit$mer
  family_name <- g$family$family

  is_negative_binomial <- gamm4_is_negative_binomial(family_name)
  if (!family_name %in% gamm4_allowed_families && !is_negative_binomial) {
    stop(
      paste0(
        "Unsupported gamm4 family for TAPS: ", family_name, ". Supported " ,
        "families are gaussian, binomial, poisson, Gamma, inverse.gaussian, " ,
        "and MASS::negative.binomial(theta = ...)."
      )
    )
  }
  if (!inherits(mer, "merMod")) {
    stop("fit$mer must be an lme4 mixed-model object.")
  }

  y <- as.numeric(g$y)
  mu <- as.numeric(lme4::getME(mer, "mu"))
  eta <- g$family$linkfun(mu)
  g_prime_mu <- 1 / g$family$mu.eta(eta)
  pseudo_response <- eta + (y - mu) * g_prime_mu
  W_diag <- as.numeric(g$weights)
  phi0 <- as.numeric(g$sig2)

  if (length(y) != length(mu) || length(y) != length(W_diag)) {
    stop("gamm4 response, fitted mean, and working weights have incompatible lengths.")
  }
  if (length(phi0) != 1L || !is.finite(phi0) || phi0 <= 0) {
    stop("gamm4 fitted scale must be a positive finite scalar.")
  }
  if (any(!is.finite(pseudo_response)) || any(!is.finite(W_diag)) ||
      any(W_diag <= 0)) {
    stop("gamm4 working response or weights are not finite and positive.")
  }

  list(
    pseudo_response = pseudo_response,
    V_phi = phi0 / W_diag,
    phi0 = phi0
  )
}

extract_random_block_gamm4 <- function(fit) {
  g <- fit$gam
  mer <- fit$mer
  Zt <- lme4::getME(mer, "Zt")
  Lambdat <- lme4::getME(mer, "Lambdat")
  Gp <- lme4::getME(mer, "Gp")
  cnms <- lme4::getME(mer, "cnms")
  smooth_names <- unique(unlist(
    lapply(g$smooth, function(s) s$lmer.name),
    use.names = FALSE
  ))
  re_names <- setdiff(names(cnms), smooth_names)

  if (length(re_names) == 0L) {
    return(list(
      Zt = Matrix::Matrix(0, 0, ncol(Zt), sparse = TRUE),
      Lambdat = Matrix::Matrix(0, 0, 0, sparse = TRUE),
      scale = as.numeric(g$sig2)
    ))
  }

  re_group <- match(re_names, names(cnms))
  re_index <- unlist(lapply(
    re_group,
    function(i) seq.int(Gp[i] + 1L, Gp[i + 1L])
  ), use.names = FALSE)

  list(
    Zt = Zt[re_index, , drop = FALSE],
    Lambdat = Lambdat[re_index, re_index, drop = FALSE],
    scale = as.numeric(g$sig2)
  )
}

gamm4_v0_inv_apply <- function(V_phi, Zt, Lambdat, scale) {
  if (nrow(Zt) == 0L) {
    return(function(v) {
      v_is_matrix <- is.matrix(v)
      v_mat <- if (v_is_matrix) v else matrix(v, ncol = 1)
      out <- v_mat / V_phi
      if (v_is_matrix) out else as.vector(out)
    })
  }

  F <- sqrt(scale) * Matrix::t(Lambdat %*% Zt)
  R_inv <- Matrix::Diagonal(x = 1 / V_phi)
  R_inv_F <- R_inv %*% F
  K <- Matrix::Diagonal(n = ncol(F)) + Matrix::crossprod(F, R_inv_F)
  K <- (K + Matrix::t(K)) / 2

  function(v) {
    v_is_matrix <- is.matrix(v)
    v_mat <- if (v_is_matrix) v else matrix(v, ncol = 1)
    R_inv_v <- v_mat / V_phi
    F_R_inv_v <- Matrix::crossprod(F, R_inv_v)
    K_inv_F_R_inv_v <- Matrix::solve(K, F_R_inv_v)
    out <- R_inv_v - R_inv_F %*% K_inv_F_R_inv_v
    out <- as.matrix(out)
    if (v_is_matrix) out else as.vector(out)
  }
}

expand_penalty_to_coefficients <- function(S_matrix, n_coef, smooth_label) {
  S_matrix <- as.matrix(S_matrix)
  if (nrow(S_matrix) != ncol(S_matrix)) {
    stop("Penalty matrix for '", smooth_label, "' must be square.")
  }
  if (nrow(S_matrix) == n_coef) return(S_matrix)
  if (nrow(S_matrix) == 0L || n_coef %% nrow(S_matrix) != 0L) {
    stop(
      "Penalty dimension for '", smooth_label,
      "' is incompatible with its model-matrix columns."
    )
  }

  n_block <- n_coef %/% nrow(S_matrix)
  kronecker(diag(n_block), S_matrix)
}

gamm4_scaled_penalty <- function(s, sp, phi0, n_coef) {
  if (is.null(s$first.sp) || is.null(s$last.sp)) {
    stop("Penalized smooth has no smoothing-parameter index.")
  }
  sp_idx <- seq.int(s$first.sp, s$last.sp)
  if (length(sp_idx) != length(s$S)) {
    stop("Penalized smooth has inconsistent penalty and smoothing-parameter indexes.")
  }
  sp_value <- sp[sp_idx]
  if (anyNA(sp_value)) {
    stop("Penalized smooth has missing smoothing-parameter values.")
  }
  Reduce(`+`, Map(
    function(S_matrix, value) {
      expand_penalty_to_coefficients(S_matrix, n_coef, s$label) *
        value / phi0
    },
    s$S, sp_value
  ))
}

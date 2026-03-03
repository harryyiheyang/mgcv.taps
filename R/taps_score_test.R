#' Score Test for a Specified Smooth Term in an mgcv GAM/BAM Model
#'
#' Computes a score test for the penalized (smooth) component of a selected
#' smooth term in a fitted \code{gam} or \code{bam} object from the \pkg{mgcv}
#' package. The remaining smooth terms are treated as nuisance components and
#' are profiled out via a structured covariance approximation. P-values can be
#' computed by several methods for quadratic forms in normal variables.
#'
#' @param fit A fitted \code{gam} or \code{bam} model object from \pkg{mgcv}.
#' @param test.component Integer. Index of the smooth term to be tested. Default is \code{1}.
#' @param null.tol Numeric. Row-norm threshold used to detect null-space basis
#'   columns of the penalty matrix when \code{getA} is unavailable. Default is \code{1e-10}.
#' @param method Character. Method for computing the p-value of the score
#'   statistic. One of \code{"davies"} (default), \code{"liu"}, \code{"hbe"},
#'   \code{"wood"}, \code{"satterthwaite"}, or \code{"saddlepoint"}.
#' @param max_eps Numeric. Absolute error tolerance passed to
#'   \code{CompQuadForm::davies}. Default is \code{1e-8}.
#' @param max_iter Integer. Maximum number of integration steps passed to
#'   \code{CompQuadForm::davies}. Default is \code{1e5}.
#' @param K_per_strata Integer. Passed to \code{extract_pseudo_response}. Default is \code{20}.
#' @param eps_delta Numeric. Tolerance passed to \code{extract_pseudo_response}. Default is \code{1e-12}.
#' @param eps_prob Numeric. Tolerance passed to \code{extract_pseudo_response}. Default is \code{1e-12}.
#' @param eps_muprime Numeric. Tolerance passed to \code{extract_pseudo_response}. Default is \code{1e-12}.
#' @param eps_mu Numeric. Tolerance passed to \code{extract_pseudo_response}. Default is \code{1e-12}.
#' @param n_threads Integer. Number of threads passed to \code{extract_pseudo_response}. Default is \code{1}.
#'
#' @return A \code{data.table} with three columns:
#'   \describe{
#'     \item{smooth.term}{Label of the tested smooth term.}
#'     \item{smooth.pvalue}{Score test p-value.}
#'     \item{method}{The p-value method used.}
#'   }
#'
#' @importFrom mgcv vcov.gam
#' @importFrom data.table data.table rbindlist
#' @importFrom ordgam testStat
#' @importFrom Matrix bdiag
#' @import CppMatrix
#' @importFrom survey pchisqsum
#' @import CompQuadForm
#'
#' @examples
#' \dontrun{
#' library(mgcv)
#' fit <- gam(y ~ s(x0) + s(x1), data = dat, method = "REML")
#' taps_score_test(fit, test.component = 1)
#' }
#'
#' @export
taps_score_test <- function(fit, test.component = 1, null.tol = 1e-10,
                            method = "davies", max_eps = 1e-8, max_iter = 1e5,
                            K_per_strata = 20, eps_delta = 1e-12,
                            eps_prob = 1e-12, eps_muprime = 1e-12,
                            eps_mu = 1e-12, n_threads = 1) {
  if (!inherits(fit, "gam")) stop("fit must be a 'gam' or 'bam' object.")

  res <- extract_pseudo_response(fit,
                                 K_per_strata = K_per_strata, eps_delta = eps_delta,
                                 eps_prob = eps_prob, eps_muprime = eps_muprime,
                                 eps_mu = eps_mu, n_threads = n_threads)

  pseudo_response <- res$pseudo_response
  V_phi           <- res$V_phi
  phi0            <- res$phi0

  beta         <- fit$coefficients
  X            <- predict(fit, newdata = fit$model, type = "lpmatrix")
  smooth_terms <- fit$smooth
  p            <- length(smooth_terms)
  phivec       <- fit$sp

  if (!is.null(res$valid_idx)) {
    idx             <- res$valid_idx
    if (sum(idx) == 0) stop("No valid observations for testing")
    pseudo_response <- pseudo_response[idx]
    V_phi           <- V_phi[idx]
    X               <- X[idx, , drop = FALSE]
  }

  smooth_index_list <- list()
  random_index_list <- list()
  S_list            <- list()

  for (i in 1:p) {
    s       <- smooth_terms[[i]]
    indices <- s$first.para:s$last.para
    reported_null_dim <- s$null.space.dim
    S_matrix <- s$S[[1]]

    if (is.null(s$getA) == 0) {
      detected_null_indices <- indices[1:reported_null_dim]
    } else {
      col_norms <- apply(S_matrix, 2, function(x) sqrt(sum(x^2)))
      detected_null_indices <- which(col_norms < null.tol)
    }

    smooth_indices <- setdiff(indices, detected_null_indices)

    if (i != test.component) {
      smooth_index_list[[i]] <- indices
      random_index_list[[i]] <- smooth_indices
      S_list[[i]]            <- S_matrix * phivec[i] / phi0
    }

    if (i == test.component) {
      Bj       <- X[, indices]
      Thetaj   <- matrixGeneralizedInverse(S_matrix / norm(S_matrix, "f"))
      Gj_apply <- function(v) {
        matrixVectorMultiply(Bj, matrixVectorMultiply(Thetaj, matrixVectorMultiply(t(Bj), v)))
      }
      random_index_list[[i]] <- smooth_indices
    }
  }

  S_list            <- Filter(Negate(is.null), S_list)
  smooth_index_list <- Filter(Negate(is.null), smooth_index_list)
  random_index_list <- Filter(Negate(is.null), random_index_list)
  smooth_index_vec  <- do.call(c, smooth_index_list)
  random_index_vec  <- do.call(c, random_index_list)
  fixed_index_vec   <- setdiff(seq_len(ncol(X)), random_index_vec)

  A        <- as.matrix(X[, fixed_index_vec])
  alpha    <- beta[fixed_index_vec]
  B_extend <- X[, smooth_index_vec]
  S_All    <- as.matrix(bdiag(S_list))
  XtX      <- matrixMultiply(t(B_extend), B_extend * (1 / V_phi))
  C        <- matrixInverse(XtX + S_All)

  Vinv_apply <- function(v) {
    if (is.matrix(v)) {
      return(apply(v, 2, function(col) Vinv_apply(col)))
    }
    part1   <- v / V_phi
    Bt_v    <- matrixVectorMultiply(t(B_extend), part1)
    C_Bt_v  <- matrixVectorMultiply(C, Bt_v)
    part2   <- matrixVectorMultiply(B_extend, C_Bt_v)
    part1 - part2 / V_phi
  }

  Vinv_A      <- sapply(1:ncol(A), function(j) Vinv_apply(A[, j]))
  XtVinvX     <- matrixMultiply(t(A), Vinv_A)
  XtVinvX_inv <- matrixGeneralizedInverse(XtVinvX)

  P_apply <- function(v) {
    Vinv_v       <- Vinv_apply(v)
    AVinv_v      <- matrixVectorMultiply(t(A), Vinv_v)
    solve_middle <- matrixVectorMultiply(XtVinvX_inv, AVinv_v)
    Vinv_v - matrixVectorMultiply(Vinv_A, solve_middle)
  }

  if (method == "satterthwaite") {
    error  <- pseudo_response - matrixVectorMultiply(A, alpha)
    r      <- Vinv_apply(error)
    Gj_r   <- Gj_apply(r)
    u      <- max(0,sum(r * Gj_r) / 2)

    q  <- ncol(Bj)
    Cj <- matrix(0, q, q)
    for (i in 1:q) {
      Pi <- P_apply(Bj[, i])
      for (j in i:q) {
        Cj[i, j] <- sum(Bj[, j] * Pi)
        if (i != j) Cj[j, i] <- Cj[i, j]
      }
    }
    PGj  <- matrixMultiply(Cj, Thetaj)
    PGj2 <- matrixMultiply(PGj, PGj)
    e    <- sum(diag(PGj)) / 2
    h    <- sum(diag(PGj2)) / 2
    kappa <- h / (2 * e)
    nu    <- 2 * e^2 / h
    pv    <- pchisq(u / kappa, df = nu, lower.tail = FALSE)
  } else {
    error  <- pseudo_response - matrixVectorMultiply(A, alpha)
    r      <- P_apply(error)
    Gj_r   <- Gj_apply(r)
    u      <- max(0,sum(r * Gj_r))
    q      <- ncol(Bj)
    eig_theta <- matrixsqrt(Thetaj)
    Theta_sqrt <- eig_theta$w
    N      <- sapply(1:q, function(i) P_apply(Bj[, i]))
    BtPB   <- crossprod(Bj, N)
    Q_small <- matrixListProduct(list(Theta_sqrt, BtPB, Theta_sqrt))
    lambda  <- eigen(Q_small, symmetric = TRUE, only.values = TRUE)$values
    lambda  <- lambda[lambda > 1e-15]

    pv <- switch(method,
                 liu        = compute_liu_pvalue(q = u, lambda = lambda),
                 hbe        = compute_hbe_pvalue(q = u, lambda = lambda),
                 wood       = compute_wood_pvalue(q = u, lambda = lambda),
                 saddlepoint = survey::pchisqsum(
                   x = u, df = rep(1, length(lambda)), a = lambda,
                   lower.tail = FALSE, method = "saddlepoint"
                 ),
                 davies = {
                   davies_res <- CompQuadForm::davies(q = u, lambda = lambda,
                                                      lim = max_iter, acc = max_eps)
                   if (davies_res$Qq <= 0 || davies_res$Qq > 1.0) {
                     compute_liu_pvalue(u, lambda)
                   } else {
                     davies_res$Qq
                   }
                 },
                 stop("method must be one of 'satterthwaite', 'davies', 'liu', 'hbe', 'wood', 'saddlepoint'")
    )
  }

  data.table(
    smooth.term    = smooth_terms[[test.component]]$label,
    smooth.pvalue  = pv,
    method         = method
  )
}

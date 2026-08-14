gammfast_laplace_penalty_setup <- function(G0) {
  p <- ncol(G0$X)
  components <- vector("list", length(G0$S))
  for (j in seq_along(G0$S)) {
    S <- matrix(0, p, p)
    ii <- G0$off[j] + seq_len(nrow(G0$S[[j]])) - 1L
    S[ii, ii] <- G0$S[[j]]
    components[[j]] <- (S + t(S)) / 2
  }
  if (!length(components)) {
    return(list(
      fixed_vectors = diag(p), penalized_vectors = matrix(0, p, 0),
      components = list()
    ))
  }
  S <- matrix(0, p, p)
  for (Sj in components) S <- S + Sj / max(1, norm(Sj, "F"))
  E <- eigen((S + t(S)) / 2, symmetric = TRUE)
  tol <- max(1, max(abs(E$values))) * 1e-9
  positive <- E$values > tol
  U1 <- E$vectors[, positive, drop = FALSE]
  U0 <- E$vectors[, !positive, drop = FALSE]
  for (j in seq_len(ncol(U0))) {
    ii <- which.max(abs(U0[, j]))
    if (U0[ii, j] < 0) U0[, j] <- -U0[, j]
  }
  for (j in seq_len(ncol(U1))) {
    ii <- which.max(abs(U1[, j]))
    if (U1[ii, j] < 0) U1[, j] <- -U1[, j]
  }
  R <- lapply(components, function(Sj) {
    A <- crossprod(U1, Sj %*% U1)
    (A + t(A)) / 2
  })
  list(fixed_vectors = U0, penalized_vectors = U1, components = R)
}

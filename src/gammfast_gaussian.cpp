// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]

#include <RcppArmadillo.h>
#ifdef _OPENMP
#include <omp.h>
#endif

using namespace Rcpp;

// [[Rcpp::export]]
List gammfast_gaussian_cache(const arma::mat& A,
                             const arma::mat& B,
                             const arma::ivec& id,
                             int n_threads = 1) {
  const int n = A.n_rows;
  const int d = A.n_cols;
  const int q = B.n_cols;
  if (B.n_rows != static_cast<arma::uword>(n) ||
      id.n_elem != static_cast<arma::uword>(n)) {
    stop("A, B, and id must have matching rows.");
  }
  if (id.min() < 1) stop("id must contain positive consecutive integers.");
  const int ng = id.max();
  std::vector<std::vector<arma::uword>> rows(ng);
  for (int i = 0; i < n; ++i) rows[id[i] - 1].push_back(i);
  for (int g = 0; g < ng; ++g) {
    if (rows[g].empty()) stop("id must contain positive consecutive integers.");
  }

  arma::mat AtA = A.t() * A;
  arma::cube BtB(q, q, ng, arma::fill::zeros);
  arma::cube BtA(q, d, ng, arma::fill::zeros);
  if (n_threads < 1) n_threads = 1;
  const int nt = std::min(n_threads, ng);

#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(nt)
#endif
  for (int g = 0; g < ng; ++g) {
    arma::uvec ii(rows[g]);
    arma::mat Bg = B.rows(ii);
    BtB.slice(g) = Bg.t() * Bg;
    BtA.slice(g) = Bg.t() * A.rows(ii);
  }

  return List::create(
    _["AtA"] = AtA,
    _["BtB"] = BtB,
    _["BtA"] = BtA,
    _["n_group"] = ng
  );
}

// [[Rcpp::export]]
List gammfast_gaussian_crossprod_cached(const arma::mat& AtA,
                                        const arma::cube& BtB,
                                        const arma::cube& BtA,
                                        const arma::mat& H,
                                        int n_threads = 1) {
  const int d = AtA.n_rows;
  const int q = H.n_rows;
  const int ng = BtB.n_slices;
  if (AtA.n_cols != static_cast<arma::uword>(d)) {
    stop("AtA must be square.");
  }
  if (H.n_cols != static_cast<arma::uword>(q) ||
      BtB.n_rows != static_cast<arma::uword>(q) ||
      BtB.n_cols != static_cast<arma::uword>(q) ||
      BtA.n_rows != static_cast<arma::uword>(q) ||
      BtA.n_cols != static_cast<arma::uword>(d) ||
      BtA.n_slices != static_cast<arma::uword>(ng)) {
    stop("Cached cross-product dimensions are inconsistent.");
  }

  arma::mat H_inv;
  if (!arma::inv_sympd(H_inv, arma::symmatu(H))) {
    stop("H must be positive definite.");
  }
  if (n_threads < 1) n_threads = 1;
  const int nt = std::min(n_threads, ng);
  std::vector<arma::mat> correction_thread(
    nt, arma::mat(d, d, arma::fill::zeros)
  );
  std::vector<double> logdet_thread(nt, 0.0);
  arma::ivec ok(ng, arma::fill::ones);

#ifdef _OPENMP
#pragma omp parallel num_threads(nt)
#endif
  {
    int thread_id = 0;
#ifdef _OPENMP
    thread_id = omp_get_thread_num();
#pragma omp for schedule(static)
#endif
    for (int g = 0; g < ng; ++g) {
      arma::mat M = arma::symmatu(H_inv + BtB.slice(g));
      arma::mat M_inv;
      double logdet_M = 0.0;
      double sign_M = 0.0;
      if (!arma::inv_sympd(M_inv, M) ||
          !arma::log_det(logdet_M, sign_M, M) || sign_M <= 0) {
        ok[g] = 0;
        continue;
      }
      arma::mat E = BtA.slice(g);
      correction_thread[thread_id] += E.t() * M_inv * E;
      logdet_thread[thread_id] += logdet_M;
    }
  }

  if (arma::any(ok == 0)) {
    stop("A cached subject-level Woodbury system was not positive definite.");
  }
  double logdet_H = 0.0;
  double sign_H = 0.0;
  if (!arma::log_det(logdet_H, sign_H, arma::symmatu(H)) || sign_H <= 0) {
    stop("H must have a positive determinant.");
  }
  arma::mat total = AtA;
  double logdet_total = ng * logdet_H;
  for (int thread_id = 0; thread_id < nt; ++thread_id) {
    total -= correction_thread[thread_id];
    logdet_total += logdet_thread[thread_id];
  }

  return List::create(
    _["crossprod"] = total,
    _["logdet"] = logdet_total
  );
}

// [[Rcpp::export]]
List gammfast_gaussian_crossprod(const arma::mat& A,
                                 const arma::mat& B,
                                 const arma::ivec& id,
                                 const arma::mat& H,
                                 int n_threads = 1) {
  const int n = A.n_rows;
  const int d = A.n_cols;
  const int q = B.n_cols;
  if (B.n_rows != static_cast<arma::uword>(n) || id.n_elem != static_cast<arma::uword>(n)) {
    stop("A, B, and id must have matching rows.");
  }
  if (H.n_rows != static_cast<arma::uword>(q) || H.n_cols != static_cast<arma::uword>(q)) {
    stop("H must be a square matrix matching ncol(B).");
  }
  if (id.min() < 1) stop("id must contain positive consecutive integers.");
  const int ng = id.max();
  std::vector<std::vector<arma::uword>> rows(ng);
  for (int i = 0; i < n; ++i) rows[id[i] - 1].push_back(i);
  for (int g = 0; g < ng; ++g) {
    if (rows[g].empty()) stop("id must contain positive consecutive integers.");
  }

  arma::mat H_inv;
  if (!arma::inv_sympd(H_inv, arma::symmatu(H))) stop("H must be positive definite.");
  if (n_threads < 1) n_threads = 1;
  const int nt = std::min(n_threads, ng);
  std::vector<arma::mat> crossprod_thread(
    nt, arma::mat(d, d, arma::fill::zeros)
  );
  std::vector<double> logdet_thread(nt, 0.0);
  arma::ivec ok(ng, arma::fill::ones);

#ifdef _OPENMP
#pragma omp parallel num_threads(nt)
#endif
  {
    int thread_id = 0;
#ifdef _OPENMP
    thread_id = omp_get_thread_num();
#pragma omp for schedule(static)
#endif
    for (int g = 0; g < ng; ++g) {
    arma::uvec ii(rows[g]);
    arma::mat Ag = A.rows(ii);
    arma::mat Bg = B.rows(ii);
    arma::mat BtB = Bg.t() * Bg;
    arma::mat M = arma::symmatu(H_inv + BtB);
    arma::mat M_inv;
    double logdet_M = 0.0;
    double sign_M = 0.0;
    if (!arma::inv_sympd(M_inv, M) || !arma::log_det(logdet_M, sign_M, M) || sign_M <= 0) {
      ok[g] = 0;
      continue;
    }
    arma::mat E = Bg.t() * Ag;
      crossprod_thread[thread_id] += Ag.t() * Ag - E.t() * M_inv * E;
      logdet_thread[thread_id] += logdet_M;
    }
  }

  if (arma::any(ok == 0)) stop("A subject-level Woodbury system was not positive definite.");
  double logdet_H = 0.0;
  double sign_H = 0.0;
  if (!arma::log_det(logdet_H, sign_H, arma::symmatu(H)) || sign_H <= 0) {
    stop("H must have a positive determinant.");
  }
  arma::mat total(d, d, arma::fill::zeros);
  double logdet_total = ng * logdet_H;
  for (int thread_id = 0; thread_id < nt; ++thread_id) {
    total += crossprod_thread[thread_id];
    logdet_total += logdet_thread[thread_id];
  }

  return List::create(
    _["crossprod"] = total,
    _["logdet"] = logdet_total
  );
}

// [[Rcpp::export]]
List gammfast_gaussian_moments(const arma::vec& residual,
                               const arma::mat& B,
                               const arma::ivec& id,
                               const arma::mat& G,
                               double sigma2,
                               int n_threads = 1) {
  const int n = residual.n_elem;
  const int q = B.n_cols;
  if (B.n_rows != static_cast<arma::uword>(n) || id.n_elem != static_cast<arma::uword>(n)) {
    stop("residual, B, and id must have matching rows.");
  }
  if (G.n_rows != static_cast<arma::uword>(q) || G.n_cols != static_cast<arma::uword>(q)) {
    stop("G must be a square matrix matching ncol(B).");
  }
  if (!R_finite(sigma2) || sigma2 <= 0) stop("sigma2 must be positive and finite.");
  if (id.min() < 1) stop("id must contain positive consecutive integers.");
  const int ng = id.max();
  std::vector<std::vector<arma::uword>> rows(ng);
  for (int i = 0; i < n; ++i) rows[id[i] - 1].push_back(i);
  for (int g = 0; g < ng; ++g) {
    if (rows[g].empty()) stop("id must contain positive consecutive integers.");
  }

  arma::mat G_inv;
  if (!arma::inv_sympd(G_inv, arma::symmatu(G))) stop("G must be positive definite.");
  arma::mat u(ng, q, arma::fill::zeros);
  arma::cube delta(q, q, ng, arma::fill::zeros);
  if (n_threads < 1) n_threads = 1;
  const int nt = std::min(n_threads, ng);
  std::vector<arma::mat> moment_thread(
    nt, arma::mat(q, q, arma::fill::zeros)
  );
  std::vector<double> rss_thread(nt, 0.0);
  arma::ivec ok(ng, arma::fill::ones);

#ifdef _OPENMP
#pragma omp parallel num_threads(nt)
#endif
  {
    int thread_id = 0;
#ifdef _OPENMP
    thread_id = omp_get_thread_num();
#pragma omp for schedule(static)
#endif
    for (int g = 0; g < ng; ++g) {
    arma::uvec ii(rows[g]);
    arma::mat Bg = B.rows(ii);
    arma::vec rg = residual.elem(ii);
    arma::mat BtB = Bg.t() * Bg;
    arma::mat D;
    if (!arma::inv_sympd(D, arma::symmatu(G_inv + BtB / sigma2))) {
      ok[g] = 0;
      continue;
    }
    arma::vec ug = D * (Bg.t() * rg) / sigma2;
    arma::vec eg = rg - Bg * ug;
    u.row(g) = ug.t();
    delta.slice(g) = D;
      moment_thread[thread_id] += ug * ug.t() + D;
      rss_thread[thread_id] += arma::dot(eg, eg) + arma::trace(D * BtB);
    }
  }

  if (arma::any(ok == 0)) stop("A subject-level posterior precision was not positive definite.");
  arma::mat moment_sum(q, q, arma::fill::zeros);
  double rss_sum = 0.0;
  for (int thread_id = 0; thread_id < nt; ++thread_id) {
    moment_sum += moment_thread[thread_id];
    rss_sum += rss_thread[thread_id];
  }

  return List::create(
    _["u"] = u,
    _["delta"] = delta,
    _["moment_sum"] = moment_sum,
    _["rss_sum"] = rss_sum
  );
}

// [[Rcpp::export]]
arma::mat gammfast_vinv_apply(const arma::mat& A,
                              const arma::mat& B,
                              const arma::ivec& id,
                              const arma::mat& G,
                              const arma::vec& R_diag,
                              int n_threads = 1) {
  const int n = A.n_rows;
  const int d = A.n_cols;
  const int q = B.n_cols;
  if (B.n_rows != static_cast<arma::uword>(n) ||
      id.n_elem != static_cast<arma::uword>(n) ||
      R_diag.n_elem != static_cast<arma::uword>(n)) {
    stop("A, B, id, and R_diag must have matching rows.");
  }
  if (G.n_rows != static_cast<arma::uword>(q) ||
      G.n_cols != static_cast<arma::uword>(q)) {
    stop("G must be a square matrix matching ncol(B).");
  }
  if (arma::any(R_diag <= 0) || !R_diag.is_finite()) {
    stop("R_diag must be positive and finite.");
  }
  if (id.min() < 1) stop("id must contain positive consecutive integers.");
  const int ng = id.max();
  std::vector<std::vector<arma::uword>> rows(ng);
  for (int i = 0; i < n; ++i) rows[id[i] - 1].push_back(i);
  for (int g = 0; g < ng; ++g) {
    if (rows[g].empty()) stop("id must contain positive consecutive integers.");
  }

  arma::mat G_inv;
  if (!arma::inv_sympd(G_inv, arma::symmatu(G))) {
    stop("G must be positive definite.");
  }
  arma::mat out(n, d, arma::fill::zeros);
  arma::ivec ok(ng, arma::fill::ones);
  if (n_threads < 1) n_threads = 1;
  const int nt = std::min(n_threads, ng);

#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(nt)
#endif
  for (int g = 0; g < ng; ++g) {
    arma::uvec ii(rows[g]);
    arma::vec root_R = arma::sqrt(R_diag.elem(ii));
    arma::mat Aw = A.rows(ii);
    arma::mat Bw = B.rows(ii);
    Aw.each_col() /= root_R;
    Bw.each_col() /= root_R;
    arma::mat D;
    if (!arma::inv_sympd(D, arma::symmatu(G_inv + Bw.t() * Bw))) {
      ok[g] = 0;
      continue;
    }
    arma::mat VwA = Aw - Bw * D * (Bw.t() * Aw);
    VwA.each_col() /= root_R;
    out.rows(ii) = VwA;
  }
  if (arma::any(ok == 0)) {
    stop("A subject-level Woodbury system was not positive definite.");
  }
  return out;
}

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

// [[Rcpp::export]]
List gammfast_re_quadratic_spectrum(const arma::mat& X,
                                    const arma::mat& B,
                                    const arma::ivec& id,
                                    const arma::mat& G,
                                    const arma::vec& R_diag,
                                    const arma::mat& penalty,
                                    const arma::mat& random_effects,
                                    int n_threads = 1) {
  const int n = X.n_rows;
  const int p = X.n_cols;
  const int k = B.n_cols;
  if (B.n_rows != static_cast<arma::uword>(n) ||
      id.n_elem != static_cast<arma::uword>(n) ||
      R_diag.n_elem != static_cast<arma::uword>(n)) {
    stop("X, B, id, and R_diag must have matching rows.");
  }
  if (G.n_rows != static_cast<arma::uword>(k) ||
      G.n_cols != static_cast<arma::uword>(k)) {
    stop("G must be a square matrix matching ncol(B).");
  }
  if (penalty.n_rows != static_cast<arma::uword>(p) ||
      penalty.n_cols != static_cast<arma::uword>(p)) {
    stop("penalty must be a square matrix matching ncol(X).");
  }
  if (!X.is_finite() || !B.is_finite() || !G.is_finite() ||
      !penalty.is_finite() || !random_effects.is_finite() ||
      !R_diag.is_finite() ||
      arma::any(R_diag <= 0)) {
    stop("The quadratic-spectrum inputs must be finite with positive R_diag.");
  }
  if (id.min() < 1) stop("id must contain positive consecutive integers.");
  const int ng = id.max();
  if (random_effects.n_rows != static_cast<arma::uword>(ng) ||
      random_effects.n_cols != static_cast<arma::uword>(k)) {
    stop("random_effects must have one row per subject and ncol(B) columns.");
  }
  std::vector<std::vector<arma::uword>> rows(ng);
  for (int i = 0; i < n; ++i) rows[id[i] - 1].push_back(i);
  for (int g = 0; g < ng; ++g) {
    if (rows[g].empty()) stop("id must contain positive consecutive integers.");
  }

  arma::mat G_inv;
  if (!arma::inv_sympd(G_inv, arma::symmatu(G))) {
    stop("G must be positive definite.");
  }
  const int q = ng * k;
  arma::mat XtWZ(p, q, arma::fill::zeros);
  arma::cube ZtWZ(k, k, ng, arma::fill::zeros);
  if (n_threads < 1) n_threads = 1;
  const int nt = std::min(n_threads, ng);
  std::vector<arma::mat> XtWX_thread(
    nt, arma::mat(p, p, arma::fill::zeros)
  );

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
      arma::vec root_R = arma::sqrt(R_diag.elem(ii));
      arma::mat Xw = X.rows(ii);
      arma::mat Bw = B.rows(ii);
      Xw.each_col() /= root_R;
      Bw.each_col() /= root_R;
      XtWX_thread[thread_id] += Xw.t() * Xw;
      XtWZ.cols(g * k, (g + 1) * k - 1) = Xw.t() * Bw;
      ZtWZ.slice(g) = Bw.t() * Bw;
    }
  }

  arma::mat XtWX(p, p, arma::fill::zeros);
  for (int thread_id = 0; thread_id < nt; ++thread_id) {
    XtWX += XtWX_thread[thread_id];
  }
  arma::mat mean_precision = arma::symmatu(XtWX + arma::symmatu(penalty));
  arma::mat mean_covariance;
  if (!arma::inv_sympd(mean_covariance, mean_precision)) {
    stop("The fitted global coefficient precision is not positive definite.");
  }

  arma::mat ZtWZ_big(q, q, arma::fill::zeros);
  for (int g = 0; g < ng; ++g) {
    ZtWZ_big.submat(g * k, g * k, (g + 1) * k - 1, (g + 1) * k - 1) =
      ZtWZ.slice(g);
  }
  arma::mat tested_information = arma::symmatu(
    ZtWZ_big - XtWZ.t() * mean_covariance * XtWZ
  );

  arma::mat normal(p + q, p + q, arma::fill::zeros);
  normal.submat(0, 0, p - 1, p - 1) = XtWX;
  normal.submat(0, p, p - 1, p + q - 1) = XtWZ;
  normal.submat(p, 0, p + q - 1, p - 1) = XtWZ.t();
  normal.submat(p, p, p + q - 1, p + q - 1) = ZtWZ_big;

  arma::mat full_penalty(p + q, p + q, arma::fill::zeros);
  full_penalty.submat(0, 0, p - 1, p - 1) = arma::symmatu(penalty);
  for (int g = 0; g < ng; ++g) {
    full_penalty.submat(
      p + g * k, p + g * k,
      p + (g + 1) * k - 1, p + (g + 1) * k - 1
    ) = G_inv;
  }
  arma::mat full_covariance;
  if (!arma::inv_sympd(
        full_covariance, arma::symmatu(normal + full_penalty))) {
    stop("The fitted joint coefficient precision is not positive definite.");
  }
  arma::mat target_rows = full_covariance.rows(p, p + q - 1);
  arma::mat target_covariance = arma::symmatu(
    target_rows * normal * target_rows.t()
  );

  arma::vec information_values;
  arma::mat information_vectors;
  if (!arma::eig_sym(
        information_values, information_vectors, tested_information)) {
    stop("The tested random-effect information eigendecomposition failed.");
  }
  const double information_scale = std::max(1.0, information_values.max());
  if (information_values.min() < -1e-8 * information_scale) {
    stop("The tested random-effect information is not positive semidefinite.");
  }
  information_values = arma::clamp(
    information_values, 0.0, arma::datum::inf
  );
  arma::mat information_root = information_vectors *
    arma::diagmat(arma::sqrt(information_values)) * information_vectors.t();
  arma::mat reference_covariance = arma::symmatu(
    information_root * target_covariance * information_root
  );
  arma::vec lambda;
  if (!arma::eig_sym(lambda, reference_covariance)) {
    stop("The random-effect reference eigendecomposition failed.");
  }
  const double lambda_scale = std::max(1.0, lambda.max());
  if (lambda.min() < -1e-8 * lambda_scale) {
    stop("The random-effect reference covariance is not positive semidefinite.");
  }
  lambda = arma::clamp(lambda, 0.0, arma::datum::inf);

  arma::vec u = arma::vectorise(random_effects.t());
  const double statistic = arma::as_scalar(u.t() * tested_information * u);

  return List::create(
    _["statistic"] = std::max(0.0, statistic),
    _["lambda"] = lambda,
    _["information_rank"] = arma::sum(
      information_values > information_scale * 1e-10
    ),
    _["n_subject"] = ng,
    _["basis_dimension"] = k
  );
}

// [[Rcpp::export]]
List gammfast_re_quadratic_random_moments(
    const arma::mat& X,
    const arma::mat& B,
    const arma::ivec& id,
    const arma::mat& G,
    const arma::vec& R_diag,
    const arma::mat& penalty,
    const arma::mat& random_effects,
    const arma::mat& probes,
    int n_threads = 1) {
  const int n = X.n_rows;
  const int p = X.n_cols;
  const int k = B.n_cols;
  if (B.n_rows != static_cast<arma::uword>(n) ||
      id.n_elem != static_cast<arma::uword>(n) ||
      R_diag.n_elem != static_cast<arma::uword>(n)) {
    stop("X, B, id, and R_diag must have matching rows.");
  }
  if (G.n_rows != static_cast<arma::uword>(k) ||
      G.n_cols != static_cast<arma::uword>(k)) {
    stop("G must be a square matrix matching ncol(B).");
  }
  if (penalty.n_rows != static_cast<arma::uword>(p) ||
      penalty.n_cols != static_cast<arma::uword>(p)) {
    stop("penalty must be a square matrix matching ncol(X).");
  }
  if (!X.is_finite() || !B.is_finite() || !G.is_finite() ||
      !penalty.is_finite() || !random_effects.is_finite() ||
      !probes.is_finite() || !R_diag.is_finite() ||
      arma::any(R_diag <= 0)) {
    stop("The random-moment inputs must be finite with positive R_diag.");
  }
  if (id.n_elem == 0 || id.min() < 1) {
    stop("id must contain positive consecutive integers.");
  }
  const int ng = id.max();
  const int q = ng * k;
  if (random_effects.n_rows != static_cast<arma::uword>(ng) ||
      random_effects.n_cols != static_cast<arma::uword>(k)) {
    stop("random_effects must have one row per subject and ncol(B) columns.");
  }
  if (probes.n_rows != static_cast<arma::uword>(q) || probes.n_cols < 1) {
    stop("probes must have n_subject * ncol(B) rows and at least one column.");
  }
  std::vector<std::vector<arma::uword>> rows(ng);
  for (int i = 0; i < n; ++i) rows[id[i] - 1].push_back(i);
  for (int g = 0; g < ng; ++g) {
    if (rows[g].empty()) stop("id must contain positive consecutive integers.");
  }

  arma::mat G_inv;
  if (!arma::inv_sympd(G_inv, arma::symmatu(G))) {
    stop("G must be positive definite.");
  }
  arma::mat XtWZ(p, q, arma::fill::zeros);
  arma::cube ZtWZ(k, k, ng, arma::fill::zeros);
  if (n_threads < 1) n_threads = 1;
  const int nt = std::min(n_threads, ng);
  std::vector<arma::mat> XtWX_thread(
    nt, arma::mat(p, p, arma::fill::zeros)
  );

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
      arma::vec root_R = arma::sqrt(R_diag.elem(ii));
      arma::mat Xw = X.rows(ii);
      arma::mat Bw = B.rows(ii);
      Xw.each_col() /= root_R;
      Bw.each_col() /= root_R;
      XtWX_thread[thread_id] += Xw.t() * Xw;
      XtWZ.cols(g * k, (g + 1) * k - 1) = Xw.t() * Bw;
      ZtWZ.slice(g) = Bw.t() * Bw;
    }
  }

  arma::mat XtWX(p, p, arma::fill::zeros);
  for (int thread_id = 0; thread_id < nt; ++thread_id) {
    XtWX += XtWX_thread[thread_id];
  }
  arma::mat A = arma::symmatu(XtWX + arma::symmatu(penalty));
  arma::mat A_inv;
  if (!arma::inv_sympd(A_inv, A)) {
    stop("The fitted global coefficient precision is not positive definite.");
  }

  arma::cube E_inv(k, k, ng, arma::fill::zeros);
  arma::mat M = A;
  for (int g = 0; g < ng; ++g) {
    arma::mat Ei_inv;
    if (!arma::inv_sympd(
          Ei_inv, arma::symmatu(ZtWZ.slice(g) + G_inv))) {
      stop("A random-effect block precision is not positive definite.");
    }
    E_inv.slice(g) = Ei_inv;
    arma::mat Bi = XtWZ.cols(g * k, (g + 1) * k - 1);
    M -= Bi * Ei_inv * Bi.t();
  }
  arma::mat M_inv;
  if (!arma::inv_sympd(M_inv, arma::symmatu(M))) {
    stop("The marginalized mean precision is not positive definite.");
  }

  auto apply_block = [&](const arma::vec& v, const arma::cube& blocks) {
    arma::vec out(q, arma::fill::zeros);
    for (int g = 0; g < ng; ++g) {
      out.subvec(g * k, (g + 1) * k - 1) =
        blocks.slice(g) * v.subvec(g * k, (g + 1) * k - 1);
    }
    return out;
  };

  auto apply_hinv = [&](const arma::vec& rb, const arma::vec& ru) {
    arma::vec Eru = apply_block(ru, E_inv);
    arma::vec beta = M_inv * (rb - XtWZ * Eru);
    arma::vec u = apply_block(ru - XtWZ.t() * beta, E_inv);
    return std::make_pair(beta, u);
  };

  auto apply_J = [&](const arma::vec& v) {
    arma::vec Dv = apply_block(v, ZtWZ);
    arma::vec out = Dv - XtWZ.t() * (A_inv * (XtWZ * v));
    return out;
  };

  auto apply_Ve = [&](const arma::vec& v) {
    arma::vec zero(p, arma::fill::zeros);
    std::pair<arma::vec, arma::vec> c1 = apply_hinv(zero, v);
    arma::vec nb = XtWX * c1.first + XtWZ * c1.second;
    arma::vec nu = XtWZ.t() * c1.first + apply_block(c1.second, ZtWZ);
    std::pair<arma::vec, arma::vec> c2 = apply_hinv(nb, nu);
    arma::vec out = c2.second;
    return out;
  };

  arma::vec u = arma::vectorise(random_effects.t());
  const double statistic = arma::dot(u, apply_J(u));
  const int n_probe = probes.n_cols;
  arma::mat probe_moments(n_probe, 4, arma::fill::zeros);

#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(nt)
#endif
  for (int j = 0; j < n_probe; ++j) {
    arma::vec z = probes.col(j);
    arma::vec v = z;
    for (int r = 0; r < 4; ++r) {
      v = apply_J(apply_Ve(v));
      probe_moments(j, r) = arma::dot(z, v);
    }
  }

  return List::create(
    _["statistic"] = std::max(0.0, statistic),
    _["moments"] = arma::mean(probe_moments, 0).t(),
    _["probe_moments"] = probe_moments,
    _["n_probe"] = n_probe,
    _["n_subject"] = ng,
    _["basis_dimension"] = k
  );
}

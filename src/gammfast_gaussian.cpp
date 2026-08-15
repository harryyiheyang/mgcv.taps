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
arma::mat gammfast_gaussian_blup_cached(const arma::cube& BtB,
                                        const arma::cube& BtA,
                                        const arma::mat& G,
                                        const arma::vec& beta,
                                        double sigma2) {
  const int q = G.n_rows;
  const int ng = BtB.n_slices;
  const int p = beta.n_elem;
  if (G.n_cols != static_cast<arma::uword>(q) ||
      BtB.n_rows != static_cast<arma::uword>(q) ||
      BtB.n_cols != static_cast<arma::uword>(q) ||
      BtA.n_rows != static_cast<arma::uword>(q) ||
      BtA.n_cols != static_cast<arma::uword>(p + 1) ||
      BtA.n_slices != static_cast<arma::uword>(ng)) {
    stop("Cached BLUP dimensions are inconsistent.");
  }
  if (!R_finite(sigma2) || sigma2 <= 0) {
    stop("sigma2 must be positive and finite.");
  }
  arma::mat Ginv;
  if (!arma::inv_sympd(Ginv, arma::symmatu(G))) {
    stop("G must be positive definite.");
  }
  const double root_sigma = std::sqrt(sigma2);
  arma::mat u(ng, q, arma::fill::zeros);
  for (int i = 0; i < ng; ++i) {
    const arma::mat Q = arma::symmatu(BtB.slice(i));
    const arma::mat M = BtA.slice(i).cols(0, p - 1) / root_sigma;
    const arma::vec b = BtA.slice(i).col(p) / root_sigma;
    arma::mat D;
    if (!arma::inv_sympd(D, arma::symmatu(Ginv + Q))) {
      stop("A cached subject posterior precision was not positive definite.");
    }
    const arma::vec t = b - Q * D * b - (M - Q * D * M) * beta;
    u.row(i) = (G * t).t();
  }
  return u;
}

// [[Rcpp::export]]
List gammfast_gaussian_projected_cached(const arma::mat& AtA,
                                        const arma::cube& BtB,
                                        const arma::cube& BtA,
                                        const arma::mat& G,
                                        double sigma2,
                                        const arma::ivec& covariance_group,
                                        bool fisher = false) {
  const int q = G.n_rows;
  const int ng = BtB.n_slices;
  const int p = AtA.n_rows - 1;
  if (p < 1 || AtA.n_cols != AtA.n_rows) {
    stop("AtA must be a square augmented cross-product matrix.");
  }
  if (G.n_cols != static_cast<arma::uword>(q) ||
      BtB.n_rows != static_cast<arma::uword>(q) ||
      BtB.n_cols != static_cast<arma::uword>(q) ||
      BtA.n_rows != static_cast<arma::uword>(q) ||
      BtA.n_cols != static_cast<arma::uword>(p + 1) ||
      BtA.n_slices != static_cast<arma::uword>(ng) ||
      covariance_group.n_elem != static_cast<arma::uword>(q)) {
    stop("Cached projected-moment dimensions are inconsistent.");
  }
  if (!R_finite(sigma2) || sigma2 <= 0) {
    stop("sigma2 must be positive and finite.");
  }
  if (covariance_group.min() < 1) {
    stop("covariance_group must contain positive consecutive integers.");
  }
  const int n_covariance_group = covariance_group.max();
  for (int j = 1; j <= n_covariance_group; ++j) {
    arma::uvec jj = arma::find(covariance_group == j);
    if (jj.n_elem == 0 || jj.max() - jj.min() + 1 != jj.n_elem) {
      stop("covariance_group must define non-empty contiguous blocks.");
    }
  }

  arma::mat Ginv;
  if (!arma::inv_sympd(Ginv, arma::symmatu(G))) {
    stop("G must be positive definite.");
  }
  const double root_sigma = std::sqrt(sigma2);
  arma::mat H = AtA.submat(0, 0, p - 1, p - 1) / sigma2;
  arma::vec h = AtA.submat(0, p, p - 1, p) / sigma2;
  arma::cube D(q, q, ng, arma::fill::zeros);
  arma::cube L(q, p, ng, arma::fill::zeros);
  arma::cube R(q, q, ng, arma::fill::zeros);
  arma::mat v(q, ng, arma::fill::zeros);

  for (int i = 0; i < ng; ++i) {
    const arma::mat Q = arma::symmatu(BtB.slice(i));
    const arma::mat M = BtA.slice(i).cols(0, p - 1) / root_sigma;
    const arma::vec b = BtA.slice(i).col(p) / root_sigma;
    if (!arma::inv_sympd(D.slice(i), arma::symmatu(Ginv + Q))) {
      stop("A cached subject posterior precision was not positive definite.");
    }
    L.slice(i) = M - Q * D.slice(i) * M;
    v.col(i) = b - Q * D.slice(i) * b;
    R.slice(i) = arma::symmatu(Q - Q * D.slice(i) * Q);
    H -= M.t() * D.slice(i) * M;
    h -= M.t() * D.slice(i) * b;
  }
  arma::mat Hinv;
  if (!arma::inv_sympd(Hinv, arma::symmatu(H))) {
    stop("The cached fixed-effect projected precision was not positive definite.");
  }
  const arma::vec beta = Hinv * h;
  arma::mat Hroot;
  if (!arma::chol(Hroot, Hinv, "lower")) {
    stop("The cached fixed-effect covariance was not positive definite.");
  }
  arma::mat u(ng, q, arma::fill::zeros);
  arma::mat tmat(ng, q, arma::fill::zeros);
  arma::mat moment_sum(q, q, arma::fill::zeros);
  arma::mat Csum(q, q, arma::fill::zeros);
  arma::cube C(q, q, ng, arma::fill::zeros);
  arma::cube U(q, p, ng, arma::fill::zeros);
  for (int i = 0; i < ng; ++i) {
    const arma::vec ti = v.col(i) - L.slice(i) * beta;
    U.slice(i) = L.slice(i) * Hroot;
    const arma::mat Ki = U.slice(i) * U.slice(i).t();
    C.slice(i) = arma::symmatu(R.slice(i) - Ki);
    const arma::vec ui = G * ti;
    const arma::mat posterior = arma::symmatu(G - G * C.slice(i) * G);
    tmat.row(i) = ti.t();
    u.row(i) = ui.t();
    moment_sum += ui * ui.t() + posterior;
    Csum += C.slice(i);
  }
  arma::vec theta;
  arma::vec score;
  arma::mat information;
  arma::ivec information_group;
  if (fisher) {
    struct Derivative {
      arma::uvec coordinates;
      arma::mat value;
      int information_group;
    };
    std::vector<Derivative> derivative;
    std::vector<double> theta_value;
    std::vector<int> information_value;
    int fs_information_group = 1;
    for (int g = 1; g <= n_covariance_group; ++g) {
      arma::uvec jj = arma::find(covariance_group == g);
      const int d = jj.n_elem;
      arma::mat Gg = G.submat(jj, jj);
      arma::mat Lg;
      if (!arma::chol(Lg, arma::symmatu(Gg), "lower")) {
        stop("A covariance block was not positive definite.");
      }
      const int info_g = d == 1 ? 0 : fs_information_group++;
      for (int a = 0; a < d; ++a) {
        for (int b = 0; b <= a; ++b) {
          arma::mat dL(d, d, arma::fill::zeros);
          if (a == b) {
            theta_value.push_back(std::log(Lg(a, a)));
            dL(a, a) = Lg(a, a);
          } else {
            theta_value.push_back(Lg(a, b));
            dL(a, b) = 1.0;
          }
          Derivative item;
          item.coordinates = jj;
          item.value = dL * Lg.t() + Lg * dL.t();
          item.information_group = info_g;
          derivative.push_back(item);
          information_value.push_back(info_g + 1);
        }
      }
    }
    const int m = derivative.size();
    theta.set_size(m);
    score.zeros(m);
    information.zeros(m, m);
    information_group.set_size(m);
    for (int j = 0; j < m; ++j) {
      theta[j] = theta_value[j];
      information_group[j] = information_value[j];
      const arma::uvec& jj = derivative[j].coordinates;
      const arma::mat& E = derivative[j].value;
      for (int i = 0; i < ng; ++i) {
        const arma::vec ti = tmat.row(i).t();
        const arma::vec tb = ti.elem(jj);
        score[j] += 0.5 * (
          arma::as_scalar(tb.t() * E * tb) -
          arma::trace(C.slice(i).submat(jj, jj) * E)
        );
      }
    }

    for (int info_g = 0; info_g < fs_information_group; ++info_g) {
      std::vector<int> pars;
      for (int j = 0; j < m; ++j) {
        if (derivative[j].information_group == info_g) pars.push_back(j);
      }
      if (pars.empty()) continue;
      arma::uvec coords;
      if (info_g == 0) {
        std::vector<arma::uword> scalar_coords;
        for (int g = 1; g <= n_covariance_group; ++g) {
          arma::uvec jj = arma::find(covariance_group == g);
          if (jj.n_elem == 1) scalar_coords.push_back(jj[0]);
        }
        coords = arma::uvec(scalar_coords);
      } else {
        coords = derivative[pars[0]].coordinates;
      }
      const int d = coords.n_elem;
      const int mb = pars.size();
      std::vector<arma::mat> Eb(mb, arma::mat(d, d, arma::fill::zeros));
      for (int a = 0; a < mb; ++a) {
        const int j = pars[a];
        if (info_g == 0) {
          const arma::uword coordinate = derivative[j].coordinates[0];
          const arma::uvec location = arma::find(coords == coordinate);
          Eb[a](location[0], location[0]) = derivative[j].value(0, 0);
        } else {
          Eb[a] = derivative[j].value;
        }
      }
      std::vector<arma::mat> Atilde(
        mb, arma::mat(p, p, arma::fill::zeros)
      );
      arma::mat local(mb, mb, arma::fill::zeros);
      for (int i = 0; i < ng; ++i) {
        const arma::mat Rb = R.slice(i).submat(coords, coords);
        const arma::mat Ub = U.slice(i).rows(coords);
        const arma::mat Kb = Ub * Ub.t();
        for (int a = 0; a < mb; ++a) {
          Atilde[a] += Ub.t() * Eb[a] * Ub;
        }
        if (info_g == 0) {
          for (int a = 0; a < mb; ++a) {
            const double ea = Eb[a](a, a);
            for (int b = 0; b <= a; ++b) {
              const double eb = Eb[b](b, b);
              local(a, b) += ea * eb * (
                Rb(a, b) * Rb(a, b) -
                2.0 * Rb(a, b) * Kb(a, b)
              );
            }
          }
        } else {
          for (int a = 0; a < mb; ++a) {
            for (int b = 0; b <= a; ++b) {
              local(a, b) +=
                arma::trace(Rb * Eb[b] * Rb * Eb[a]) -
                arma::trace(Rb * Eb[b] * Kb * Eb[a]) -
                arma::trace(Kb * Eb[b] * Rb * Eb[a]);
            }
          }
        }
      }
      for (int a = 0; a < mb; ++a) {
        for (int b = 0; b <= a; ++b) {
          const double value = 0.5 * (
            local(a, b) + arma::accu(Atilde[a] % Atilde[b])
          );
          information(pars[a], pars[b]) = value;
          information(pars[b], pars[a]) = value;
        }
      }
    }
  }
  return List::create(
    _["beta"] = beta, _["u"] = u, _["t"] = tmat,
    _["moment_sum"] = arma::symmatu(moment_sum),
    _["C_sum"] = arma::symmatu(Csum),
    _["mean_covariance"] = Hinv,
    _["theta"] = theta, _["score"] = score,
    _["information"] = information,
    _["information_group"] = information_group
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

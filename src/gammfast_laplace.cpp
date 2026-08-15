#include <RcppArmadillo.h>
#ifdef _OPENMP
#include <omp.h>
#endif

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]
using namespace Rcpp;

namespace {

arma::mat inverse_spd(const arma::mat& A, const char* message) {
  arma::mat out;
  const arma::mat S = arma::symmatu(A);
  if (!arma::inv_sympd(out, S)) stop(message);
  return out;
}

}  // namespace

// [[Rcpp::export]]
List gammfast_laplace_influence_cached(
    const arma::mat& X,
    const arma::mat& B,
    const arma::ivec& id,
    const arma::mat& G,
    const arma::cube& working_BtB,
    const arma::cube& working_BtA,
    const arma::mat& working_mean_covariance,
    const arma::cube& determinant_BtB,
    const arma::cube& determinant_BtA,
    const arma::mat& determinant_mean_covariance,
    const arma::vec& determinant_derivative,
    const arma::mat& u,
    double scale = 1.0,
    int n_threads = 1) {
  const int n = X.n_rows;
  const int p = X.n_cols;
  const int q = B.n_cols;
  if (B.n_rows != static_cast<arma::uword>(n) ||
      id.n_elem != static_cast<arma::uword>(n) ||
      determinant_derivative.n_elem != static_cast<arma::uword>(n)) {
    stop("Cached Laplace influence inputs must have matching rows.");
  }
  if (G.n_rows != static_cast<arma::uword>(q) ||
      G.n_cols != static_cast<arma::uword>(q) ||
      working_mean_covariance.n_rows != static_cast<arma::uword>(p) ||
      working_mean_covariance.n_cols != static_cast<arma::uword>(p) ||
      determinant_mean_covariance.n_rows != static_cast<arma::uword>(p) ||
      determinant_mean_covariance.n_cols != static_cast<arma::uword>(p) ||
      working_BtB.n_rows != static_cast<arma::uword>(q) ||
      working_BtB.n_cols != static_cast<arma::uword>(q) ||
      determinant_BtB.n_rows != static_cast<arma::uword>(q) ||
      determinant_BtB.n_cols != static_cast<arma::uword>(q) ||
      working_BtA.n_rows != static_cast<arma::uword>(q) ||
      working_BtA.n_cols < static_cast<arma::uword>(p) ||
      determinant_BtA.n_rows != static_cast<arma::uword>(q) ||
      determinant_BtA.n_cols < static_cast<arma::uword>(p)) {
    stop("Cached Laplace influence matrix dimensions are inconsistent.");
  }
  if (!R_finite(scale) || scale <= 0) {
    stop("The cached Laplace influence scale must be positive and finite.");
  }
  if (id.n_elem == 0 || id.min() < 1) {
    stop("Cached Laplace influence IDs must be positive consecutive integers.");
  }
  const int ng = id.max();
  if (working_BtB.n_slices != static_cast<arma::uword>(ng) ||
      working_BtA.n_slices != static_cast<arma::uword>(ng) ||
      determinant_BtB.n_slices != static_cast<arma::uword>(ng) ||
      determinant_BtA.n_slices != static_cast<arma::uword>(ng) ||
      u.n_rows != static_cast<arma::uword>(ng) ||
      u.n_cols != static_cast<arma::uword>(q)) {
    stop("Cached Laplace influence subject dimensions are inconsistent.");
  }
  if (!X.is_finite() || !B.is_finite() || !G.is_finite() ||
      !working_BtB.is_finite() || !working_BtA.is_finite() ||
      !working_mean_covariance.is_finite() ||
      !determinant_BtB.is_finite() || !determinant_BtA.is_finite() ||
      !determinant_mean_covariance.is_finite() ||
      !determinant_derivative.is_finite() ||
      !u.is_finite()) {
    stop("Cached Laplace influence inputs must be finite.");
  }

  std::vector<std::vector<arma::uword>> rows(ng);
  for (int i = 0; i < n; ++i) rows[id[i] - 1].push_back(i);
  for (int g = 0; g < ng; ++g) {
    if (rows[g].empty()) {
      stop("Cached Laplace influence IDs must be consecutive.");
    }
  }
  arma::mat Ginv = inverse_spd(
    G, "The cached Laplace influence covariance must be positive definite."
  );
  arma::cube working_D(q, q, ng, arma::fill::zeros);
  arma::cube determinant_D(q, q, ng, arma::fill::zeros);
  arma::mat effective(n, p, arma::fill::zeros);
  arma::vec leverage(n, arma::fill::zeros);
  if (n_threads < 1) n_threads = 1;
  const int nt = std::min(n_threads, ng);
  arma::ivec ok(ng, arma::fill::ones);

#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(nt)
#endif
  for (int g = 0; g < ng; ++g) {
    arma::uvec ii(rows[g]);
    const arma::mat Bg = B.rows(ii);
    const arma::mat determinant_M =
      determinant_BtA.slice(g).cols(0, p - 1);
    if (!arma::inv_sympd(
          working_D.slice(g),
          arma::symmatu(Ginv + working_BtB.slice(g)
        )) || !arma::inv_sympd(
          determinant_D.slice(g),
          arma::symmatu(Ginv + determinant_BtB.slice(g))
        )) {
      ok[g] = 0;
      continue;
    }
    const arma::mat BD = Bg * determinant_D.slice(g);
    effective.rows(ii) = X.rows(ii) - BD * determinant_M;
    leverage.elem(ii) = scale * arma::sum(BD % Bg, 1);
  }
  if (arma::any(ok == 0)) {
    stop("A cached Laplace influence subject Hessian was not positive definite.");
  }
  const arma::mat effective_covariance =
    effective * determinant_mean_covariance;
  leverage += arma::sum(effective_covariance % effective, 1);
  const arma::vec a = 0.5 * determinant_derivative % leverage / scale;

  arma::vec smooth_rhs = X.t() * a;
  arma::mat Ba(q, ng, arma::fill::zeros);
  std::vector<arma::vec> rhs_thread(
    nt, arma::vec(p, arma::fill::zeros)
  );
#ifdef _OPENMP
#pragma omp parallel num_threads(nt)
#endif
  {
    int tid = 0;
#ifdef _OPENMP
    tid = omp_get_thread_num();
#pragma omp for schedule(static)
#endif
    for (int g = 0; g < ng; ++g) {
      arma::uvec ii(rows[g]);
      Ba.col(g) = B.rows(ii).t() * a.elem(ii);
      const arma::mat M = working_BtA.slice(g).cols(0, p - 1);
      rhs_thread[tid] -= M.t() * working_D.slice(g) * Ba.col(g);
    }
  }
  for (int t = 0; t < nt; ++t) smooth_rhs += rhs_thread[t];
  const arma::vec smooth_influence =
    working_mean_covariance * smooth_rhs;

  std::vector<arma::mat> influence_thread(
    nt, arma::mat(q, q, arma::fill::zeros)
  );
#ifdef _OPENMP
#pragma omp parallel num_threads(nt)
#endif
  {
    int tid = 0;
#ifdef _OPENMP
    tid = omp_get_thread_num();
#pragma omp for schedule(static)
#endif
    for (int g = 0; g < ng; ++g) {
      const arma::mat M = working_BtA.slice(g).cols(0, p - 1);
      const arma::vec subject_influence = working_D.slice(g) * (
        scale * Ba.col(g) - M * smooth_influence
      );
      const arma::vec ug = u.row(g).t();
      influence_thread[tid] +=
        subject_influence * ug.t() + ug * subject_influence.t();
    }
  }
  arma::mat influence_sum(q, q, arma::fill::zeros);
  for (int t = 0; t < nt; ++t) influence_sum += influence_thread[t];
  influence_sum = arma::symmatu(influence_sum);
  return List::create(
    _["influence_sum"] = influence_sum,
    _["a"] = a,
    _["leverage"] = leverage
  );
}

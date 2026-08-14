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
List gammfast_laplace_variance_step(
    const arma::mat& X_penalized,
    const arma::mat& B,
    const arma::ivec& id,
    const arma::mat& G,
    const arma::mat& smooth_precision,
    const arma::vec& working_weight,
    const arma::vec& determinant_weight,
    const arma::vec& determinant_derivative,
    const arma::mat& u,
    int n_threads = 1) {
  const int n = B.n_rows;
  const int ns = X_penalized.n_cols;
  const int q = B.n_cols;
  if (X_penalized.n_rows != static_cast<arma::uword>(n) ||
      id.n_elem != static_cast<arma::uword>(n) ||
      working_weight.n_elem != static_cast<arma::uword>(n) ||
      determinant_weight.n_elem != static_cast<arma::uword>(n) ||
      determinant_derivative.n_elem != static_cast<arma::uword>(n)) {
    stop("Laplace variance-step inputs must have matching rows.");
  }
  if (G.n_rows != static_cast<arma::uword>(q) ||
      G.n_cols != static_cast<arma::uword>(q)) {
    stop("The Laplace variance-step covariance must match the random design.");
  }
  if (smooth_precision.n_rows != static_cast<arma::uword>(ns) ||
      smooth_precision.n_cols != static_cast<arma::uword>(ns)) {
    stop("The Laplace variance-step smooth precision is inconsistent.");
  }
  if (id.n_elem == 0 || id.min() < 1) {
    stop("Laplace variance-step IDs must be positive consecutive integers.");
  }
  const int ng = id.max();
  if (u.n_rows != static_cast<arma::uword>(ng) ||
      u.n_cols != static_cast<arma::uword>(q)) {
    stop("The Laplace variance-step BLUP matrix must have one row per UID.");
  }
  if (!X_penalized.is_finite() || !B.is_finite() || !G.is_finite() ||
      !smooth_precision.is_finite() || !working_weight.is_finite() ||
      !determinant_weight.is_finite() ||
      !determinant_derivative.is_finite() || !u.is_finite() ||
      arma::any(working_weight < 0.0) ||
      arma::any(determinant_weight < 0.0)) {
    stop("Laplace variance-step inputs must be finite with nonnegative weights.");
  }

  std::vector<std::vector<arma::uword>> rows(ng);
  for (int i = 0; i < n; ++i) rows[id[i] - 1].push_back(i);
  for (int g = 0; g < ng; ++g) {
    if (rows[g].empty()) stop("Laplace variance-step IDs must be consecutive.");
  }
  if (n_threads < 1) n_threads = 1;
  const int nt = std::min(n_threads, ng);
  const arma::mat precision = inverse_spd(
    G, "The Laplace variance-step covariance must be positive definite."
  );

  arma::mat schur = smooth_precision +
    X_penalized.t() * (X_penalized.each_col() % working_weight);
  std::vector<arma::mat> schur_thread(
    nt, arma::mat(ns, ns, arma::fill::zeros)
  );
  arma::ivec ok(ng, arma::fill::ones);
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
      const arma::mat Bg = B.rows(ii);
      const arma::mat Xpg = X_penalized.rows(ii);
      const arma::mat weighted_B = Bg.each_col() % working_weight.elem(ii);
      arma::mat C;
      const arma::mat subject_precision = arma::symmatu(
        precision + Bg.t() * weighted_B
      );
      if (!arma::inv_sympd(C, subject_precision)) {
        ok[g] = 0;
        continue;
      }
      if (ns > 0) {
        const arma::mat cross = Xpg.t() * weighted_B;
        schur_thread[tid] -= cross * C * cross.t();
      }
    }
  }
  if (arma::any(ok == 0)) {
    stop("A Laplace variance-step UID Hessian was not positive definite.");
  }
  for (int t = 0; t < nt; ++t) schur += schur_thread[t];
  arma::mat schur_inverse(ns, ns, arma::fill::zeros);
  if (ns > 0) {
    schur_inverse = inverse_spd(
      schur,
      "The Laplace variance-step mean Schur complement was not positive definite."
    );
  }

  const bool shared_determinant = arma::approx_equal(
    working_weight, determinant_weight, "absdiff", 0.0
  );
  arma::mat determinant_schur_inverse = schur_inverse;
  if (!shared_determinant) {
    arma::mat determinant_schur = smooth_precision +
      X_penalized.t() * (X_penalized.each_col() % determinant_weight);
    std::vector<arma::mat> determinant_schur_thread(
      nt, arma::mat(ns, ns, arma::fill::zeros)
    );
    ok.ones();
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
        const arma::mat Bg = B.rows(ii);
        const arma::mat Xpg = X_penalized.rows(ii);
        const arma::mat weighted_B =
          Bg.each_col() % determinant_weight.elem(ii);
        arma::mat C;
        const arma::mat subject_precision = arma::symmatu(
          precision + Bg.t() * weighted_B
        );
        if (!arma::inv_sympd(C, subject_precision)) {
          ok[g] = 0;
          continue;
        }
        if (ns > 0) {
          const arma::mat cross = Xpg.t() * weighted_B;
          determinant_schur_thread[tid] -= cross * C * cross.t();
        }
      }
    }
    if (arma::any(ok == 0)) {
      stop("A Laplace determinant UID Hessian was not positive definite.");
    }
    for (int t = 0; t < nt; ++t) {
      determinant_schur += determinant_schur_thread[t];
    }
    if (ns > 0) {
      determinant_schur_inverse = inverse_spd(
        determinant_schur,
        "The Laplace determinant mean Schur complement was not positive definite."
      );
    }
  }

  arma::vec a(n, arma::fill::zeros);
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(nt)
#endif
  for (int g = 0; g < ng; ++g) {
    arma::uvec ii(rows[g]);
    const arma::mat Bg = B.rows(ii);
    const arma::mat Xpg = X_penalized.rows(ii);
    const arma::mat weighted_B =
      Bg.each_col() % determinant_weight.elem(ii);
    const arma::mat C = inverse_spd(
      precision + Bg.t() * weighted_B,
      "A Laplace variance-step leverage Hessian was not positive definite."
    );
    const arma::mat BC = Bg * C;
    arma::vec leverage = arma::sum(BC % Bg, 1);
    if (ns > 0) {
      const arma::mat cross = Xpg.t() * weighted_B;
      const arma::mat effective = Xpg - BC * cross.t();
      leverage += arma::sum(
        (effective * determinant_schur_inverse) % effective, 1
      );
    }
    a.elem(ii) = 0.5 * determinant_derivative.elem(ii) % leverage;
  }

  arma::vec smooth_rhs = X_penalized.t() * a;
  std::vector<arma::vec> smooth_rhs_thread(
    nt, arma::vec(ns, arma::fill::zeros)
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
      const arma::mat Bg = B.rows(ii);
      if (ns > 0) {
        const arma::mat Xpg = X_penalized.rows(ii);
        const arma::mat weighted_B = Bg.each_col() % working_weight.elem(ii);
        const arma::mat C = inverse_spd(
          precision + Bg.t() * weighted_B,
          "A Laplace variance-step influence Hessian was not positive definite."
        );
        const arma::mat cross = Xpg.t() * weighted_B;
        smooth_rhs_thread[tid] -= cross * C * (Bg.t() * a.elem(ii));
      }
    }
  }
  for (int t = 0; t < nt; ++t) smooth_rhs += smooth_rhs_thread[t];
  arma::vec smooth_influence(ns, arma::fill::zeros);
  if (ns > 0) smooth_influence = schur_inverse * smooth_rhs;

  std::vector<arma::mat> moment_thread(
    nt, arma::mat(q, q, arma::fill::zeros)
  );
  std::vector<arma::mat> conditional_thread(
    nt, arma::mat(q, q, arma::fill::zeros)
  );
  std::vector<arma::mat> mean_thread(
    nt, arma::mat(q, q, arma::fill::zeros)
  );
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
      arma::uvec ii(rows[g]);
      const arma::mat Bg = B.rows(ii);
      const arma::mat Xpg = X_penalized.rows(ii);
      const arma::mat weighted_B = Bg.each_col() % working_weight.elem(ii);
      const arma::mat C = inverse_spd(
        precision + Bg.t() * weighted_B,
        "A Laplace variance-step moment Hessian was not positive definite."
      );
      const arma::mat cross = Xpg.t() * weighted_B;
      arma::vec rhs = Bg.t() * a.elem(ii);
      if (ns > 0) rhs -= cross.t() * smooth_influence;
      const arma::vec subject_influence = C * rhs;
      arma::mat mean_correction(q, q, arma::fill::zeros);
      if (ns > 0) {
        mean_correction = C * cross.t() * schur_inverse * cross * C;
      }
      const arma::vec ug = u.row(g).t();
      const arma::mat influence_correction =
        subject_influence * ug.t() + ug * subject_influence.t();
      conditional_thread[tid] += C;
      mean_thread[tid] += mean_correction;
      influence_thread[tid] += influence_correction;
      moment_thread[tid] += ug * ug.t() + C + mean_correction -
        influence_correction;
    }
  }
  arma::mat moment_sum(q, q, arma::fill::zeros);
  arma::mat conditional_sum(q, q, arma::fill::zeros);
  arma::mat mean_sum(q, q, arma::fill::zeros);
  arma::mat influence_sum(q, q, arma::fill::zeros);
  for (int t = 0; t < nt; ++t) {
    moment_sum += moment_thread[t];
    conditional_sum += conditional_thread[t];
    mean_sum += mean_thread[t];
    influence_sum += influence_thread[t];
  }
  moment_sum = arma::symmatu(moment_sum);
  conditional_sum = arma::symmatu(conditional_sum);
  mean_sum = arma::symmatu(mean_sum);
  influence_sum = arma::symmatu(influence_sum);
  return List::create(
    _["moment_sum"] = moment_sum,
    _["conditional_sum"] = conditional_sum,
    _["mean_sum"] = mean_sum,
    _["influence_sum"] = influence_sum,
    _["precision_score"] = arma::symmatu(
      moment_sum - static_cast<double>(ng) * G
    ),
    _["schur_inverse"] = schur_inverse,
    _["a"] = a
  );
}

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace Rcpp;

namespace {

arma::mat sympd_inverse(const arma::mat& A, const char* message) {
  arma::mat out;
  if (!arma::inv_sympd(out, arma::symmatu(A))) stop(message);
  return out;
}

arma::vec block_apply_vec(const arma::cube& blocks, const arma::vec& x) {
  const int d = blocks.n_rows;
  const int ng = blocks.n_slices;
  arma::vec out(x.n_elem, arma::fill::zeros);
  for (int g = 0; g < ng; ++g) {
    const int lo = g * d;
    const int hi = lo + d - 1;
    out.subvec(lo, hi) = blocks.slice(g) * x.subvec(lo, hi);
  }
  return out;
}

} // namespace

// [[Rcpp::export]]
List gammfast_variance_quadratic(const arma::vec& response,
                                 const arma::mat& X,
                                 const arma::vec& beta,
                                 const arma::mat& B,
                                 const arma::ivec& id,
                                 const arma::mat& G,
                                 const arma::vec& R_diag,
                                 const arma::mat& penalty,
                                 const arma::mat& probes,
                                 bool exact,
                                 double eigen_tol = 1e-8,
                                 int n_threads = 1) {
  const int n = response.n_elem;
  const int p = X.n_cols;
  const int k = B.n_cols;
  if (X.n_rows != static_cast<arma::uword>(n) ||
      B.n_rows != static_cast<arma::uword>(n) ||
      id.n_elem != static_cast<arma::uword>(n) ||
      R_diag.n_elem != static_cast<arma::uword>(n)) {
    stop("response, X, B, id, and R_diag must have matching rows.");
  }
  if (beta.n_elem != static_cast<arma::uword>(p) ||
      penalty.n_rows != static_cast<arma::uword>(p) ||
      penalty.n_cols != static_cast<arma::uword>(p)) {
    stop("beta and penalty must match ncol(X).");
  }
  if (G.n_rows != static_cast<arma::uword>(k) ||
      G.n_cols != static_cast<arma::uword>(k)) {
    stop("G must be square and match ncol(B).");
  }
  if (!response.is_finite() || !X.is_finite() || !beta.is_finite() ||
      !B.is_finite() || !G.is_finite() || !R_diag.is_finite() ||
      !penalty.is_finite() || arma::any(R_diag <= 0)) {
    stop("Variance-quadratic inputs must be finite and R_diag must be positive.");
  }
  if (id.n_elem == 0 || id.min() < 1) {
    stop("id must contain positive consecutive integers.");
  }
  const int ng = id.max();
  std::vector<std::vector<arma::uword>> rows(ng);
  for (int i = 0; i < n; ++i) rows[id[i] - 1].push_back(i);
  for (int g = 0; g < ng; ++g) {
    if (rows[g].empty()) stop("id must contain positive consecutive integers.");
  }

  arma::vec ge;
  arma::mat gu;
  if (!arma::eig_sym(ge, gu, arma::symmatu(G))) {
    stop("The fitted random covariance eigendecomposition failed.");
  }
  const double gscale = std::max(1.0, arma::abs(ge).max());
  arma::uvec gkeep = arma::find(ge > eigen_tol * gscale);
  const int g_negative = arma::sum(ge < -eigen_tol * gscale);
  const double g_min = ge.min();
  if (gkeep.n_elem == 0) {
    return List::create(
      _["statistic"] = 0.0, _["lambda"] = NumericVector(0),
      _["moments"] = NumericVector(0), _["probe_moments"] = NumericMatrix(0),
      _["n_subject"] = ng, _["basis_dimension"] = 0,
      _["fixed_rank"] = 0, _["smooth_rank"] = 0,
      _["g_negative"] = g_negative, _["g_min_eigen"] = g_min,
      _["penalty_negative"] = 0, _["penalty_min_eigen"] = 0.0,
      _["reference_negative"] = 0, _["reference_min_eigen"] = 0.0
    );
  }
  arma::mat Lg = gu.cols(gkeep) * arma::diagmat(arma::sqrt(ge.elem(gkeep)));
  const int d = Lg.n_cols;
  const int q = ng * d;

  arma::vec se;
  arma::mat su;
  if (!arma::eig_sym(se, su, arma::symmatu(penalty))) {
    stop("The global penalty eigendecomposition failed.");
  }
  const double sscale = std::max(1.0, arma::abs(se).max());
  arma::uvec spos = arma::find(se > eigen_tol * sscale);
  arma::uvec szero = arma::find(se <= eigen_tol * sscale);
  const int penalty_negative = arma::sum(se < -eigen_tol * sscale);
  const double penalty_min = se.min();
  const int ns = spos.n_elem;
  const int nf = szero.n_elem;
  arma::mat Us(n, ns, arma::fill::zeros);
  if (ns > 0) {
    Us = X * su.cols(spos);
    Us.each_row() %= (1.0 / arma::sqrt(se.elem(spos))).t();
  }
  arma::mat Xf = X * su.cols(szero);
  arma::vec beta_f = su.cols(szero).t() * beta;
  arma::vec residual = response - Xf * beta_f;

  arma::cube T(d, d, ng, arma::fill::zeros);
  arma::mat CU(q, ns, arma::fill::zeros);
  arma::mat CX(q, nf, arma::fill::zeros);
  arma::vec Cr(q, arma::fill::zeros);
  if (n_threads < 1) n_threads = 1;
  const int nt = std::min(n_threads, ng);
  std::vector<arma::mat> UtAU_thread(nt, arma::mat(ns, ns, arma::fill::zeros));
  std::vector<arma::mat> UtAX_thread(nt, arma::mat(ns, nf, arma::fill::zeros));
  std::vector<arma::mat> XtAX_thread(nt, arma::mat(nf, nf, arma::fill::zeros));
  std::vector<arma::vec> UtAr_thread(nt, arma::vec(ns, arma::fill::zeros));
  std::vector<arma::vec> XtAr_thread(nt, arma::vec(nf, arma::fill::zeros));
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
      arma::mat Li = B.rows(ii) * Lg;
      arma::mat Ui = Us.rows(ii);
      arma::mat Xi = Xf.rows(ii);
      arma::vec ri = residual.elem(ii);
      arma::vec rinv = 1.0 / R_diag.elem(ii);
      arma::mat M = arma::join_rows(Ui, Xi);
      M = arma::join_rows(M, ri);
      M = arma::join_rows(M, Li);
      arma::mat AiM = M.each_col() % rinv;
      int pos = 0;
      arma::mat AiU(ii.n_elem, ns, arma::fill::zeros);
      arma::mat AiX(ii.n_elem, nf, arma::fill::zeros);
      if (ns > 0) { AiU = AiM.cols(pos, pos + ns - 1); pos += ns; }
      if (nf > 0) { AiX = AiM.cols(pos, pos + nf - 1); pos += nf; }
      arma::vec Air = AiM.col(pos); ++pos;
      arma::mat AiL = AiM.cols(pos, pos + d - 1);
      const int lo = g * d;
      const int hi = lo + d - 1;
      T.slice(g) = Li.t() * AiL;
      if (ns > 0) CU.rows(lo, hi) = Li.t() * AiU;
      if (nf > 0) CX.rows(lo, hi) = Li.t() * AiX;
      Cr.subvec(lo, hi) = Li.t() * Air;
      if (ns > 0) {
        UtAU_thread[tid] += Ui.t() * AiU;
        UtAr_thread[tid] += Ui.t() * Air;
      }
      if (nf > 0) {
        XtAX_thread[tid] += Xi.t() * AiX;
        XtAr_thread[tid] += Xi.t() * Air;
      }
      if (ns > 0 && nf > 0) UtAX_thread[tid] += Ui.t() * AiX;
    }
  }

  arma::mat UtAU(ns, ns, arma::fill::zeros), UtAX(ns, nf, arma::fill::zeros);
  arma::mat XtAX(nf, nf, arma::fill::zeros);
  arma::vec UtAr(ns, arma::fill::zeros), XtAr(nf, arma::fill::zeros);
  for (int t = 0; t < nt; ++t) {
    UtAU += UtAU_thread[t]; UtAX += UtAX_thread[t]; XtAX += XtAX_thread[t];
    UtAr += UtAr_thread[t]; XtAr += XtAr_thread[t];
  }
  arma::mat Cinv(ns, ns, arma::fill::zeros);
  if (ns > 0) Cinv = sympd_inverse(arma::eye(ns, ns) + UtAU,
                                   "The global smooth Woodbury system was not positive definite.");
  arma::mat LPX = CX;
  arma::vec LPr = Cr;
  arma::mat H = XtAX;
  arma::vec XPr = XtAr;
  if (ns > 0) {
    LPX -= CU * Cinv * UtAX;
    LPr -= CU * Cinv * UtAr;
    H -= UtAX.t() * Cinv * UtAX;
    XPr -= UtAX.t() * Cinv * UtAr;
  }
  arma::mat Hinv(nf, nf, arma::fill::zeros);
  if (nf > 0) {
    Hinv = sympd_inverse(H, "The fixed-effect projected precision was not positive definite.");
    LPr -= LPX * Hinv * XPr;
  }
  const double statistic = arma::dot(LPr, LPr);

  auto apply_A = [&](const arma::vec& v) {
    arma::vec out = block_apply_vec(T, v);
    if (ns > 0) out -= CU * (Cinv * (CU.t() * v));
    if (nf > 0) out -= LPX * (Hinv * (LPX.t() * v));
    return out;
  };

  if (exact) {
    arma::mat A(q, q, arma::fill::zeros);
    for (int g = 0; g < ng; ++g) {
      const int lo = g * d;
      const int hi = lo + d - 1;
      A.submat(lo, lo, hi, hi) = T.slice(g);
    }
    if (ns > 0) A -= CU * Cinv * CU.t();
    if (nf > 0) A -= LPX * Hinv * LPX.t();
    arma::vec lambda;
    if (!arma::eig_sym(lambda, arma::symmatu(A))) {
      stop("The post-estimation quadratic reference eigendecomposition failed.");
    }
    const double lscale = std::max(1.0, arma::abs(lambda).max());
    const int reference_negative = arma::sum(lambda < -eigen_tol * lscale);
    const double reference_min = lambda.min();
    lambda.elem(arma::find(lambda <= eigen_tol * lscale)).zeros();
    return List::create(
      _["statistic"] = std::max(0.0, statistic), _["lambda"] = lambda,
      _["moments"] = NumericVector(0), _["probe_moments"] = NumericMatrix(0),
      _["n_subject"] = ng, _["basis_dimension"] = d,
      _["fixed_rank"] = nf, _["smooth_rank"] = ns,
      _["g_negative"] = g_negative, _["g_min_eigen"] = g_min,
      _["penalty_negative"] = penalty_negative,
      _["penalty_min_eigen"] = penalty_min,
      _["reference_negative"] = reference_negative,
      _["reference_min_eigen"] = reference_min
    );
  }

  if (probes.n_rows < static_cast<arma::uword>(q) || probes.n_cols < 2) {
    stop("Randomized probes must cover the positive subject covariance rank.");
  }
  const int np = probes.n_cols;
  arma::mat probe_moments(np, 4, arma::fill::zeros);
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(nt)
#endif
  for (int j = 0; j < np; ++j) {
    arma::vec z = probes.col(j).head(q);
    arma::vec v = z;
    for (int r = 0; r < 4; ++r) {
      v = apply_A(v);
      probe_moments(j, r) = std::max(0.0, arma::dot(z, v));
    }
  }
  return List::create(
    _["statistic"] = std::max(0.0, statistic), _["lambda"] = NumericVector(0),
    _["moments"] = arma::mean(probe_moments, 0).t(),
    _["probe_moments"] = probe_moments,
    _["n_subject"] = ng, _["basis_dimension"] = d,
    _["fixed_rank"] = nf, _["smooth_rank"] = ns,
    _["g_negative"] = g_negative, _["g_min_eigen"] = g_min,
    _["penalty_negative"] = penalty_negative,
    _["penalty_min_eigen"] = penalty_min,
    _["reference_negative"] = 0, _["reference_min_eigen"] = 0.0
  );
}

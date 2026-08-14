#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace Rcpp;

namespace {

struct PenaltyParts {
  arma::mat Us;
  arma::mat Xf;
  arma::mat positive_vectors;
  arma::mat zero_vectors;
  arma::vec sqrt_positive;
  int negative;
  double minimum;
};

arma::mat sympd_inverse(const arma::mat& A, const char* message) {
  arma::mat out;
  const arma::mat S = arma::symmatu(A);
  if (!arma::inv_sympd(out, S)) stop(message);
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

PenaltyParts decompose_penalty(const arma::mat& X,
                               const arma::mat& penalty,
                               double eigen_tol) {
  arma::vec values;
  arma::mat vectors;
  if (!arma::eig_sym(values, vectors, arma::symmatu(penalty))) {
    stop("The global penalty eigendecomposition failed.");
  }
  const double scale = std::max(1.0, arma::abs(values).max());
  arma::uvec positive = arma::find(values > eigen_tol * scale);
  arma::uvec zero = arma::find(values <= eigen_tol * scale);
  PenaltyParts out;
  out.negative = arma::sum(values < -eigen_tol * scale);
  out.minimum = values.min();
  out.positive_vectors = vectors.cols(positive);
  out.zero_vectors = vectors.cols(zero);
  out.sqrt_positive = arma::sqrt(values.elem(positive));
  out.Us = X * out.positive_vectors;
  if (positive.n_elem > 0) {
    out.Us.each_row() %= (1.0 / out.sqrt_positive).t();
  }
  out.Xf = X * vectors.cols(zero);
  return out;
}

} // namespace

// [[Rcpp::export]]
List gammfast_projected_moments(const arma::vec& response,
                                const arma::mat& X,
                                const arma::mat& B,
                                const arma::ivec& id,
                                const arma::mat& G,
                                const arma::mat& penalty,
                                bool return_projection = false,
                                double eigen_tol = 1e-10,
                                int n_threads = 1) {
  const int n = response.n_elem;
  const int p = X.n_cols;
  const int q = B.n_cols;
  if (X.n_rows != static_cast<arma::uword>(n) ||
      B.n_rows != static_cast<arma::uword>(n) ||
      id.n_elem != static_cast<arma::uword>(n)) {
    stop("response, X, B, and id must have matching rows.");
  }
  if (G.n_rows != static_cast<arma::uword>(q) ||
      G.n_cols != static_cast<arma::uword>(q)) {
    stop("G must be square and match ncol(B).");
  }
  if (penalty.n_rows != static_cast<arma::uword>(p) ||
      penalty.n_cols != static_cast<arma::uword>(p)) {
    stop("penalty must be square and match ncol(X).");
  }
  if (!response.is_finite() || !X.is_finite() || !B.is_finite() ||
      !G.is_finite() || !penalty.is_finite()) {
    stop("Projected-moment inputs must be finite.");
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

  arma::mat Ginv = sympd_inverse(G, "G must be positive definite.");
  PenaltyParts parts = decompose_penalty(X, penalty, eigen_tol);
  const int ns = parts.Us.n_cols;
  const int nf = parts.Xf.n_cols;

  arma::mat UtAU(ns, ns, arma::fill::zeros);
  arma::mat UtAX(ns, nf, arma::fill::zeros);
  arma::mat XtAX(nf, nf, arma::fill::zeros);
  arma::vec UtAr(ns, arma::fill::zeros);
  arma::vec XtAr(nf, arma::fill::zeros);
  arma::mat U2(ns, ns, arma::fill::zeros);
  arma::mat U2X(ns, nf, arma::fill::zeros);
  arma::mat X2(nf, nf, arma::fill::zeros);
  double trace_Ainv = 0.0;

  for (int g = 0; g < ng; ++g) {
    arma::uvec ii(rows[g]);
    arma::mat Bg = B.rows(ii);
    arma::mat D = sympd_inverse(
      Ginv + Bg.t() * Bg,
      "A subject-level posterior precision was not positive definite."
    );
    arma::mat AiU = parts.Us.rows(ii);
    arma::mat AiX = parts.Xf.rows(ii);
    arma::vec Air = response.elem(ii);
    if (ns > 0) AiU -= Bg * D * (Bg.t() * AiU);
    if (nf > 0) AiX -= Bg * D * (Bg.t() * AiX);
    Air -= Bg * D * (Bg.t() * Air);
    if (ns > 0) {
      UtAU += parts.Us.rows(ii).t() * AiU;
      UtAr += parts.Us.rows(ii).t() * Air;
      U2 += AiU.t() * AiU;
    }
    if (nf > 0) {
      XtAX += parts.Xf.rows(ii).t() * AiX;
      XtAr += parts.Xf.rows(ii).t() * Air;
      X2 += AiX.t() * AiX;
    }
    if (ns > 0 && nf > 0) {
      UtAX += parts.Us.rows(ii).t() * AiX;
      U2X += AiU.t() * AiX;
    }
    trace_Ainv += ii.n_elem - arma::trace(D * (Bg.t() * Bg));
  }

  arma::mat Cinv(ns, ns, arma::fill::zeros);
  arma::vec smooth_response(ns, arma::fill::zeros);
  arma::mat smooth_fixed(ns, nf, arma::fill::zeros);
  double trace_global = trace_Ainv;
  if (ns > 0) {
    Cinv = sympd_inverse(
      arma::eye(ns, ns) + UtAU,
      "The global smooth Woodbury system was not positive definite."
    );
    smooth_response = Cinv * UtAr;
    if (nf > 0) smooth_fixed = Cinv * UtAX;
    trace_global -= arma::trace(Cinv * U2);
  }

  arma::mat H = XtAX;
  arma::vec XPr = XtAr;
  if (ns > 0 && nf > 0) {
    H -= UtAX.t() * smooth_fixed;
    XPr -= UtAX.t() * smooth_response;
  }
  arma::mat Hinv(nf, nf, arma::fill::zeros);
  arma::vec fixed_response(nf, arma::fill::zeros);
  double trace_P = trace_global;
  if (nf > 0) {
    Hinv = sympd_inverse(
      H, "The fixed-effect projected precision was not positive definite."
    );
    fixed_response = Hinv * XPr;
    arma::mat VgX2 = X2;
    if (ns > 0) {
      VgX2 -= U2X.t() * smooth_fixed;
      VgX2 -= smooth_fixed.t() * U2X;
      VgX2 += smooth_fixed.t() * U2 * smooth_fixed;
    }
    trace_P -= arma::trace(Hinv * VgX2);
  }

  arma::vec beta(p, arma::fill::zeros);
  arma::mat mean_covariance(p, p, arma::fill::zeros);
  arma::mat coefficient_transform(p, ns, arma::fill::zeros);
  if (ns > 0) {
    arma::vec smooth_coefficient = smooth_response;
    if (nf > 0) smooth_coefficient -= smooth_fixed * fixed_response;
    coefficient_transform = parts.positive_vectors;
    coefficient_transform.each_row() %= (1.0 / parts.sqrt_positive).t();
    beta += coefficient_transform * smooth_coefficient;
    arma::mat smooth_covariance = Cinv;
    if (nf > 0) {
      smooth_covariance += smooth_fixed * Hinv * smooth_fixed.t();
    }
    mean_covariance += coefficient_transform * smooth_covariance *
      coefficient_transform.t();
  }
  if (nf > 0) {
    beta += parts.zero_vectors * fixed_response;
    mean_covariance += parts.zero_vectors * Hinv * parts.zero_vectors.t();
    if (ns > 0) {
      arma::mat cross_covariance = -smooth_fixed * Hinv;
      arma::mat full_cross = coefficient_transform * cross_covariance *
        parts.zero_vectors.t();
      mean_covariance += full_cross + full_cross.t();
    }
  }

  arma::mat u(ng, q, arma::fill::zeros);
  arma::mat bpr(ng, q, arma::fill::zeros);
  arma::mat moment_sum(q, q, arma::fill::zeros);
  arma::mat btpb_sum(q, q, arma::fill::zeros);
  arma::vec Py(n, arma::fill::zeros);
  arma::mat PX;
  if (return_projection) PX.zeros(n, p);
  if (n_threads < 1) n_threads = 1;
  const int nt = std::min(n_threads, ng);
  std::vector<arma::mat> moment_thread(
    nt, arma::mat(q, q, arma::fill::zeros)
  );
  std::vector<arma::mat> btpb_thread(
    nt, arma::mat(q, q, arma::fill::zeros)
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
      arma::mat Bg = B.rows(ii);
      arma::mat D;
      const arma::mat subject_precision = arma::symmatu(Ginv + Bg.t() * Bg);
      if (!arma::inv_sympd(D, subject_precision)) {
        ok[g] = 0;
        continue;
      }
      arma::mat AiU = parts.Us.rows(ii);
      arma::mat AiX = parts.Xf.rows(ii);
      arma::vec Air = response.elem(ii);
      arma::mat AiB = Bg;
      if (ns > 0) AiU -= Bg * D * (Bg.t() * AiU);
      if (nf > 0) AiX -= Bg * D * (Bg.t() * AiX);
      Air -= Bg * D * (Bg.t() * Air);
      AiB -= Bg * D * (Bg.t() * AiB);

      arma::mat CU(q, ns, arma::fill::zeros);
      arma::mat LPX(q, nf, arma::fill::zeros);
      if (ns > 0) CU = Bg.t() * AiU;
      if (nf > 0) {
        LPX = Bg.t() * AiX;
        if (ns > 0) LPX -= CU * smooth_fixed;
      }
      arma::vec BPr = Bg.t() * Air;
      if (ns > 0) BPr -= CU * smooth_response;
      if (nf > 0) BPr -= LPX * fixed_response;
      arma::mat BtPB = Bg.t() * AiB;
      if (ns > 0) BtPB -= CU * Cinv * CU.t();
      if (nf > 0) BtPB -= LPX * Hinv * LPX.t();
      BtPB = arma::symmatu(BtPB);
      arma::vec ug = G * BPr;
      arma::mat posterior = arma::symmatu(G - G * BtPB * G);
      u.row(g) = ug.t();
      bpr.row(g) = BPr.t();
      moment_thread[tid] += ug * ug.t() + posterior;
      btpb_thread[tid] += BtPB;

      arma::vec Pyr = Air;
      if (ns > 0) Pyr -= AiU * smooth_response;
      arma::mat VgX = AiX;
      if (ns > 0 && nf > 0) VgX -= AiU * smooth_fixed;
      if (nf > 0) Pyr -= VgX * fixed_response;
      Py.elem(ii) = Pyr;

      if (return_projection && ns > 0) {
        arma::mat PUs = AiU * Cinv;
        if (nf > 0) {
          PUs -= VgX * Hinv * (UtAX.t() * Cinv);
        }
        PUs.each_row() %= parts.sqrt_positive.t();
        PX.rows(ii) = PUs * parts.positive_vectors.t();
      }
    }
  }
  if (arma::any(ok == 0)) {
    stop("A subject-level posterior precision was not positive definite.");
  }
  for (int t = 0; t < nt; ++t) {
    moment_sum += moment_thread[t];
    btpb_sum += btpb_thread[t];
  }
  const double rss_sum = std::max(
    0.0, arma::dot(Py, Py) + static_cast<double>(n) - trace_P
  );
  return List::create(
    _["beta"] = beta,
    _["mean_covariance"] = arma::symmatu(mean_covariance),
    _["u"] = u, _["bpr"] = bpr,
    _["moment_sum"] = arma::symmatu(moment_sum),
    _["btpb_sum"] = arma::symmatu(btpb_sum),
    _["Py"] = Py, _["PX"] = PX,
    _["trace_P"] = trace_P, _["rss_sum"] = rss_sum,
    _["n_subject"] = ng, _["fixed_rank"] = nf,
    _["smooth_rank"] = ns,
    _["penalty_negative"] = parts.negative,
    _["penalty_min_eigen"] = parts.minimum
  );
}

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

  PenaltyParts parts = decompose_penalty(X, penalty, eigen_tol);
  const int penalty_negative = parts.negative;
  const double penalty_min = parts.minimum;
  const int ns = parts.Us.n_cols;
  const int nf = parts.Xf.n_cols;
  arma::mat Us = parts.Us;
  arma::mat Xf = parts.Xf;
  arma::vec beta_f = parts.zero_vectors.t() * beta;
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

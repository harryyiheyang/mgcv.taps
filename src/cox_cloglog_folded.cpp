// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]

#include <RcppArmadillo.h>
#include <algorithm>
#include <vector>
#include <cmath>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace Rcpp;

inline double clip(double x, double lo, double hi){
  return x < lo ? lo : (x > hi ? hi : x);
}

// [[Rcpp::export]]
List cox_cloglog_folded(
  NumericVector time,
  IntegerVector status,
  IntegerVector strata,
  NumericVector eta,
  List bh_time_list,
  List bh_hazard_list,
  IntegerVector K_per_strata,
  double eps_delta   = 1e-12,   // lower bound for ΔH
  double eps_prob    = 1e-12,   // lower bound for probabilities (used only on q)
  double eps_muprime = 1e-12,   // lower bound for dmu/deta
  double clip_eta_min = -20.0,
  double clip_eta_max =  20.0,
  int n_threads = 1
){
  const int n = time.size();
  if(eta.size()!=n || status.size()!=n || strata.size()!=n)
    stop("time/status/strata/eta must have same length");

  IntegerVector strata_in = clone(strata);
  int S = 0;
  for(int i=0;i<n;++i) S = std::max(S, strata_in[i]);
  if(S<=0) stop("Invalid strata labeling");

  // resolve K per strata
  IntegerVector Ks(S);
  if(K_per_strata.size()==1){
    for(int s=0;s<S;++s) Ks[s] = std::max(1, K_per_strata[0]);
  } else if(K_per_strata.size()==S){
    for(int s=0;s<S;++s) Ks[s] = std::max(1, K_per_strata[s]);
  } else {
    stop("K_per_strata must have length 1 or S");
  }

  // per-strata baseline series and selected indices
  std::vector< std::vector<double> > bh_time(S), bh_haz(S), bh_delta(S);
  std::vector< std::vector<int> >    sel_idx(S);

  for(int s=0;s<S;++s){
    if(s >= bh_time_list.size() || s >= bh_hazard_list.size())
      stop("bh_time_list/bh_hazard_list must have length S");

    NumericVector tvec = bh_time_list[s];
    NumericVector hvec = bh_hazard_list[s];
    if(tvec.size() != hvec.size())
      stop("bh_time and bh_hazard length mismatch in strata %d", s+1);

    const int J = tvec.size();
    if(J==0){
      bh_time[s].clear();
      bh_haz[s].clear();
      bh_delta[s].clear();
      sel_idx[s].clear();
      continue;
    }

    bh_time[s]  = std::vector<double>(tvec.begin(), tvec.end());
    bh_haz[s]   = std::vector<double>(hvec.begin(), hvec.end());
    bh_delta[s] = std::vector<double>(J);

    for(int j=0;j<J;++j){
      const double prev = (j==0 ? 0.0 : bh_haz[s][j-1]);
      bh_delta[s][j] = bh_haz[s][j] - prev;
    }

    const double Hmax = bh_haz[s][J-1];
    int K_s = std::min(Ks[s], J);
    if(Hmax <= 0.0 || K_s<=0){
      sel_idx[s].clear();
      continue;
    }

    std::vector<int> chosen;
    chosen.reserve(K_s);
    for(int m=1; m<=K_s; ++m){
      double target = (Hmax * m) / K_s;
      int j = std::lower_bound(bh_haz[s].begin(), bh_haz[s].end(), target) - bh_haz[s].begin();
      if(j >= J) j = J-1;
      while(!chosen.empty() && j == chosen.back() && j+1 < J) j = j+1;
      if(chosen.empty() || j != chosen.back()) chosen.push_back(j);
    }
    if(chosen.empty() || chosen.back() != J-1) chosen.push_back(J-1);
    sel_idx[s] = chosen;
  }

  // outputs
  NumericVector w_star(n, 0.0);
  NumericVector z_star(n, NA_REAL);
  NumericVector Lambda_i(n, 0.0);
  IntegerVector N_i(n, 0);

  #ifdef _OPENMP
  if(n_threads < 1) n_threads = 1;
  #pragma omp parallel for num_threads(n_threads) schedule(dynamic, 1000)
  #endif
  for(int i=0;i<n;++i){
    const int s = strata_in[i] - 1;
    if(s < 0 || s >= S) continue;

    const std::vector<double>& tJ = bh_time[s];
    const std::vector<double>& HJ = bh_haz[s];
    const std::vector<double>& dJ = bh_delta[s];
    const std::vector<int>&    SJ = sel_idx[s];

    const double ti   = time[i];
    const int    delt = status[i] ? 1 : 0;
    double etai = clip(eta[i], clip_eta_min, clip_eta_max);

    if(SJ.empty() || tJ.empty()){
      w_star[i]   = 0.0;
      z_star[i]   = NA_REAL;
      Lambda_i[i] = 0.0;
      N_i[i]      = delt;
      continue;
    }

    int idx_ti = std::upper_bound(tJ.begin(), tJ.end(), ti) - tJ.begin() - 1;
    if(idx_ti < 0){
      Lambda_i[i] = 0.0;
      w_star[i]   = 0.0;
      z_star[i]   = NA_REAL;
      N_i[i]      = delt;
      continue;
    }

    if(idx_ti >= (int)tJ.size() - 1 && ti > tJ.back()){
      double H_last = HJ.back();
      double t_last = tJ.back();
      double lambda_last = 0.0;
      if(HJ.size() > 1){
        lambda_last = (HJ.back() - HJ[HJ.size()-2]) / (tJ.back() - tJ[tJ.size()-2]);
      }
      Lambda_i[i] = H_last + lambda_last * (ti - t_last);
      idx_ti = (int)tJ.size() - 1;
    } else {
      if(idx_ti >= (int)tJ.size()) idx_ti = (int)tJ.size() - 1;
      Lambda_i[i] = HJ[idx_ti];
    }

    double Wi_sum = 0.0;
    double Wz_sum = 0.0;

    for(size_t k=0; k<SJ.size(); ++k){
      const int j_start    = (k == 0 ? -1 : SJ[k-1]);
      const int j_end_full = SJ[k];
      const int j_end_eff  = std::min(j_end_full, idx_ti);
      if(j_end_eff < j_start + 1) continue;

      double Delta_bin = 0.0;
      for(int j = j_start + 1; j <= j_end_eff; ++j)
        Delta_bin += dJ[j];
        if(Delta_bin < eps_delta) Delta_bin = eps_delta;

        const double t_start = (j_start >= 0 ? tJ[j_start] : 0.0);
        const double t_end   = tJ[j_end_eff];

        int y_ij = 0;
        if(delt == 1){
          bool event_in_bin = (ti >= t_start && ti <= t_end);
          bool is_last_bin  = (j_end_eff == idx_ti);
          if(event_in_bin && is_last_bin) y_ij = 1;
        }

        // ------- 核心改动：只裁剪 q，mu = 1 - q -------
          double a = Delta_bin * std::exp(etai);
          if(a < 0.0) a = 0.0;

          double q = std::exp(-a);
          q = std::max(q, eps_prob);           // 仅下界裁剪 q
          double mu = 1.0 - q;                 // 保持 mu + q = 1 的解析一致性

          // dmu/deta = a * q
          double dmu_deta = a * q;
          if(dmu_deta < eps_muprime) dmu_deta = eps_muprime;

          // Var(y) = mu * q
          double var_y = mu * q;
          if(var_y < eps_prob*eps_prob) var_y = eps_prob*eps_prob;

          double w_ij = (dmu_deta * dmu_deta) / var_y; // = (a^2 q)/mu
          double z_ij = etai + ((double)y_ij - mu) / dmu_deta;

          Wi_sum += w_ij;
          Wz_sum += w_ij * z_ij;
    }

    if(Wi_sum > 0.0){
      double zi = Wz_sum / Wi_sum;
      z_star[i] = clip(zi, etai - 20.0, etai + 20.0);
      w_star[i] = Wi_sum;
    } else {
      z_star[i] = NA_REAL;
      w_star[i] = 0.0;
    }

    N_i[i] = delt;
  }

  return List::create(
    _["w_star"]   = w_star,
    _["z_star"]   = z_star,
    _["Lambda_i"] = Lambda_i,
    _["N_i"]      = N_i,
    _["sel_idx"]  = sel_idx
  );
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]

#include <RcppArmadillo.h>
#include <algorithm>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace {

arma::mat row_weight(const arma::mat& x, const arma::vec& w) {
  arma::mat out = x;
  out.each_col() %= w;
  return out;
}

arma::mat chol_solve(const arma::mat& chol, const arma::mat& rhs) {
  if (rhs.n_cols == 0L) return arma::mat(chol.n_rows, 0L);
  arma::mat y = arma::solve(
    arma::trimatl(chol), rhs, arma::solve_opts::fast
  );
  return arma::solve(
    arma::trimatu(chol.t()), y, arma::solve_opts::fast
  );
}

class NativeReUidProfile {
public:
  arma::mat dense;
  arma::mat weighted_dense;
  arma::vec weight;
  arma::mat local_penalty;
  std::vector<arma::uvec> rows;
  std::vector<arma::mat> local;
  std::vector<arma::mat> weighted_local;
  std::vector<arma::mat> chol;
  std::vector<arma::mat> cross;
  arma::mat schur;
  int n_threads;

  NativeReUidProfile(
    const arma::mat& dense_,
    const arma::mat& local_,
    const arma::ivec& group,
    const arma::vec& weight_,
    const arma::mat& dense_penalty,
    const arma::mat& local_penalty_,
    int n_threads_
  ) :
    dense(dense_),
    weighted_dense(row_weight(dense_, weight_)),
    weight(weight_),
    local_penalty((local_penalty_ + local_penalty_.t()) / 2.0),
    n_threads(std::max(1, n_threads_)) {

    const int n_group = group.max();
    n_threads = std::min(n_threads, n_group);
    rows.resize(n_group);
    local.resize(n_group);
    weighted_local.resize(n_group);
    chol.resize(n_group);
    cross.resize(n_group);

    std::vector<std::vector<arma::uword> > row_buffer(n_group);
    for (arma::uword i = 0; i < group.n_elem; ++i) {
      row_buffer[group[i] - 1].push_back(i);
    }
    for (int g = 0; g < n_group; ++g) {
      rows[g] = arma::conv_to<arma::uvec>::from(row_buffer[g]);
      local[g] = local_.rows(rows[g]);
      weighted_local[g] = row_weight(local[g], weight.elem(rows[g]));
    }

    arma::mat dense_cross = dense.t() * weighted_dense + dense_penalty;
    schur = (dense_cross + dense_cross.t()) / 2.0;
    std::vector<arma::mat> schur_part(n_group);
    std::vector<int> failed(n_group, 0);

    #ifdef _OPENMP
    #pragma omp parallel for num_threads(n_threads) schedule(static)
    #endif
    for (int g = 0; g < n_group; ++g) {
      arma::mat cg = local[g].t() * weighted_local[g] + local_penalty;
      cg = (cg + cg.t()) / 2.0;
      if (!arma::chol(chol[g], cg, "lower")) {
        failed[g] = 1;
        continue;
      }
      cross[g] = local[g].t() * weighted_dense.rows(rows[g]);
      const arma::mat solved_cross = chol_solve(chol[g], cross[g]);
      schur_part[g] = cross[g].t() * solved_cross;
    }

    for (int g = 0; g < n_group; ++g) {
      if (failed[g]) {
        Rcpp::stop("A uid-local random-effect block is not positive definite.");
      }
      schur -= schur_part[g];
    }
    schur = (schur + schur.t()) / 2.0;
  }

  arma::mat apply(const arma::mat& value, const arma::mat& schur_inverse) const {
    const int n_group = rows.size();
    const arma::mat weighted_value = row_weight(value, weight);
    const arma::mat rhs_dense = dense.t() * weighted_value;
    std::vector<arma::mat> rhs_local(n_group);
    std::vector<arma::mat> dense_adjustment(n_group);

    #ifdef _OPENMP
    #pragma omp parallel for num_threads(n_threads) schedule(static)
    #endif
    for (int g = 0; g < n_group; ++g) {
      rhs_local[g] = local[g].t() * weighted_value.rows(rows[g]);
      const arma::mat local_zero = chol_solve(chol[g], rhs_local[g]);
      dense_adjustment[g] = cross[g].t() * local_zero;
    }

    arma::mat adjusted_rhs = rhs_dense;
    for (int g = 0; g < n_group; ++g) {
      adjusted_rhs -= dense_adjustment[g];
    }
    const arma::mat dense_coef = schur_inverse * adjusted_rhs;
    arma::mat out = weighted_value;

    #ifdef _OPENMP
    #pragma omp parallel for num_threads(n_threads) schedule(static)
    #endif
    for (int g = 0; g < n_group; ++g) {
      const arma::mat local_rhs = rhs_local[g] - cross[g] * dense_coef;
      const arma::mat local_coef = chol_solve(chol[g], local_rhs);
      arma::mat correction = weighted_local[g] * local_coef;
      if (dense.n_cols > 0L) {
        correction += weighted_dense.rows(rows[g]) * dense_coef;
      }
      out.rows(rows[g]) -= correction;
    }
    return out;
  }
};

}

// [[Rcpp::export]]
Rcpp::List native_re_uid_profile_create_cpp(
  const arma::mat& dense,
  const arma::mat& local,
  const arma::ivec& group,
  const arma::vec& weight,
  const arma::mat& dense_penalty,
  const arma::mat& local_penalty,
  int n_threads = 1
) {
  const arma::uword n = weight.n_elem;
  if (n == 0L) {
    Rcpp::stop("Native random-effect design has no observations.");
  }
  if (dense.n_rows != n || local.n_rows != n || group.n_elem != n) {
    Rcpp::stop("Native random-effect design inputs have inconsistent row counts.");
  }
  if (local.n_cols == 0L) {
    Rcpp::stop("Native random-effect design has no uid-local columns.");
  }
  if (dense_penalty.n_rows != dense.n_cols ||
      dense_penalty.n_cols != dense.n_cols) {
    Rcpp::stop("Dense nuisance penalty has inconsistent dimensions.");
  }
  if (local_penalty.n_rows != local.n_cols ||
      local_penalty.n_cols != local.n_cols) {
    Rcpp::stop("Uid-local nuisance penalty has inconsistent dimensions.");
  }
  if (!dense.is_finite() || !local.is_finite() || !weight.is_finite() ||
      !dense_penalty.is_finite() || !local_penalty.is_finite()) {
    Rcpp::stop("Native random-effect design inputs must be finite.");
  }
  if (arma::any(weight <= 0.0)) {
    Rcpp::stop("Native random-effect weights must be positive.");
  }
  if (group.min() != 1 || group.max() < 1) {
    Rcpp::stop("Uid group indexes must start at one.");
  }
  std::vector<int> group_count(group.max(), 0);
  for (arma::uword i = 0; i < group.n_elem; ++i) {
    if (group[i] < 1 || group[i] > group.max()) {
      Rcpp::stop("Uid group indexes must be consecutive positive integers.");
    }
    group_count[group[i] - 1]++;
  }
  if (std::find(group_count.begin(), group_count.end(), 0) != group_count.end()) {
    Rcpp::stop("Every uid group must contain at least one observation.");
  }

  Rcpp::XPtr<NativeReUidProfile> profile(
    new NativeReUidProfile(
      dense, local, group, weight, dense_penalty, local_penalty, n_threads
    ),
    true
  );
  return Rcpp::List::create(
    Rcpp::_ ["pointer"] = profile,
    Rcpp::_ ["schur"] = profile->schur,
    Rcpp::_ ["n_group"] = group.max(),
    Rcpp::_ ["block_size"] = local.n_cols
  );
}

// [[Rcpp::export]]
arma::mat native_re_uid_profile_apply_cpp(
  SEXP pointer,
  const arma::mat& value,
  const arma::mat& schur_inverse
) {
  Rcpp::XPtr<NativeReUidProfile> profile(pointer);
  if (profile.get() == NULL) {
    Rcpp::stop("Native random-effect profile pointer is no longer valid.");
  }
  if (value.n_rows != profile->weight.n_elem) {
    Rcpp::stop("Profile input has an inconsistent number of rows.");
  }
  if (schur_inverse.n_rows != profile->dense.n_cols ||
      schur_inverse.n_cols != profile->dense.n_cols) {
    Rcpp::stop("Dense Schur inverse has inconsistent dimensions.");
  }
  return profile->apply(value, schur_inverse);
}

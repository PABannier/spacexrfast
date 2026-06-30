#include <Rcpp.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <set>
#include <string>
#include <thread>
#include <vector>

using Rcpp::IntegerVector;
using Rcpp::List;
using Rcpp::CharacterVector;
using Rcpp::LogicalVector;
using Rcpp::NumericMatrix;
using Rcpp::NumericVector;

namespace {

struct LikelihoodTables {
  std::vector<double> q;
  std::vector<double> sq;
  std::vector<double> x;
  int n_row = 0;
  int n_col = 0;
  int k_val = 0;
};

struct QDerivatives {
  double d0 = 0.0;
  double d1 = 0.0;
  double d2 = 0.0;
};

struct FitResult {
  std::vector<double> weights;
  bool converged = false;
};

struct TypePresence {
  bool all_pairs = true;
  bool all_pairs_class = false;
  double singlet_score = 0.0;
};

inline double clamp_double(double x, double lo, double hi) {
  return std::min(std::max(x, lo), hi);
}

inline double q_at(const LikelihoodTables& tables, int y, int m_1based) {
  const int row = std::min(std::max(y, 0), tables.k_val);
  const int col = std::min(std::max(m_1based - 1, 0), tables.n_col - 1);
  return tables.q[row + col * tables.n_row];
}

inline double sq_at(const LikelihoodTables& tables, int y, int m_1based) {
  const int row = std::min(std::max(y, 0), tables.k_val);
  const int col = std::min(std::max(m_1based - 1, 0), tables.n_col - 1);
  return tables.sq[row + col * tables.n_row];
}

inline QDerivatives calc_q_all_one(double y_value, double lambda, const LikelihoodTables& tables) {
  const double epsilon = 1e-4;
  const double delta = 1e-6;
  const double x_max = tables.x.back();
  lambda = clamp_double(lambda, epsilon, x_max - epsilon);

  const int y = std::min(static_cast<int>(y_value), tables.k_val);
  const double l = std::floor(std::sqrt(lambda / delta));
  const double m_left = std::min(l - 9.0, 40.0);
  const double m_right = std::max(std::ceil(std::sqrt(std::max(l - 48.7499, 0.0) * 4.0)) - 2.0, 0.0);
  int m = static_cast<int>(m_left + m_right);
  m = std::min(std::max(m, 1), tables.n_col - 1);

  const double ti1 = tables.x[m - 1];
  const double ti = tables.x[m];
  const double hi = ti - ti1;
  const double fti1 = q_at(tables, y, m);
  const double fti = q_at(tables, y, m + 1);
  const double zi1 = sq_at(tables, y, m);
  const double zi = sq_at(tables, y, m + 1);

  const double diff1 = lambda - ti1;
  const double diff2 = ti - lambda;
  const double diff3 = fti / hi - zi * hi / 6.0;
  const double diff4 = fti1 / hi - zi1 * hi / 6.0;
  const double zdi = zi / hi;
  const double zdi1 = zi1 / hi;

  QDerivatives out;
  out.d0 = zdi * std::pow(diff1, 3.0) / 6.0 +
           zdi1 * std::pow(diff2, 3.0) / 6.0 +
           diff3 * diff1 + diff4 * diff2;
  out.d1 = zdi * diff1 * diff1 / 2.0 -
           zdi1 * diff2 * diff2 / 2.0 +
           diff3 - diff4;
  out.d2 = zdi * diff1 + zdi1 * diff2;
  return out;
}

inline double log_l_score(const std::vector<double>& y,
                          const std::vector<double>& prediction,
                          const LikelihoodTables& tables) {
  double score = 0.0;
  for (size_t g = 0; g < y.size(); ++g) {
    score += -calc_q_all_one(y[g], prediction[g], tables).d0;
  }
  return score;
}

inline void psd_project_2x2(double& a, double& b, double& c, double epsilon = 0.001) {
  const double trace_half = 0.5 * (a + c);
  const double diff_half = 0.5 * (a - c);
  const double radius = std::sqrt(diff_half * diff_half + b * b);
  double lambda1 = trace_half + radius;
  double lambda2 = trace_half - radius;
  lambda1 = std::max(lambda1, epsilon);
  lambda2 = std::max(lambda2, epsilon);

  if (std::abs(b) < 1e-30 && std::abs(a - c) < 1e-30) {
    a = lambda1;
    b = 0.0;
    c = lambda1;
    return;
  }

  double v1x = b;
  double v1y = lambda1 - a;
  const double n1 = std::sqrt(v1x * v1x + v1y * v1y);
  if (n1 < 1e-30) {
    v1x = lambda1 - c;
    v1y = b;
  }
  const double n2 = std::sqrt(v1x * v1x + v1y * v1y);
  v1x /= n2;
  v1y /= n2;
  const double v2x = -v1y;
  const double v2y = v1x;

  a = lambda1 * v1x * v1x + lambda2 * v2x * v2x;
  b = lambda1 * v1x * v1y + lambda2 * v2x * v2y;
  c = lambda1 * v1y * v1y + lambda2 * v2y * v2y;
}

inline double objective_1d(double d, double D, double u) {
  return 0.5 * D * u * u - d * u;
}

inline double objective_2d(double d0, double d1, double D00, double D01, double D11,
                           double u0, double u1) {
  return 0.5 * (D00 * u0 * u0 + 2.0 * D01 * u0 * u1 + D11 * u1 * u1) -
         d0 * u0 - d1 * u1;
}

inline std::vector<double> solve_qp_unconstrained_bounds(const std::vector<double>& d,
                                                         double D00,
                                                         double D01,
                                                         double D11,
                                                         const std::vector<double>& lb) {
  const int n = static_cast<int>(d.size());
  if (n == 1) {
    return {std::max(d[0] / D00, lb[0])};
  }

  const double det = D00 * D11 - D01 * D01;
  std::vector<std::vector<double>> candidates;
  if (std::abs(det) > 1e-30) {
    candidates.push_back({
      (D11 * d[0] - D01 * d[1]) / det,
      (-D01 * d[0] + D00 * d[1]) / det
    });
  }
  candidates.push_back({lb[0], (d[1] - D01 * lb[0]) / D11});
  candidates.push_back({(d[0] - D01 * lb[1]) / D00, lb[1]});
  candidates.push_back({lb[0], lb[1]});

  double best_obj = std::numeric_limits<double>::infinity();
  std::vector<double> best{lb[0], lb[1]};
  for (auto cand : candidates) {
    cand[0] = std::max(cand[0], lb[0]);
    cand[1] = std::max(cand[1], lb[1]);
    const double obj = objective_2d(d[0], d[1], D00, D01, D11, cand[0], cand[1]);
    if (obj < best_obj) {
      best_obj = obj;
      best = cand;
    }
  }
  return best;
}

inline std::vector<double> solve_qp_equality_bounds(const std::vector<double>& d,
                                                    double D00,
                                                    double D01,
                                                    double D11,
                                                    const std::vector<double>& lb,
                                                    double rhs) {
  const int n = static_cast<int>(d.size());
  if (n == 1) {
    return {rhs};
  }

  const double lower = lb[0];
  const double upper = rhs - lb[1];
  double u0 = lower;
  const double quad = D00 - 2.0 * D01 + D11;
  const double lin = D01 * rhs - D11 * rhs - d[0] + d[1];
  if (std::abs(quad) > 1e-30) {
    u0 = -lin / quad;
  }
  u0 = clamp_double(u0, lower, upper);
  return {u0, rhs - u0};
}

FitResult solve_irwls(const std::vector<double>& y,
                      const std::vector<int>& types,
                      const std::vector<double>& ref,
                      int n_genes,
                      int n_types,
                      double nUMI,
                      const LikelihoodTables& tables,
                      int n_iter,
                      double min_change,
                      bool constrain) {
  const int k = static_cast<int>(types.size());
  std::vector<double> solution(k, 1.0 / static_cast<double>(k));
  double change = 1.0;
  int iterations = 0;

  std::vector<double> s0(n_genes), s1(n_genes), prediction(n_genes);
  for (int g = 0; g < n_genes; ++g) {
    s0[g] = nUMI * ref[g + types[0] * n_genes];
    if (k == 2) {
      s1[g] = nUMI * ref[g + types[1] * n_genes];
    }
  }

  while (change > min_change && iterations < n_iter) {
    std::vector<double> current(k);
    for (int i = 0; i < k; ++i) {
      current[i] = std::max(solution[i], 0.0);
    }

    const double threshold = std::max(1e-4, nUMI * 1e-7);
    for (int g = 0; g < n_genes; ++g) {
      double pred = s0[g] * current[0];
      if (k == 2) {
        pred += s1[g] * current[1];
      }
      pred = std::abs(pred);
      prediction[g] = pred < threshold ? threshold : pred;
    }

    double grad0 = 0.0, grad1 = 0.0;
    double h00 = 0.0, h01 = 0.0, h11 = 0.0;
    for (int g = 0; g < n_genes; ++g) {
      const QDerivatives qd = calc_q_all_one(y[g], prediction[g], tables);
      grad0 += -qd.d1 * s0[g];
      h00 += -qd.d2 * s0[g] * s0[g];
      if (k == 2) {
        grad1 += -qd.d1 * s1[g];
        h01 += -qd.d2 * s0[g] * s1[g];
        h11 += -qd.d2 * s1[g] * s1[g];
      }
    }

    std::vector<double> d(k);
    d[0] = -grad0;
    if (k == 2) {
      d[1] = -grad1;
      psd_project_2x2(h00, h01, h11);
      const double trace_half = 0.5 * (h00 + h11);
      const double diff_half = 0.5 * (h00 - h11);
      const double norm_factor = trace_half + std::sqrt(diff_half * diff_half + h01 * h01);
      h00 /= norm_factor;
      h01 /= norm_factor;
      h11 /= norm_factor;
      d[0] /= norm_factor;
      d[1] /= norm_factor;
      h00 += 1e-7;
      h11 += 1e-7;
    } else {
      h00 = std::max(h00, 0.001);
      const double norm_factor = h00;
      h00 = 1.0 + 1e-7;
      d[0] /= norm_factor;
    }

    std::vector<double> lb(k);
    for (int i = 0; i < k; ++i) {
      lb[i] = -current[i];
    }

    std::vector<double> step;
    if (constrain) {
      const double rhs = 1.0 - std::accumulate(current.begin(), current.end(), 0.0);
      step = solve_qp_equality_bounds(d, h00, h01, h11, lb, rhs);
    } else {
      step = solve_qp_unconstrained_bounds(d, h00, h01, h11, lb);
    }

    const double alpha = 0.3;
    std::vector<double> new_solution(k);
    change = 0.0;
    for (int i = 0; i < k; ++i) {
      new_solution[i] = current[i] + alpha * step[i];
      const double diff = new_solution[i] - solution[i];
      change += diff * diff;
    }
    change = std::sqrt(change);
    solution.swap(new_solution);
    ++iterations;
  }

  FitResult result;
  result.weights = solution;
  result.converged = change <= min_change;
  return result;
}

double fit_score(const std::vector<double>& y,
                 const std::vector<int>& types,
                 const std::vector<double>& ref,
                 int n_genes,
                 int n_types,
                 double nUMI,
                 const LikelihoodTables& tables,
                 int n_iter,
                 double min_change,
                 bool constrain) {
  const FitResult fit = solve_irwls(y, types, ref, n_genes, n_types, nUMI, tables,
                                    n_iter, min_change, constrain);
  std::vector<double> prediction(n_genes, 0.0);
  for (int g = 0; g < n_genes; ++g) {
    for (size_t k = 0; k < types.size(); ++k) {
      prediction[g] += nUMI * ref[g + types[k] * n_genes] * fit.weights[k];
    }
  }
  return log_l_score(y, prediction, tables);
}

double singlet_score(const std::vector<double>& y,
                     int type,
                     const std::vector<double>& ref,
                     int n_genes,
                     int n_types,
                     double nUMI,
                     const LikelihoodTables& tables,
                     double min_change,
                     bool constrain) {
  if (!constrain) {
    return fit_score(y, {type}, ref, n_genes, n_types, nUMI, tables, 25, min_change, false);
  }
  std::vector<double> prediction(n_genes);
  for (int g = 0; g < n_genes; ++g) {
    prediction[g] = nUMI * ref[g + type * n_genes];
  }
  return log_l_score(y, prediction, tables);
}

TypePresence check_type_presence(const std::vector<int>& candidates,
                                 const std::vector<double>& score_mat,
                                 double min_score,
                                 int my_type,
                                 const std::vector<int>& class_ids,
                                 double ql_score_cutoff,
                                 const std::vector<double>& singlet_scores) {
  const int n = static_cast<int>(candidates.size());
  const bool has_class = std::any_of(class_ids.begin(), class_ids.end(), [](int x) { return x >= 0; });
  TypePresence out;
  out.singlet_score = singlet_scores[my_type];
  out.all_pairs = true;
  out.all_pairs_class = has_class;

  std::set<int> other_class;
  other_class.insert(my_type);
  for (int i = 0; i < n - 1; ++i) {
    const int type1 = candidates[i];
    for (int j = i + 1; j < n; ++j) {
      const int type2 = candidates[j];
      if (score_mat[i + j * n] < min_score + ql_score_cutoff) {
        if (type1 != my_type && type2 != my_type) {
          out.all_pairs = false;
        }
        if (has_class) {
          const bool first_class = class_ids[my_type] == class_ids[type1];
          const bool second_class = class_ids[my_type] == class_ids[type2];
          if (!first_class && !second_class) {
            out.all_pairs_class = false;
          }
          if (first_class) {
            other_class.insert(type1);
          }
          if (second_class) {
            other_class.insert(type2);
          }
        }
      }
    }
  }

  if (!has_class) {
    out.all_pairs_class = out.all_pairs;
  }
  if (out.all_pairs_class && !out.all_pairs && other_class.size() > 1) {
    for (int type : other_class) {
      if (type != my_type) {
        out.singlet_score = std::min(out.singlet_score, singlet_scores[type]);
      }
    }
  }
  return out;
}

struct FastDoubletWorker {
  const std::vector<double>& beads;
  const std::vector<double>& nUMI;
  const std::vector<double>& ref;
  const std::vector<double>& all_weights;
  const std::vector<unsigned char>& conv_all_in;
  const std::vector<int>& class_ids;
  const LikelihoodTables& tables;
  const int n_cells;
  const int n_genes;
  const int n_types;
  const bool constrain;
  const double min_change;
  const double confidence_threshold;
  const double doublet_threshold;
  const bool return_diagnostics;

  std::vector<int>& spot_class;
  std::vector<int>& first_type;
  std::vector<int>& second_type;
  std::vector<unsigned char>& first_class;
  std::vector<unsigned char>& second_class;
  std::vector<double>& first_weight;
  std::vector<double>& second_weight;
  std::vector<double>& min_score_out;
  std::vector<double>& singlet_score_out;
  std::vector<unsigned char>& conv_doublet;
  std::vector<std::vector<double>>& score_mats;
  std::vector<std::vector<double>>& singlet_scores_out;
  std::vector<std::vector<int>>& diagnostic_candidates;

  FastDoubletWorker(const std::vector<double>& beads,
                    const std::vector<double>& nUMI,
                    const std::vector<double>& ref,
                    const std::vector<double>& all_weights,
                    const std::vector<unsigned char>& conv_all_in,
                    const std::vector<int>& class_ids,
                    const LikelihoodTables& tables,
                    int n_cells,
                    int n_genes,
                    int n_types,
                    bool constrain,
                    double min_change,
                    double confidence_threshold,
                    double doublet_threshold,
                    bool return_diagnostics,
                    std::vector<int>& spot_class,
                    std::vector<int>& first_type,
                    std::vector<int>& second_type,
                    std::vector<unsigned char>& first_class,
                    std::vector<unsigned char>& second_class,
                    std::vector<double>& first_weight,
                    std::vector<double>& second_weight,
                    std::vector<double>& min_score_out,
                    std::vector<double>& singlet_score_out,
                    std::vector<unsigned char>& conv_doublet,
                    std::vector<std::vector<double>>& score_mats,
                    std::vector<std::vector<double>>& singlet_scores_out,
                    std::vector<std::vector<int>>& diagnostic_candidates)
      : beads(beads), nUMI(nUMI), ref(ref), all_weights(all_weights),
        conv_all_in(conv_all_in), class_ids(class_ids), tables(tables),
        n_cells(n_cells), n_genes(n_genes), n_types(n_types), constrain(constrain),
        min_change(min_change), confidence_threshold(confidence_threshold),
        doublet_threshold(doublet_threshold), return_diagnostics(return_diagnostics),
        spot_class(spot_class), first_type(first_type), second_type(second_type),
        first_class(first_class), second_class(second_class), first_weight(first_weight),
        second_weight(second_weight), min_score_out(min_score_out),
        singlet_score_out(singlet_score_out), conv_doublet(conv_doublet),
        score_mats(score_mats), singlet_scores_out(singlet_scores_out),
        diagnostic_candidates(diagnostic_candidates) {}

  void run(std::size_t begin, std::size_t end) {
    for (std::size_t cell = begin; cell < end; ++cell) {
      std::vector<double> y(n_genes);
      for (int g = 0; g < n_genes; ++g) {
        y[g] = beads[cell + static_cast<std::size_t>(g) * n_cells];
      }

      std::vector<int> candidates;
      candidates.reserve(n_types);
      for (int t = 0; t < n_types; ++t) {
        if (all_weights[cell + static_cast<std::size_t>(t) * n_cells] > 0.01) {
          candidates.push_back(t);
        }
      }
      if (candidates.empty()) {
        const int fallback = std::min(3, n_types);
        for (int t = 0; t < fallback; ++t) {
          candidates.push_back(t);
        }
      }
      if (candidates.size() == 1) {
        if (candidates[0] == 0) {
          candidates.push_back(1);
        } else {
          candidates.push_back(0);
        }
      }

      std::vector<double> singlet_scores(n_types, std::numeric_limits<double>::infinity());
      for (int type : candidates) {
        singlet_scores[type] = singlet_score(y, type, ref, n_genes, n_types, nUMI[cell],
                                             tables, min_change, constrain);
      }

      const int n_cand = static_cast<int>(candidates.size());
      std::vector<double> score_mat(static_cast<std::size_t>(n_cand) * n_cand, 0.0);
      double min_score = 0.0;
      int first = -1;
      int second = -1;
      for (int i = 0; i < n_cand - 1; ++i) {
        const int type1 = candidates[i];
        for (int j = i + 1; j < n_cand; ++j) {
          const int type2 = candidates[j];
          const double score = fit_score(y, {type1, type2}, ref, n_genes, n_types,
                                         nUMI[cell], tables, 25, min_change, constrain);
          score_mat[i + j * n_cand] = score;
          score_mat[j + i * n_cand] = score;
          if (second < 0 || score < min_score) {
            first = type1;
            second = type2;
            min_score = score;
          }
        }
      }

      TypePresence type1_pres = check_type_presence(candidates, score_mat, min_score, first,
                                                    class_ids, confidence_threshold, singlet_scores);
      TypePresence type2_pres = check_type_presence(candidates, score_mat, min_score, second,
                                                    class_ids, confidence_threshold, singlet_scores);

      int spot = 0;
      bool first_class_value = false;
      bool second_class_value = false;
      double singlet_score_value = 0.0;
      if (!type1_pres.all_pairs_class && !type2_pres.all_pairs_class) {
        spot = 0;
        singlet_score_value = min_score + 2.0 * doublet_threshold;
      } else if (type1_pres.all_pairs_class && !type2_pres.all_pairs_class) {
        first_class_value = !type1_pres.all_pairs;
        singlet_score_value = type1_pres.singlet_score;
        spot = 3;
      } else if (!type1_pres.all_pairs_class && type2_pres.all_pairs_class) {
        first_class_value = !type2_pres.all_pairs;
        singlet_score_value = type2_pres.singlet_score;
        std::swap(first, second);
        spot = 3;
      } else {
        spot = 2;
        singlet_score_value = std::min(type1_pres.singlet_score, type2_pres.singlet_score);
        first_class_value = !type1_pres.all_pairs;
        second_class_value = !type2_pres.all_pairs;
        if (type2_pres.singlet_score < type1_pres.singlet_score) {
          std::swap(first, second);
          first_class_value = !type2_pres.all_pairs;
          second_class_value = !type1_pres.all_pairs;
        }
      }

      if (singlet_score_value - min_score < doublet_threshold) {
        spot = 1;
      }

      FitResult final_fit = solve_irwls(y, {first, second}, ref, n_genes, n_types, nUMI[cell],
                                        tables, 50, min_change, constrain);
      const double weight_sum = std::accumulate(final_fit.weights.begin(), final_fit.weights.end(), 0.0);
      double w0 = final_fit.weights[0] / weight_sum;
      double w1 = final_fit.weights[1] / weight_sum;

      spot_class[cell] = spot;
      first_type[cell] = first + 1;
      second_type[cell] = second + 1;
      first_class[cell] = first_class_value;
      second_class[cell] = second_class_value;
      first_weight[cell] = w0;
      second_weight[cell] = w1;
      min_score_out[cell] = min_score;
      singlet_score_out[cell] = singlet_score_value;
      conv_doublet[cell] = final_fit.converged;

      if (return_diagnostics) {
        score_mats[cell] = score_mat;
        singlet_scores_out[cell] = singlet_scores;
        diagnostic_candidates[cell] = candidates;
      }
    }
  }
};

template <typename T>
std::vector<T> matrix_to_vector(const Rcpp::Matrix<REALSXP>& matrix) {
  return std::vector<T>(matrix.begin(), matrix.end());
}

std::vector<double> numeric_matrix_to_vector(const NumericMatrix& matrix) {
  return std::vector<double>(matrix.begin(), matrix.end());
}

} // namespace

// [[Rcpp::export]]
List fast_doublet_chunk_cpp(NumericMatrix beads,
                            NumericVector nUMI,
                            NumericMatrix ref_profiles,
                            NumericMatrix all_weights,
                            LogicalVector conv_all,
                            IntegerVector class_ids,
                            NumericMatrix Q_mat,
                            NumericMatrix SQ_mat,
                            NumericVector X_vals,
                            int K_val,
                            bool constrain,
                            double min_change,
                            double confidence_threshold,
                            double doublet_threshold,
                            bool return_diagnostics,
                            int max_cores) {
  const int n_cells = beads.nrow();
  const int n_genes = beads.ncol();
  const int ref_genes = ref_profiles.nrow();
  const int n_types = ref_profiles.ncol();
  if (ref_genes != n_genes) {
    Rcpp::stop("fast_doublet_chunk_cpp: ref_profiles rows must match bead columns.");
  }
  if (all_weights.nrow() != n_cells || all_weights.ncol() != n_types) {
    Rcpp::stop("fast_doublet_chunk_cpp: all_weights dimensions must be cells x types.");
  }
  if (n_types < 2) {
    Rcpp::stop("fast_doublet_chunk_cpp: at least two cell types are required.");
  }

  LikelihoodTables tables;
  tables.q = numeric_matrix_to_vector(Q_mat);
  tables.sq = numeric_matrix_to_vector(SQ_mat);
  tables.x = std::vector<double>(X_vals.begin(), X_vals.end());
  tables.n_row = Q_mat.nrow();
  tables.n_col = Q_mat.ncol();
  tables.k_val = K_val;

  std::vector<double> beads_vec = numeric_matrix_to_vector(beads);
  std::vector<double> nUMI_vec(nUMI.begin(), nUMI.end());
  std::vector<double> ref_vec = numeric_matrix_to_vector(ref_profiles);
  std::vector<double> all_weights_vec = numeric_matrix_to_vector(all_weights);
  std::vector<int> class_vec(class_ids.begin(), class_ids.end());
  std::vector<unsigned char> conv_all_vec(conv_all.size());
  for (R_xlen_t i = 0; i < conv_all.size(); ++i) {
    conv_all_vec[i] = static_cast<unsigned char>(conv_all[i] == TRUE);
  }

  std::vector<int> spot_class(n_cells);
  std::vector<int> first_type(n_cells);
  std::vector<int> second_type(n_cells);
  std::vector<unsigned char> first_class(n_cells);
  std::vector<unsigned char> second_class(n_cells);
  std::vector<double> first_weight(n_cells);
  std::vector<double> second_weight(n_cells);
  std::vector<double> min_score_out(n_cells);
  std::vector<double> singlet_score_out(n_cells);
  std::vector<unsigned char> conv_doublet(n_cells);
  std::vector<std::vector<double>> score_mats(n_cells);
  std::vector<std::vector<double>> singlet_scores(n_cells);
  std::vector<std::vector<int>> diagnostic_candidates(n_cells);

  FastDoubletWorker worker(beads_vec, nUMI_vec, ref_vec, all_weights_vec, conv_all_vec,
                           class_vec, tables, n_cells, n_genes, n_types, constrain,
                           min_change, confidence_threshold, doublet_threshold,
                           return_diagnostics, spot_class, first_type, second_type,
                           first_class, second_class, first_weight, second_weight,
                           min_score_out, singlet_score_out, conv_doublet, score_mats,
                           singlet_scores, diagnostic_candidates);

  const int n_threads = std::max(1, std::min(max_cores, n_cells));
  std::vector<std::thread> threads;
  threads.reserve(n_threads);
  for (int thread = 0; thread < n_threads; ++thread) {
    const std::size_t begin = static_cast<std::size_t>(thread) * n_cells / n_threads;
    const std::size_t end = static_cast<std::size_t>(thread + 1) * n_cells / n_threads;
    threads.emplace_back([&worker, begin, end]() {
      worker.run(begin, end);
    });
  }
  for (auto& thread : threads) {
    thread.join();
  }

  CharacterVector spot_labels(n_cells);
  const char* labels[] = {"reject", "singlet", "doublet_certain", "doublet_uncertain"};
  for (int i = 0; i < n_cells; ++i) {
    spot_labels[i] = labels[spot_class[i]];
  }

  LogicalVector first_class_out(n_cells), second_class_out(n_cells), conv_doublet_out(n_cells);
  for (int i = 0; i < n_cells; ++i) {
    first_class_out[i] = first_class[i] != 0;
    second_class_out[i] = second_class[i] != 0;
    conv_doublet_out[i] = conv_doublet[i] != 0;
  }

  List out = List::create(
    Rcpp::Named("spot_class") = spot_labels,
    Rcpp::Named("first_type") = IntegerVector(first_type.begin(), first_type.end()),
    Rcpp::Named("second_type") = IntegerVector(second_type.begin(), second_type.end()),
    Rcpp::Named("first_class") = first_class_out,
    Rcpp::Named("second_class") = second_class_out,
    Rcpp::Named("first_weight") = NumericVector(first_weight.begin(), first_weight.end()),
    Rcpp::Named("second_weight") = NumericVector(second_weight.begin(), second_weight.end()),
    Rcpp::Named("min_score") = NumericVector(min_score_out.begin(), min_score_out.end()),
    Rcpp::Named("singlet_score") = NumericVector(singlet_score_out.begin(), singlet_score_out.end()),
    Rcpp::Named("conv_doublet") = conv_doublet_out
  );

  if (return_diagnostics) {
    List score_list(n_cells), singlet_list(n_cells);
    for (int cell = 0; cell < n_cells; ++cell) {
      const int n_cand = static_cast<int>(diagnostic_candidates[cell].size());
      NumericMatrix score_mat(n_cand, n_cand);
      std::copy(score_mats[cell].begin(), score_mats[cell].end(), score_mat.begin());
      IntegerVector candidate_types(n_cand);
      for (int j = 0; j < n_cand; ++j) {
        candidate_types[j] = diagnostic_candidates[cell][j] + 1;
      }
      score_mat.attr("candidate_type_index") = candidate_types;
      score_list[cell] = score_mat;
      singlet_list[cell] = NumericVector(singlet_scores[cell].begin(), singlet_scores[cell].end());
    }
    out["score_mat"] = score_list;
    out["singlet_scores"] = singlet_list;
  }

  return out;
}

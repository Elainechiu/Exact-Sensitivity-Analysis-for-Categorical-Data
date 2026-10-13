/*
 ================================================================================
 sensitivityIxJ.cpp
 
 Code for "Exact, Nonparametric Sensitivity Analysis for Observational Studies
 of Contingency Tables"
 Authors: Elaine K. Chiu and Hyunseung Kang
 
 C++ routines of sensitivityIxJ: exact sensitivity analysis for I x J
 contingency tables under the generic bias sensitivity model. The file contains
 
 - the exact computation of the worst-case null distribution for any I x J
 table and any count-based test statistic, and
 - the exact first and second moments of the cell counts, used by the
 normal approximation.
 ================================================================================
 
 SETTING
 Rows are treatment levels with margins N_i., columns are outcome levels with
 margins N_.j, and N is the total count. Treatment assignment follows
 
 P(Z_s = i | F)  proportional to  exp(xi_i(x_s) + gamma * delta_i * u_s),
 gamma >= 0,  delta_i in {0, 1},  u_s in [0, 1].
 
 The worst case is attained at a binary u, and the null distribution depends
 on u only through
 u-bar_j = #{s : r_s = j, u_s = 1},   j = 1, ..., J.
 Every argument named u below is this vector u-bar.
 
 NOTATION
 I_1      = {i : delta_i = 1}
 N-bar    = sum_{i in I_1} N_i.
 t-bar_j  = sum_{i in I_1} t_ij
 u-bar    = sum_j u-bar_j, the number of subjects with u_s = 1
 q-bar    = number of subjects with u_s = 1 assigned to a treatment in I_1
 n-bar_j  = number of subjects with outcome j and u_s = 1 assigned to a
 treatment in I_1, so that q-bar = sum_j n-bar_j
 Omega_t  = set of tables with the observed margins whose statistic is at or
 beyond the observed value
 
 EXACT P-VALUE
 The exact p-value at a candidate u is
 
 alpha(u) = sum_{t-bar} W(t-bar) S(t-bar, u) / (kappa_0 S_0(u)).
 
 Each term is defined below with the function that computes it. All
 functions work on the log scale and return the log of the term.
 
 kappa(t)      Count of the assignments of one table t, apart from the
 factor that depends on u.
 kappa(t) = prod_j t-bar_j! (N_.j - t-bar_j)! / prod_ij t_ij!
 Computed by  log_c_t
 
 W(t-bar)      Sum of kappa(t) over the tables in Omega_t that share the
 same t-bar. Does not depend on u or gamma.
 W(t-bar) = sum_{t in Omega_t : t-bar(t) = t-bar} kappa(t)
 Computed by  build_W_cpp (enumerates Omega_t), or
 build_W_from_tables (from a given list of tables)
 
 S(t-bar, u)   Factor of the numerator that depends on u, shared by all
 tables with the same t-bar.
 S(t-bar, u) = sum_{q-bar} exp(gamma q-bar)
 sum_{n-bar_1 + ... + n-bar_J = q-bar}
 prod_j C(u-bar_j, n-bar_j)
 C(N_.j - u-bar_j, t-bar_j - n-bar_j)
 Computed by  log_S
 
 kappa_0       Factor of the denominator that does not depend on u.
 kappa_0 = N-bar! (N - N-bar)! / prod_i N_i.!
 Computed by  log_c_0
 
 S_0(u)        Factor of the denominator that depends on u.
 S_0(u) = sum_{q-bar} exp(gamma q-bar)
 C(u-bar, q-bar) C(N - u-bar, N-bar - q-bar)
 Computed by  log_den_bd
 
 alpha(u)      The exact p-value, combining the terms above.
 Computed by  alpha_bd, at every candidate u
 
 The probability of a single table t at u is
 P_u(N = t) = kappa(t) S(t-bar, u) / (kappa_0 S_0(u)),
 computed by point_prob_bd.
 
 The computation has two stages:
 Stage 1  build_W_cpp enumerates Omega_t once and returns W(t-bar).
 Stage 2  alpha_bd evaluates alpha(u) at every candidate u.
 Because W does not depend on u or gamma, stage 1 is run only once.
 
 BUILT-IN TEST STATISTICS (evaluated in C++ inside the enumeration)
 ordinal         T = sum_ij w_i v_j N_ij                  (score matrix)
 chi-square      T = sum_ij (N_ij - E_ij)^2 / E_ij,   E_ij = N_i. N_.j / N
 G^2             T = 2 sum_ij N_ij log(N_ij / E_ij)
 Freeman-Halton  tables ordered by their null probability
 P(t) = prod_i N_i.! prod_j N_.j! / (N! prod_ij t_ij!)
 Any other count-based statistic can be supplied as an R function.
 ================================================================================
 */

#include <Rcpp.h>
#include <cmath>
#include <vector>
#include <map>
#include <limits>
#include <algorithm>
using namespace Rcpp;

static const double NEG_INF = -std::numeric_limits<double>::infinity();

// log C(n, k), with C = 0 (log = -Inf) outside 0 <= k <= n
// =============================================================================
// FUNCTION  lchoose_                                            [internal helper]
// PURPOSE   log of the binomial coefficient C(n, k), returning -Inf (log 0) when
//           k < 0, n < 0 or k > n.  Every counting kernel in Part A is a product
//           of binomial coefficients, evaluated in log space to avoid overflow.
// PAPER     Building block of the kernels in Theorem 4 / Corollary C.3.
// INPUT     n, k  integers.
// OUTPUT    double, log C(n, k).
// EXAMPLE   (C++)  lchoose_(5, 2)  ==  std::log(10.0);   lchoose_(3, 5)  ==  -Inf
// =============================================================================
static inline double lchoose_(int n, int k) {
  if (n < 0 || k < 0 || k > n) return NEG_INF;
  return std::lgamma(n + 1.0) - std::lgamma(k + 1.0) - std::lgamma(n - k + 1.0);
}
// log(e^a + e^b), stable
// =============================================================================
// FUNCTION  lse                                                 [internal helper]
// PURPOSE   log(exp(a) + exp(b)) computed stably (log-sum-exp of two terms), with
//           -Inf treated as log 0.  Used to accumulate sums of kernels in log space.
// INPUT     a, b  log-scale values.
// OUTPUT    double.
// EXAMPLE   (C++)  lse(std::log(2.0), std::log(3.0))  ==  std::log(5.0)
// =============================================================================
static inline double lse(double a, double b) {
  if (a == NEG_INF) return b;
  if (b == NEG_INF) return a;
  double m = std::max(a, b);
  return m + std::log(std::exp(a - m) + std::exp(b - m));
}

/*
 --------------------------------------------------------------------------------
 log_S(u, m, TD, J, gamma)                                       internal helper
 
 Log of S(t-bar, u), the factor of the numerator that depends on u. It is
 shared by all tables with the same t-bar:
 
 S(t-bar, u) = sum_{q-bar} exp(gamma q-bar)
 sum_{n-bar_1 + ... + n-bar_J = q-bar}
 prod_j C(u-bar_j, n-bar_j) C(N_.j - u-bar_j, t-bar_j - n-bar_j)
 
 For each outcome level j, n-bar_j ranges over
 max(0, t-bar_j - (N_.j - u-bar_j)) <= n-bar_j <= min(u-bar_j, t-bar_j),
 so that no more subjects with u_s = 1, or with u_s = 0, are used than exist.
 The J factors are combined by a convolution over q-bar on the log scale, and
 each q-bar is then weighted by exp(gamma q-bar).
 
 Input   u      u-bar_j, length J
 m      outcome margins N_.j, length J
 TD     t-bar_j, length J
 J      number of outcome levels
 gamma  sensitivity parameter
 Output  log S(t-bar, u); -Inf if no admissible n-bar exists.
 Note    At gamma = 0, S(t-bar, u) = prod_j C(N_.j, t-bar_j).
 --------------------------------------------------------------------------------
 */

static double log_S(const int* u, const int* m, const int* TD, int J, double gamma) {
  std::vector<double> g(1, 0.0);                          // log-convolution so far, index = q-bar
  for (int j = 0; j < J; ++j) {
    int lo = std::max(0, TD[j] - (m[j] - u[j]));
    int hi = std::min(u[j], TD[j]);
    if (lo > hi) return NEG_INF;                          // no admissible split: contributes 0
    std::vector<double> f(hi - lo + 1);
    for (int a = lo; a <= hi; ++a)
      f[a - lo] = lchoose_(u[j], a) + lchoose_(m[j] - u[j], TD[j] - a);
    std::vector<double> h(g.size() + hi, NEG_INF);
    for (size_t q = 0; q < g.size(); ++q) {
      if (g[q] == NEG_INF) continue;
      for (int a = lo; a <= hi; ++a) h[q + a] = lse(h[q + a], g[q] + f[a - lo]);
    }
    g.swap(h);
  }
  double out = NEG_INF;
  for (size_t Q = 0; Q < g.size(); ++Q)
    if (g[Q] != NEG_INF) out = lse(out, g[Q] + gamma * (double)Q);
  return out;
}

/*
 --------------------------------------------------------------------------------
 log_den_bd(ubar, N, ND, gamma)                                  
 
 Log of S_0(u-bar), the factor of the denominator that depends on u:
 
 S_0(u-bar) = sum_{q-bar} exp(gamma q-bar) C(u-bar, q-bar) C(N - u-bar, N-bar - q-bar),
 
 where q-bar ranges over
 max(0, N-bar - (N - u-bar)) <= q-bar <= min(u-bar, N-bar).
 S_0 depends on u only through the total u-bar = sum_j u-bar_j.
 
 Input   ubar   u-bar, the number of subjects with u_s = 1
 N      total number of subjects
 ND     N-bar, the number of subjects assigned to a treatment in I_1
 gamma  sensitivity parameter
 Output  log S_0(u-bar).
 Note    At gamma = 0, S_0(u-bar) = C(N, N-bar) for every u-bar.
 --------------------------------------------------------------------------------
 */
// [[Rcpp::export]]
double log_den_bd(int ubar, int N, int ND, double gamma) {
  double out = NEG_INF;
  int lo = std::max(0, ND - (N - ubar)), hi = std::min(ubar, ND);
  for (int Q = lo; Q <= hi; ++Q)
    out = lse(out, lchoose_(ubar, Q) + lchoose_(N - ubar, ND - Q) + gamma * (double)Q);
  return out;
}

/*
 --------------------------------------------------------------------------------
 log_c_t(tab, delta)                                             
 
 Log of kappa(t), the factor of the numerator for one table t that does not
 depend on u:
 
 kappa(t) = prod_j t-bar_j! (N_.j - t-bar_j)! / prod_i t_ij!,
 
 where t-bar_j = sum_{i in I_1} t_ij. Within outcome level j, kappa(t) counts
 the ways to divide the t-bar_j subjects assigned to I_1 among the treatments
 in I_1, and the other N_.j - t-bar_j subjects among the remaining treatments,
 so that the cell counts equal t_ij.
 
 Input   tab    I x J integer table; rows are treatment levels
 delta  binary vector of length I marking the treatments in I_1
 Output  log kappa(t).
 Note    If delta = 0, then t-bar_j = 0 and kappa(t) = prod_j N_.j! / prod_i t_ij!.
 --------------------------------------------------------------------------------
 */

// [[Rcpp::export]]
double log_c_t(IntegerMatrix tab, IntegerVector delta) {
  int I = tab.nrow(), J = tab.ncol();
  double out = 0.0;
  for (int j = 0; j < J; ++j) {
    int TD = 0, m = 0;
    for (int i = 0; i < I; ++i) {
      m += tab(i, j); if (delta[i] == 1) TD += tab(i, j);
      out -= std::lgamma(tab(i, j) + 1.0);
    }
    out += std::lgamma(TD + 1.0) + std::lgamma(m - TD + 1.0);
  }
  return out;
}

/*
 --------------------------------------------------------------------------------
 log_c_0(treatment_margins, delta)
 
 Log of kappa_0, the factor of the denominator that does not depend on u:
 
 kappa_0 = N-bar! (N - N-bar)! / prod_i N_i.!,
 
 where N-bar = sum_{i in I_1} N_i.. Given which N-bar subjects are assigned to
 I_1, kappa_0 counts the ways to divide them among the treatments in I_1, and
 the other N - N-bar subjects among the remaining treatments, so that the
 treatment counts equal N_i..
 
 Input   treatment_margins  N_i., length I
 delta              binary vector of length I marking the treatments in I_1
 Output  log kappa_0.
 Note    If delta = 0, then N-bar = 0 and kappa_0 = N! / prod_i N_i.!.
 --------------------------------------------------------------------------------
 */
// [[Rcpp::export]]
double log_c_0(IntegerVector treatment_margins, IntegerVector delta) {
  int I = treatment_margins.size(), N = 0, ND = 0;
  double out = 0.0;
  for (int i = 0; i < I; ++i) {
    N += treatment_margins[i]; if (delta[i] == 1) ND += treatment_margins[i];
    out -= std::lgamma(treatment_margins[i] + 1.0);
  }
  return out + std::lgamma(ND + 1.0) + std::lgamma(N - ND + 1.0);
}

/*
 --------------------------------------------------------------------------------
 point_prob_bd(tab, delta, u, gamma)
 
 Probability of one single table t under the sensitivity model at a binary u:
 
 P_u(N = t) = kappa(t) S(t-bar, u) / (kappa_0 S_0(u-bar)),
 
 computed on the log scale from log_c_t, log_S, log_c_0 and log_den_bd. The
 margins N_i., N_.j, N, N-bar, t-bar_j and u-bar are taken from tab, delta
 and u.
 
 Input   tab    I x J integer table; rows are treatment levels
 delta  binary vector of length I marking the treatments in I_1
 u      u-bar_j, length J, with 0 <= u-bar_j <= N_.j
 gamma  sensitivity parameter
 Output  P_u(N = t), a probability in [0, 1].
 Errors  if u does not have length J, if delta does not have length I, or if
 some u-bar_j lies outside [0, N_.j].
 Note    At gamma = 0, or when u-bar_j = 0 for all j, or u-bar_j = N_.j for all
 j, P_u(N = t) is the multivariate hypergeometric probability
 prod_i N_i.! prod_j N_.j! / (N! prod_ij t_ij!).
 --------------------------------------------------------------------------------
 */
// [[Rcpp::export]]
double point_prob_bd(IntegerMatrix tab, IntegerVector delta, IntegerVector u, double gamma) {
  int I = tab.nrow(), J = tab.ncol();
  if (u.size() != J) stop("u must have length J (number of outcome levels).");
  if (delta.size() != I) stop("delta must have length I (number of treatment levels).");
  std::vector<int> m(J, 0), TD(J, 0), uu(J);
  IntegerVector n(I);
  int N = 0, ND = 0, ubar = 0;
  for (int j = 0; j < J; ++j) {
    uu[j] = u[j]; ubar += u[j];
    for (int i = 0; i < I; ++i) { m[j] += tab(i, j); if (delta[i] == 1) TD[j] += tab(i, j); }
    if (u[j] < 0 || u[j] > m[j]) stop("u_j must lie in [0, N_.j].");
  }
  for (int i = 0; i < I; ++i) {
    n[i] = 0; for (int j = 0; j < J; ++j) n[i] += tab(i, j);
    N += n[i]; if (delta[i] == 1) ND += n[i];
  }
  double lnum = log_c_t(tab, delta) + log_S(uu.data(), m.data(), TD.data(), J, gamma);
  double lden = log_c_0(n, delta) + log_den_bd(ubar, N, ND, gamma);
  return std::exp(lnum - lden);
}

/*
 --------------------------------------------------------------------------------
 alpha_bd(TDmat, logW, u_space, m, gamma, N, ND, logc0)
 
 Exact p-value at every candidate u-bar, given W(t-bar) from build_W_cpp or
 build_W_from_tables:
 
 alpha(u) = sum_r W(t-bar_r) S(t-bar_r, u) / (kappa_0 S_0(u-bar)),
 
 where the sum runs over the distinct t-bar vectors (rows of TDmat). W does not
 depend on u or gamma, so one W serves every candidate and every gamma.
 
 Input   TDmat    distinct t-bar vectors, one per row (J columns)
 logW     log W(t-bar), one entry per row of TDmat
 u_space  candidate u-bar vectors, one per row (J columns), with
 0 <= u-bar_j <= N_.j
 m        outcome margins N_.j, length J
 gamma    sensitivity parameter
 N        total number of subjects
 ND       N-bar, the number of subjects assigned to a treatment in I_1
 logc0    log kappa_0, from log_c_0
 Output  vector of alpha(u), one entry per row of u_space.
 Errors  if TDmat or u_space does not have J columns, if logW does not have one
 entry per row of TDmat, or if some u-bar_j lies outside [0, N_.j].
 Note    The row u-bar = 0 gives the p-value under random assignment.
 --------------------------------------------------------------------------------
 */
// [[Rcpp::export]]
NumericVector alpha_bd(IntegerMatrix TDmat, NumericVector logW, IntegerMatrix u_space,
                       IntegerVector m, double gamma, int N, int ND, double logc0) {
  int J = m.size(), R = TDmat.nrow(), U = u_space.nrow();
  if (TDmat.ncol() != J || u_space.ncol() != J) stop("TDmat and u_space must have J columns.");
  if (logW.size() != R) stop("logW must have one entry per row of TDmat.");
  std::vector<int> mm(J), TD(J), uu(J);
  for (int j = 0; j < J; ++j) mm[j] = m[j];
  NumericVector out(U);
  for (int k = 0; k < U; ++k) {
    int ubar = 0;
    for (int j = 0; j < J; ++j) {
      uu[j] = u_space(k, j); ubar += uu[j];
      if (uu[j] < 0 || uu[j] > mm[j]) stop("u_space row %d: u_j outside [0, N_.j].", k + 1);
    }
    double lden = logc0 + log_den_bd(ubar, N, ND, gamma);
    double acc = NEG_INF;
    for (int r = 0; r < R; ++r) {
      for (int j = 0; j < J; ++j) TD[j] = TDmat(r, j);
      double ls = log_S(uu.data(), mm.data(), TD.data(), J, gamma);
      if (ls != NEG_INF) acc = lse(acc, logW[r] + ls);
    }
    out[k] = (acc == NEG_INF) ? 0.0 : std::exp(acc - lden);
    if ((k & 1023) == 0) Rcpp::checkUserInterrupt();
  }
  return out;
}

/*
 --------------------------------------------------------------------------------
 Enumerator
 
 Stage 1 enumeration used by build_W_cpp. Visits every I x J table of
 nonnegative integers with row margins N_i. and column margins N_.j, and builds
 
 W(t-bar) = sum_{t in Omega_t : t-bar(t) = t-bar} kappa(t),
 
 where Omega_t is the set of tables with T(t) >= c (or T(t) <= c when greater
 is FALSE).
 
 Cells are filled row by row, each within the capacity left in its row and
 column. The last cell of each row and the whole last row are fixed by the
 margins, so each table is visited exactly once. For each complete table, the
 statistic T(t) is evaluated, and if t is in Omega_t, kappa(t) is added to
 W(t-bar) on the log scale. Tables are stored only when keep_tables is TRUE.
 
 Members
 linear_stat, chisq_stat, g2_stat, freeman_halton_stat
 built-in statistics in their reduced forms (see below)
 passes     TRUE if T(t) >= c - tol, or T(t) <= c + tol when greater is FALSE;
 uses the R function stat_fun when no built-in statistic is set
 leaf       processes one complete table: counts it, tests membership in
 Omega_t, computes t-bar and log kappa(t), and adds kappa(t) to
 W(t-bar); stores the table if keep_tables is TRUE
 rec        fills cell (i, j) with every feasible value and moves to the
 next cell
 Counters
 n_visited  number of tables with the given margins
 n_tables   number of tables in Omega_t
 --------------------------------------------------------------------------------
 */
struct Enumerator {
  int I, J;
  std::vector<int> n, m, delta;        // row margins, col margins, delta (marks I_1)
  bool use_score; std::vector<double> score;   // row-major I x J score matrix w_i v_j
  int stat_type = 0;                   // 0 = R callback, 1 = linear (score), 2 = chi-square,
                                       // 3 = G^2 (likelihood ratio), 4 = Freeman-Halton (probability ordering)
  std::vector<double> inv_E;           // row-major 1/E_ij for the chi-square statistic
  std::vector<double> klogk;           // k log k, k = 0..N, for G^2
  std::vector<double> lfact;           // log k!,  k = 0..N, for the Freeman-Halton statistic
  Function stat_fun;                   // R callback used when no score matrix is given
  double threshold; bool greater; double tol;
  bool keep;
  std::vector<int> t;                  // current table, row-major
  std::vector<int> crem;               // remaining column capacity
  std::map<std::vector<int>, double> W;   // t-bar -> log W(t-bar)
  long long n_tables = 0, n_visited = 0;
  std::vector<int> kept;               // kept tables, flattened row-major

  Enumerator(Function f) : stat_fun(f) {}

  double linear_stat() const {
    double s = 0.0;
    for (int k = 0; k < I * J; ++k) s += score[k] * t[k];
    return s;
  }
  /*
   Built-in statistics in reduced form. With the margins fixed, each reduced
   form is a monotone function of the original statistic, so it selects the
   same tables. The threshold c must be given in the same reduced form.
   chi-square      sum_ij (t_ij - E_ij)^2 / E_ij = sum_ij t_ij^2 / E_ij - N
   reduced form: sum_ij t_ij^2 / E_ij
   G^2             2 sum_ij t_ij log(t_ij / E_ij) = 2 (sum_ij t_ij log t_ij - C),
   C = sum_i N_i. log N_i. + sum_j N_.j log N_.j - N log N
   reduced form: sum_ij t_ij log t_ij
   Freeman-Halton  tables ordered by the null probability
   P(t) = prod_i N_i.! prod_j N_.j! / (N! prod_ij t_ij!);
   P(t) <= P(t_obs) exactly when
   sum_ij log t_ij! >= sum_ij log t_obs,ij!
   reduced form: sum_ij log t_ij!
   */
  
  double chisq_stat() const {
    double s = 0.0;
    for (int k = 0; k < I * J; ++k) s += inv_E[k] * (double)t[k] * (double)t[k];
    return s;
  }
  double g2_stat() const {
    double s = 0.0;
    for (int k = 0; k < I * J; ++k) s += klogk[t[k]];
    return s;
  }
  double freeman_halton_stat() const {
    double s = 0.0;
    for (int k = 0; k < I * J; ++k) s += lfact[t[k]];
    return s;
  }
  // membership in Omega_t; tol guards against floating-point ties (T(t) == c)
  bool passes() {
    double val;
    if (stat_type == 1) val = linear_stat();
    else if (stat_type == 2) val = chisq_stat();
    else if (stat_type == 3) val = g2_stat();
    else if (stat_type == 4) val = freeman_halton_stat();
    else {
      IntegerMatrix tab(I, J);
      for (int i = 0; i < I; ++i) for (int j = 0; j < J; ++j) tab(i, j) = t[i * J + j];
      val = as<double>(stat_fun(tab));
    }
    return greater ? (val >= threshold - tol) : (val <= threshold + tol);
  }
  // a complete table: test membership, then fold kappa(t) into W(t-bar)
  void leaf() {
    ++n_visited;
    if (!passes()) return;
    ++n_tables;
    std::vector<int> TD(J, 0);
    double lc = 0.0;
    for (int j = 0; j < J; ++j) {
      for (int i = 0; i < I; ++i) {
        int v = t[i * J + j];
        if (delta[i] == 1) TD[j] += v;
        lc -= std::lgamma(v + 1.0);
      }
      lc += std::lgamma(TD[j] + 1.0) + std::lgamma(m[j] - TD[j] + 1.0);
    }
    auto it = W.find(TD);
    if (it == W.end()) W.emplace(TD, lc); else it->second = lse(it->second, lc);
    if (keep) kept.insert(kept.end(), t.begin(), t.end());
    if ((n_visited & 65535) == 0) Rcpp::checkUserInterrupt();
  }
  // fill cell (i, j); rrem = remaining capacity of row i
  void rec(int i, int j, int rrem) {
    if (i == I - 1) {                              // last row forced by column remainders
      int sum = 0;
      for (int jj = 0; jj < J; ++jj) { t[i * J + jj] = crem[jj]; sum += crem[jj]; }
      if (sum == n[i]) leaf();                     // always true when margins are consistent
      return;
    }
    if (j == J - 1) {                              // last column of row i forced
      if (rrem > crem[j]) return;
      t[i * J + j] = rrem; crem[j] -= rrem;
      rec(i + 1, 0, n[i + 1]);
      crem[j] += rrem;
      return;
    }
    int cap_after = 0;                             // capacity left in columns > j
    for (int jj = j + 1; jj < J; ++jj) cap_after += crem[jj];
    int lo = std::max(0, rrem - cap_after), hi = std::min(rrem, crem[j]);
    for (int v = lo; v <= hi; ++v) {
      t[i * J + j] = v; crem[j] -= v;
      rec(i, j + 1, rrem - v);
      crem[j] += v;
    }
  }
};

/*
 --------------------------------------------------------------------------------
 build_W_cpp(row_margins, col_margins, delta, threshold, greater, score,
 stat_fun, tol, keep_tables, statistic)
 
 Stage 1: enumerates Omega_t and returns W(t-bar) for every distinct t-bar,
 
 W(t-bar) = sum_{t in Omega_t : t-bar(t) = t-bar} kappa(t).
 
 W does not depend on u or gamma, so one call serves every candidate u and
 every gamma in alpha_bd.
 
 The statistic is chosen as follows:
 - if score is supplied, the linear statistic T(t) = sum_ij score_ij t_ij;
 - otherwise, if statistic is 2, 3 or 4, the corresponding built-in statistic;
 - otherwise, the R function stat_fun.
 Built-in statistics are evaluated in reduced form, and threshold must be given
 in the same form:
 chi-square      sum_ij t_ij^2 / E_ij,   E_ij = N_i. N_.j / N
 G^2             sum_ij t_ij log t_ij
 Freeman-Halton  sum_ij log t_ij!   (larger value = smaller null probability)
 
 Input   row_margins  N_i., length I
 col_margins  N_.j, length J
 delta        binary vector of length I marking the treatments in I_1
 threshold    critical value c
 greater      TRUE keeps T(t) >= c, FALSE keeps T(t) <= c
 (default TRUE)
 score        optional I x J score matrix
 stat_fun     optional R function of an I x J table
 tol          tolerance in the comparison with c (default 1e-9)
 keep_tables  if TRUE, also return the tables of Omega_t
 (default FALSE)
 statistic    0 = score or stat_fun, 2 = chi-square, 3 = G^2,
 4 = Freeman-Halton (default 0)
 Output  list with
 TDmat      distinct t-bar vectors, one per row
 logW       log W(t-bar), one entry per row of TDmat
 n_tables   number of tables in Omega_t
 n_visited  number of tables with the given margins
 tables     (if keep_tables is TRUE) tables of Omega_t, one per row,
 with cells stored row by row
 Errors  if delta does not have length I, if the row and column margins have
 different sums, if no statistic is supplied, or if score is not I x J.
 --------------------------------------------------------------------------------
 */
// [[Rcpp::export]]
List build_W_cpp(IntegerVector row_margins, IntegerVector col_margins, IntegerVector delta,
                 double threshold, bool greater = true,
                 Nullable<NumericMatrix> score = R_NilValue,
                 Nullable<Function> stat_fun = R_NilValue,
                 double tol = 1e-9, bool keep_tables = false, int statistic = 0) {
  int I = row_margins.size(), J = col_margins.size();
  if (delta.size() != I) stop("delta must have length I.");
  long long sr = 0, sc = 0;
  for (int i = 0; i < I; ++i) sr += row_margins[i];
  for (int j = 0; j < J; ++j) sc += col_margins[j];
  if (sr != sc) stop("row and column margins must sum to the same N.");
  if (score.isNull() && stat_fun.isNull() && statistic == 0) stop("supply score, a built-in statistic code, or stat_fun.");

  Function f = stat_fun.isNull() ? Function("identity") : Function(stat_fun.get());
  Enumerator E(f);
  E.I = I; E.J = J;
  E.n.assign(row_margins.begin(), row_margins.end());
  E.m.assign(col_margins.begin(), col_margins.end());
  E.delta.assign(delta.begin(), delta.end());
  E.use_score = !score.isNull();
  E.stat_type = E.use_score ? 1 : statistic;
  if (E.stat_type == 2) {
    E.inv_E.resize(I * J);
    for (int i = 0; i < I; ++i) for (int j = 0; j < J; ++j)
      E.inv_E[i * J + j] = (double)sr / ((double)row_margins[i] * (double)col_margins[j]);
  }
  if (E.stat_type == 3) { E.klogk.resize(sr + 1); E.klogk[0] = 0.0; for (long long k = 1; k <= sr; ++k) E.klogk[k] = k * std::log((double)k); }
  if (E.stat_type == 4) { E.lfact.resize(sr + 1); for (long long k = 0; k <= sr; ++k) E.lfact[k] = std::lgamma(k + 1.0); }
  if (E.use_score) {
    NumericMatrix S(score.get());
    if (S.nrow() != I || S.ncol() != J) stop("score must be I x J.");
    E.score.resize(I * J);
    for (int i = 0; i < I; ++i) for (int j = 0; j < J; ++j) E.score[i * J + j] = S(i, j);
  }
  E.threshold = threshold; E.greater = greater; E.tol = tol; E.keep = keep_tables;
  E.t.assign(I * J, 0); E.crem = E.m;
  E.rec(0, 0, E.n[0]);

  int R = E.W.size();
  IntegerMatrix TDmat(R, J); NumericVector logW(R);
  int r = 0;
  for (auto& kv : E.W) { for (int j = 0; j < J; ++j) TDmat(r, j) = kv.first[j]; logW[r] = kv.second; ++r; }
  List out = List::create(_["TDmat"] = TDmat, _["logW"] = logW,
                          _["n_tables"] = (double)E.n_tables, _["n_visited"] = (double)E.n_visited);
  if (keep_tables) {
    IntegerMatrix tabs(E.n_tables, I * J);
    for (long long k = 0; k < E.n_tables; ++k) for (int c = 0; c < I * J; ++c) tabs(k, c) = E.kept[k * I * J + c];
    out["tables"] = tabs;
  }
  return out;
}

/*
 --------------------------------------------------------------------------------
 build_W_from_tables(tabs, I, J, delta)
 
 Stage 1 from a given list of the tables in Omega_t. Computes t-bar and
 kappa(t) for each table and returns W(t-bar) for every distinct t-bar,
 
 W(t-bar) = sum_{t in Omega_t : t-bar(t) = t-bar} kappa(t).
 
 Gives the same TDmat and logW as build_W_cpp when tabs holds the tables of
 Omega_t, for example the tables returned by build_W_cpp with keep_tables = TRUE.
 
 Input   tabs   integer matrix with one table per row and I * J columns, cells
 stored row by row
 I      number of treatment levels
 J      number of outcome levels
 delta  binary vector of length I marking the treatments in I_1
 Output  list with
 TDmat     distinct t-bar vectors, one per row
 logW      log W(t-bar), one entry per row of TDmat
 n_tables  number of tables in tabs
 Errors  if tabs does not have I * J columns.
 --------------------------------------------------------------------------------
 */
// [[Rcpp::export]]
List build_W_from_tables(IntegerMatrix tabs, int I, int J, IntegerVector delta) {
  if (tabs.ncol() != I * J) stop("each row of tabs must hold I*J cells (row-major).");
  std::map<std::vector<int>, double> W;
  std::vector<int> TD(J), m(J);
  for (int k = 0; k < tabs.nrow(); ++k) {
    double lc = 0.0;
    for (int j = 0; j < J; ++j) {
      TD[j] = 0; m[j] = 0;
      for (int i = 0; i < I; ++i) {
        int v = tabs(k, i * J + j);
        m[j] += v; if (delta[i] == 1) TD[j] += v;
        lc -= std::lgamma(v + 1.0);
      }
      lc += std::lgamma(TD[j] + 1.0) + std::lgamma(m[j] - TD[j] + 1.0);
    }
    auto it = W.find(TD);
    if (it == W.end()) W.emplace(TD, lc); else it->second = lse(it->second, lc);
    if ((k & 65535) == 0) Rcpp::checkUserInterrupt();
  }
  int R = W.size(); IntegerMatrix TDmat(R, J); NumericVector logW(R); int r = 0;
  for (auto& kv : W) { for (int j = 0; j < J; ++j) TDmat(r, j) = kv.first[j]; logW[r] = kv.second; ++r; }
  return List::create(_["TDmat"] = TDmat, _["logW"] = logW, _["n_tables"] = (double)tabs.nrow());
}

/*
 ================================================================================
 EXACT MOMENTS OF THE CELL COUNTS
 
 First and second moments of the cell counts under the sensitivity model, used
 by the normal approximation. The helpers Factorial, choose_C, log_choose and
 compute_log_sum_exp are defined first.
 ================================================================================
 */

#include <stdio.h>
#include <limits.h>
#include <float.h>




/*
 --------------------------------------------------------------------------------
 Factorial(nValue)
 
 n! by direct multiplication in double precision. Used only by choose_C.
 
 Input   nValue  integer n >= 2, passed as a double
 Output  n!
 Note    For n = 0 or n = 1 the result is 0, not 1; choose_C never calls
 Factorial with these values. The result overflows for n > 170.
 --------------------------------------------------------------------------------
 */

double Factorial(double nValue)
{
  double result = nValue;
  double result_next;
  double pc = nValue;
  do
  {
    result_next = result*(pc-1);
    result = result_next;
    pc--;
  }while(pc>2);
  nValue = result;
  return nValue;
}

/*
 --------------------------------------------------------------------------------
 choose_C(n, k)
 
 Binomial coefficient C(n, k) in double precision. Used by the moment functions
 zu_moments_two_treatment, zu_moments_three_treatment and
 zu_moments_four_treatment.
 
 The cases k < 0, n < 0, k > n, k = 0, k = n, k = 1 and k = n - 1 are handled
 directly; all other cases use Factorial.
 
 Input   n, k  nonnegative integers, passed as doubles
 Output  C(n, k); 0 if k < 0, n < 0 or k > n.
 Note    In the general case (2 <= k <= n - 2), Factorial overflows for
 n > 170 and the result is not finite.
 --------------------------------------------------------------------------------
 */

// [[Rcpp::export]]
double choose_C(double n, double k)
{

  if(n<-0.5){
    return(0.0);
  }

  if(k<-0.5){
    return(0.0);
  }

  if (fabs(n - k) < 1e-7 || fabs(k)  < 1e-7) {
    return 1.0;
  }


  if(k>n){
    return(0.0);
  }



  if( fabs(k-1.0) < 1e-7 || fabs(k - (n-1)) < 1e-7)
  {
    return n;
  }


  return Factorial(n) /(Factorial(k)*Factorial((n - k)));
}



/*
 --------------------------------------------------------------------------------
 log_choose(n, k)
 
 log C(n, k); -Inf if k < 0 or k > n. Same result as lchoose_. Used by
 zu_moments_two_treatment_log, zu_moments_three_treatment_log and
 zu_moments_four_treatment_log.
 
 Input   n, k  integers
 Output  log C(n, k).
 --------------------------------------------------------------------------------
 */

double log_choose(int n, int k) {
  if (k < 0 || k > n) return -INFINITY; // log(0) = -Inf
  return lgamma(n + 1.0) - lgamma(k + 1.0) - lgamma(n - k + 1.0);
}


/*
 --------------------------------------------------------------------------------
 compute_log_sum_exp(log_terms)
 
 log(sum_k exp(log_terms[k])), computed by factoring out the largest term.
 Terms equal to -Inf contribute 0. Used by the zu_moments_*_log functions to
 normalize the distribution of (Q_1, ..., Q_I) without overflow.
 
 Input   log_terms  vector of log-scale values
 Output  log of the sum; -Inf if log_terms is empty or every term is -Inf.
 --------------------------------------------------------------------------------
 */

double compute_log_sum_exp(const std::vector<double>& log_terms) {
  if (log_terms.empty()) return -INFINITY; // log(0) = -Inf
  double max_log = -INFINITY;
  for (const auto& lt : log_terms) {
    if (lt > max_log) max_log = lt;
  }
  if (max_log == -INFINITY) return -INFINITY; // All terms are -Inf
  double sum = 0.0;
  for (const auto& lt : log_terms) {
    sum += std::exp(lt - max_log);
  }
  return max_log + std::log(sum);
}






/*
 * This function calculates the denominator of any table with two treatments, thus, will
 * be the denominator of two by two table, generic bias
 * the shared_divisor command is meant for scaling the denominator and the numerator at the same time
 * We always call the lower bound of something with the notation lb_XX, and upper bound ub_XX
 * We follow the convention in Rosenbaum's that the treatment margins is n_i, outcome margins are m_j
 */





/*The following function assists the computation for normal approximation
 *
 */

/*
 * This function calculates the first two noncentral moment in a two treatment
 * setup, namely, E(Z_1^Tu) and E((Z_1^TU)^2)
 * The first returned value is the expectation, the second returned value is the
 * second non-central moment of Z1^Tu
 */





/*
 --------------------------------------------------------------------------------
 zu_moments_two_treatment_log(n1, n2, N, gamma_delta, Us, shared_divisor)
 
 First and second moments of Q_1, the number of subjects with u_s = 1 assigned
 to treatment 1, for a table with two treatment levels. At a binary u with
 u-bar = Us,
 
 P(Q_1 = q_1)  proportional to  exp(gamma_1 q_1) C(Us, q_1) C(N - Us, n1 - q_1),
 
 where gamma_1 = gamma_delta[1] - gamma_delta[2] and q_1 ranges over
 max(0, n1 + Us - N) <= q_1 <= min(Us, n1).
 Sums are computed on the log scale.
 
 Input   n1, n2          treatment margins N_1., N_2.; n2 is not used, since
 N_2. = N - N_1.
 N               total number of subjects
 gamma_delta     gamma * delta, length 2
 Us              u-bar, the number of subjects with u_s = 1
 shared_divisor  positive scaling constant applied to every term; it
 cancels and does not affect the result
 Output  c(E[Q_1], E[Q_1^2]).
 Errors  if every term of the distribution is zero.
 Note    If Us = 0, the result is (0, 0); if Us = N, it is (n1, n1^2). At
 gamma = 0, Q_1 is hypergeometric, with E[Q_1] = n1 Us / N.
 --------------------------------------------------------------------------------
 */

// [[Rcpp::export]]
NumericVector zu_moments_two_treatment_log(int n1, int n2, int N, NumericVector gamma_delta, int Us, double shared_divisor) {
  double gamma = gamma_delta[0] - gamma_delta[1];
  double z1u_first = 0.0;
  double z1u_second = 0.0;

  // Sanity check
  if (Us == 0) {   
    z1u_first = 0.0;
    z1u_second = 0.0;
    return NumericVector::create(z1u_first, z1u_second);
  }

  if (Us == N) {
    z1u_first = static_cast<double>(n1);
    z1u_second = std::pow(static_cast<double>(n1), 2);
    return NumericVector::create(z1u_first, z1u_second);
  }

  // Define bounds for q1
  int lb_q1 = std::max(0, n1 + Us - N);
  int ub_q1 = std::min(Us, n1);

  // Vectors to store log-weighted terms for convolution and each moment
  std::vector<double> log_weights;
  std::vector<double> log_z1u_terms;
  std::vector<double> log_z1u_sq_terms;

  // Iterate through all valid q1 to compute log-weighted terms
  for (int q1 = lb_q1; q1 <= ub_q1; ++q1) {
    // Compute the exponent term
    double exponent = gamma * q1;

    // Compute log(kernel)
    double log_kernel = log_choose(Us, q1) + log_choose(N - Us, n1 - q1) - std::log(shared_divisor);

    // Compute log(weighted_kernel) = exponent + log_kernel
    double log_weighted = exponent + log_kernel;

    // Store the log-weighted term for the convolution sum
    log_weights.push_back(log_weighted);

    // Store log(q1 * w) = log(q1) + log_weighted, handle q1=0
    if (q1 > 0) {
      log_z1u_terms.push_back(std::log(static_cast<double>(q1)) + log_weighted);
      log_z1u_sq_terms.push_back(std::log(static_cast<double>(q1) * static_cast<double>(q1)) + log_weighted);
    }
  }

  // Compute log-sum-exp for convolution sum and each moment
  double log_convolution_sum = compute_log_sum_exp(log_weights);
  double log_z1u_sum = compute_log_sum_exp(log_z1u_terms);
  double log_z1u_sq_sum = compute_log_sum_exp(log_z1u_sq_terms);

  // Check if convolution sum is zero (log_convolution_sum == -Inf)
  if (log_convolution_sum == -INFINITY) {
    Rcpp::stop("Convolution sum is zero; division by zero detected.");
  }

  // Compute the moments by exponentiating the differences
  // Handle cases where log_z1u_sum or log_z1u_sq_sum might be -Inf
  z1u_first = (log_z1u_sum == -INFINITY) ? 0.0 : std::exp(log_z1u_sum - log_convolution_sum);
  z1u_second = (log_z1u_sq_sum == -INFINITY) ? 0.0 : std::exp(log_z1u_sq_sum - log_convolution_sum);

  return NumericVector::create(z1u_first, z1u_second);
}



/*
 --------------------------------------------------------------------------------
 zu_moments_three_treatment_log(n1, n2, n3, N, gamma_delta, Us, shared_divisor)
 
 First and second moments of (Q_1, Q_2), where Q_i is the number of subjects
 with u_s = 1 assigned to treatment i, for a table with three treatment levels.
 At a binary u with u-bar = Us, and with Q_3 = Us - Q_1 - Q_2,
 
 P(Q_1 = q_1, Q_2 = q_2)  proportional to
 exp(gamma_1 q_1 + gamma_2 q_2)
 C(Us, q_1) C(Us - q_1, q_2) C(N - Us, n1 - q_1) C(N - Us - n1 + q_1, n2 - q_2),
 
 where gamma_i = gamma_delta[i] - gamma_delta[3]. Because delta is binary,
 gamma_1 and gamma_2 are equal whenever both are nonzero, which the code uses
 to write the exponent in terms of q = q_1 + q_2. Sums are computed on the log
 scale.
 
 Input   n1, n2, n3      treatment margins N_1., N_2., N_3.; n3 is not used,
 since N_3. = N - N_1. - N_2.
 N               total number of subjects
 gamma_delta     gamma * delta, length 3, with delta binary
 Us              u-bar, the number of subjects with u_s = 1
 shared_divisor  positive scaling constant applied to every term; it
 cancels and does not affect the result
 Output  c(E[Q_1], E[Q_2], E[Q_1^2], E[Q_2^2], E[Q_1 Q_2]).
 Errors  if every term of the distribution is zero.
 Note    If Us = 0, every moment is 0; if Us = N, the result is
 (n1, n2, n1^2, n2^2, n1 n2). At gamma = 0, (Q_1, Q_2, Q_3) is
 multivariate hypergeometric, with E[Q_i] = N_i. Us / N.
 --------------------------------------------------------------------------------
 */

// [[Rcpp::export]]
NumericVector zu_moments_three_treatment_log(int n1, int n2, int n3, int N, NumericVector gamma_delta, int Us, double shared_divisor) {
  double gamma1 = gamma_delta[0] - gamma_delta[2];
  double gamma2 = gamma_delta[1] - gamma_delta[2];
  double gamma = (gamma1 != 0.0) ? gamma1 : (gamma2 != 0.0 ? gamma2 : 0.0);

  // Sanity check
  if (Us == 0) {   
    return NumericVector::create(0.0, 0.0, 0.0, 0.0, 0.0);
  }

  if (Us == N) {
    return NumericVector::create(
      static_cast<double>(n1),
      static_cast<double>(n2),
      static_cast<double>(n1 * n1),
      static_cast<double>(n2 * n2),
      static_cast<double>(n1 * n2)
    );
  }

  int lb_q = std::max(0, n1 + n2 + Us - N);
  int ub_q = std::min(n1 + n2, Us);
  int lb_q1 = std::max(0, n1 + Us - N);
  int ub_q1 = std::min(Us, n1);
  int lb_q2 = std::max(0, n2 + Us - N);
  int ub_q2 = std::min(n2, Us);

  // Vectors to store log-weighted terms for convolution and each moment
  std::vector<double> log_weights;
  std::vector<double> log_z1u_terms;
  std::vector<double> log_z1u_sq_terms;
  std::vector<double> log_z2u_terms;
  std::vector<double> log_z2u_sq_terms;
  std::vector<double> log_z1u_z2u_terms;

  // Iterate through all valid q and q1 to compute log-weighted terms
  for (int q = lb_q; q <= ub_q; ++q) {
    for (int q1 = lb_q1; q1 <= ub_q1; ++q1) {
      int q2 = q - q1;
      if (q2 < lb_q2 || q2 > ub_q2) continue;

      // Compute the exponent term based on gamma1 and gamma2
      double exponent;
      if (gamma1 == 0.0) {
        exponent = gamma * (q - q1);
      } else if (gamma2 == 0.0) {
        exponent = gamma * (q - q2);
      } else {
        exponent = gamma * q;
      }

      // Compute log(exp_term * kernel) = exponent + log(kernel)
      double log_kernel = log_choose(Us, q1) + log_choose(Us - q1, q2) +
        log_choose(N - Us, n1 - q1) + log_choose((N - Us) - (n1 - q1), n2 - q2) -
        std::log(shared_divisor);

      double log_weighted = exponent + log_kernel;

      // Store the log-weighted term for the convolution sum
      log_weights.push_back(log_weighted);

      // Store log(q1 * w) = log(q1) + log_weighted, handle q1=0
      if (q1 > 0) {
        log_z1u_terms.push_back(std::log(static_cast<double>(q1)) + log_weighted);
        log_z1u_sq_terms.push_back(std::log(static_cast<double>(q1 * q1)) + log_weighted);
      }

      // Store log(q2 * w) = log(q2) + log_weighted, handle q2=0
      if (q2 > 0) {
        log_z2u_terms.push_back(std::log(static_cast<double>(q2)) + log_weighted);
        log_z2u_sq_terms.push_back(std::log(static_cast<double>(q2 * q2)) + log_weighted);
      }

      // Store log(q1 * q2 * w) = log(q1) + log(q2) + log_weighted, handle q1=0 or q2=0
      if (q1 > 0 && q2 > 0) {
        log_z1u_z2u_terms.push_back(std::log(static_cast<double>(q1)) + std::log(static_cast<double>(q2)) + log_weighted);
      }
    }
  }

  // Compute log-sum-exp for convolution sum and each moment
  double log_convolution_sum = compute_log_sum_exp(log_weights);
  double log_z1u_sum = compute_log_sum_exp(log_z1u_terms);
  double log_z1u_sq_sum = compute_log_sum_exp(log_z1u_sq_terms);
  double log_z2u_sum = compute_log_sum_exp(log_z2u_terms);
  double log_z2u_sq_sum = compute_log_sum_exp(log_z2u_sq_terms);
  double log_z1u_z2u_sum = compute_log_sum_exp(log_z1u_z2u_terms);

  // Check if convolution sum is zero (log_convolution_sum == -Inf)
  if (log_convolution_sum == -INFINITY) {
    Rcpp::stop("Convolution sum is zero; division by zero detected.");
  }

  // Compute the moments by exponentiating the differences
  double E_z1u = (log_z1u_sum == -INFINITY) ? 0.0 : std::exp(log_z1u_sum - log_convolution_sum);
  double E_z2u = (log_z2u_sum == -INFINITY) ? 0.0 : std::exp(log_z2u_sum - log_convolution_sum);
  double E_z1u_sq = (log_z1u_sq_sum == -INFINITY) ? 0.0 : std::exp(log_z1u_sq_sum - log_convolution_sum);
  double E_z2u_sq = (log_z2u_sq_sum == -INFINITY) ? 0.0 : std::exp(log_z2u_sq_sum - log_convolution_sum);
  double E_z1u_z2u = (log_z1u_z2u_sum == -INFINITY) ? 0.0 : std::exp(log_z1u_z2u_sum - log_convolution_sum);

  return NumericVector::create(E_z1u, E_z2u, E_z1u_sq, E_z2u_sq, E_z1u_z2u);
}





/*
 * This function obtains the moments for z1u, z2u, z3u, and their pairwise
 * products when there are four treatment levels, where z4u = Us - z1u - z2u - z3u.
 * The joint distribution is the same as in zu_moments_three_treatment_log with
 * one additional factor, and the accumulation is done in log space.
 *
 * The inputs are n1, n2, n3, n4 (the number of units in each treatment level),
 * N (the total number of units), gamma_delta (a vector of length four), Us (the
 * number of units with u = 1), and shared_divisor.
 *
 * The returned vector contains:
 * (E(z1u), E(z2u), E(z3u),
 *  E((z1u)^2), E((z2u)^2), E((z3u)^2),
 *  E(z1u*z2u), E(z1u*z3u), E(z2u*z3u))
 *
 * The output is the same as zu_moments_four_treatment,
 * but this version is more numerically stable.
 */

// [[Rcpp::export]]
NumericVector zu_moments_four_treatment_log(int n1, int n2, int n3, int n4, int N, NumericVector gamma_delta, int Us, double shared_divisor) {

  // Calculate gamma values
  double gamma1 = gamma_delta[0] - gamma_delta[3];
  double gamma2 = gamma_delta[1] - gamma_delta[3];
  double gamma3 = gamma_delta[2] - gamma_delta[3];

  
  double gamma = (gamma1 != 0.0) ? gamma1 : ((gamma2 != 0.0) ? gamma2 : gamma3);

  // Initialize cumulants for first and second moments, and product moments
  // First moments
  std::vector<double> log_z1u_terms;
  std::vector<double> log_z2u_terms;
  std::vector<double> log_z3u_terms;

  // Second moments
  std::vector<double> log_z1u_sq_terms;
  std::vector<double> log_z2u_sq_terms;
  std::vector<double> log_z3u_sq_terms;

  // Product moments
  std::vector<double> log_z1u_z2u_terms;
  std::vector<double> log_z1u_z3u_terms;
  std::vector<double> log_z2u_z3u_terms;

  // Vector to store log-weighted terms for convolution sum
  std::vector<double> log_weights;

  // Sanity check
  if (Us == 0) {  
    return NumericVector::create(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
  }

  if (Us == N) {
    // When Us == N, all units are treated; moments are deterministic
    double z1u_first = static_cast<double>(n1);
    double z2u_first = static_cast<double>(n2);
    double z3u_first = static_cast<double>(n3);

    double z1u_second = std::pow(static_cast<double>(n1), 2);
    double z2u_second = std::pow(static_cast<double>(n2), 2);
    double z3u_second = std::pow(static_cast<double>(n3), 2);

    double z1u_z2u = static_cast<double>(n1) * static_cast<double>(n2);
    double z1u_z3u = static_cast<double>(n1) * static_cast<double>(n3);
    double z2u_z3u = static_cast<double>(n2) * static_cast<double>(n3);

    return NumericVector::create(z1u_first, z2u_first, z3u_first,
                                 z1u_second, z2u_second, z3u_second,
                                 z1u_z2u, z1u_z3u, z2u_z3u);
  }

  // Define bounds for q1, q2, q3
  int lb_q1 = std::max(0, n1 + Us - N);
  int ub_q1 = std::min(Us, n1);

  for (int q1 = lb_q1; q1 <= ub_q1; ++q1) {
    // After selecting q1, define bounds for q2
    int lb_q2 = std::max(0, n2 + Us - N - q1);
    int ub_q2 = std::min(Us - q1, n2);

    for (int q2 = lb_q2; q2 <= ub_q2; ++q2) {
      // After selecting q1 and q2, define bounds for q3
      int lb_q3 = std::max(0, n3 + Us - N - q1 - q2);
      int ub_q3 = std::min(Us - q1 - q2, n3);

      for (int q3 = lb_q3; q3 <= ub_q3; ++q3) {
        // Compute q4 based on q1, q2, q3
        int q4 = Us - q1 - q2 - q3;

        // Ensure q4 is within valid bounds
        if (q4 < 0 || q4 > n4) continue;

        // Compute the exponent term: gamma1*q1 + gamma2*q2 + gamma3*q3
        double exponent = gamma1 * q1 + gamma2 * q2 + gamma3 * q3;

        // Compute log(kernel)
        double log_kernel = log_choose(Us, q1) +
          log_choose(Us - q1, q2) +
          log_choose(Us - q1 - q2, q3) +
          log_choose(N - Us, n1 - q1) +
          log_choose(N - Us - (n1 - q1), n2 - q2) +
          log_choose(N - Us - (n1 - q1) - (n2 - q2), n3 - q3) -
          std::log(shared_divisor);

        // Compute log(weighted_kernel) = exponent + log_kernel
        double log_weighted = exponent + log_kernel;

        // Store the log-weighted term for the convolution sum
        log_weights.push_back(log_weighted);

        // Store log terms for first moments
        if (q1 > 0) {
          log_z1u_terms.push_back(std::log(static_cast<double>(q1)) + log_weighted);
        }
        if (q2 > 0) {
          log_z2u_terms.push_back(std::log(static_cast<double>(q2)) + log_weighted);
        }
        if (q3 > 0) {
          log_z3u_terms.push_back(std::log(static_cast<double>(q3)) + log_weighted);
        }

        // Store log terms for second moments
        if (q1 > 0) {
          log_z1u_sq_terms.push_back(std::log(static_cast<double>(q1) * static_cast<double>(q1)) + log_weighted);
        }
        if (q2 > 0) {
          log_z2u_sq_terms.push_back(std::log(static_cast<double>(q2) * static_cast<double>(q2)) + log_weighted);
        }
        if (q3 > 0) {
          log_z3u_sq_terms.push_back(std::log(static_cast<double>(q3) * static_cast<double>(q3)) + log_weighted);
        }

        // Store log terms for product moments
        if (q1 > 0 && q2 > 0) {
          log_z1u_z2u_terms.push_back(std::log(static_cast<double>(q1)) + std::log(static_cast<double>(q2)) + log_weighted);
        }
        if (q1 > 0 && q3 > 0) {
          log_z1u_z3u_terms.push_back(std::log(static_cast<double>(q1)) + std::log(static_cast<double>(q3)) + log_weighted);
        }
        if (q2 > 0 && q3 > 0) {
          log_z2u_z3u_terms.push_back(std::log(static_cast<double>(q2)) + std::log(static_cast<double>(q3)) + log_weighted);
        }
      }
    }
  }

  // Compute log-sum-exp for convolution sum and each moment
  double log_convolution_sum = compute_log_sum_exp(log_weights);
  double log_z1u_sum = compute_log_sum_exp(log_z1u_terms);
  double log_z2u_sum = compute_log_sum_exp(log_z2u_terms);
  double log_z3u_sum = compute_log_sum_exp(log_z3u_terms);

  double log_z1u_sq_sum = compute_log_sum_exp(log_z1u_sq_terms);
  double log_z2u_sq_sum = compute_log_sum_exp(log_z2u_sq_terms);
  double log_z3u_sq_sum = compute_log_sum_exp(log_z3u_sq_terms);

  double log_z1u_z2u_sum = compute_log_sum_exp(log_z1u_z2u_terms);
  double log_z1u_z3u_sum = compute_log_sum_exp(log_z1u_z3u_terms);
  double log_z2u_z3u_sum = compute_log_sum_exp(log_z2u_z3u_terms);

  // Check if convolution sum is zero (log_convolution_sum == -Inf)
  if (log_convolution_sum == -INFINITY) {
    Rcpp::stop("Convolution sum is zero; division by zero detected.");
  }

  // Compute the moments by exponentiating the differences
  // Handle cases where log_sums might be -Inf
  double E_z1u = (log_z1u_sum == -INFINITY) ? 0.0 : std::exp(log_z1u_sum - log_convolution_sum);
  double E_z2u = (log_z2u_sum == -INFINITY) ? 0.0 : std::exp(log_z2u_sum - log_convolution_sum);
  double E_z3u = (log_z3u_sum == -INFINITY) ? 0.0 : std::exp(log_z3u_sum - log_convolution_sum);

  double E_z1u_sq = (log_z1u_sq_sum == -INFINITY) ? 0.0 : std::exp(log_z1u_sq_sum - log_convolution_sum);
  double E_z2u_sq = (log_z2u_sq_sum == -INFINITY) ? 0.0 : std::exp(log_z2u_sq_sum - log_convolution_sum);
  double E_z3u_sq = (log_z3u_sq_sum == -INFINITY) ? 0.0 : std::exp(log_z3u_sq_sum - log_convolution_sum);

  double E_z1u_z2u = (log_z1u_z2u_sum == -INFINITY) ? 0.0 : std::exp(log_z1u_z2u_sum - log_convolution_sum);
  double E_z1u_z3u = (log_z1u_z3u_sum == -INFINITY) ? 0.0 : std::exp(log_z1u_z3u_sum - log_convolution_sum);
  double E_z2u_z3u = (log_z2u_z3u_sum == -INFINITY) ? 0.0 : std::exp(log_z2u_z3u_sum - log_convolution_sum);

  return NumericVector::create(E_z1u, E_z2u, E_z3u,
                               E_z1u_sq, E_z2u_sq, E_z3u_sq,
                               E_z1u_z2u, E_z1u_z3u, E_z2u_z3u);
}




/*
 * This function computes first two central moments of z1^TR_1 for two by two table
 */


// =============================================================================
// FUNCTION  zr_moments_two_by_two                            [exported to R]
// PURPOSE   Exact mean and variance of the free cell N_11 of a 2 x 2 table under the
//           generic bias model at u-bar = (u1, u2), from the moments of Q_1:
//             E[N_11]   = rho_{1,1} E[Q_1] + rho_{1,0} (N_1. - E[Q_1]),
//             Var[N_11] = (w_{1,1} - w_{1,0}) E[Q_1]
//                         - (w_{1,1}/u-bar + w_{1,0}/(N - u-bar)) E[Q_1^2]
//                         + N_1. (N - u-bar - N_1. + 2 E[Q_1]) w_{1,0} / (N - u-bar)
//                         + Var[Q_1] (rho_{1,1} - rho_{1,0})^2,
//           with rho_{j,1} = u-bar_j / u-bar, rho_{j,0} = (N_.j - u-bar_j)/(N - u-bar)
//           (q11, q10 in the code) and w_{j,k} the within-group variances of the
//           outcome indicators (w11, w10).  When u-bar is 0 or N the hypergeometric
//           moments are returned.
// PAPER     Theorem D.1 (mean and variance of N_ij), Appendix D.
// INPUT     n1, n2 = N_i.;  m1, m2 = N_.j;  N;  gamma_delta (length 2);  u1, u2;
//           shared_divisor.
// OUTPUT    c(E[N_11], Var[N_11]).
// EXAMPLE   (R)
//   zr_moments_two_by_two(20L, 30L, 15L, 35L, 50L, gamma_delta = c(0, 1), u1 = 0L, u2 = 35L, shared_divisor = 1)
//   ## with u1 = u2 = 0 (no confounding) the result is c(15*20/50, 20*30*15*35/50^2/49)
// =============================================================================
// [[Rcpp::export]]
NumericVector zr_moments_two_by_two(int n1, int n2, int m1, int m2, int N, NumericVector gamma_delta, int u1, int u2, double shared_divisor) {
  NumericVector z1u_moments = zu_moments_two_treatment_log(n1,n2,N,gamma_delta,u1+u2,shared_divisor);
  double z1u_first= z1u_moments[0];
  double z1u_second= z1u_moments[1];
  double Us = u1 + u2;
  //if Us = 0 or N, becomes the usual RCT case
  if(Us==0 || Us==N)
    return(NumericVector::create((double) m1*n1/N,(double) n1*n2*m1*m2/N/N/(N-1)));

  double q11 = u1/Us;
  double q10 = (m1-u1)/(N-Us);
  double z1r1_first = q11*z1u_first+q10*(n1-z1u_first);
  double w11, w10;
  if(Us==1)
    w11 = 0;
  else
    w11  = (double) u1/(Us-1)*(1-2*q11) + (double) Us/(Us-1)*q11*q11;

  if(Us==N-1)
    w10 = 0;
  else
    w10   = (double)(m1 - u1)/(N-Us-1)*(1-2*q10) + (double) (N-Us)/(N-Us-1)*pow(q10,2);
  double z1r1_var = (w11-w10)*z1u_first-(z1u_second)*(w11/Us+w10/(N-Us))
    +n1*(N-Us-n1+2*z1u_first)*w10/(N-Us)+(z1u_second-z1u_first*z1u_first)*(q11-q10)*(q11-q10);

  return(NumericVector::create(z1r1_first,z1r1_var));

}



// [[Rcpp::export]]
NumericVector zr_moments_two_by_three(int n1, int n2, int m1, int m2, int m3, int N, NumericVector gamma_delta, int u1, int u2, int u3, double shared_divisor) {
  NumericVector z1u_moments = zu_moments_two_treatment_log(n1,n2,N,gamma_delta,u1+u2+u3,shared_divisor);
  double z1u_first= z1u_moments[0];
  double z1u_second= z1u_moments[1];
  double Us = u1 + u2+u3;
  //if Us = 0 or N, becomes the usual RCT case
  if(Us==0 || Us==N)
    return(NumericVector::create((double) m1*n1/N,(double) m2*n1/N, (double) n1*(N-n1)*m1*(N-m1)/N/N/(N-1),(double) n1*(N-n1)*m2*(N-m2)/N/N/(N-1),(double) m1*m2*n1*(n1-N)/N/N/(N-1)));

  double q11 = u1/Us;
  double q10 = (m1-u1)/(N-Us);
  double q21 = u2/Us;
  double q20 = (m2-u2)/(N-Us);
  double z1r1_first = q11*z1u_first+q10*(n1-z1u_first);
  double z1r2_first = q21*z1u_first+q20*(n1-z1u_first);
  double w11, w10, w21, w20;
  if(Us==1)
  { w11 = 0;
    w21 = 0; }
  else
  {w11  = (double) u1/(Us-1)*(1-2*q11) + (double) Us/(Us-1)*q11*q11;
    w21  = (double) u2/(Us-1)*(1-2*q21) + (double) Us/(Us-1)*q21*q21;}
  if(Us==N-1)
  {w10 = 0;
    w20 = 0;
  }
  else
  {
    w10   = (double)(m1 - u1)/(N-Us-1)*(1-2*q10) + (double) (N-Us)/(N-Us-1)*pow(q10,2);
    w20   = (double)(m2 - u2)/(N-Us-1)*(1-2*q20) + (double) (N-Us)/(N-Us-1)*pow(q20,2);

  }
  double z1r1_var = (w11-w10)*z1u_first-(z1u_second)*(w11/Us+w10/(N-Us))
    +n1*(N-Us-n1+2*z1u_first)*w10/(N-Us)+(z1u_second-z1u_first*z1u_first)*(q11-q10)*(q11-q10);
  double z1r2_var = (w21-w20)*z1u_first-(z1u_second)*(w21/Us+w20/(N-Us))
    +n1*(N-Us-n1+2*z1u_first)*w20/(N-Us)+(z1u_second-z1u_first*z1u_first)*(q21-q20)*(q21-q20);
  double z1r1_z1r2_cov;
  if(Us>1 && Us <N-1){
    z1r1_z1r2_cov = q11*q21/(Us-1)*(z1u_second-Us*z1u_first)+q10*q20/(N-Us-1)*n1*(n1-N+Us)
    +q10*q20/(N-Us-1)*(z1u_second-(2*n1-N+Us)*z1u_first)
    +(q11*q21+q10*q20-q10*q21-q20*q11)*(z1u_second-z1u_first*z1u_first);
  }
  if(Us==1){
    z1r1_z1r2_cov = q10*q20/(N-Us-1)*n1*(n1-N+Us)
    +q10*q20/(N-Us-1)*(z1u_second-(2*n1-N+Us)*z1u_first)
    +(q11*q21+q10*q20-q10*q21-q20*q11)*(z1u_second-z1u_first*z1u_first);

  }
  if(Us==N-1){
    z1r1_z1r2_cov = q11*q21/(Us-1)*(z1u_second-Us*z1u_first)
    +(q11*q21+q10*q20-q10*q21-q20*q11)*(z1u_second-z1u_first*z1u_first);
  }


  return(NumericVector::create(z1r1_first,z1r2_first, z1r1_var, z1r2_var, z1r1_z1r2_cov));

}




/*
 FUNCTION  zr_moments_two_by_four
 PURPOSE   Exact means, variances and pairwise covariances of the free cells
 N_11, N_12, N_13 of a 2 x 4 table at u-bar = (u1, ..., u4), from the
 moments of Q_1.  Same formulas as zr_moments_two_by_three, with one
 more outcome level.
 INPUT     n1, n2 = N_i.;  m1..m4 = N_.j;  N;  gamma_delta (length 2);  u1..u4;
 shared_divisor.
 OUTPUT    c(E[N_11], E[N_12], E[N_13], Var[N_11], Var[N_12], Var[N_13],
 Cov[N_11,N_12], Cov[N_11,N_13], Cov[N_12,N_13]).
 EXAMPLE   (R)
 zr_moments_two_by_four(30L, 30L, 15L, 15L, 15L, 15L, 60L, c(0, 1), 0L, 0L, 15L, 15L, 1)
 */

// [[Rcpp::export]]
NumericVector zr_moments_two_by_four(int n1, int n2, int m1, int m2, int m3, int m4, int N, NumericVector gamma_delta, int u1, int u2, int u3, int u4, double shared_divisor) {
  NumericVector z1u_moments = zu_moments_two_treatment_log(n1,n2,N,gamma_delta,u1+u2+u3+u4,shared_divisor);
  double z1u_first= z1u_moments[0];
  double z1u_second= z1u_moments[1];
  double Us = u1 + u2+u3+u4;
  //if Us = 0 or N, becomes the usual RCT case
  if(Us==0 || Us==N)
    return(NumericVector::create((double) m1*n1/N,(double) m2*n1/N, (double) m3*n1/N, (double) n1*(N-n1)*m1*(N-m1)/N/N/(N-1), (double) n1*(N-n1)*m2*(N-m2)/N/N/(N-1),
                                 (double) n1*(N-n1)*m3*(N-m3)/N/N/(N-1), (double) m1*m2*n1*(n1-N)/N/N/(N-1),(double) m1*m3*n1*(n1-N)/N/N/(N-1),(double) m2*m3*n1*(n1-N)/N/N/(N-1)));

  double q11 = u1/Us;
  double q10 = (m1-u1)/(N-Us);
  double q21 = u2/Us;
  double q20 = (m2-u2)/(N-Us);
  double q31 = u3/Us;
  double q30 = (m3-u3)/(N-Us);
  double z1r1_first = q11*z1u_first+q10*(n1-z1u_first);
  double z1r2_first = q21*z1u_first+q20*(n1-z1u_first);
  double z1r3_first = q31*z1u_first+q30*(n1-z1u_first);
  double w11, w10, w21, w20, w31, w30;
  if(Us==1)
  { w11 = 0;
    w21 = 0;
    w31 = 0; }
  else
  {w11  = (double) u1/(Us-1)*(1-2*q11) + (double) Us/(Us-1)*q11*q11;
    w21  = (double) u2/(Us-1)*(1-2*q21) + (double) Us/(Us-1)*q21*q21;
    w31  = (double) u3/(Us-1)*(1-2*q31) + (double) Us/(Us-1)*q31*q31;}
  if(Us==N-1)
  {w10 = 0;
    w20 = 0;
    w30 = 0;
  }
  else
  {
    w10   = (double)(m1 - u1)/(N-Us-1)*(1-2*q10) + (double) (N-Us)/(N-Us-1)*pow(q10,2);
    w20   = (double)(m2 - u2)/(N-Us-1)*(1-2*q20) + (double) (N-Us)/(N-Us-1)*pow(q20,2);
    w30   = (double)(m3 - u3)/(N-Us-1)*(1-2*q30) + (double) (N-Us)/(N-Us-1)*pow(q30,2);


  }
  double z1r1_var = (w11-w10)*z1u_first-(z1u_second)*(w11/Us+w10/(N-Us))
    +n1*(N-Us-n1+2*z1u_first)*w10/(N-Us)+(z1u_second-z1u_first*z1u_first)*(q11-q10)*(q11-q10);
  double z1r2_var = (w21-w20)*z1u_first-(z1u_second)*(w21/Us+w20/(N-Us))
    +n1*(N-Us-n1+2*z1u_first)*w20/(N-Us)+(z1u_second-z1u_first*z1u_first)*(q21-q20)*(q21-q20);
  double z1r3_var = (w31-w30)*z1u_first-(z1u_second)*(w31/Us+w30/(N-Us))
    +n1*(N-Us-n1+2*z1u_first)*w30/(N-Us)+(z1u_second-z1u_first*z1u_first)*(q31-q30)*(q31-q30);
  double z1r1_z1r2_cov;
  double z1r1_z1r3_cov;
  double z1r2_z1r3_cov;
  if(Us>1 && Us <N-1){
    z1r1_z1r2_cov = q11*q21/(Us-1)*(z1u_second-Us*z1u_first)+q10*q20/(N-Us-1)*n1*(n1-N+Us)
    +q10*q20/(N-Us-1)*(z1u_second-(2*n1-N+Us)*z1u_first)
    +(q11*q21+q10*q20-q10*q21-q20*q11)*(z1u_second-z1u_first*z1u_first);

    z1r1_z1r3_cov = q11*q31/(Us-1)*(z1u_second-Us*z1u_first)+q10*q30/(N-Us-1)*n1*(n1-N+Us)
      +q10*q30/(N-Us-1)*(z1u_second-(2*n1-N+Us)*z1u_first)
      +(q11*q31+q10*q30-q10*q31-q30*q11)*(z1u_second-z1u_first*z1u_first);

      z1r2_z1r3_cov = q21*q31/(Us-1)*(z1u_second-Us*z1u_first)+q20*q30/(N-Us-1)*n1*(n1-N+Us)
        +q20*q30/(N-Us-1)*(z1u_second-(2*n1-N+Us)*z1u_first)
        +(q21*q31+q20*q30-q20*q31-q30*q21)*(z1u_second-z1u_first*z1u_first);
  }
  if(Us==1){
    z1r1_z1r2_cov = q10*q20/(N-Us-1)*n1*(n1-N+Us)
    +q10*q20/(N-Us-1)*(z1u_second-(2*n1-N+Us)*z1u_first)
    +(q11*q21+q10*q20-q10*q21-q20*q11)*(z1u_second-z1u_first*z1u_first);

    z1r1_z1r3_cov = q10*q30/(N-Us-1)*n1*(n1-N+Us)
      +q10*q30/(N-Us-1)*(z1u_second-(2*n1-N+Us)*z1u_first)
      +(q11*q31+q10*q30-q10*q31-q30*q11)*(z1u_second-z1u_first*z1u_first);

      z1r2_z1r3_cov = q20*q30/(N-Us-1)*n1*(n1-N+Us)
        +q20*q30/(N-Us-1)*(z1u_second-(2*n1-N+Us)*z1u_first)
        +(q21*q31+q20*q30-q20*q31-q30*q21)*(z1u_second-z1u_first*z1u_first);

  }
  if(Us==N-1){
    z1r1_z1r2_cov = q11*q21/(Us-1)*(z1u_second-Us*z1u_first)
    +(q11*q21+q10*q20-q10*q21-q20*q11)*(z1u_second-z1u_first*z1u_first);

    z1r1_z1r3_cov = q11*q31/(Us-1)*(z1u_second-Us*z1u_first)
      +(q11*q31+q10*q30-q10*q31-q30*q11)*(z1u_second-z1u_first*z1u_first);

    z1r2_z1r3_cov = q21*q31/(Us-1)*(z1u_second-Us*z1u_first)
      +(q21*q31+q20*q30-q20*q31-q30*q21)*(z1u_second-z1u_first*z1u_first);
  }


  return(NumericVector::create(z1r1_first,z1r2_first, z1r3_first, z1r1_var, z1r2_var, z1r3_var, z1r1_z1r2_cov, z1r1_z1r3_cov, z1r2_z1r3_cov));

}




/*
 FUNCTION  zr_moments_two_by_five
 PURPOSE   Exact means, variances and pairwise covariances of the free cells
 N_11, ..., N_14 of a 2 x 5 table at u-bar = (u1, ..., u5), from the
 moments of Q_1.  Same formulas as zr_moments_two_by_three, with two
 more outcome levels.
 INPUT     n1, n2 = N_i.;  m1..m5 = N_.j;  N;  gamma_delta (length 2);  u1..u5;
 shared_divisor.
 OUTPUT    length 14: c(E[N_11], ..., E[N_14], Var[N_11], ..., Var[N_14],
 Cov[N_11,N_12], Cov[N_11,N_13], Cov[N_11,N_14],
 Cov[N_12,N_13], Cov[N_12,N_14], Cov[N_13,N_14]).
 EXAMPLE   (R)
 zr_moments_two_by_five(40L, 40L, 16L, 16L, 16L, 16L, 16L, 80L, c(0, 1), 0L, 0L, 0L, 16L, 16L, 1)
 */

// [[Rcpp::export]]
NumericVector zr_moments_two_by_five(int n1, int n2, int m1, int m2, int m3, int m4, int m5, int N, NumericVector gamma_delta, int u1, int u2, int u3, int u4,int u5, double shared_divisor) {
  NumericVector z1u_moments = zu_moments_two_treatment_log(n1,n2,N,gamma_delta,u1+u2+u3+u4+u5,shared_divisor);
  double z1u_first= z1u_moments[0];
  double z1u_second= z1u_moments[1];
  double Us = u1 + u2+u3+u4+u5;
  //if Us = 0 or N, becomes the usual RCT case
  if(Us==0 || Us==N)
    return(NumericVector::create((double) m1*n1/N,(double) m2*n1/N, (double) m3*n1/N,(double) m4*n1/N,(double) n1*(N-n1)*m1*(N-m1)/N/N/(N-1),(double) n1*(N-n1)*m2*(N-m2)/N/N/(N-1), (double) n1*(N-n1)*m3*(N-m3)/N/N/(N-1), (double) n1*(N-n1)*m4*(N-m4)/N/N/(N-1),
                                 (double) m1*m2*n1*(n1-N)/N/N/(N-1),(double) m1*m3*n1*(n1-N)/N/N/(N-1),(double) m1*m4*n1*(n1-N)/N/N/(N-1), (double) m2*m3*n1*(n1-N)/N/N/(N-1), (double) m2*m4*n1*(n1-N)/N/N/(N-1), (double) m3*m4*n1*(n1-N)/N/N/(N-1)));

  double q11 = u1/Us;
  double q10 = (m1-u1)/(N-Us);
  double q21 = u2/Us;
  double q20 = (m2-u2)/(N-Us);
  double q31 = u3/Us;
  double q30 = (m3-u3)/(N-Us);
  double q41 = u4/Us;
  double q40 = (m4-u4)/(N-Us);
  double z1r1_first = q11*z1u_first+q10*(n1-z1u_first);
  double z1r2_first = q21*z1u_first+q20*(n1-z1u_first);
  double z1r3_first = q31*z1u_first+q30*(n1-z1u_first);
  double z1r4_first = q41*z1u_first+q40*(n1-z1u_first);
  double w11, w10, w21, w20, w31, w30, w41, w40;
  if(Us==1)
  { w11 = 0;
    w21 = 0;
    w31 = 0;
    w41 = 0; }
  else
  {w11  = (double) u1/(Us-1)*(1-2*q11) + (double) Us/(Us-1)*q11*q11;
    w21  = (double) u2/(Us-1)*(1-2*q21) + (double) Us/(Us-1)*q21*q21;
    w31  = (double) u3/(Us-1)*(1-2*q31) + (double) Us/(Us-1)*q31*q31;
    w41  = (double) u4/(Us-1)*(1-2*q41) + (double) Us/(Us-1)*q41*q41;

  }
  if(Us==N-1)
  {w10 = 0;
    w20 = 0;
    w30 = 0;
    w40 = 0;
  }
  else
  {
    w10   = (double)(m1 - u1)/(N-Us-1)*(1-2*q10) + (double) (N-Us)/(N-Us-1)*pow(q10,2);
    w20   = (double)(m2 - u2)/(N-Us-1)*(1-2*q20) + (double) (N-Us)/(N-Us-1)*pow(q20,2);
    w30   = (double)(m3 - u3)/(N-Us-1)*(1-2*q30) + (double) (N-Us)/(N-Us-1)*pow(q30,2);
    w40   = (double)(m4 - u4)/(N-Us-1)*(1-2*q40) + (double) (N-Us)/(N-Us-1)*pow(q40,2);


  }
  double z1r1_var = (w11-w10)*z1u_first-(z1u_second)*(w11/Us+w10/(N-Us))
    +n1*(N-Us-n1+2*z1u_first)*w10/(N-Us)+(z1u_second-z1u_first*z1u_first)*(q11-q10)*(q11-q10);
  double z1r2_var = (w21-w20)*z1u_first-(z1u_second)*(w21/Us+w20/(N-Us))
    +n1*(N-Us-n1+2*z1u_first)*w20/(N-Us)+(z1u_second-z1u_first*z1u_first)*(q21-q20)*(q21-q20);
  double z1r3_var = (w31-w30)*z1u_first-(z1u_second)*(w31/Us+w30/(N-Us))
    +n1*(N-Us-n1+2*z1u_first)*w30/(N-Us)+(z1u_second-z1u_first*z1u_first)*(q31-q30)*(q31-q30);

  double z1r4_var = (w41-w40)*z1u_first-(z1u_second)*(w41/Us+w40/(N-Us))
    +n1*(N-Us-n1+2*z1u_first)*w40/(N-Us)+(z1u_second-z1u_first*z1u_first)*(q41-q40)*(q41-q40);

  double z1r1_z1r2_cov;
  double z1r1_z1r3_cov;
  double z1r1_z1r4_cov;
  double z1r2_z1r3_cov;
  double z1r2_z1r4_cov;
  double z1r3_z1r4_cov;

  if(Us>1 && Us <N-1){
    z1r1_z1r2_cov = q11*q21/(Us-1)*(z1u_second-Us*z1u_first)+q10*q20/(N-Us-1)*n1*(n1-N+Us)
    +q10*q20/(N-Us-1)*(z1u_second-(2*n1-N+Us)*z1u_first)
    +(q11*q21+q10*q20-q10*q21-q20*q11)*(z1u_second-z1u_first*z1u_first);

    z1r1_z1r3_cov = q11*q31/(Us-1)*(z1u_second-Us*z1u_first)+q10*q30/(N-Us-1)*n1*(n1-N+Us)
      +q10*q30/(N-Us-1)*(z1u_second-(2*n1-N+Us)*z1u_first)
      +(q11*q31+q10*q30-q10*q31-q30*q11)*(z1u_second-z1u_first*z1u_first);

      z1r1_z1r4_cov = q11*q41/(Us-1)*(z1u_second-Us*z1u_first)+q10*q40/(N-Us-1)*n1*(n1-N+Us)
        +q10*q40/(N-Us-1)*(z1u_second-(2*n1-N+Us)*z1u_first)
        +(q11*q41+q10*q40-q10*q41-q40*q11)*(z1u_second-z1u_first*z1u_first);

        z1r2_z1r3_cov = q21*q31/(Us-1)*(z1u_second-Us*z1u_first)+q20*q30/(N-Us-1)*n1*(n1-N+Us)
          +q20*q30/(N-Us-1)*(z1u_second-(2*n1-N+Us)*z1u_first)
          +(q21*q31+q20*q30-q20*q31-q30*q21)*(z1u_second-z1u_first*z1u_first);

          z1r2_z1r4_cov = q21*q41/(Us-1)*(z1u_second-Us*z1u_first)+q20*q40/(N-Us-1)*n1*(n1-N+Us)
            +q20*q40/(N-Us-1)*(z1u_second-(2*n1-N+Us)*z1u_first)
            +(q21*q41+q20*q40-q20*q41-q40*q21)*(z1u_second-z1u_first*z1u_first);

            z1r3_z1r4_cov = q31*q41/(Us-1)*(z1u_second-Us*z1u_first)+q30*q40/(N-Us-1)*n1*(n1-N+Us)
              +q30*q40/(N-Us-1)*(z1u_second-(2*n1-N+Us)*z1u_first)
              +(q31*q41+q30*q40-q30*q41-q40*q31)*(z1u_second-z1u_first*z1u_first);
  }
  if(Us==1){
    z1r1_z1r2_cov = q10*q20/(N-Us-1)*n1*(n1-N+Us)
    +q10*q20/(N-Us-1)*(z1u_second-(2*n1-N+Us)*z1u_first)
    +(q11*q21+q10*q20-q10*q21-q20*q11)*(z1u_second-z1u_first*z1u_first);

    z1r1_z1r3_cov = q10*q30/(N-Us-1)*n1*(n1-N+Us)
      +q10*q30/(N-Us-1)*(z1u_second-(2*n1-N+Us)*z1u_first)
      +(q11*q31+q10*q30-q10*q31-q30*q11)*(z1u_second-z1u_first*z1u_first);

      z1r1_z1r4_cov = q10*q40/(N-Us-1)*n1*(n1-N+Us)
        +q10*q40/(N-Us-1)*(z1u_second-(2*n1-N+Us)*z1u_first)
        +(q11*q41+q10*q40-q10*q41-q40*q11)*(z1u_second-z1u_first*z1u_first);

        z1r2_z1r3_cov = q20*q30/(N-Us-1)*n1*(n1-N+Us)
          +q20*q30/(N-Us-1)*(z1u_second-(2*n1-N+Us)*z1u_first)
          +(q21*q31+q20*q30-q20*q31-q30*q21)*(z1u_second-z1u_first*z1u_first);

          z1r2_z1r4_cov = q20*q40/(N-Us-1)*n1*(n1-N+Us)
            +q20*q40/(N-Us-1)*(z1u_second-(2*n1-N+Us)*z1u_first)
            +(q21*q41+q20*q40-q20*q41-q40*q21)*(z1u_second-z1u_first*z1u_first);

            z1r3_z1r4_cov = q30*q40/(N-Us-1)*n1*(n1-N+Us)
              +q30*q40/(N-Us-1)*(z1u_second-(2*n1-N+Us)*z1u_first)
              +(q31*q41+q30*q40-q30*q41-q40*q31)*(z1u_second-z1u_first*z1u_first);

  }
  if(Us==N-1){
    z1r1_z1r2_cov = q11*q21/(Us-1)*(z1u_second-Us*z1u_first)
    +(q11*q21+q10*q20-q10*q21-q20*q11)*(z1u_second-z1u_first*z1u_first);

    z1r1_z1r3_cov = q11*q31/(Us-1)*(z1u_second-Us*z1u_first)
      +(q11*q31+q10*q30-q10*q31-q30*q11)*(z1u_second-z1u_first*z1u_first);

    z1r1_z1r4_cov = q11*q41/(Us-1)*(z1u_second-Us*z1u_first)
      +(q11*q41+q10*q40-q10*q41-q40*q11)*(z1u_second-z1u_first*z1u_first);

    z1r2_z1r3_cov = q21*q31/(Us-1)*(z1u_second-Us*z1u_first)
      +(q21*q31+q20*q30-q20*q31-q30*q21)*(z1u_second-z1u_first*z1u_first);

    z1r2_z1r4_cov = q21*q41/(Us-1)*(z1u_second-Us*z1u_first)
      +(q21*q41+q20*q40-q20*q41-q40*q21)*(z1u_second-z1u_first*z1u_first);

    z1r3_z1r4_cov = q31*q41/(Us-1)*(z1u_second-Us*z1u_first)
      +(q31*q41+q30*q40-q30*q41-q40*q31)*(z1u_second-z1u_first*z1u_first);
  }


  return(NumericVector::create(z1r1_first,z1r2_first, z1r3_first, z1r4_first, z1r1_var, z1r2_var, z1r3_var, z1r4_var, z1r1_z1r2_cov, z1r1_z1r3_cov, z1r1_z1r4_cov, z1r2_z1r3_cov, z1r2_z1r4_cov, z1r3_z1r4_cov));

}



/*
 * This function computes the moments of z1^TR_1 and z2^Tr1 for three by two table
 */


/*
 FUNCTION  zr_moments_three_by_two
 PURPOSE   Exact means, variances and covariance of the free cells N_11, N_21 of a
 3 x 2 table at u-bar = (u1, u2), from the joint moments of (Q_1, Q_2)
 returned by zu_moments_three_treatment_log.  The covariance of two cells
 in the same outcome level but different treatments involves E[Q_1 Q_2].
 INPUT     n1, n2, n3 = N_i.;  m1, m2 = N_.j;  N;  gamma_delta (length 3);  u1, u2;
 shared_divisor.
 OUTPUT    c(E[N_11], E[N_21], Var[N_11], Var[N_21], Cov[N_11,N_21]).
 EXAMPLE   (R)
 zr_moments_three_by_two(60L, 10L, 20L, 15L, 75L, 90L, gamma_delta = c(0, 0, 1), u1 = 0L, u2 = 75L, shared_divisor = 1)
## u-bar at a corner of the allocation space (u1 = 0, u2 = N_.2)
 */

// [[Rcpp::export]]
NumericVector zr_moments_three_by_two(int n1, int n2, int n3, int m1, int m2, int N, NumericVector gamma_delta, int u1, int u2, double shared_divisor) {
  NumericVector zu_moments = zu_moments_three_treatment_log(n1,n2,n3,N,gamma_delta,u1+u2,shared_divisor);
  double z1u_first= zu_moments[0];
  double z2u_first= zu_moments[1];
  double z1u_second=zu_moments[2];
  double z2u_second=zu_moments[3];
  double z1u_z2u_product = zu_moments[4];
  double Us = u1 + u2;
  //if Us = 0 or N, becomes the usual RCT case
  if(Us==0 || Us==N)
    return(NumericVector::create((double) m1*n1/N,(double) m1*n2/N,
                                 (double) n1*(N-n1)*m1*(N-m1)/N/N/(N-1),(double) n2*(N-n2)*m1*(N-m1)/N/N/(N-1),
                                 (double) m1*(m1-N)*n1*n2/N/N/(N-1)));

  double q11 = u1/Us;
  double q10 = (m1-u1)/(N-Us);
  double z1r1_first = q11*z1u_first+q10*(n1-z1u_first);
  double z2r1_first = q11*z2u_first+q10*(n2-z2u_first);
  double w11, w10;
  if(Us==1)
    w11 = 0;
  else
    w11  = (double) u1/(Us-1)*(1-2*q11) + (double) Us/(Us-1)*q11*q11;

  if(Us==N-1)
    w10 = 0;
  else
    w10   = (double)(m1 - u1)/(N-Us-1)*(1-2*q10) + (double) (N-Us)/(N-Us-1)*pow(q10,2);
  double z1r1_var = (w11-w10)*z1u_first-(z1u_second)*(w11/Us+w10/(N-Us))
    +n1*(N-Us-n1+2*z1u_first)*w10/(N-Us)+(z1u_second-z1u_first*z1u_first)*(q11-q10)*(q11-q10);
  double z2r1_var = (w11-w10)*z2u_first-(z2u_second)*(w11/Us+w10/(N-Us))
    +n2*(N-Us-n2+2*z2u_first)*w10/(N-Us)+(z2u_second-z2u_first*z2u_first)*(q11-q10)*(q11-q10);

  double z1r1_z2r1_cov;
  if(Us>1 && Us<N-1){
    z1r1_z2r1_cov = (q11*q11+q10*q10-2*q11*q10)*(z1u_z2u_product-z1u_first*z2u_first)
    +u1*(u1-Us)/Us/Us/(Us-1)*z1u_z2u_product
    +(m1-u1)*(m1-u1-N+Us)/(N-Us)/(N-Us)/(N-Us-1)*(z1u_z2u_product+n1*n2)
    -(m1-u1)*(m1-u1-N+Us)/(N-Us)/(N-Us)/(N-Us-1)*(n1*z2u_first+n2*z1u_first);

  }
  if(Us==1){
    z1r1_z2r1_cov = (q11*q11+q10*q10-2*q11*q10)*(z1u_z2u_product-z1u_first*z2u_first)
    +(m1-u1)*(m1-u1-N+Us)/(N-Us)/(N-Us)/(N-Us-1)*(z1u_z2u_product+n1*n2)
    -(m1-u1)*(m1-u1-N+Us)/(N-Us)/(N-Us)/(N-Us-1)*(n1*z2u_first+n2*z1u_first);


  }
  if(Us==N-1){
    z1r1_z2r1_cov = (q11*q11+q10*q10-2*q11*q10)*(z1u_z2u_product-z1u_first*z2u_first)
    +u1*(u1-Us)/Us/Us/(Us-1)*z1u_z2u_product;

  }


  return(NumericVector::create(z1r1_first,z2r1_first,z1r1_var,z2r1_var,z1r1_z2r1_cov));

}





/*
 FUNCTION  zr_moments_three_by_three
 PURPOSE   Exact means, variances and all six pairwise covariances of the free
 cells N_11, N_12, N_21, N_22 of a 3 x 3 table at u-bar = (u1, u2, u3),
 from the joint moments of (Q_1, Q_2).  Covers all three covariance
 cases: same treatment and different outcome, different treatment and
 same outcome, and different treatment and different outcome.
 INPUT     n1, n2, n3 = N_i.;  m1, m2, m3 = N_.j;  N;  gamma_delta (length 3);
 u1, u2, u3;  shared_divisor.
 OUTPUT    length 14: c(E[N_11], E[N_12], E[N_21], E[N_22],
 Var[N_11], Var[N_12], Var[N_21], Var[N_22],
 Cov[N_11,N_12], Cov[N_11,N_21], Cov[N_11,N_22],
 Cov[N_12,N_21], Cov[N_12,N_22], Cov[N_21,N_22]).
 EXAMPLE   (R)
 zr_moments_three_by_three(35L, 24L, 27L, 24L, 36L, 26L, 86L, c(0, 1, 1), 0L, 10L, 20L, 1)
 */

// [[Rcpp::export]]
NumericVector zr_moments_three_by_three(int n1, int n2, int n3, int m1, int m2, int m3, int N, NumericVector gamma_delta, int u1, int u2, int u3, double shared_divisor) {
  NumericVector zu_moments = zu_moments_three_treatment_log(n1,n2,n3, N,gamma_delta,u1+u2+u3,shared_divisor);
  double z1u_first= zu_moments[0];
  double z2u_first = zu_moments[1];
  double z1u_second = zu_moments[2];
  double z2u_second = zu_moments[3];
  double z1u_z2u_product = zu_moments[4];
  double Us = u1 + u2+u3;
  //if Us = 0 or N, becomes the usual RCT case
  if(Us==0 || Us==N)
    return(NumericVector::create((double) m1*n1/N,(double)m2*n1/N,(double) m1*n2/N, (double) m2*n2/N,
                                 (double) n1*(N-n1)*m1*(N-m1)/N/N/(N-1),(double) n1*(N-n1)*m2*(N-m2)/N/N/(N-1),(double) n2*(N-n2)*m1*(N-m1)/N/N/(N-1),(double) n2*(N-n2)*m2*(N-m2)/N/N/(N-1),
                                 (double) m1*m2*n1*(n1-N)/N/N/(N-1), (double) n1*n2*m1*(m1-N)/N/N/(N-1),(double) m1*m2*n1*n2/N/N/(N-1),(double) m1*m2*n1*n2/N/N/(N-1), (double) n1*n2*m2*(m2-N)/N/N/(N-1),(double) m1*m2*n2*(n2-N)/N/N/(N-1)));

  double q11 = u1/Us;
  double q10 = (m1-u1)/(N-Us);
  double q21 = u2/Us;
  double q20 = (m2-u2)/(N-Us);
  double z1r1_first = q11*z1u_first+q10*(n1-z1u_first);
  double z1r2_first = q21*z1u_first+q20*(n1-z1u_first);
  double z2r1_first = q11*z2u_first+q10*(n2-z2u_first);
  double z2r2_first = q21*z2u_first+q20*(n2-z2u_first);
  double w11, w10, w21, w20;
  if(Us==1)
  { w11 = 0;
    w21 = 0;
  }
  else
  {w11  = (double) u1/(Us-1)*(1-2*q11) + (double) Us/(Us-1)*q11*q11;
    w21  = (double) u2/(Us-1)*(1-2*q21) + (double) Us/(Us-1)*q21*q21;}
  if(Us==N-1)
  {w10 = 0;
    w20 = 0;
  }
  else
  {
    w10   = (double)(m1 - u1)/(N-Us-1)*(1-2*q10) + (double) (N-Us)/(N-Us-1)*pow(q10,2);
    w20   = (double)(m2 - u2)/(N-Us-1)*(1-2*q20) + (double) (N-Us)/(N-Us-1)*pow(q20,2);

  }
  double z1r1_var = (w11-w10)*z1u_first-(z1u_second)*(w11/Us+w10/(N-Us))
    +n1*(N-Us-n1+2*z1u_first)*w10/(N-Us)+(z1u_second-z1u_first*z1u_first)*(q11-q10)*(q11-q10);
  double z1r2_var = (w21-w20)*z1u_first-(z1u_second)*(w21/Us+w20/(N-Us))
    +n1*(N-Us-n1+2*z1u_first)*w20/(N-Us)+(z1u_second-z1u_first*z1u_first)*(q21-q20)*(q21-q20);
  double z2r1_var = (w11-w10)*z2u_first-(z2u_second)*(w11/Us+w10/(N-Us))
    +n2*(N-Us-n2+2*z2u_first)*w10/(N-Us)+(z2u_second-z2u_first*z2u_first)*(q11-q10)*(q11-q10);
  double z2r2_var = (w21-w20)*z2u_first-(z2u_second)*(w21/Us+w20/(N-Us))
    +n2*(N-Us-n2+2*z2u_first)*w20/(N-Us)+(z2u_second-z2u_first*z2u_first)*(q21-q20)*(q21-q20);
  double z1r1_z1r2_cov, z1r1_z2r1_cov, z1r1_z2r2_cov, z1r2_z2r1_cov, z1r2_z2r2_cov, z2r1_z2r2_cov;
  if(Us>1 && Us <N-1){
    z1r1_z1r2_cov = q11*q21/(Us-1)*(z1u_second-Us*z1u_first)+q10*q20/(N-Us-1)*n1*(n1-N+Us)
    +q10*q20/(N-Us-1)*(z1u_second-(2*n1-N+Us)*z1u_first)
    +(q11*q21+q10*q20-q10*q21-q20*q11)*(z1u_second-z1u_first*z1u_first);

    z2r1_z2r2_cov = q11*q21/(Us-1)*(z2u_second-Us*z2u_first)+q10*q20/(N-Us-1)*n2*(n2-N+Us)
      +q10*q20/(N-Us-1)*(z2u_second-(2*n2-N+Us)*z2u_first)
      +(q11*q21+q10*q20-q10*q21-q20*q11)*(z2u_second-z2u_first*z2u_first);


      z1r1_z2r1_cov = (q11*q11+q10*q10-2*q11*q10)*(z1u_z2u_product-z1u_first*z2u_first)
        +u1*(u1-Us)/Us/Us/(Us-1)*z1u_z2u_product
      +(m1-u1)*(m1-u1-N+Us)/(N-Us)/(N-Us)/(N-Us-1)*(z1u_z2u_product+n1*n2)
      -(m1-u1)*(m1-u1-N+Us)/(N-Us)/(N-Us)/(N-Us-1)*(n1*z2u_first+n2*z1u_first);

      z1r2_z2r2_cov = (q21*q21+q20*q20-2*q21*q20)*(z1u_z2u_product-z1u_first*z2u_first)
        +u2*(u2-Us)/Us/Us/(Us-1)*z1u_z2u_product
      +(m2-u2)*(m2-u2-N+Us)/(N-Us)/(N-Us)/(N-Us-1)*(z1u_z2u_product+n1*n2)
      -(m2-u2)*(m2-u2-N+Us)/(N-Us)/(N-Us)/(N-Us-1)*(n1*z2u_first+n2*z1u_first);

      z1r1_z2r2_cov = z1r2_z2r1_cov = (q11*q21+q10*q20-q21*q10-q11*q20)*(z1u_z2u_product-z1u_first*z2u_first)
        +q11*q21/(Us-1)*z1u_z2u_product
      +(q10*q20)/(N-Us-1)*(z1u_z2u_product+n1*n2)
      -q20*q10/(N-Us-1)*(n1*z2u_first+n2*z1u_first);



  }
  if(Us==1){
    z1r1_z1r2_cov = q10*q20/(N-Us-1)*n1*(n1-N+Us)
    +q10*q20/(N-Us-1)*(z1u_second-(2*n1-N+Us)*z1u_first)
    +(q11*q21+q10*q20-q10*q21-q20*q11)*(z1u_second-z1u_first*z1u_first);

    z2r1_z2r2_cov = q10*q20/(N-Us-1)*n2*(n2-N+Us)
      +q10*q20/(N-Us-1)*(z2u_second-(2*n2-N+Us)*z2u_first)
      +(q11*q21+q10*q20-q10*q21-q20*q11)*(z2u_second-z2u_first*z2u_first);


      z1r1_z2r1_cov = (q11*q11+q10*q10-2*q11*q10)*(z1u_z2u_product-z1u_first*z2u_first)
        +(m1-u1)*(m1-u1-N+Us)/(N-Us)/(N-Us)/(N-Us-1)*(z1u_z2u_product+n1*n2)
        -(m1-u1)*(m1-u1-N+Us)/(N-Us)/(N-Us)/(N-Us-1)*(n1*z2u_first+n2*z1u_first);

        z1r2_z2r2_cov = (q21*q21+q20*q20-2*q21*q20)*(z1u_z2u_product-z1u_first*z2u_first)
          +(m2-u2)*(m2-u2-N+Us)/(N-Us)/(N-Us)/(N-Us-1)*(z1u_z2u_product+n1*n2)
          -(m2-u2)*(m2-u2-N+Us)/(N-Us)/(N-Us)/(N-Us-1)*(n1*z2u_first+n2*z1u_first);

          z1r1_z2r2_cov = z1r2_z2r1_cov = (q11*q21+q10*q20-q21*q10-q11*q20)*(z1u_z2u_product-z1u_first*z2u_first)
            +(q10*q20)/(N-Us-1)*(z1u_z2u_product+n1*n2)
            -q20*q10/(N-Us-1)*(n1*z2u_first+n2*z1u_first);



  }
  if(Us==N-1){
    z1r1_z1r2_cov = q11*q21/(Us-1)*(z1u_second-Us*z1u_first)
    +(q11*q21+q10*q20-q10*q21-q20*q11)*(z1u_second-z1u_first*z1u_first);

    z2r1_z2r2_cov = q11*q21/(Us-1)*(z2u_second-Us*z2u_first)
      +(q11*q21+q10*q20-q10*q21-q20*q11)*(z2u_second-z2u_first*z2u_first);


    z1r1_z2r1_cov = (q11*q11+q10*q10-2*q11*q10)*(z1u_z2u_product-z1u_first*z2u_first)
      +u1*(u1-Us)/Us/Us/(Us-1)*z1u_z2u_product;

    z1r2_z2r2_cov = (q21*q21+q20*q20-2*q21*q20)*(z1u_z2u_product-z1u_first*z2u_first)
      +u2*(u2-Us)/Us/Us/(Us-1)*z1u_z2u_product;


    z1r1_z2r2_cov = z1r2_z2r1_cov = (q11*q21+q10*q20-q21*q10-q11*q20)*(z1u_z2u_product-z1u_first*z2u_first)
      +q11*q21/(Us-1)*z1u_z2u_product;


  }


  return(NumericVector::create(z1r1_first,z1r2_first, z2r1_first, z2r2_first,
                               z1r1_var, z1r2_var, z2r1_var,z2r2_var,
                               z1r1_z1r2_cov, z1r1_z2r1_cov,
                               z1r1_z2r2_cov, z1r2_z2r1_cov,
                               z1r2_z2r2_cov, z2r1_z2r2_cov));

}

/*
 FUNCTION  zr_moments_four_by_two
 PURPOSE   Exact means, variances and pairwise covariances of the free cells
 N_11, N_21, N_31 of a 4 x 2 table at u-bar = (u1, u2), from the joint
 moments of (Q_1, Q_2, Q_3) returned by zu_moments_four_treatment_log.
 All three covariances are between cells in the same outcome level and
 different treatments.
 INPUT     n1..n4 = N_i.;  m1, m2 = N_.j;  N;  gamma_delta (length 4);  u1, u2;
 shared_divisor.
 OUTPUT    c(E[N_11], E[N_21], E[N_31], Var[N_11], Var[N_21], Var[N_31],
 Cov[N_11,N_21], Cov[N_11,N_31], Cov[N_21,N_31]).
 EXAMPLE   (R)
 zr_moments_four_by_two(10L, 10L, 10L, 10L, 12L, 28L, 40L, 0.5 * c(0, 0, 1, 1), 0L, 28L, 1)
 */

// [[Rcpp::export]]
NumericVector zr_moments_four_by_two(int n1, int n2, int n3, int n4, int m1, int m2, int N, NumericVector gamma_delta, int u1, int u2, double shared_divisor) {
  NumericVector zu_moments = zu_moments_four_treatment_log(n1,n2,n3,n4,N,gamma_delta,u1+u2,shared_divisor);
  double z1u_first= zu_moments[0];
  double z2u_first= zu_moments[1];
  double z3u_first= zu_moments[2];
  double z1u_second=zu_moments[3];
  double z2u_second=zu_moments[4];
  double z3u_second=zu_moments[5];
  double z1u_z2u_product=zu_moments[6];
  double z1u_z3u_product=zu_moments[7];
  double z2u_z3u_product=zu_moments[8];
  double Us = u1 + u2;
  //if Us = 0 or N, becomes the usual RCT case
  if(Us==0 || Us==N)
    return(NumericVector::create((double) m1*n1/N,(double) m1*n2/N,(double) m1*n3/N,
                                 (double) n1*(N-n1)*m1*(N-m1)/N/N/(N-1),(double) n2*(N-n2)*m1*(N-m1)/N/N/(N-1), (double) n3*(N-n3)*m1*(N-m1)/N/N/(N-1),
                                 (double) m1*(m1-N)*n1*n2/N/N/(N-1),(double) m1*(m1-N)*n1*n3/N/N/(N-1),(double) m1*(m1-N)*n2*n3/N/N/(N-1)));

  double q11 = u1/Us;
  double q10 = (m1-u1)/(N-Us);
  double z1r1_first = q11*z1u_first+q10*(n1-z1u_first);
  double z2r1_first = q11*z2u_first+q10*(n2-z2u_first);
  double z3r1_first = q11*z3u_first+q10*(n3-z3u_first);
  double w11, w10;
  if(Us==1)
    w11 = 0;
  else
    w11  = (double) u1/(Us-1)*(1-2*q11) + (double) Us/(Us-1)*q11*q11;

  if(Us==N-1)
    w10 = 0;
  else
    w10   = (double)(m1 - u1)/(N-Us-1)*(1-2*q10) + (double) (N-Us)/(N-Us-1)*pow(q10,2);
  double z1r1_var = (w11-w10)*z1u_first-(z1u_second)*(w11/Us+w10/(N-Us))
    +n1*(N-Us-n1+2*z1u_first)*w10/(N-Us)+(z1u_second-z1u_first*z1u_first)*(q11-q10)*(q11-q10);
  double z2r1_var = (w11-w10)*z2u_first-(z2u_second)*(w11/Us+w10/(N-Us))
    +n2*(N-Us-n2+2*z2u_first)*w10/(N-Us)+(z2u_second-z2u_first*z2u_first)*(q11-q10)*(q11-q10);
  double z3r1_var = (w11-w10)*z3u_first-(z3u_second)*(w11/Us+w10/(N-Us))
    +n3*(N-Us-n3+2*z3u_first)*w10/(N-Us)+(z3u_second-z3u_first*z3u_first)*(q11-q10)*(q11-q10);

  double z1r1_z2r1_cov, z1r1_z3r1_cov, z2r1_z3r1_cov;
  if(Us>1 && Us<N-1){
    z1r1_z2r1_cov = (q11*q11+q10*q10-2*q11*q10)*(z1u_z2u_product-z1u_first*z2u_first)
    +u1*(u1-Us)/Us/Us/(Us-1)*z1u_z2u_product
    +(m1-u1)*(m1-u1-N+Us)/(N-Us)/(N-Us)/(N-Us-1)*(z1u_z2u_product+n1*n2)
    -(m1-u1)*(m1-u1-N+Us)/(N-Us)/(N-Us)/(N-Us-1)*(n1*z2u_first+n2*z1u_first);

    z1r1_z3r1_cov = (q11*q11+q10*q10-2*q11*q10)*(z1u_z3u_product-z1u_first*z3u_first)
      +u1*(u1-Us)/Us/Us/(Us-1)*z1u_z3u_product
    +(m1-u1)*(m1-u1-N+Us)/(N-Us)/(N-Us)/(N-Us-1)*(z1u_z3u_product+n1*n3)
    -(m1-u1)*(m1-u1-N+Us)/(N-Us)/(N-Us)/(N-Us-1)*(n1*z3u_first+n3*z1u_first);

    z2r1_z3r1_cov = (q11*q11+q10*q10-2*q11*q10)*(z2u_z3u_product-z2u_first*z3u_first)
      +u1*(u1-Us)/Us/Us/(Us-1)*z2u_z3u_product
    +(m1-u1)*(m1-u1-N+Us)/(N-Us)/(N-Us)/(N-Us-1)*(z2u_z3u_product+n2*n3)
    -(m1-u1)*(m1-u1-N+Us)/(N-Us)/(N-Us)/(N-Us-1)*(n2*z3u_first+n3*z2u_first);

  }
  if(Us==1){
    z1r1_z2r1_cov = (q11*q11+q10*q10-2*q11*q10)*(z1u_z2u_product-z1u_first*z2u_first)
    +(m1-u1)*(m1-u1-N+Us)/(N-Us)/(N-Us)/(N-Us-1)*(z1u_z2u_product+n1*n2)
    -(m1-u1)*(m1-u1-N+Us)/(N-Us)/(N-Us)/(N-Us-1)*(n1*z2u_first+n2*z1u_first);

    z1r1_z3r1_cov = (q11*q11+q10*q10-2*q11*q10)*(z1u_z3u_product-z1u_first*z3u_first)
      +(m1-u1)*(m1-u1-N+Us)/(N-Us)/(N-Us)/(N-Us-1)*(z1u_z3u_product+n1*n3)
      -(m1-u1)*(m1-u1-N+Us)/(N-Us)/(N-Us)/(N-Us-1)*(n1*z3u_first+n3*z1u_first);

      z2r1_z3r1_cov = (q11*q11+q10*q10-2*q11*q10)*(z2u_z3u_product-z2u_first*z3u_first)
        +(m1-u1)*(m1-u1-N+Us)/(N-Us)/(N-Us)/(N-Us-1)*(z2u_z3u_product+n2*n3)
        -(m1-u1)*(m1-u1-N+Us)/(N-Us)/(N-Us)/(N-Us-1)*(n2*z3u_first+n3*z2u_first);




  }
  if(Us==N-1){
    z1r1_z2r1_cov = (q11*q11+q10*q10-2*q11*q10)*(z1u_z2u_product-z1u_first*z2u_first)
    +u1*(u1-Us)/Us/Us/(Us-1)*z1u_z2u_product;

    z1r1_z3r1_cov = (q11*q11+q10*q10-2*q11*q10)*(z1u_z3u_product-z1u_first*z3u_first)
      +u1*(u1-Us)/Us/Us/(Us-1)*z1u_z3u_product;

    z2r1_z3r1_cov = (q11*q11+q10*q10-2*q11*q10)*(z2u_z3u_product-z2u_first*z3u_first)
      +u1*(u1-Us)/Us/Us/(Us-1)*z2u_z3u_product;

  }


  return(NumericVector::create(z1r1_first,z2r1_first,z3r1_first,z1r1_var,z2r1_var,z3r1_var,z1r1_z2r1_cov,z1r1_z3r1_cov, z2r1_z3r1_cov));

}








/*
 * The following functions are kept for designer's numerically checking.
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 */





/*
 FUNCTION  zu_moments_two_treatment
 PURPOSE   Same quantity as zu_moments_two_treatment_log (E[Q_1], E[Q_1^2] for
 I = 2), computed by direct summation with choose_C instead of in log
 space; faster, but can overflow for large N or gamma.  shared_divisor
 rescales numerator and denominator alike.  The R wrappers call the
 _log version.
 INPUT     n1, n2 = N_i.;  N;  gamma_delta (length 2);  Us = u-bar;
 shared_divisor.
 OUTPUT    c(E[Q_1], E[Q_1^2]).
 EXAMPLE   (R)
 zu_moments_two_treatment(20L, 30L, 50L, c(0, 1), 10L, 1e6)
## equals zu_moments_two_treatment_log(20L, 30L, 50L, c(0, 1), 10L, 1e6) up to rounding
 */

// [[Rcpp::export]]
NumericVector zu_moments_two_treatment(int n1, int n2, int N, NumericVector gamma_delta, int Us, double shared_divisor) {
  double gamma = gamma_delta[0] - gamma_delta[1];
  double z1u_first = 0.0;
  double z1u_second = 0.0;
  double convolution_sum = 0.0;

  // Sanity check
  if (Us == 0) {   
    z1u_first = 0.0;
    z1u_second = 0.0;
    return NumericVector::create(z1u_first, z1u_second);
  }

  if (Us == N) {
    z1u_first = n1;
    z1u_second = std::pow(n1, 2);
    return NumericVector::create(z1u_first, z1u_second);
  }

  // General case
  int lb_q1 = std::max(0, n1 + Us - N);
  int ub_q1 = std::min(Us, n1);
  double z1u_first_cumulant = 0.0;
  double z1u_second_cumulant = 0.0;

  for (int q1 = lb_q1; q1 <= ub_q1; ++q1) {
    // Precompute the kernel
    double weighted_kernel = choose_C(Us, q1) * choose_C(N - Us, n1 - q1) / shared_divisor * exp(gamma * q1);

    convolution_sum += weighted_kernel;
    z1u_first_cumulant += q1 * weighted_kernel;
    z1u_second_cumulant += std::pow(q1, 2) * weighted_kernel;
  }

  // Check for division by zero
  if (convolution_sum == 0.0) {
    stop("Convolution sum is zero; division by zero detected.");
  }

  // Normalize cumulants by the convolution sum
  z1u_first = z1u_first_cumulant / convolution_sum;
  z1u_second = z1u_second_cumulant / convolution_sum;

  return NumericVector::create(z1u_first, z1u_second);
}


/*
 FUNCTION  zu_moments_three_treatment
 PURPOSE   Same quantity as zu_moments_three_treatment_log (the first two moments
 of Q_1 and Q_2 and the cross moment E[Q_1 Q_2], for I = 3), computed by
 direct summation instead of in log space; faster, but may overflow for
 large tables.  These moments feed the covariance computations.  The R
 wrappers call the _log version.
 INPUT     n1, n2, n3 = N_i.;  N;  gamma_delta (length 3);  Us = u-bar;
 shared_divisor.
 OUTPUT    c(E[Q_1], E[Q_2], E[Q_1^2], E[Q_2^2], E[Q_1 Q_2]).
 EXAMPLE   (R)
 zu_moments_three_treatment(35L, 24L, 27L, 86L, c(0, 1, 1), 30L, 1e6)
## equals zu_moments_three_treatment_log(35L, 24L, 27L, 86L, c(0, 1, 1), 30L, 1e6) up to rounding
 */

// [[Rcpp::export]]
NumericVector zu_moments_three_treatment(int n1, int n2, int n3, int N, NumericVector gamma_delta, int Us, double shared_divisor) {
  double gamma1 = gamma_delta[0] - gamma_delta[2];
  double gamma2 = gamma_delta[1] - gamma_delta[2];
  double gamma = (gamma1 != 0) ? gamma1 : (gamma2 != 0 ? gamma2 : 0);

  double convolution_sum = 0.0;

  // Sanity check
  if (Us == 0) {   
    return NumericVector::create(0.0, 0.0, 0.0, 0.0, 0.0);
  }

  if (Us == N) {
    return NumericVector::create(n1, n2, std::pow(n1, 2), std::pow(n2, 2), n1 * n2);
  }

  int lb_q = std::max(0, n1 + n2 + Us - N);
  int ub_q = std::min(n1 + n2, Us);
  int lb_q1 = std::max(0, n1 + Us - N);
  int ub_q1 = std::min(Us, n1);
  int lb_q2 = std::max(0, n2 + Us - N);
  int ub_q2 = std::min(n2, Us);

  double z1u_first_cumulant = 0.0, z1u_second_cumulant = 0.0, z2u_first_cumulant = 0.0;
  double z2u_second_cumulant = 0.0, z1u_z2u_product_cumulant = 0.0;

  for (int q = lb_q; q <= ub_q; ++q) {
    for (int q1 = lb_q1; q1 <= ub_q1; ++q1) {
      int q2 = q - q1;
      if (q2 > ub_q2 || q2 < lb_q2) {
        continue;
      }

      double exp_term = (gamma1 == 0) ? exp(gamma * (q - q1)) :
        (gamma2 == 0) ? exp(gamma * (q - q2)) : exp(gamma * q);

      double kernel = choose_C(Us, q1) * choose_C(Us - q1, q2) *
        choose_C(N - Us, n1 - q1) * choose_C((N - Us) - (n1 - q1), n2 - q2) / shared_divisor;

      double weighted_kernel = exp_term * kernel;

      convolution_sum += weighted_kernel;
      z1u_first_cumulant += q1 * weighted_kernel;
      z1u_second_cumulant += std::pow(q1, 2) * weighted_kernel;
      z2u_first_cumulant += q2 * weighted_kernel;
      z2u_second_cumulant += std::pow(q2, 2) * weighted_kernel;
      z1u_z2u_product_cumulant += q1 * q2 * weighted_kernel;
    }
  }

  if (convolution_sum == 0.0) stop("Convolution sum is zero; division by zero detected.");

  return NumericVector::create(
    z1u_first_cumulant / convolution_sum,
    z2u_first_cumulant / convolution_sum,
    z1u_second_cumulant / convolution_sum,
    z2u_second_cumulant / convolution_sum,
    z1u_z2u_product_cumulant / convolution_sum
  );
}



/*
 FUNCTION  zu_moments_four_treatment
 PURPOSE   Same quantity as zu_moments_four_treatment_log (the first two moments
 of Q_1, Q_2, Q_3 and their pairwise cross moments, for I = 4), computed
 by direct summation instead of in log space; faster, but may overflow
 for large tables.  These moments feed the covariance computations.  The
 R wrappers call the _log version.
 INPUT     n1..n4 = N_i.;  N;  gamma_delta (length 4);  Us = u-bar;
 shared_divisor.
 OUTPUT    length 9: c(E[Q_1], E[Q_2], E[Q_3], E[Q_1^2], E[Q_2^2], E[Q_3^2],
 E[Q_1 Q_2], E[Q_1 Q_3], E[Q_2 Q_3]).
 EXAMPLE   (R)
 zu_moments_four_treatment(10L, 10L, 10L, 10L, 40L, 0.5 * c(0, 0, 1, 1), 12L, 1e6)
## equals zu_moments_four_treatment_log(10L, 10L, 10L, 10L, 40L, 0.5 * c(0, 0, 1, 1), 12L, 1e6) up to rounding
 */

// [[Rcpp::export]]
NumericVector zu_moments_four_treatment(int n1, int n2, int n3, int n4, int N, NumericVector gamma_delta, int Us, double shared_divisor) {
  // Calculate gamma values
  double gamma1 = gamma_delta[0] - gamma_delta[3];
  double gamma2 = gamma_delta[1] - gamma_delta[3];
  double gamma3 = gamma_delta[2] - gamma_delta[3];
  double gamma = (gamma1 != 0) ? gamma1 : (gamma2 != 0 ? gamma2 : gamma3);

  double convolution_sum = 0.0;
  // Sanity check
  if (Us == 0) {   
    return NumericVector::create(0.0, 0.0, 0.0, 0.0, 0.0,0.0,0.0,0.0,0.0);
  }

  if (Us == N) {
    return NumericVector::create(n1, n2,n3,std::pow(n1, 2), std::pow(n2, 2),std::pow(n3,2), n1 * n2, n1*n3, n2*n3);
  }

  //general case, define the cumulant

  double z1u_first_cumulant = 0.0, z1u_second_cumulant = 0.0, z2u_first_cumulant = 0.0, z2u_second_cumulant = 0.0,
    z3u_first_cumulant = 0.0, z3u_second_cumulant = 0.0;
  double z1u_z2u_product_cumulant = 0.0, z1u_z3u_product_cumulant = 0.0, z2u_z3u_product_cumulant = 0.0;

  // Define bounds for q, q1, q2, and q3
  int lb_q = std::max(0, n1 + n2 + n3 + Us - N);
  int ub_q = std::min(n1 + n2 + n3, Us);
  int lb_q1 = std::max(0, n1 + Us - N);
  int ub_q1 = std::min(Us, n1);

  for (int q = lb_q; q <= ub_q; ++q) {
    for (int q1 = lb_q1; q1 <= ub_q1; ++q1) {
      int lb_q2 = std::max(0, n1 + n2 - q1 + Us - N);
      int ub_q2 = std::min(n2, Us - q1);

      for (int q2 = lb_q2; q2 <= ub_q2; ++q2) {
        int lb_q3 = std::max(0, n1 + n2 + n3 - q1 - q2 + Us - N);
        int ub_q3 = std::min(Us - q1 - q2, n3);
        int q3 = q - q1 - q2;

        if (q3 > ub_q3 || q3 < lb_q3) {
          continue;
        } else {
          // Determine the exponential term based on which gamma values are zero
          double exp_term;
          if (gamma1 == 0 && gamma2 == 0 && gamma3 == 0) {
            exp_term = exp(gamma * (q - q1 - q2 - q3));
          } else if (gamma1 == 0 && gamma2 == 0) {
            exp_term = exp(gamma * (q - q1 - q2));
          } else if (gamma1 == 0 && gamma3 == 0) {
            exp_term = exp(gamma * (q - q1 - q3));
          } else if (gamma2 == 0 && gamma3 == 0) {
            exp_term = exp(gamma * (q - q2 - q3));
          } else if (gamma1 == 0) {
            exp_term = exp(gamma * (q - q1));
          } else if (gamma2 == 0) {
            exp_term = exp(gamma * (q - q2));
          } else if (gamma3 == 0) {
            exp_term = exp(gamma * (q - q3));
          } else {
            exp_term = exp(gamma * q);
          }
          double weighted_kernel = exp_term / shared_divisor *
            choose_C(Us, q1) * choose_C(Us - q1, q2) * choose_C(Us - q1 - q2, q3) *
            choose_C(N - Us, n1 - q1) * choose_C(N - Us - n1 + q1, n2 - q2) *
            choose_C(N - Us - n1 + q1 - n2 + q2, n3 - q3);
          // Accumulate the convolution sum
          convolution_sum += weighted_kernel;
          z1u_first_cumulant += q1 * weighted_kernel;
          z1u_second_cumulant += std::pow(q1, 2) * weighted_kernel;
          z2u_first_cumulant += q2 * weighted_kernel;
          z2u_second_cumulant += std::pow(q2, 2) * weighted_kernel;
          z3u_first_cumulant += q3 * weighted_kernel;
          z3u_second_cumulant += std::pow(q3, 2) * weighted_kernel;
          z1u_z2u_product_cumulant += q1 * q2 * weighted_kernel;
          z1u_z3u_product_cumulant += q1 * q3 * weighted_kernel;
          z2u_z3u_product_cumulant += q2 * q3 * weighted_kernel;



        }
      }
    }
  }

  if (convolution_sum == 0.0) stop("Convolution sum is zero; division by zero detected.");

  return NumericVector::create(
    z1u_first_cumulant / convolution_sum,
    z2u_first_cumulant / convolution_sum,
    z3u_first_cumulant / convolution_sum,
    z1u_second_cumulant / convolution_sum,
    z2u_second_cumulant / convolution_sum,
    z3u_second_cumulant / convolution_sum,
    z1u_z2u_product_cumulant / convolution_sum,
    z1u_z3u_product_cumulant / convolution_sum,
    z2u_z3u_product_cumulant / convolution_sum

  );
}

















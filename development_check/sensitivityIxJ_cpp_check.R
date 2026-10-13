# =============================================================================
# sensitivityIxJ_cpp_check.R
#
# Code for "Exact, Nonparametric Sensitivity Analysis for Observational Studies
# of Contingency Tables"
# Authors: Elaine K. Chiu and Hyunseung Kang
#
# Unit tests for library/sensitivityIxJ.cpp.
# This file runs the building block functions in sensitivityIxJ.cpp to check
# they work as intended by comparing them with some special cases/manual
# computation results.
#
# Run from the top-level folder (the one with README.txt):
#   source("development_check/sensitivityIxJ_cpp_check.R")
# Each check prints PASS or FAIL; a summary is printed at the end.
# =============================================================================


library(Rcpp)
library(BiasedUrn)

# ---- Setup ------------------------------------------------------------------
# Wrappers that expose the internal helpers of library/sensitivityIxJ.cpp to R.
wrapper_code <- '
// [[Rcpp::export]]
double lchoose_test_(int n, int k) { return lchoose_(n, k); }

// [[Rcpp::export]]
double lse_test_(double a, double b) { return lse(a, b); }

// [[Rcpp::export]]
double log_S_test_(IntegerVector u, IntegerVector m, IntegerVector TD, double gamma) {
  std::vector<int> uu(u.begin(), u.end()), mm(m.begin(), m.end()), tt(TD.begin(), TD.end());
  return log_S(uu.data(), mm.data(), tt.data(), u.size(), gamma);
}

// [[Rcpp::export]]
double Factorial_test_(double n) { return Factorial(n); }

// [[Rcpp::export]]
double log_choose_test_(int n, int k) { return log_choose(n, k); }

// [[Rcpp::export]]
double compute_log_sum_exp_test_(NumericVector x) {
  std::vector<double> v(x.begin(), x.end());
  return compute_log_sum_exp(v);
}
'

# Compile the library code together with the wrappers.
lib_code <- readLines("library/sensitivityIxJ.cpp")
sourceCpp(code = paste(c(lib_code, wrapper_code), collapse = "\n"))

# Record and print the result of one check.
results <- logical(0)
check <- function(label, expr) {
  ok <- tryCatch(isTRUE(expr), error = function(e) FALSE)
  results[label] <<- ok
  cat(if (ok) "PASS  " else "FAIL  ", label, "\n", sep = "")
}

# TRUE if x and y agree within numerical tolerance.
near <- function(x, y, tol = 1e-10) isTRUE(all.equal(x, y, tolerance = tol))


# =============================================================================
# lchoose_(n, k): log of the binomial coefficient C(n, k).
#                 Returns -Inf when k < 0, n < 0 or k > n.
# =============================================================================

# Purpose:  small values match hand calculation
# Input:    (5, 2), (10, 3), (6, 0), (6, 6)
# Expected: log(10), log(120), 0, 0
# Paper:    building block of the counting kernels
check("lchoose_: small values",
      near(c(lchoose_test_(5L, 2L), lchoose_test_(10L, 3L),
             lchoose_test_(6L, 0L), lchoose_test_(6L, 6L)),
           c(log(10), log(120), 0, 0)))

# Purpose:  agrees with base R for all n up to 60
# Input:    every (n, k) with 0 <= k <= n <= 60
# Expected: same as base R lchoose(n, k)
# Paper:    building block of the counting kernels
grid <- do.call(rbind, lapply(0:60, function(n) cbind(n = n, k = 0:n)))
check("lchoose_: matches base R",
      near(mapply(lchoose_test_, grid[, "n"], grid[, "k"]),
           lchoose(grid[, "n"], grid[, "k"])))

# Purpose:  impossible combinations give -Inf (count of 0)
# Input:    (3, 5), (3, -1), (-1, 0)
# Expected: -Inf
# Paper:    kernels are zero outside the feasible range
check("lchoose_: out of range gives -Inf",
      identical(c(lchoose_test_(3L, 5L), lchoose_test_(3L, -1L),
                  lchoose_test_(-1L, 0L)),
                rep(-Inf, 3)))

# Purpose:  no overflow for large n (C(1000, 500) is about 1e299)
# Input:    (1000, 500), (5000, 2500)
# Expected: finite, same as base R
# Paper:    log scale keeps the exact method feasible for large N
check("lchoose_: large n",
      near(c(lchoose_test_(1000L, 500L), lchoose_test_(5000L, 2500L)),
           c(lchoose(1000, 500), lchoose(5000, 2500))))



# =============================================================================
# lse(a, b): log(exp(a) + exp(b)), computed without overflow.
#            -Inf is treated as log 0.
# =============================================================================

# Purpose:  small values match direct calculation
# Input:    (log 2, log 3), (0, 0), (log 1, log 9)
# Expected: log(5), log(2), log(10)
# Paper:    accumulates sums of counting kernels on the log scale
check("lse: small values",
      near(c(lse_test_(log(2), log(3)), lse_test_(0, 0), lse_test_(log(1), log(9))),
           c(log(5), log(2), log(10))))

# Purpose:  -Inf acts as log 0, so adding it leaves the other term unchanged
# Input:    (-Inf, log 7), (log 7, -Inf), (-Inf, -Inf)
# Expected: log(7), log(7), -Inf
# This definition is for simplicity that if there is any impossible
# combination it contributes zero probability mass 
check("lse: -Inf is the identity",
      near(lse_test_(-Inf, log(7)), log(7)) &&
        near(lse_test_(log(7), -Inf), log(7)) &&
        identical(lse_test_(-Inf, -Inf), -Inf))

# Purpose:  order of the arguments does not matter
# Input:    100 random pairs (a, b) in [-50, 50]
# Expected: lse(a, b) equal to lse(b, a)
# Paper:    accumulates sums of counting kernels on the log scale
set.seed(1)
a <- runif(100, -50, 50)
b <- runif(100, -50, 50)
check("lse: symmetric",
      near(mapply(lse_test_, a, b), mapply(lse_test_, b, a)))

# Purpose:  agrees with a reference log-sum-exp on random pairs
# Input:    the same 100 random pairs
# Expected: max(a, b) + log(1 + exp(-|a - b|))
# Paper:    accumulates sums of counting kernels on the log scale
check("lse: matches reference formula",
      near(mapply(lse_test_, a, b), pmax(a, b) + log1p(exp(-abs(a - b)))))

# Purpose:  no overflow or underflow at extreme magnitudes
#           (exp(1000) overflows and exp(-1000) underflows in double precision)
# Input:    (1000, 1000), (-1000, -1000), (0, -800)
# Expected: 1000 + log(2), -1000 + log(2), 0
# Paper:    log scale keeps the exact method feasible for large N
check("lse: extreme magnitudes",
      near(c(lse_test_(1000, 1000), lse_test_(-1000, -1000), lse_test_(0, -800)),
           c(1000 + log(2), -1000 + log(2), 0)))



# =============================================================================
# log_S(u, m, TD, J, gamma): log of S(t-bar, u), the factor of the numerator
#                            that depends on u.
# =============================================================================

# Purpose:  one outcome level, value worked out by hand
# Input:    u = 2, m = 4, TD = 2, gamma = 0.7
# Expected: log(1 + 4 e^gamma + e^(2 gamma))
# Paper:    factor of the numerator that depends on u
check("log_S: one level, hand calculation",
      near(log_S_test_(2L, 4L, 2L, 0.7), log(1 + 4 * exp(0.7) + exp(1.4))))

# Purpose:  several outcome levels, value worked out by hand
# Input:    u = (1, 1, 1), m = (2, 2, 2), TD = (1, 1, 1), gamma = 1.5
# Expected: 3 * log(1 + e^gamma)
# Paper:    factor of the numerator that depends on u
check("log_S: three levels, hand calculation",
      near(log_S_test_(c(1L, 1L, 1L), c(2L, 2L, 2L), c(1L, 1L, 1L), 1.5),
           3 * log(1 + exp(1.5))))

# Purpose:  at gamma = 0, S(t-bar, u) = prod_j C(N_.j, t-bar_j) for any u
# Input:    m = (24, 36, 26), TD = (12, 18, 21), u = (0, 10, 20) and (5, 5, 5)
# Expected: sum(lchoose(m, TD))
# Paper:    random assignment as the special case gamma = 0
m_S  <- c(24L, 36L, 26L)
TD_S <- c(12L, 18L, 21L)
check("log_S: gamma = 0 identity",
      near(log_S_test_(c(0L, 10L, 20L), m_S, TD_S, 0), sum(lchoose(m_S, TD_S))) &&
        near(log_S_test_(c(5L, 5L, 5L),   m_S, TD_S, 0), sum(lchoose(m_S, TD_S))))

# Purpose:  with no subject having u_s = 1, gamma has no effect
# Input:    u = (0, 0, 0), same m and TD, gamma = 2
# Expected: sum(lchoose(m, TD))
# Paper:    factor of the numerator that depends on u
check("log_S: u = 0 gives the gamma = 0 value",
      near(log_S_test_(c(0L, 0L, 0L), m_S, TD_S, 2), sum(lchoose(m_S, TD_S))))

# Purpose:  with every subject having u_s = 1, n-bar_j = t-bar_j is forced
# Input:    u = m, same TD, gamma = 0.5
# Expected: 0.5 * sum(TD) + sum(lchoose(m, TD))
# Paper:    factor of the numerator that depends on u
check("log_S: u = N_.j forces n-bar = t-bar",
      near(log_S_test_(m_S, m_S, TD_S, 0.5), 0.5 * sum(TD_S) + sum(lchoose(m_S, TD_S))))

# Purpose:  no admissible n-bar gives -Inf (contribution 0)
# Input:    u = (1, 1), m = (2, 2), TD = (3, 0), where t-bar_1 > N_.1
# Expected: -Inf
# Paper:    kernels are zero outside the feasible range
check("log_S: infeasible t-bar gives -Inf",
      identical(log_S_test_(c(1L, 1L), c(2L, 2L), c(3L, 0L), 1), -Inf))

# Purpose:  S(t-bar, u) is a sum of exp(gamma q-bar) terms, so it cannot
#           decrease as gamma grows
# Input:    u = (0, 10, 20), same m and TD, gamma from 0 to 5 in steps of 0.25
# Expected: all successive differences >= 0
# Paper:    larger gamma means stronger unmeasured confounding
v_S <- vapply(seq(0, 5, by = 0.25),
              function(g) log_S_test_(c(0L, 10L, 20L), m_S, TD_S, g), numeric(1))
check("log_S: nondecreasing in gamma", all(diff(v_S) >= -1e-12))

# Purpose:  no overflow for large margins
# Input:    m = (200, 300, 250), TD = (100, 150, 120), u = (50, 100, 80)
# Expected: finite at gamma = 1; equal to sum(lchoose(m, TD)) at gamma = 0
# Paper:    log scale keeps the exact method feasible for large N
mL <- c(200L, 300L, 250L); TL <- c(100L, 150L, 120L); uL <- c(50L, 100L, 80L)
check("log_S: large margins",
      is.finite(log_S_test_(uL, mL, TL, 1)) &&
        near(log_S_test_(uL, mL, TL, 0), sum(lchoose(mL, TL))))


# =============================================================================
# log_den_bd(ubar, N, ND, gamma): log of S_0(u-bar), the factor of the
#                                 denominator that depends on u.
# =============================================================================

# Purpose:  value worked out by hand
# Input:    ubar = 2, N = 4, ND = 2, gamma = 0.7
# Expected: log(1 + 4 e^gamma + e^(2 gamma))
# Paper:    factor of the denominator that depends on u
check("log_den_bd: hand calculation (ubar = 2)",
      near(log_den_bd(2L, 4L, 2L, 0.7), log(1 + 4 * exp(0.7) + exp(1.4))))

# Purpose:  value worked out by hand, unequal group sizes
# Input:    ubar = 1, N = 3, ND = 1, gamma = 1.2
# Expected: log(2 + e^gamma)
# Paper:    factor of the denominator that depends on u
check("log_den_bd: hand calculation (ubar = 1)",
      near(log_den_bd(1L, 3L, 1L, 1.2), log(2 + exp(1.2))))

# Purpose:  at gamma = 0, S_0 = C(N, ND) for any ubar
# Input:    N = 86, ND = 51, ubar = 0, 30, 86
# Expected: lchoose(86, 51)
# Paper:    random assignment as the special case gamma = 0
check("log_den_bd: gamma = 0 identity",
      near(c(log_den_bd(0L, 86L, 51L, 0), log_den_bd(30L, 86L, 51L, 0),
             log_den_bd(86L, 86L, 51L, 0)),
           rep(lchoose(86, 51), 3)))

# Purpose:  with no subject having u_s = 1, gamma has no effect
# Input:    ubar = 0, N = 20, ND = 8, gamma = 1.7
# Expected: lchoose(20, 8)
# Paper:    factor of the denominator that depends on u
check("log_den_bd: ubar = 0",
      near(log_den_bd(0L, 20L, 8L, 1.7), lchoose(20, 8)))

# Purpose:  with every subject having u_s = 1, q-bar = ND is forced
# Input:    ubar = 20, N = 20, ND = 8, gamma = 1.7
# Expected: 1.7 * 8 + lchoose(20, 8)
# Paper:    factor of the denominator that depends on u
check("log_den_bd: ubar = N",
      near(log_den_bd(20L, 20L, 8L, 1.7), 1.7 * 8 + lchoose(20, 8)))

# Purpose:  no subject assigned to I_1, or every subject assigned to I_1
# Input:    ND = 0 and ND = N, with ubar = 6, N = 15, gamma = 0.9
# Expected: 0 and 0.9 * 6
# Paper:    factor of the denominator that depends on u
check("log_den_bd: ND = 0 and ND = N",
      near(log_den_bd(6L, 15L, 0L, 0.9), 0) &&
        near(log_den_bd(6L, 15L, 15L, 0.9), 0.9 * 6))

# Purpose:  S_0(u-bar) is S(t-bar, u) with a single outcome level
# Input:    (ubar, N, ND) = (30, 86, 51), (5, 12, 7); gamma = 0.4, 2.5
# Expected: log_S with u = ubar, m = N, TD = ND
# Paper:    numerator and denominator share the same counting structure
check("log_den_bd: matches log_S with one outcome level",
      near(log_den_bd(30L, 86L, 51L, 0.4), log_S_test_(30L, 86L, 51L, 0.4)) &&
        near(log_den_bd(5L, 12L, 7L, 2.5),   log_S_test_(5L, 12L, 7L, 2.5)))

# Purpose:  S_0 is a sum of exp(gamma q-bar) terms, so it cannot decrease
#           as gamma grows
# Input:    ubar = 30, N = 86, ND = 51, gamma from 0 to 5 in steps of 0.25
# Expected: all successive differences >= 0
# Paper:    larger gamma means stronger unmeasured confounding
v_D <- vapply(seq(0, 5, by = 0.25),
              function(g) log_den_bd(30L, 86L, 51L, g), numeric(1))
check("log_den_bd: nondecreasing in gamma", all(diff(v_D) >= -1e-12))

# Purpose:  no overflow for large N
# Input:    ubar = 300, N = 1000, ND = 400
# Expected: finite at gamma = 2; equal to lchoose(1000, 400) at gamma = 0
# Paper:    log scale keeps the exact method feasible for large N
check("log_den_bd: large N",
      is.finite(log_den_bd(300L, 1000L, 400L, 2)) &&
        near(log_den_bd(300L, 1000L, 400L, 0), lchoose(1000, 400)))




# =============================================================================
# log_c_t(tab, delta): log of kappa(t), the factor of the numerator for one
#                      table t that does not depend on u.
# =============================================================================

tab_C <- matrix(c(1L, 2L,
                  1L, 1L,
                  2L, 0L), nrow = 3, byrow = TRUE)
tab_B <- matrix(c(12L, 18L,  5L,
                  6L, 12L,  6L,
                  6L,  6L, 15L), nrow = 3, byrow = TRUE)
d_B   <- c(0L, 1L, 1L)

# Purpose:  value worked out by hand
# Input:    tab_C (3 x 2), delta = (0, 1, 1)
# Expected: log(3)
# Paper:    factor of the numerator that does not depend on u
check("log_c_t: hand calculation",
      near(log_c_t(tab_C, c(0L, 1L, 1L)), log(3)))

# Purpose:  with one treatment in I_1 and one outside, there is only one way
#           to divide the subjects
# Input:    2 x 3 table, delta = (0, 1)
# Expected: 0
check("log_c_t: two treatment levels give 0",
      near(log_c_t(matrix(c(3L, 1L, 4L, 1L, 5L, 9L), nrow = 2, byrow = TRUE),
                   c(0L, 1L)), 0))

# Purpose:  with delta = 0, kappa(t) = prod_j N_.j! / prod_i t_ij!
# Input:    tab_B, delta = (0, 0, 0)
# Expected: sum(lfactorial(colSums(tab_B))) - sum(lfactorial(tab_B))
check("log_c_t: delta = 0",
      near(log_c_t(tab_B, c(0L, 0L, 0L)),
           sum(lfactorial(colSums(tab_B))) - sum(lfactorial(tab_B))))

# Purpose:  swapping delta with 1 - delta leaves kappa(t) unchanged
# Input:    tab_B, delta = (0, 1, 1) and (1, 0, 0)
# Expected: equal values
check("log_c_t: invariant to complementing delta",
      near(log_c_t(tab_B, d_B), log_c_t(tab_B, 1L - d_B)))

# Purpose:  reordering the treatments within I_1 leaves kappa(t) unchanged
# Input:    tab_B with rows 2 and 3 swapped, delta = (0, 1, 1)
# Expected: same value as for tab_B
check("log_c_t: invariant to reordering rows within I_1",
      near(log_c_t(tab_B[c(1, 3, 2), ], d_B), log_c_t(tab_B, d_B)))

# Purpose:  kappa(t) prod_j C(N_.j, t-bar_j) is the number of treatment
#           assignments that produce table t, for any delta
# Input:    tab_B, delta = (0, 1, 1), (1, 0, 1), (0, 0, 1)
# Expected: sum(lfactorial(colSums(tab_B))) - sum(lfactorial(tab_B))
n_assign <- sum(lfactorial(colSums(tab_B))) - sum(lfactorial(tab_B))
check("log_c_t: recovers the number of assignments",
      all(vapply(list(c(0L, 1L, 1L), c(1L, 0L, 1L), c(0L, 0L, 1L)), function(d) {
        tbar <- colSums(tab_B[d == 1L, , drop = FALSE])
        near(log_c_t(tab_B, d) + sum(lchoose(colSums(tab_B), tbar)), n_assign)
      }, logical(1))))

# Purpose:  no overflow for large counts
# Input:    3 x 3 table with cell counts in the hundreds, delta = (0, 1, 1)
# Expected: finite value
tab_L <- matrix(c(150L, 200L, 120L,
                  100L, 180L,  90L,
                  130L, 110L, 160L), nrow = 3, byrow = TRUE)
check("log_c_t: large counts", is.finite(log_c_t(tab_L, d_B)))



# =============================================================================
# log_c_0(treatment_margins, delta): log of kappa_0, the factor of the
#                                    denominator that does not depend on u.
# =============================================================================

n_0 <- c(35, 24, 27)
d_0 <- c(0, 1, 1)

# Purpose:  value worked out by hand
# Input:    treatment_margins = (1, 1, 2), delta = (0, 1, 1)
# Expected: log(3)
check("log_c_0: hand calculation",
      near(log_c_0(c(1, 1, 2), c(0, 1, 1)), log(3)))

# Purpose:  with one treatment outside I_1, kappa_0 is a binomial coefficient
# Input:    treatment_margins = (35, 24, 27), delta = (0, 1, 1)
# Expected: lchoose(51, 24)
check("log_c_0: one treatment outside I_1",
      near(log_c_0(n_0, d_0), lchoose(51, 24)))

# Purpose:  with one treatment in I_1 and one outside, there is only one way
#           to divide the subjects
# Input:    treatment_margins = (20, 30), delta = (0, 1)
# Expected: 0
check("log_c_0: two treatment levels give 0",
      near(log_c_0(c(20, 30), c(0, 1)), 0))

# Purpose:  with delta = 0, kappa_0 = N! / prod_i N_i.!
# Input:    treatment_margins = (35, 24, 27), delta = (0, 0, 0)
# Expected: lfactorial(86) - sum(lfactorial(n))
check("log_c_0: delta = 0",
      near(log_c_0(n_0, c(0, 0, 0)), lfactorial(86) - sum(lfactorial(n_0))))

# Purpose:  swapping delta with 1 - delta leaves kappa_0 unchanged
# Input:    treatment_margins = (35, 24, 27), delta = (0, 1, 1) and (1, 0, 0)
# Expected: equal values
check("log_c_0: invariant to complementing delta",
      near(log_c_0(n_0, d_0), log_c_0(n_0, 1 - d_0)))

# Purpose:  kappa_0 C(N, N-bar) is the total number of treatment assignments,
#           for any delta
# Input:    treatment_margins = (35, 24, 27), delta = (0, 1, 1), (1, 0, 1), (0, 0, 1)
# Expected: lfactorial(86) - sum(lfactorial(n))
check("log_c_0: recovers the number of assignments",
      all(vapply(list(c(0, 1, 1), c(1, 0, 1), c(0, 0, 1)), function(d)
        near(log_c_0(n_0, d) + lchoose(sum(n_0), sum(n_0[d == 1])),
             lfactorial(86) - sum(lfactorial(n_0))),
        logical(1))))

# Purpose:  kappa_0 equals kappa(t) for the one-column table of treatment margins
# Input:    treatment_margins = (35, 24, 27) as a 3 x 1 table, delta = (0, 1, 1)
# Expected: same value as log_c_t
check("log_c_0: matches log_c_t for a one-column table",
      near(log_c_0(n_0, d_0), log_c_t(matrix(n_0, ncol = 1), d_0)))

# Purpose:  no overflow for large margins
# Input:    treatment_margins = (400, 350, 250), delta = (0, 1, 1)
# Expected: lchoose(600, 350)
check("log_c_0: large margins",
      near(log_c_0(c(400, 350, 250), d_0), lchoose(600, 350)))





# =============================================================================
# point_prob_bd(tab, delta, u, gamma): probability of one table t under the
#                                      sensitivity model at a binary u.
# In this section we mainly apply our understanding for hypergeometric
# to check the correctness. First, when either gamma = 0 or all u entries are 
# 1 or 0, then the probability mass collapse to the one from either multivariate
# hypergeometric or hypergeometric case
# =============================================================================

# Multivariate hypergeometric probability of drawing counts x from an urn with
# n[j] balls of colour j.
dmvhyper <- function(x, n) prod(choose(n, x)) / choose(sum(n), sum(x))

# Probability of a table at gamma = 0, drawing the rows one after another.
dtable0 <- function(tab) {
  remaining <- colSums(tab)
  p <- 1
  for (i in seq_len(nrow(tab))) {
    p <- p * dmvhyper(tab[i, ], remaining)
    remaining <- remaining - tab[i, ]
  }
  p
}

tab_2 <- matrix(c(3, 2,
                  1, 4), nrow = 2, byrow = TRUE)   # N_i. = (5, 5), N_.j = (4, 6)

# Purpose:  at gamma = 0, a 2 x 2 table follows the hypergeometric distribution
# Input:    tab_2, delta = (0, 1), gamma = 0, u = (0, 0) and (2, 3)
# Expected: dhyper(3, 4, 6, 5)
check("point_prob_bd: 2 x 2 at gamma = 0 matches dhyper",
      near(point_prob_bd(tab_2, c(0, 1), c(0, 0), 0), dhyper(3, 4, 6, 5)) &&
        near(point_prob_bd(tab_2, c(0, 1), c(2, 3), 0), dhyper(3, 4, 6, 5)))

# Purpose:  when every subject has the same u, gamma has no effect
#           so the probability mass collapse to the one from 
#           hypergeometric again
# Input:    tab_2, delta = (0, 1), gamma = 2, u = (0, 0) and u = (4, 6)
# Expected: dhyper(3, 4, 6, 5)
check("point_prob_bd: constant u matches dhyper for any gamma",
      near(point_prob_bd(tab_2, c(0, 1), c(0, 0), 2), dhyper(3, 4, 6, 5)) &&
        near(point_prob_bd(tab_2, c(0, 1), c(4, 6), 2), dhyper(3, 4, 6, 5)))

# Purpose:  at gamma = 0, an I x J table follows the multivariate
#           hypergeometric distribution
# Input:    tab_B (3 x 3), delta = (0, 1, 1), gamma = 0, u = (0, 10, 20)
# Expected: dtable0(tab_B)
check("point_prob_bd: 3 x 3 at gamma = 0 matches multivariate hypergeometric",
      near(point_prob_bd(tab_B, c(0, 1, 1), c(0, 10, 20), 0), dtable0(tab_B)))

# Purpose:  with u_s = 1 exactly for the subjects with outcome 2, a 2 x 2
#           table follows Fisher's noncentral hypergeometric distribution
#           with odds ratio exp(gamma)
#           This is based on our Corollary to Theorem 3
# Input:    tab_2, u = (0, 6), gamma = 1.2, delta = (0, 1) and (1, 0)
# Expected: dFNCHypergeo(t_22, 6, 4, 5, exp(1.2)) and
#           dFNCHypergeo(t_12, 6, 4, 5, exp(1.2))
check("point_prob_bd: matches Fisher's noncentral hypergeometric",
      near(point_prob_bd(tab_2, c(0, 1), c(0, 6), 1.2),
           BiasedUrn::dFNCHypergeo(4, 6, 4, 5, exp(1.2)), tol = 1e-6) &&
        near(point_prob_bd(tab_2, c(1, 0), c(0, 6), 1.2),
             BiasedUrn::dFNCHypergeo(2, 6, 4, 5, exp(1.2)), tol = 1e-6))

# Purpose:  value worked out by hand with u not tied to the outcome
# Input:    table with all cells 1, delta = (0, 1), u = (1, 1), gamma = 0.9
# Expected: (1 + e^gamma)^2 / (1 + 4 e^gamma + e^(2 gamma))
check("point_prob_bd: hand calculation",
      near(point_prob_bd(matrix(1, 2, 2), c(0, 1), c(1, 1), 0.9),
           (1 + exp(0.9))^2 / (1 + 4 * exp(0.9) + exp(1.8))))

# Purpose:  probabilities of all 2 x 3 tables with fixed margins sum to 1
# Input:    N_i. = (3, 4), N_.j = (2, 3, 2), delta = (0, 1), u = (1, 2, 1),
#           gamma = 1.5
# Expected: 1
cm_23 <- c(2, 3, 2)
p_23 <- 0
for (a in 0:2) for (b in 0:3) {
  c3 <- 3 - a - b
  if (c3 < 0 || c3 > 2) next
  row1 <- c(a, b, c3)
  p_23 <- p_23 + point_prob_bd(rbind(row1, cm_23 - row1), c(0, 1), c(1, 2, 1), 1.5)
}
check("point_prob_bd: 2 x 3 probabilities sum to 1", near(p_23, 1))

# Purpose:  probabilities of all 3 x 2 tables with fixed margins sum to 1
# Input:    N_i. = (2, 2, 2), N_.j = (3, 3), delta = (0, 1, 1), u = (2, 1),
#           gamma = 0.8
# Expected: 1
rm_32 <- c(2, 2, 2)
p_32 <- 0
for (a in 0:2) for (b in 0:2) {
  c3 <- 3 - a - b
  if (c3 < 0 || c3 > 2) next
  col1 <- c(a, b, c3)
  p_32 <- p_32 + point_prob_bd(cbind(col1, rm_32 - col1), c(0, 1, 1), c(2, 1), 0.8)
}
check("point_prob_bd: 3 x 2 probabilities sum to 1", near(p_32, 1))

# Purpose:  invalid inputs stop with an error
# Input:    tab_2 with u of length 3; delta of length 3; u_1 = 5 > N_.1 = 4
# Expected: an error in each case
is_error <- function(expr) inherits(try(expr, silent = TRUE), "try-error")
check("point_prob_bd: invalid inputs give errors",
      is_error(point_prob_bd(tab_2, c(0, 1), c(0, 0, 0), 1)) &&
        is_error(point_prob_bd(tab_2, c(0, 1, 1), c(0, 0), 1)) &&
        is_error(point_prob_bd(tab_2, c(0, 1), c(5, 0), 1)))

# Purpose:  no overflow for large counts
# Input:    tab_L (3 x 3, counts in the hundreds), delta = (0, 1, 1),
#           u = (50, 100, 80), gamma = 1
# Expected: a finite value in [0, 1]
p_L <- point_prob_bd(tab_L, c(0, 1, 1), c(50, 100, 80), 1)
check("point_prob_bd: large counts", is.finite(p_L) && p_L >= 0 && p_L <= 1)



# =============================================================================
# alpha_bd(TDmat, logW, u_space, m, gamma, N, ND, logc0): exact p-value at
#                                                        every candidate u-bar.
# =============================================================================

# Stage 1 (build_W_cpp) followed by stage 2 (alpha_bd) for one table.
run_alpha <- function(tab, delta, score, threshold, greater, u_space, gamma) {
  W <- build_W_cpp(rowSums(tab), colSums(tab), delta, threshold, greater, score)
  alpha_bd(W$TDmat, W$logW, u_space, colSums(tab), gamma,
           sum(tab), sum(rowSums(tab)[delta == 1]), log_c_0(rowSums(tab), delta))
}

sc_2 <- matrix(c(1, 0,
                 0, 0), nrow = 2, byrow = TRUE)    # T = t_11

# Purpose:  at gamma = 0, a 2 x 2 table gives the one-sided Fisher exact test
# Input:    tab_2, T = t_11, threshold 3, delta = (0, 1), u = (0, 0), gamma = 0
# Expected: fisher.test(tab_2, alternative = "greater")$p.value,
#           equal to phyper(2, 4, 6, 5, lower.tail = FALSE)
a_2 <- run_alpha(tab_2, c(0, 1), sc_2, 3, TRUE, matrix(c(0, 0), nrow = 1), 0)
check("alpha_bd: 2 x 2 at gamma = 0 matches fisher.test",
      near(a_2, fisher.test(tab_2, alternative = "greater")$p.value) &&
        near(a_2, phyper(2, 4, 6, 5, lower.tail = FALSE)))

# Purpose:  at gamma = 0 the p-value does not depend on u
# Input:    tab_2, T = t_11, threshold 3, u = (0, 0), (2, 3), (4, 0), gamma = 0
# Expected: 66 / 252 for every row
check("alpha_bd: gamma = 0 gives the same value for every u",
      near(run_alpha(tab_2, c(0, 1), sc_2, 3, TRUE,
                     matrix(c(0, 0,  2, 3,  4, 0), ncol = 2, byrow = TRUE), 0),
           rep(66 / 252, 3)))

# Purpose:  with u_s = 1 exactly for the subjects with outcome 2, the p-value
#           follows Fisher's noncentral hypergeometric distribution
#           Please see the Corollary to Theorem 3 
# Input:    tab_2, T = t_11, threshold 3, delta = (0, 1), u = (0, 6), gamma = 1.2
# Expected: P(t_22 >= 4) = sum of dFNCHypergeo(x, 6, 4, 5, exp(1.2)) for x = 4, 5
check("alpha_bd: matches Fisher's noncentral hypergeometric",
      near(run_alpha(tab_2, c(0, 1), sc_2, 3, TRUE, matrix(c(0, 6), nrow = 1), 1.2),
           sum(sapply(4:5, function(x) BiasedUrn::dFNCHypergeo(x, 6, 4, 5, exp(1.2)))),
           tol = 1e-6))

# Small 3 x 3 table and every table with its margins.
tab_A  <- matrix(c(2, 0, 0,
                   1, 1, 1,
                   0, 1, 1), nrow = 3, byrow = TRUE)
cm_A   <- colSums(tab_A)                            # (3, 2, 2)
d_A    <- c(0, 1, 1)
sc_A   <- outer(0:2, 0:2) * 1                       # linear scores
tobs_A <- sum(sc_A * tab_A)                         # 9
tabs_A <- list()
for (a in 0:2) for (b in 0:(2 - a)) for (c1 in 0:3) for (d in 0:(3 - c1)) {
  r1 <- c(a, b, 2 - a - b)
  r2 <- c(c1, d, 3 - c1 - d)
  r3 <- cm_A - r1 - r2
  if (any(r3 < 0)) next
  tabs_A[[length(tabs_A) + 1]] <- unname(rbind(r1, r2, r3))
}
in_A <- vapply(tabs_A, function(t) sum(sc_A * t) >= tobs_A, logical(1))

# Purpose:  alpha(u) equals the sum of point probabilities over Omega_t
# Input:    tab_A, linear scores, threshold 9, gamma = 0.7, five candidate u
# Expected: sum of point_prob_bd over the tables with T >= 9, for each u
U_A <- matrix(c(0, 0, 0,
                1, 0, 0,
                1, 2, 2,
                3, 2, 2,
                2, 1, 1), ncol = 3, byrow = TRUE)
want_A <- apply(U_A, 1, function(u)
  sum(vapply(tabs_A[in_A], function(t) point_prob_bd(t, d_A, u, 0.7), numeric(1))))
check("alpha_bd: equals the sum of point_prob_bd over Omega_t",
      near(run_alpha(tab_A, d_A, sc_A, tobs_A, TRUE, U_A, 0.7), want_A))

# Purpose:  when Omega_t contains every table, the p-value is 1
# Input:    tab_A, linear scores, threshold -1, gamma = 1.2, five candidate u
# Expected: 1 for every row
check("alpha_bd: equals 1 when Omega_t contains every table",
      near(run_alpha(tab_A, d_A, sc_A, -1, TRUE, U_A, 1.2), rep(1, 5)))

# Purpose:  the upper and lower tails add to 1 for an integer-valued statistic
# Input:    tab_A, linear scores, P(T >= 9) and P(T <= 8), u = (1, 2, 2),
#           gamma = 1.3
# Expected: 1
u_A <- matrix(c(1, 2, 2), nrow = 1)
check("alpha_bd: upper and lower tails add to 1",
      near(run_alpha(tab_A, d_A, sc_A, tobs_A,     TRUE,  u_A, 1.3) +
             run_alpha(tab_A, d_A, sc_A, tobs_A - 1, FALSE, u_A, 1.3), 1))

# Purpose:  the row u = 0 gives the randomization p-value for any gamma
# Input:    tab_A, linear scores, threshold 9, u = (0, 0, 0), gamma = 0 and 2
# Expected: equal values
u_0 <- matrix(c(0, 0, 0), nrow = 1)
check("alpha_bd: u = 0 does not depend on gamma",
      near(run_alpha(tab_A, d_A, sc_A, tobs_A, TRUE, u_0, 2),
           run_alpha(tab_A, d_A, sc_A, tobs_A, TRUE, u_0, 0)))





# =============================================================================
# build_W_cpp(...): stage 1, enumerates Omega_t and returns W(t-bar).
# =============================================================================

key_of <- function(x) paste(as.vector(t(x)), collapse = ",")   # cells row by row
keys_A <- vapply(tabs_A, key_of, character(1))
rm_A   <- rowSums(tab_A)

# Purpose:  counts, t-bar and W worked out by hand for a 2 x 2 table
# Input:    tab_2, T = t_11, threshold 3, delta = (0, 1)
# Expected: n_visited = 5, n_tables = 2, TDmat rows (0, 5) and (1, 4), logW = (0, 0)
W_2 <- build_W_cpp(rowSums(tab_2), colSums(tab_2), c(0, 1), 3, TRUE, sc_2)
check("build_W_cpp: 2 x 2 hand calculation",
      W_2$n_visited == 5 && W_2$n_tables == 2 &&
        identical(unname(W_2$TDmat) * 1, matrix(c(0, 5, 1, 4), 2, byrow = TRUE)) &&
        near(W_2$logW, c(0, 0)))

# Purpose:  every table with the given margins is visited exactly once
# Input:    tab_A margins, threshold -1 (every table kept), keep_tables = TRUE
# Expected: n_visited = n_tables = number of tables listed in R; no duplicates;
#           same set of tables
W_all <- build_W_cpp(rm_A, cm_A, d_A, -1, TRUE, sc_A, keep_tables = TRUE)
kept_all <- apply(W_all$tables, 1, paste, collapse = ",")
check("build_W_cpp: visits each table exactly once",
      W_all$n_visited == length(tabs_A) && W_all$n_tables == length(tabs_A) &&
        !any(duplicated(kept_all)) && setequal(kept_all, keys_A))

# Purpose:  the kept tables are exactly those with T(t) >= c
# Input:    tab_A, linear scores, threshold 9, keep_tables = TRUE
# Expected: same set as the tables listed in R with T >= 9
W_up <- build_W_cpp(rm_A, cm_A, d_A, tobs_A, TRUE, sc_A, keep_tables = TRUE)
check("build_W_cpp: Omega_t matches the tables with T >= c",
      W_up$n_tables == sum(in_A) &&
        setequal(apply(W_up$tables, 1, paste, collapse = ","), keys_A[in_A]))

# Purpose:  W(t-bar) equals the sum of kappa(t), from log_c_t, over the kept
#           tables with that t-bar
# Input:    tab_A, linear scores, threshold 9, keep_tables = TRUE
# Expected: same logW for every row of TDmat
kept_up <- lapply(seq_len(nrow(W_up$tables)), function(k)
  matrix(W_up$tables[k, ], nrow = 3, byrow = TRUE))
tbar_up <- t(vapply(kept_up, function(x) colSums(x[d_A == 1, , drop = FALSE]), numeric(3)))
logk_up <- vapply(kept_up, function(x) log_c_t(x, d_A), numeric(1))
want_W  <- apply(W_up$TDmat, 1, function(tb) {
  same <- apply(tbar_up, 1, function(r) all(r == tb))
  log(sum(exp(logk_up[same])))
})
check("build_W_cpp: logW equals the sum of kappa(t) by t-bar",
      near(W_up$logW, want_W))

# Purpose:  the lower tail keeps T(t) <= c, and the two tails split the tables
# Input:    tab_A, linear scores, T >= 9 and T <= 8
# Expected: n_tables(T <= 8) = number of tables with T <= 8 listed in R;
#           n_tables(T >= 9) + n_tables(T <= 8) = n_visited
W_lo <- build_W_cpp(rm_A, cm_A, d_A, tobs_A - 1, FALSE, sc_A)
check("build_W_cpp: lower tail and partition",
      W_lo$n_tables == sum(!in_A) &&
        W_up$n_tables + W_lo$n_tables == W_up$n_visited)

# Purpose:  with every table kept, W recovers the total number of assignments
# Input:    tab_A margins, threshold -1
# Expected: log of sum_r W(t-bar_r) prod_j C(N_.j, t-bar_rj) equals
#           lfactorial(N) - sum(lfactorial(N_i.))
lhs_all <- log(sum(exp(W_all$logW +
                         apply(W_all$TDmat, 1, function(tb) sum(lchoose(cm_A, tb))))))
check("build_W_cpp: W recovers the total number of assignments",
      near(lhs_all, lfactorial(sum(cm_A)) - sum(lfactorial(rm_A))))

# Purpose:  an R function gives the same result as the same score matrix
# Input:    tab_A, stat_fun = function(t) sum(sc_A * t), threshold 9
# Expected: same TDmat, logW and n_tables as with score = sc_A
W_fun <- build_W_cpp(rm_A, cm_A, d_A, tobs_A, TRUE, NULL, function(t) sum(sc_A * t))
check("build_W_cpp: stat_fun matches score",
      identical(W_fun$TDmat, W_up$TDmat) && near(W_fun$logW, W_up$logW) &&
        W_fun$n_tables == W_up$n_tables)

# Purpose:  a supplied score is used even when statistic is also set
# Input:    tab_A, score = sc_A, statistic = 2, threshold 9
# Expected: same TDmat, logW and n_tables as with score = sc_A alone
W_both <- build_W_cpp(rm_A, cm_A, d_A, tobs_A, TRUE, sc_A, statistic = 2)
check("build_W_cpp: score takes priority over statistic",
      identical(W_both$TDmat, W_up$TDmat) && near(W_both$logW, W_up$logW) &&
        W_both$n_tables == W_up$n_tables)

E_A <- outer(rm_A, cm_A) / sum(cm_A)

# Purpose:  built-in chi-square keeps the tables with full chi-square >= observed
# Input:    tab_A, statistic = 2, threshold sum(tab_A^2 / E)
# Expected: kept tables equal those listed in R with sum((t - E)^2 / E) >= observed
chi  <- function(x) sum((x - E_A)^2 / E_A)
in_c <- vapply(tabs_A, function(x) chi(x) >= chi(tab_A) - 1e-9, logical(1))
W_c  <- build_W_cpp(rm_A, cm_A, d_A, sum(tab_A^2 / E_A), TRUE,
                    statistic = 2, keep_tables = TRUE)
check("build_W_cpp: chi-square matches the full statistic",
      setequal(apply(W_c$tables, 1, paste, collapse = ","), keys_A[in_c]))

# Purpose:  built-in G^2 keeps the tables with full G^2 >= observed
# Input:    tab_A, statistic = 3, threshold sum(t log t) of tab_A
# Expected: kept tables equal those listed in R with 2 sum(t log(t / E)) >= observed
g2   <- function(x) 2 * sum(ifelse(x > 0, x * log(x / E_A), 0))
xlx  <- function(x) sum(ifelse(x > 0, x * log(x), 0))
in_g <- vapply(tabs_A, function(x) g2(x) >= g2(tab_A) - 1e-9, logical(1))
W_g  <- build_W_cpp(rm_A, cm_A, d_A, xlx(tab_A), TRUE,
                    statistic = 3, keep_tables = TRUE)
check("build_W_cpp: G^2 matches the full statistic",
      setequal(apply(W_g$tables, 1, paste, collapse = ","), keys_A[in_g]))

# Purpose:  built-in Freeman-Halton keeps the tables no more probable than the
#           observed table under random assignment
# Input:    tab_A, statistic = 4, threshold sum(lfactorial(tab_A))
# Expected: kept tables equal those listed in R with dtable0(t) <= dtable0(tab_A)
p_obs <- dtable0(tab_A)
in_f  <- vapply(tabs_A, function(x) dtable0(x) <= p_obs * (1 + 1e-9), logical(1))
W_f   <- build_W_cpp(rm_A, cm_A, d_A, sum(lfactorial(tab_A)), TRUE,
                     statistic = 4, keep_tables = TRUE)
check("build_W_cpp: Freeman-Halton matches the null probability ordering",
      setequal(apply(W_f$tables, 1, paste, collapse = ","), keys_A[in_f]))

# Purpose:  invalid inputs stop with an error
# Input:    delta of length 2; margins with different sums; no statistic;
#           score of the wrong size
# Expected: an error in each case
check("build_W_cpp: invalid inputs give errors",
      is_error(build_W_cpp(rm_A, cm_A, c(0, 1), tobs_A, TRUE, sc_A)) &&
        is_error(build_W_cpp(rm_A, cm_A + 1, d_A, tobs_A, TRUE, sc_A)) &&
        is_error(build_W_cpp(rm_A, cm_A, d_A, tobs_A, TRUE)) &&
        is_error(build_W_cpp(rm_A, cm_A, d_A, tobs_A, TRUE, matrix(0, 2, 2))))




# =============================================================================
# build_W_from_tables(tabs, I, J, delta): stage 1 from a given list of tables.
# =============================================================================

# Purpose:  same W as build_W_cpp when given the tables build_W_cpp kept
# Input:    tables of Omega_t for tab_A (linear scores, threshold 9), delta = (0, 1, 1)
# Expected: same TDmat, logW and n_tables as build_W_cpp
W_tab <- build_W_from_tables(W_up$tables, 3, 3, d_A)
check("build_W_from_tables: matches build_W_cpp",
      identical(W_tab$TDmat, W_up$TDmat) && near(W_tab$logW, W_up$logW) &&
        W_tab$n_tables == W_up$n_tables)



# =============================================================================
# Factorial(nValue): n! by direct multiplication, used by choose_C.
# We note that though choose_C is only used for legacy library
# =============================================================================

# Purpose:  small values match hand calculation
# Input:    n = 2, 3, 5, 10
# Expected: 2, 6, 120, 3628800
check("Factorial: small values",
      near(c(Factorial_test_(2), Factorial_test_(3),
             Factorial_test_(5), Factorial_test_(10)),
           c(2, 6, 120, 3628800)))

# Purpose:  agrees with base R for every n from 2 to 170
# Input:    n = 2, ..., 170
# Expected: factorial(n), within relative tolerance 1e-10
check("Factorial: matches base R for 2 <= n <= 170",
      near(vapply(2:170, Factorial_test_, numeric(1)), factorial(2:170)))

# Purpose:  n = 0 and n = 1 give 0, as stated in the comment
# Input:    n = 0, 1
# Expected: 0, 0
check("Factorial: n = 0 and n = 1 give 0",
      identical(c(Factorial_test_(0), Factorial_test_(1)), c(0, 0)))

# Purpose:  the result overflows above 170, as stated in the comment
# Input:    n = 171
# Expected: Inf
check("Factorial: overflows for n > 170",
      identical(Factorial_test_(171), Inf))


# =============================================================================
# choose_C(n, k): binomial coefficient C(n, k) in double precision.
# =============================================================================

# Purpose:  negative arguments give 0
# Input:    (n, k) = (-1, 0), (3, -1)
# Expected: 0, 0
check("choose_C: negative arguments give 0",
      identical(c(choose_C(-1, 0), choose_C(3, -1)), c(0, 0)))

# Purpose:  k > n gives 0
# Input:    (n, k) = (3, 5), (0, 1)
# Expected: 0, 0
check("choose_C: k > n gives 0",
      identical(c(choose_C(3, 5), choose_C(0, 1)), c(0, 0)))

# Purpose:  k = 0 and k = n give 1
# Input:    (n, k) = (0, 0), (7, 0), (7, 7)
# Expected: 1, 1, 1
check("choose_C: k = 0 and k = n give 1",
      identical(c(choose_C(0, 0), choose_C(7, 0), choose_C(7, 7)), c(1, 1, 1)))

# Purpose:  k = 1 and k = n - 1 give n
# Input:    (n, k) = (7, 1), (7, 6), (2, 1)
# Expected: 7, 7, 2
check("choose_C: k = 1 and k = n - 1 give n",
      identical(c(choose_C(7, 1), choose_C(7, 6), choose_C(2, 1)), c(7, 7, 2)))

# Purpose:  general case matches hand calculation
# Input:    (n, k) = (10, 3), (6, 3)
# Expected: 120, 20
check("choose_C: general case, hand calculation",
      near(c(choose_C(10, 3), choose_C(6, 3)), c(120, 20)))

# Purpose:  agrees with base R for every (n, k) with 0 <= k <= n <= 170
# Input:    all such (n, k)
# Expected: choose(n, k), relative error below 1e-10 for every pair
grid_C <- do.call(rbind, lapply(0:170, function(n) cbind(n = n, k = 0:n)))
got_C  <- mapply(choose_C, grid_C[, "n"], grid_C[, "k"])
want_C <- choose(grid_C[, "n"], grid_C[, "k"])
check("choose_C: matches base R for n <= 170",
      max(abs(got_C / want_C - 1)) < 1e-10)

# Purpose:  the special cases still work for large n
# Input:    (n, k) = (500, 0), (500, 1), (500, 499), (500, 500)
# Expected: 1, 500, 500, 1
check("choose_C: special cases work for large n",
      identical(c(choose_C(500, 0), choose_C(500, 1),
                  choose_C(500, 499), choose_C(500, 500)),
                c(1, 500, 500, 1)))

# Purpose:  the general case is not finite for n > 170, as stated in the comment
# Input:    (n, k) = (171, 2), (200, 100)
# Expected: not finite
check("choose_C: general case is not finite for n > 170",
      !is.finite(choose_C(171, 2)) && !is.finite(choose_C(200, 100)))




### A simple way to check whether there is any failure in the checks

cat("\n", sum(results), " / ", length(results), " checks passed\n", sep = "")
if (any(!results)) cat("Failed:\n", paste0("  ", names(results)[!results], collapse = "\n"), "\n", sep = "")



# =============================================================================
# sensitivityIxJ_R_check.R
#
# Code for "Exact, Nonparametric Sensitivity Analysis for Observational Studies
# of Contingency Tables"
#
# Unit tests for library/sensitivityIxJ.R.
# Each function is checked against special cases with known answers, against
# values computed directly in R, and against the other functions of the
# interface.
#
# Run from the top-level folder (the one with README.txt):
#   source("development_check/sensitivityIxJ_R_check.R")
# Each check prints PASS or FAIL; a summary is printed at the end.
# =============================================================================
library(Rcpp)
# ---- Setup ------------------------------------------------------------------
sourceCpp("library/sensitivityIxJ.cpp")
sourceCpp("library/SIS_methods.cpp")
source("library/sensitivityIxJ.R")
# Record and print the result of one check.
results <- logical(0)
check <- function(label, expr) {
  ok <- tryCatch(isTRUE(expr), error = function(e) FALSE)
  results[label] <<- ok
  cat(if (ok) "PASS  " else "FAIL  ", label, "\n", sep = "")
}
# TRUE if x and y agree within numerical tolerance.
near <- function(x, y, tol = 1e-10) isTRUE(all.equal(x, y, tolerance = tol))
# TRUE if the expression raises an error.
is_error <- function(expr) inherits(try(expr, silent = TRUE), "try-error")
# 3 x 3 table used by several sections (N_i. = (35, 24, 27), N_.j = (24, 36, 26)).
tab_3 <- matrix(c(12, 18,  5,
                  6, 12,  6,
                  6,  6, 15), nrow = 3, byrow = TRUE)
# =============================================================================
# .as_treatment_rows(tb, row): put the table in treatment-rows orientation.
# =============================================================================
# Purpose:  a treatment-rows table is returned unchanged
# Input:    tab_3, row = "treatment"
# Expected: the same cell counts, same dimensions
check(".as_treatment_rows: treatment rows unchanged",
      identical(unname(.as_treatment_rows(tab_3, "treatment")) * 1, tab_3))
# Purpose:  an outcome-rows table is transposed, so both orientations agree
# Input:    tab_3 with row = "treatment", and its transpose with row = "outcome"
# Expected: identical results
check(".as_treatment_rows: outcome rows are transposed",
      identical(.as_treatment_rows(tab_3, "treatment"),
                .as_treatment_rows(t(tab_3), "outcome")))
# Purpose:  the result is an integer matrix, whatever the input type
# Input:    a double matrix and a table object
# Expected: storage mode "integer" in both cases
check(".as_treatment_rows: result is an integer matrix",
      identical(storage.mode(.as_treatment_rows(tab_3, "treatment")), "integer") &&
        identical(storage.mode(.as_treatment_rows(as.table(tab_3), "treatment")),
                  "integer"))
# =============================================================================
# .check_common(tb, gamma, delta, row): validate the arguments shared by the
#                                       exact routines.
# =============================================================================
# Purpose:  a valid call returns delta as an integer vector
# Input:    tab_3, gamma = 1, delta = (0, 1, 1), row = "treatment"
# Expected: the integer vector 0 1 1
check(".check_common: returns delta as an integer vector",
      identical(.check_common(tab_3, 1, c(0, 1, 1), "treatment"), c(0L, 1L, 1L)))
# Purpose:  gamma = 0 is allowed, and so is a delta with a single 1
# Input:    gamma = 0 with delta = (0, 1, 1); gamma = 2 with delta = (0, 0, 1)
# Expected: no error in either case
check(".check_common: accepts gamma = 0 and a single nonzero delta",
      !is_error(.check_common(tab_3, 0, c(0, 1, 1), "treatment")) &&
        !is_error(.check_common(tab_3, 2, c(0, 0, 1), "treatment")))
# Purpose:  an invalid row label is rejected
# Input:    row = "Treatment" and row = "outcomes"
# Expected: an error in each case
check(".check_common: invalid row label gives an error",
      is_error(.check_common(tab_3, 1, c(0, 1, 1), "Treatment")) &&
        is_error(.check_common(tab_3, 1, c(0, 1, 1), "outcomes")))
# Purpose:  gamma must be a single nonnegative number
# Input:    gamma = -1; gamma = c(1, 2)
# Expected: an error in each case
check(".check_common: invalid gamma gives an error",
      is_error(.check_common(tab_3, -1, c(0, 1, 1), "treatment")) &&
        is_error(.check_common(tab_3, c(1, 2), c(0, 1, 1), "treatment")))
# Purpose:  delta must be binary and have one entry per treatment level
# Input:    delta = (1, 3, 3); delta = (0, 1); delta = (0, 0.5, 1)
# Expected: an error in each case
check(".check_common: invalid delta gives an error",
      is_error(.check_common(tab_3, 1, c(1, 3, 3), "treatment")) &&
        is_error(.check_common(tab_3, 1, c(0, 1), "treatment")) &&
        is_error(.check_common(tab_3, 1, c(0, 0.5, 1), "treatment")))
# Purpose:  a constant delta is rejected, since the confounder would shift
#           every treatment level equally
# Input:    delta = (0, 0, 0) and delta = (1, 1, 1)
# Expected: an error in each case
check(".check_common: constant delta gives an error",
      is_error(.check_common(tab_3, 1, c(0, 0, 0), "treatment")) &&
        is_error(.check_common(tab_3, 1, c(1, 1, 1), "treatment")))
# =============================================================================
# .check_u_space(U, outcome_margins): validate a matrix of candidate u-bar
#                                     vectors.
# =============================================================================
m_U <- c(24, 36, 26)
# Purpose:  a valid matrix is returned as an integer matrix, values unchanged
# Input:    rows (0, 0, 10) and (5, 5, 5), outcome margins (24, 36, 26)
# Expected: the same values, storage mode "integer"
U_in <- matrix(c(0, 0, 10,
                 5, 5,  5), nrow = 2, byrow = TRUE)
U_out <- .check_u_space(U_in, m_U)
check(".check_u_space: returns an integer matrix, values unchanged",
      identical(storage.mode(U_out), "integer") && identical(U_out * 1, U_in))
# Purpose:  the two extreme candidates are accepted
# Input:    the all-zero row, and the row equal to the outcome margins
# Expected: no error
check(".check_u_space: accepts the extreme candidates",
      !is_error(.check_u_space(matrix(0, 1, 3), m_U)) &&
        !is_error(.check_u_space(matrix(m_U, 1), m_U)))
# Purpose:  the number of columns must equal the number of outcome levels
# Input:    a 2-column matrix against margins of length 3
# Expected: an error
check(".check_u_space: wrong number of columns gives an error",
      is_error(.check_u_space(matrix(c(0, 0), 1), m_U)))
# Purpose:  an entry above its outcome margin is rejected
# Input:    row (0, 0, 30), where 30 exceeds N_.3 = 26
# Expected: an error
check(".check_u_space: entry above the outcome margin gives an error",
      is_error(.check_u_space(matrix(c(0, 0, 30), 1), m_U)))
# Purpose:  a negative entry is rejected
# Input:    row (0, -1, 0)
# Expected: an error
check(".check_u_space: negative entry gives an error",
      is_error(.check_u_space(matrix(c(0, -1, 0), 1), m_U)))
# Purpose:  the input must be a numeric matrix
# Input:    a numeric vector; a character matrix
# Expected: an error in each case
check(".check_u_space: non-matrix or non-numeric input gives an error",
      is_error(.check_u_space(c(0, 0, 10), m_U)) &&
        is_error(.check_u_space(matrix("0", 1, 3), m_U)))
# Purpose:  an infeasible row is caught wherever it sits in the matrix
# Input:    three rows, with the violation in the middle row
# Expected: an error
check(".check_u_space: infeasible row in the middle is caught",
      is_error(.check_u_space(matrix(c(0, 0,  0,
                                       0, 0, 30,
                                       1, 1,  1), nrow = 3, byrow = TRUE), m_U)))
# =============================================================================
# full.candidate.set(outcome_margins): every candidate u-bar vector.
# =============================================================================
# Purpose:  the set has prod_j (N_.j + 1) rows
# Input:    outcome_margins = (24, 36, 26), (4, 4), (3, 2, 2, 1)
# Expected: 24975, 25, 72
check("full.candidate.set: size is prod(N_.j + 1)",
      all(vapply(list(c(24, 36, 26), c(4, 4), c(3, 2, 2, 1)),
                 function(m) nrow(full.candidate.set(m)) == prod(m + 1),
                 logical(1))))
# =============================================================================
# full.candidate.set.rows(outcome_margins, k): selected rows of the full
#                                              candidate set.
# =============================================================================
# Purpose:  agrees with full.candidate.set row by row
# Input:    outcome_margins = (4, 3, 2), every row index in order
# Expected: the same matrix
m_r <- c(4, 3, 2)
check("full.candidate.set.rows: matches full.candidate.set",
      identical(unname(full.candidate.set.rows(m_r, seq_len(prod(m_r + 1)))),
                unname(full.candidate.set(m_r))))
# Purpose:  scattered indices, including the first and last, give the right rows
# Input:    outcome_margins = (24, 36, 26), k = 1, 26, 1000, 12345, 24975
# Expected: the same rows as full.candidate.set at those indices
m_s <- c(24, 36, 26)
k_s <- c(1, 26, 1000, 12345, prod(m_s + 1))
check("full.candidate.set.rows: scattered indices",
      identical(unname(full.candidate.set.rows(m_s, k_s)),
                unname(full.candidate.set(m_s)[k_s, ])))
# =============================================================================
# .u_space_ordinal(outcome_margins): candidate set for an ordinal test, filling
#                                    the highest outcome level first.
# =============================================================================
# Purpose:  contents worked out by hand
# Input:    outcome_margins = (3, 2, 2), N = 7
# Expected: 8 rows; rows 1, 3, 5, 8 are (0,0,0), (0,0,2), (0,2,2), (3,2,2)
U_o <- .u_space_ordinal(c(3, 2, 2))
check(".u_space_ordinal: hand calculation",
      nrow(U_o) == 8 &&
        all(U_o[1, ] == c(0, 0, 0)) && all(U_o[3, ] == c(0, 0, 2)) &&
        all(U_o[5, ] == c(0, 2, 2)) && all(U_o[8, ] == c(3, 2, 2)))
# Purpose:  N + 1 rows, row k has k subjects with u_s = 1, and a level is used
#           only once every higher level is full
# Input:    outcome_margins = (10, 4, 6, 5)
# Expected: nrow = N + 1; row sums 0, 1, ..., N; in every row, u-bar_j > 0
#           implies u-bar_k = N_.k for all k > j
m_o <- c(10, 4, 6, 5)
U_p <- .u_space_ordinal(m_o)
check(".u_space_ordinal: N + 1 rows filling the highest level first",
      nrow(U_p) == sum(m_o) + 1 &&
        identical(as.integer(rowSums(U_p)), 0:sum(m_o)) &&
        all(apply(U_p, 1, function(u)
          all(vapply(seq_len(length(m_o) - 1), function(j)
            u[j] == 0 || all(u[(j + 1):length(m_o)] == m_o[(j + 1):length(m_o)]),
            logical(1))))))
# =============================================================================
# .u_space_sign_score(outcome_margins): candidate matrix for a binary outcome.
# =============================================================================
# Purpose:  the two rows are the reference point and the closed-form worst case
# Input:    outcome_margins = (15, 75)
# Expected: a 2 x 2 integer matrix with rows (0, 0) and (0, 75)
U_s <- .u_space_sign_score(c(15, 75))
check(".u_space_sign_score: hand calculation",
      identical(unname(U_s) * 1, matrix(c(0, 0, 0, 75), ncol = 2, byrow = TRUE)) &&
        identical(storage.mode(U_s), "integer"))
# Purpose:  both rows are feasible, and the first outcome level is never used
# Input:    outcome_margins = (4, 6) and (100, 1)
# Expected: 2 rows within [0, N_.j]; column 1 is zero in both rows
check(".u_space_sign_score: feasible, first level unused",
      all(vapply(list(c(4, 6), c(100, 1)), function(m) {
        U <- .u_space_sign_score(m)
        nrow(U) == 2 && all(U >= 0) && all(sweep(U, 2, m, "<=")) &&
          all(U[, 1] == 0)
      }, logical(1))))
# =============================================================================
# .prep_from_W(res, n, m, d, J): collect the stage-1 output and the constants
#                                stage 2 needs.
# =============================================================================
n_P   <- rowSums(tab_3)                      # 35, 24, 27
m_P   <- colSums(tab_3)                      # 24, 36, 26
d_P   <- c(0, 1, 1)
sc_P  <- outer(c(0, 1, 2), c(0, 1, 2))
res_P <- build_W_cpp(n_P, m_P, d_P, sum(sc_P * tab_3), TRUE, sc_P, NULL)
# Purpose:  the stage-1 fields are passed through unchanged, and the derived
#           fields are N, N-bar and log kappa_0
# Input:    res_P with n = (35, 24, 27), m = (24, 36, 26), delta = (0, 1, 1)
# Expected: TDmat, logW, n_tables as in res_P; N = 86; ND = 24 + 27 = 51;
#           logc0 = log_c_0(n, delta)
prep_P <- .prep_from_W(res_P, n_P, m_P, d_P, 3)
check(".prep_from_W: fields are passed through and derived correctly",
      identical(prep_P$TDmat, res_P$TDmat) &&
        identical(prep_P$logW, res_P$logW) &&
        identical(prep_P$n_tables, res_P$n_tables) &&
        all(prep_P$m == m_P) && prep_P$J == 3 &&
        prep_P$N == 86 && prep_P$ND == 51 &&
        near(prep_P$logc0, log_c_0(n_P, d_P)))
# Purpose:  ND sums the treatment margins with delta_i = 1, so delta and
#           1 - delta give values of ND that add to N
# Input:    the same margins with delta = (0, 1, 1), (1, 0, 0) and (0, 0, 1)
# Expected: ND = 51, 35 and 27; the first two add to 86
nd_of <- function(d) .prep_from_W(res_P, n_P, m_P, d, 3)$ND
check(".prep_from_W: ND sums the treatment margins where delta is 1",
      nd_of(c(0, 1, 1)) == 51 && nd_of(c(1, 0, 0)) == 35 &&
        nd_of(c(0, 0, 1)) == 27 &&
        nd_of(c(0, 1, 1)) + nd_of(c(1, 0, 0)) == 86)
# =============================================================================
# .rct(prep, gamma): randomization p-value, evaluated at u = 0.
# =============================================================================
# Purpose:  at u = 0 the p-value does not depend on gamma
# Input:    prep_P from .prep_from_W; gamma = 0 and gamma = 2
# Expected: equal values
check(".rct: invariant to gamma at u = 0",
      near(.rct(prep_P, gamma = 0), .rct(prep_P, gamma = 2)))
# =============================================================================
# possible.table(threshold, table, direction, transform.fun): enumerate Omega_t.
# =============================================================================
tb_PT   <- matrix(c(3, 2, 1,
                    0, 2, 4,
                    0, 1, 5), nrow = 3, byrow = TRUE)
stat_PT <- function(t) sum(outer(0:2, 0:2) * t)
# Purpose:  every table in Omega_t keeps the margins of the observed table,
#           meets the statistic threshold, and the observed table is included
# Input:    tb_PT, threshold = stat_PT(tb_PT), direction = "greater than"
# Expected: Omega_t is nonempty; every table has the margins of tb_PT and
#           stat_PT(t) >= threshold; tb_PT itself appears in Omega_t
Omega_PT <- possible.table(threshold = stat_PT(tb_PT), table = tb_PT,
                           direction = "greater than", transform.fun = stat_PT)
check("possible.table: margins and threshold are respected, table included",
      length(Omega_PT) > 0 &&
        all(sapply(Omega_PT, function(t)
          all(rowSums(t) == rowSums(tb_PT)) && all(colSums(t) == colSums(tb_PT)))) &&
        all(sapply(Omega_PT, stat_PT) >= stat_PT(tb_PT)) &&
        any(sapply(Omega_PT, function(t) all(as.vector(unclass(t)) == as.vector(tb_PT)))))
# Purpose:  a threshold outside the range of transform.fun gives an empty list
# Input:    tb_PT, threshold = -1, direction = "less than" (stat_PT is never
#           negative, since all scores and counts are nonnegative)
# Expected: an empty list
check("possible.table: unreachable threshold gives an empty list",
      length(possible.table(threshold = -1, table = tb_PT,
                            direction = "less than",
                            transform.fun = stat_PT)) == 0)
# =============================================================================
# generic.I.by.J.sensitivity.point.probability(table, u_allocation, gamma, delta):
#                                       exact point probability of one table.
# =============================================================================
# Multivariate hypergeometric probability of a table, computed with lgamma.
hyper_prob <- function(tb)
  exp(sum(lgamma(rowSums(tb) + 1)) + sum(lgamma(colSums(tb) + 1)) -
        lgamma(sum(tb) + 1) - sum(lgamma(tb + 1)))
# Purpose:  at gamma = 0, the point probability reduces to the multivariate
#           hypergeometric probability of the table, for any u_allocation
# Input:    tab_3, u_allocation = (0, 0, 0) and (5, 10, 8), gamma = 0,
#           delta = (0, 1, 1)
# Expected: both equal hyper_prob(tab_3)
check("generic.I.by.J.sensitivity.point.probability: gamma = 0 gives the hypergeometric probability",
      near(generic.I.by.J.sensitivity.point.probability(tab_3, u_allocation = c(0, 0, 0),
                                                        gamma = 0, delta = c(0, 1, 1)),
           hyper_prob(tab_3)) &&
        near(generic.I.by.J.sensitivity.point.probability(tab_3, u_allocation = c(5, 10, 8),
                                                          gamma = 0, delta = c(0, 1, 1)),
             hyper_prob(tab_3)))
# Purpose:  for a fixed gamma and u_allocation, the point probabilities form a
#           probability distribution over every table with the observed margins
# Input:    all tables with the margins of tb_SM (N_i. = (3, 3), N_.j = (2, 4));
#           gamma = 1, u_allocation = (1, 3), delta = (0, 1)
# Expected: the point probabilities sum to 1
tb_SM <- matrix(c(2, 1,
                  0, 3), nrow = 2, byrow = TRUE)
all_tabs_SM <- possible.table(threshold = -1, table = tb_SM,
                              direction = "greater than", transform.fun = function(t) 0)
probs_SM <- sapply(all_tabs_SM, function(t)
  generic.I.by.J.sensitivity.point.probability(t, u_allocation = c(1, 3),
                                               gamma = 1, delta = c(0, 1)))
check("generic.I.by.J.sensitivity.point.probability: sums to 1 over all tables with the margins",
      near(sum(probs_SM), 1))
# Purpose:  at u = 0, the point probability of a 2 x 2 table equals the
#           hypergeometric probability of its N_11 cell, for any gamma
# Input:    every table with the margins of tb2; gamma = 1 (no effect when
#           u_allocation = 0); delta = (0, 1)
# Expected: the point probability at u_allocation = (0, 0) equals
#           dhyper(t[1, 1], m1, m2, n1) for every such t
tb2 <- matrix(c(3, 2,
                1, 4), nrow = 2, byrow = TRUE)
m1 <- sum(tb2[, 1]); m2 <- sum(tb2[, 2]); n1 <- sum(tb2[1, ])   # column margins, first row margin
all_tabs_2x2 <- possible.table(threshold = -1, table = tb2,
                               direction = "greater than", transform.fun = function(t) 0)
check("generic.I.by.J.sensitivity.point.probability: matches dhyper for a 2 x 2 table",
      all(sapply(all_tabs_2x2, function(t)
        near(generic.I.by.J.sensitivity.point.probability(t, u_allocation = c(0, 0),
                                                          gamma = 1, delta = c(0, 1)),
             dhyper(t[1, 1], m1, m2, n1)))))
# =============================================================================
# exact.general.sen.IxJ(obs.table, gamma, delta, ...): worst-case p-value for a
#                                                      general statistic.
# =============================================================================
# Purpose:  the randomization p-value (u = 0) of an upper-tail test on N_11
#           equals the upper-tail hypergeometric p-value
# Input:    tb2, transform.fun = N_11 cell count, gamma = 1 (rct.prob is
#           evaluated at u = 0 and does not depend on gamma)
# Expected: rct.prob equals 1 - phyper(observed N_11 - 1, m1, m2, n1)
res_rct <- exact.general.sen.IxJ(obs.table = tb2, gamma = 1, delta = c(0, 1),
                                 transform.fun = function(t) t[1, 1])
check("exact.general.sen.IxJ: rct.prob matches the phyper upper tail",
      near(res_rct$rct.prob, 1 - phyper(tb2[1, 1] - 1, m1, m2, n1)))
# =============================================================================
# exact.score.sen.IxJ(obs.table, gamma, delta, treatment.scores, outcome.scores):
#                     worst-case p-value for the ordinal test.
# =============================================================================
# Purpose:  for a binary outcome (J = 2), only two allocations are evaluated,
#           and the maximizer is the closed-form u+ = (0, N_.2)
# Input:    a 3 x 2 table, gamma = 1, delta = (0, 1, 1), nondecreasing scores
# Expected: n.candidates equals 2; maximizer equals (0, N_.2)
tb_ES2 <- matrix(c(3, 2,
                   1, 4,
                   2, 2), nrow = 3, byrow = TRUE)
r_ES2 <- exact.score.sen.IxJ(obs.table = tb_ES2, gamma = 1, delta = c(0, 1, 1),
                             treatment.scores = c(0, 1, 2), outcome.scores = c(0, 1))
check("exact.score.sen.IxJ: J = 2 evaluates exactly 2 candidates at u+",
      r_ES2$n.candidates == 2 && length(r_ES2$maximizer) == 2 &&
        all(r_ES2$maximizer == c(0, sum(tb_ES2[, 2]))))
# Purpose:  for J >= 3, the small ordinal candidate set (N + 1 allocations)
#           gives the same worst-case p-value as searching the full candidate
#           set with the same ordinal statistic via exact.general.sen.IxJ
# Input:    a 3 x 3 table, gamma = 1, delta = (0, 1, 1), nondecreasing scores
# Expected: rct.prob and max.prob agree between the two functions;
#           exact.score.sen.IxJ evaluates N + 1 candidates, fewer than the
#           full candidate set
tb_ES3 <- matrix(c(2, 1, 0,
                   1, 2, 1,
                   0, 1, 2), nrow = 3, byrow = TRUE)
r_small <- exact.score.sen.IxJ(obs.table = tb_ES3, gamma = 1, delta = c(0, 1, 1),
                               treatment.scores = c(0, 1, 2), outcome.scores = c(0, 1, 2))
r_full  <- exact.general.sen.IxJ(obs.table = tb_ES3, gamma = 1, delta = c(0, 1, 1),
                                 treatment.scores = c(0, 1, 2), outcome.scores = c(0, 1, 2))
check("exact.score.sen.IxJ: small candidate set matches the full search",
      near(r_small$rct.prob, r_full$rct.prob) &&
        near(r_small$max.prob, r_full$max.prob) &&
        r_small$n.candidates == sum(tb_ES3) + 1 &&
        r_small$n.candidates < r_full$n.candidates)
# Purpose:  a nonmonotone delta is rejected, since the small candidate set is
#           then not guaranteed to contain the worst case
# Input:    tb_ES3, delta = (0, 1, 0) (not nondecreasing)
# Expected: an error
check("exact.score.sen.IxJ: nonmonotone delta gives an error",
      is_error(exact.score.sen.IxJ(obs.table = tb_ES3, gamma = 1, delta = c(0, 1, 0),
                                   treatment.scores = c(0, 1, 2), outcome.scores = c(0, 1, 2))))
# =============================================================================
# construct_A_matrix(I, J): map from the free cells of an I x J table to all
#                           cells.
# =============================================================================
# Purpose:  A %*% free cells, plus the remainder that depends only on the
#           margins, reconstructs the observed table exactly
# Input:    tb_A, a 3 x 3 table; its free cells (N_11, N_12, N_21, N_22) and
#           margins
# Expected: as.numeric(t(tb_A)) - A %*% free_cells equals N_1. at cell (1,3),
#           N_2. at cell (2,3), N_.1 at cell (3,1), N_.2 at cell (3,2),
#           N - N_1. - N_2. - N_.1 - N_.2 at cell (3,3), and 0 at the four
#           free cells
tb_A <- matrix(c(2, 3, 1,
                 4, 1, 2,
                 0, 5, 3), nrow = 3, byrow = TRUE)
n_A <- rowSums(tb_A); m_A <- colSums(tb_A); N_A <- sum(tb_A)
free_A <- c(tb_A[1, 1], tb_A[1, 2], tb_A[2, 1], tb_A[2, 2])
A33 <- construct_A_matrix(3, 3)
expected_remainder <- c(0, 0, n_A[1],
                        0, 0, n_A[2],
                        m_A[1], m_A[2], N_A - n_A[1] - n_A[2] - m_A[1] - m_A[2])
check("construct_A_matrix: A %*% free cells plus the margin remainder equals tb_A",
      identical(unname(as.numeric(t(tb_A)) - as.numeric(A33 %*% free_A)),
                unname(expected_remainder)))
# =============================================================================
# create_A_vector(I, J, treatment.scores, outcome.scores): flatten w_i v_j.
# =============================================================================
# Purpose:  create_A_vector agrees with row-major flattening of outer(w, v),
#           for a non-square case
# Input:    I = 3, J = 2, treatment.scores = (0, 1, 2), outcome.scores = (0, 1)
# Expected: identical to as.vector(t(outer(w, v)))
check("create_A_vector: matches row-major flattening of outer(w, v)",
      identical(create_A_vector(3, 2, c(0, 1, 2), c(0, 1)),
                as.vector(t(outer(c(0, 1, 2), c(0, 1))))))
# Purpose:  the dot product with a row-major flattened table equals the
#           ordinal statistic sum_ij w_i v_j N_ij
# Input:    a 3 x 2 table tb_V, treatment.scores = (0, 1, 2), outcome.scores = (0, 1)
# Expected: sum(create_A_vector(...) * as.vector(t(tb_V))) equals
#           sum(outer(w, v) * tb_V)
tb_V <- matrix(c(2, 3,
                 4, 1,
                 0, 5), nrow = 3, byrow = TRUE)
w_V <- c(0, 1, 2); v_V <- c(0, 1)
check("create_A_vector: dot product with a flattened table matches the ordinal statistic",
      near(sum(create_A_vector(3, 2, w_V, v_V) * as.vector(t(tb_V))),
           sum(outer(w_V, v_V) * tb_V)))
# ---- Summary ----------------------------------------------------------------
cat("\n", sum(results), " / ", length(results), " checks passed\n", sep = "")
if (any(!results)) cat("Failed:\n", paste0("  ", names(results)[!results], collapse = "\n"), "\n", sep = "")
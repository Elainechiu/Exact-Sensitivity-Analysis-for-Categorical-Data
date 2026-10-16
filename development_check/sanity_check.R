# =============================================================================
# sanity_check.R
#
# Code for "Exact, Nonparametric Sensitivity Analysis for Observational Studies
# of Contingency Tables"
# Authors: Elaine K. Chiu and Hyunseung Kang
#
# Checks the engine against results known independently of the paper.
# Theorem 4 gives the exact null distribution of the table at a binary u,
#     P_u(N = t) = kappa(t) S(t-bar, u) / (kappa_0 S_0(u)),
# and Theorems 1-3 give the candidate set over which the worst case is taken.
# If the implementation is right, it must reproduce three facts that do not
# depend on the paper:
#   (1) at gamma = 0 the model is random assignment, so P_u(N = t) is the
#       multivariate hypergeometric probability for every u, and every
#       worst-case p-value equals the permutation p-value, which for the
#       probability-ordering statistic is what fisher.test() computes;
#   (2) for every (u, gamma), P_u(.) is a probability distribution on the
#       tables with the observed margins, so the point probabilities sum to 1;
#   (3) the mean and variance of an ordinal statistic computed by brute force
#       from P_u(.) equal the closed-form moments of Theorem D.1 used by the
#       normal approximation.
#
# Run from the top-level folder (the one with README.txt):
#   source("development_check/sanity_check.R")
# Each check prints PASS or FAIL; a summary is printed at the end.
# =============================================================================
library(Rcpp)

# ---- Setup ------------------------------------------------------------------
sourceCpp("library/sensitivityIxJ.cpp")
source("library/sensitivityIxJ.R")
set.seed(20260922)

# Record and print the result of one check.
results <- logical(0)
check <- function(label, expr) {
  ok <- tryCatch(isTRUE(expr), error = function(e) FALSE)
  results[label] <<- ok
  cat(if (ok) "PASS  " else "FAIL  ", label, "\n", sep = "")
}

# TRUE if x and y agree within numerical tolerance.
near <- function(x, y, tol = 1e-10) isTRUE(all.equal(x, y, tolerance = tol))

# A random I x J table with N subjects and no empty row or column.
rand_table <- function(I, J, N) {
  z <- sample.int(I, N, TRUE); r <- sample.int(J, N, TRUE)
  tb <- matrix(0L, I, J)
  for (s in seq_len(N)) tb[z[s], r[s]] <- tb[z[s], r[s]] + 1L
  if (any(rowSums(tb) == 0) || any(colSums(tb) == 0)) return(rand_table(I, J, N))
  tb
}

# Log of the multivariate hypergeometric probability of a table: the
# permutation distribution under random assignment with fixed margins.
log_mvhyper <- function(t)
  sum(lgamma(rowSums(t) + 1)) + sum(lgamma(colSums(t) + 1)) -
  lgamma(sum(t) + 1) - sum(lgamma(t + 1))

# Every table with the margins of tb (threshold 0 on a constant statistic).
all_tables <- function(tb)
  possible.table(threshold = 0, table = tb, direction = "greater than",
                 transform.fun = function(t) 1)

# Point probability of one table at (u-bar, gamma, delta).
point_prob <- function(t, u, gamma, delta)
  generic.I.by.J.sensitivity.point.probability(table = as.matrix(t), u_allocation = u,
                                               gamma = gamma, delta = delta)

# Permutation p-value P(T(N) >= T(tb)) under random assignment.
perm_pvalue <- function(tb, stat) {
  tabs <- all_tables(tb); c_obs <- stat(tb)
  sum(sapply(tabs, function(t)
    if (stat(as.matrix(t)) >= c_obs - 1e-9) exp(log_mvhyper(as.matrix(t))) else 0))
}


# =============================================================================
# (1) gamma = 0: the model collapses to random assignment.
# =============================================================================

# Purpose:  at gamma = 0 the point probability equals the multivariate
#           hypergeometric probability, whatever u-bar is
# Input:    20 random tables in each of 2x2, 2x3, 3x3, 4x2, 3x4, 5x2, 2x5
#           (N between 8 and 30), delta = (0, 1, ..., 1), a random u-bar
# Expected: relative difference below 1e-10 on all 140 tables
worst <- 0; n_cases <- 0
for (d in list(c(2,2), c(2,3), c(3,3), c(4,2), c(3,4), c(5,2), c(2,5))) for (rep in 1:20) {
  tb <- rand_table(d[1], d[2], sample(8:30, 1)); delta <- c(0, rep(1, d[1] - 1))
  u  <- sapply(colSums(tb), function(x) sample(0:x, 1))
  p_engine <- point_prob(tb, u, 0, delta); p_hyper <- exp(log_mvhyper(tb))
  worst <- max(worst, abs(p_engine - p_hyper) / p_hyper); n_cases <- n_cases + 1
}
check(sprintf("gamma = 0: point probability is multivariate hypergeometric (%d random tables, worst rel. diff %.1e)",
              n_cases, worst),
      worst < 1e-10)

# Purpose:  Theorem 1 at gamma = 0: alpha(T, r, u) is the same at every
#           candidate of U_C and equals the permutation p-value
# Input:    a 3x3 table, chi-square statistic, delta = (0, 1, 1), default
#           u_space (the full candidate set U_C)
# Expected: all.prob constant across U_C; max.prob equals the permutation
#           p-value computed by enumeration
tb_1 <- matrix(c(5,3,2, 4,6,3, 2,3,6), 3, byrow = TRUE)
chisq_stat <- function(t) { E <- outer(rowSums(t), colSums(t)) / sum(t); sum((t - E)^2 / E) }
r_1 <- exact.general.sen.IxJ(obs.table = tb_1, gamma = 0, delta = c(0, 1, 1), transform.fun = chisq_stat)
p_perm_1 <- perm_pvalue(tb_1, chisq_stat)
check(sprintf("Theorem 1, chi-square, gamma = 0: alpha constant over %d candidates and equals permutation p-value %.8f",
              r_1$n.candidates, p_perm_1),
      diff(range(r_1$all.prob)) < 1e-12 && near(r_1$max.prob, p_perm_1))

# Purpose:  Theorem 2 at gamma = 0: the maximum over U_O and the p-value under
#           random assignment both equal the permutation p-value
# Input:    a 3x3 table, scores (0, 1, 2) x (0, 1, 2), delta = (0, 1, 1)
# Expected: max.prob == rct.prob == permutation p-value
tb_2 <- matrix(c(6,3,1, 3,4,3, 1,3,6), 3, byrow = TRUE)
w_2 <- c(0, 1, 2); v_2 <- c(0, 1, 2)
r_2 <- exact.score.sen.IxJ(obs.table = tb_2, gamma = 0, delta = c(0, 1, 1),
                           treatment.scores = w_2, outcome.scores = v_2)
p_perm_2 <- perm_pvalue(tb_2, function(t) sum(outer(w_2, v_2) * t))
check(sprintf("Theorem 2, ordinal, gamma = 0: max over U_O = rct.prob = permutation p-value %.8f", p_perm_2),
      near(r_2$max.prob, p_perm_2) && near(r_2$rct.prob, p_perm_2))

# Purpose:  Theorem 3 at gamma = 0: the sign-score worst case on a 2x2 table
#           is the one-sided Fisher exact test
# Input:    a 2x2 table, scores (0, 1) x (0, 1), delta = (0, 1)
# Expected: max.prob == fisher.test(alternative = "greater")$p.value
tb_3 <- matrix(c(9,4, 5,10), 2, byrow = TRUE)
r_3 <- exact.score.sen.IxJ(obs.table = tb_3, gamma = 0, delta = c(0, 1),
                           treatment.scores = c(0, 1), outcome.scores = c(0, 1))
p_fisher_3 <- fisher.test(tb_3, alternative = "greater")$p.value
check(sprintf("Theorem 3, 2x2 sign-score, gamma = 0: %.10f = one-sided fisher.test %.10f", r_3$max.prob, p_fisher_3),
      near(r_3$max.prob, p_fisher_3))

# =============================================================================
# (2) Total probability mass is 1 for every (u, gamma).
# =============================================================================

# Purpose:  for each table, the point probabilities over all tables with the
#           observed margins sum to 1 at every (u-bar, gamma) tried; this
#           exercises the denominator kappa_0 S_0(u) against the numerator
#           over the whole space
# Input:    a 3x3, a 4x2, a 2x4 and a 4x3 table (the last with the
#           nonmonotone delta = (0, 1, 0, 1)); u-bar = 0, u-bar = N_.j, and
#           two random u-bar; gamma = 0, 0.7, 2.5
# Expected: |mass - 1| below 1e-10 in all 12 combinations per table
mass_cases <- list(
  list(tb = matrix(c(4,3,1, 2,5,3, 1,2,6), 3, byrow = TRUE),           delta = c(0,1,1)),
  list(tb = matrix(c(6,4, 5,5, 3,7, 2,8), 4, byrow = TRUE),            delta = c(0,0,1,1)),
  list(tb = matrix(c(5,2,4,3, 2,6,3,5), 2, byrow = TRUE),              delta = c(0,1)),
  list(tb = matrix(c(3,2,2, 1,4,2, 2,1,3, 2,2,3), 4, byrow = TRUE),    delta = c(0,1,0,1)))
for (cs in mass_cases) {
  tabs <- all_tables(cs$tb); m <- colSums(cs$tb)
  us <- rbind(rep(0, length(m)), m,
              sapply(m, function(x) sample(0:x, 1)), sapply(m, function(x) sample(0:x, 1)))
  worst <- 0
  for (gam in c(0, 0.7, 2.5)) for (k in seq_len(nrow(us))) {
    mass <- sum(sapply(tabs, function(t) point_prob(t, us[k, ], gam, cs$delta)))
    worst <- max(worst, abs(mass - 1))
  }
  check(sprintf("total mass = 1: %dx%d, delta = (%s), %d tables, 4 u-bar x 3 gamma (worst |mass - 1| %.1e)",
                nrow(cs$tb), ncol(cs$tb), paste(cs$delta, collapse = ","), length(tabs), worst),
        worst < 1e-10)
}


# =============================================================================
# (3) Exact moments of an ordinal statistic equal the Theorem D.1 moments
#     used by the normal approximation.
# =============================================================================

# Purpose:  the mean and variance of T = sum_ij w_i v_j N_ij computed by brute
#           force from P_u(.) over all tables equal mu_T and var_T from
#           norm_single_u_allocation_p_value (Theorem D.1); the kernel
#           counting of Theorem 4 and the moment algebra of Appendix D are two
#           independent derivations
# Input:    a 3x3 table, the 3x3 table of Appendix F.2, a 4x2 and a 2x4 table,
#           each with nondecreasing scores and delta; u-bar = 0, u-bar = N_.j,
#           and one random u-bar; gamma = 0 and 1
# Expected: |mean diff| below 1e-8 and |var diff| below 1e-6 in all 6
#           combinations per table
mom_cases <- list(
  list(tb = matrix(c(6,3,1, 3,4,3, 1,3,6), 3, byrow = TRUE),      delta = c(0,1,1),   w = c(0,1,2),   v = c(0,1,2)),
  list(tb = matrix(c(12,18,5, 6,12,6, 6,6,15), 3, byrow = TRUE),  delta = c(0,1,1),   w = c(0,1,2.5), v = c(0,1,2)),
  list(tb = matrix(c(6,4, 5,5, 3,7, 2,8), 4, byrow = TRUE),       delta = c(0,1,1,1), w = c(1,2,3,4), v = c(0,1)),
  list(tb = matrix(c(5,2,4,3, 2,6,3,5), 2, byrow = TRUE),         delta = c(0,1),     w = c(0,1),     v = c(1,2,3,4)))
for (cs in mom_cases) {
  tabs <- all_tables(cs$tb); m <- colSums(cs$tb); S <- outer(cs$w, cs$v)
  Tvals <- sapply(tabs, function(t) sum(S * as.matrix(t)))
  us <- rbind(rep(0, length(m)), m, sapply(m, function(x) sample(0:x, 1)))
  worst_mean <- 0; worst_var <- 0
  for (gam in c(0, 1)) for (k in seq_len(nrow(us))) {
    p <- sapply(tabs, function(t) point_prob(t, us[k, ], gam, cs$delta))
    mean_exact <- sum(Tvals * p); var_exact <- sum(Tvals^2 * p) - mean_exact^2
    nm <- norm_single_u_allocation_p_value(obs.table = cs$tb, gamma = gam, delta = cs$delta,
                                           u_allocation = us[k, ],
                                           treatment.scores = cs$w, outcome.scores = cs$v)
    worst_mean <- max(worst_mean, abs(mean_exact - nm$mu_T))
    worst_var  <- max(worst_var,  abs(var_exact  - nm$var_T))
  }
  check(sprintf("moments: %dx%d, delta = (%s), 3 u-bar x 2 gamma (worst diff: mean %.1e, var %.1e)",
                nrow(cs$tb), ncol(cs$tb), paste(cs$delta, collapse = ","), worst_mean, worst_var),
        worst_mean < 1e-8 && worst_var < 1e-6)
}


# =============================================================================
# Summary
# =============================================================================
cat(sprintf("\n%d passed, %d failed\n", sum(results), sum(!results)))
# =============================================================================
# verify_same_names.R
#
# Code for "Exact, Nonparametric Sensitivity Analysis for Observational Studies
# of Contingency Tables"
#
# Runs the previous engine (legacy/sensitivityIxJ_orig.cpp and .R: one
# numerator per table dimension, enumeration in R) and the current engine
# (library/sensitivityIxJ.cpp and .R) under the same public function names,
# on the tables that appear in the paper, and checks that the reference sets
# and p-values agree. Timings of both engines are printed in a table.
#
# Run from the top-level folder (the one with README.txt):
#   source("development_check/verify_same_names.R")
# Each check prints PASS or FAIL; a summary is printed at the end. The legacy
# enumeration is done in R, so the larger tables take minutes on that side.
# =============================================================================
library(Rcpp)

# ---- Setup ------------------------------------------------------------------
# The two engines are loaded into separate environments so that the same
# function names can be called side by side.
env_old <- new.env()
sourceCpp("legacy/sensitivityIxJ_orig.cpp", env = env_old)
sys.source("legacy/sensitivityIxJ_orig.R", envir = env_old)
env_new <- new.env()
sourceCpp("library/sensitivityIxJ.cpp", env = env_new)
sys.source("library/sensitivityIxJ.R", envir = env_new)
sourceCpp("library/SIS_methods.cpp")     # table sampler, shared by both engines

# Record and print the result of one check.
results <- logical(0)
check <- function(label, expr) {
  ok <- tryCatch(isTRUE(expr), error = function(e) FALSE)
  results[label] <<- ok
  cat(if (ok) "PASS  " else "FAIL  ", label, "\n", sep = "")
}

# TRUE if x and y agree within numerical tolerance.
near <- function(x, y, tol = 1e-8) isTRUE(all.equal(x, y, tolerance = tol))

# Evaluate expr and return its value with the elapsed seconds.
tm <- function(expr) {
  t0 <- proc.time()[["elapsed"]]; v <- expr
  list(v = v, s = proc.time()[["elapsed"]] - t0)
}

chisq_stat <- function(tb) { E <- outer(rowSums(tb), colSums(tb)) / sum(tb); sum((tb - E)^2 / E) }

# Tables from the paper: Appendix F.1 (Table F.1), Appendix F.2, and the
# data analysis of Section 6 (girls and boys, prior and profile scores). The
# last case uses the chi-square statistic through the R-callback path.
cases <- list(
  list(src = "F.1 block 1", tb = matrix(c(2,3,0, 0,1,4, 0,1,4), 3, byrow = TRUE), u = c(0,0,3),  gamma = 1, delta = c(0,1,1), w = c(0,1,2), v = c(0,1,2)),
  list(src = "F.1 block 2", tb = matrix(c(2,2,1, 1,2,2, 1,2,2), 3, byrow = TRUE), u = c(0,3,5),  gamma = 1, delta = c(0,1,1), w = c(0,1,2), v = c(0,1,2)),
  list(src = "F.1 block 3", tb = matrix(c(3,2,1, 0,2,4, 0,1,5), 3, byrow = TRUE), u = c(0,4,10), gamma = 1, delta = c(0,1,1), w = c(0,1,2), v = c(0,1,2)),
  list(src = "F.1 block 4", tb = matrix(c(3,3,0, 1,2,3, 2,3,1), 3, byrow = TRUE), u = c(0,8,4),  gamma = 1, delta = c(0,1,1), w = c(0,1,2), v = c(0,1,2)),
  list(src = "F.1 2x3 b1",  tb = matrix(c(4,6,0, 1,3,6), 2, byrow = TRUE),        u = c(0,1,6),  gamma = 1, delta = c(0,1),   w = c(0,1),   v = c(0,1,2)),
  list(src = "F.1 2x3 b2",  tb = matrix(c(2,4,4, 2,2,6), 2, byrow = TRUE),        u = c(3,6,10), gamma = 1, delta = c(0,1),   w = c(0,1),   v = c(0,1,2)),
  list(src = "F.2 exp 1.1 (N=86)",  tb = matrix(c(12,18,5, 6,12,6, 6,6,15), 3, byrow = TRUE),  u = c(0,10,20), gamma = 1,   delta = c(0,1,1), w = c(0,1,2.5), v = c(0,1,2)),
  list(src = "F.2 girls (N=94)",    tb = matrix(c(12,3,0, 18,12,3, 17,25,4), 3, byrow = TRUE), u = c(0,40,7),  gamma = 0.5, delta = c(0,1,1), w = c(0,1,2.5), v = c(0,1,2)),
  list(src = "F.2 boys (N=112)",    tb = matrix(c(10,8,1, 29,11,3, 20,24,6), 3, byrow = TRUE), u = c(0,20,10), gamma = 0.5, delta = c(0,1,1), w = c(0,1,2.5), v = c(0,1,2)),
  list(src = "Data girls, prior",   tb = matrix(c(12,3,0, 18,12,3, 17,25,4), 3, byrow = TRUE), u = c(0,15,7),  gamma = log(2), delta = c(0,1,1), w = c(0,0.25,1.5), v = c(0,1,1.5)),
  list(src = "Data boys, prior",    tb = matrix(c(10,8,1, 29,11,3, 20,24,6), 3, byrow = TRUE), u = c(0,20,10), gamma = log(2), delta = c(0,1,1), w = c(0,0.25,1.5), v = c(0,1,1.5)),
  list(src = "Data girls, profile", tb = matrix(c(12,3,0, 18,12,3, 17,25,4), 3, byrow = TRUE), u = c(0,15,7),  gamma = log(2), delta = c(0,1,1), w = c(0,0.2,1.6),  v = c(0,0.67,0.7)),
  list(src = "Data boys, profile",  tb = matrix(c(10,8,1, 29,11,3, 20,24,6), 3, byrow = TRUE), u = c(0,20,10), gamma = log(2), delta = c(0,1,1), w = c(0,0.2,1.6),  v = c(0,0.67,0.7)),
  list(src = "F.2 (N=86) chi-sq",   tb = matrix(c(12,18,5, 6,12,6, 6,6,15), 3, byrow = TRUE),  u = c(0,10,20), gamma = 1,   delta = c(0,1,1), chisq = TRUE)
)


# =============================================================================
# exact.general.sen.IxJ at one u-bar: reference set size and p-value, old
# engine (possible.table in R, then table_space) vs new engine (internal
# enumeration), with timings.
# =============================================================================

# Purpose:  for each table of the paper, the two engines find the same number
#           of tables in Omega_t and the same p-value at the given u-bar
# Input:    the 14 cases above, each at its u-bar, gamma and delta; the
#           ordinal statistic with the stated scores, or chi-square
# Expected: |Omega_t| equal; p-values equal within 1e-8
rows <- list()
for (cs in cases) {
  if (isTRUE(cs$chisq)) { stat <- chisq_stat } else { sc <- outer(cs$w, cs$v); stat <- function(tb) sum(sc * tb) }
  o_enum <- tm(env_old$possible.table(threshold = stat(cs$tb), table = cs$tb, direction = "greater than", transform.fun = stat))
  o_p    <- tm(env_old$exact.general.sen.IxJ(u_space = t(as.matrix(cs$u)), obs.table = cs$tb, table_space = o_enum$v,
                                             gamma = cs$gamma, delta = cs$delta))
  n_p    <- tm(env_new$exact.general.sen.IxJ(u_space = t(as.matrix(cs$u)), obs.table = cs$tb, table_space = NULL,
                                             gamma = cs$gamma, delta = cs$delta, transform.fun = stat))
  same <- near(o_p$v$max.prob, n_p$v$max.prob) && length(o_enum$v) == n_p$v$n.tables
  rows[[length(rows) + 1]] <- data.frame(
    case = cs$src, N = sum(cs$tb),
    tables_old = length(o_enum$v), tables_new = n_p$v$n.tables, distinct_TD = n_p$v$n.distinct.tbar,
    p_old = signif(o_p$v$max.prob, 7), p_new = signif(n_p$v$max.prob, 7), same = same,
    enum_old_s = round(o_enum$s, 2), enum_new_s = round(n_p$v$time[["enumerate.W"]], 3),
    p_old_s = round(o_p$s, 2), p_new_s = round(n_p$v$time[["sweep"]], 3))
  check(sprintf("%-22s |Omega_t| old %d / new %d, p old %.6f / new %.6f", cs$src,
                length(o_enum$v), n_p$v$n.tables, o_p$v$max.prob, n_p$v$max.prob),
        same)
}
timings <- do.call(rbind, rows); rownames(timings) <- NULL
cat("\nTimings (seconds): enumeration and single-u p-value, old vs new engine\n")
print(timings, row.names = FALSE); cat("\n")


# =============================================================================
# exact.score.sen.IxJ: old vs new on the data-analysis table.
# =============================================================================

# Purpose:  the ordinal worst-case p-value and its maximizer agree between the
#           engines on the girls' table with the profile-likelihood scores
# Input:    girls' 3x3 table, gamma = log 2, delta = (0,1,1),
#           scores (0, 0.2, 1.6) x (0, 0.67, 0.7)
# Expected: max.prob equal within 1e-8 and the same maximizer
girls <- matrix(c(12,3,0, 18,12,3, 17,25,4), 3, byrow = TRUE)
o <- tm(env_old$exact.score.sen.IxJ(obs.table = girls, gamma = log(2), delta = c(0,1,1), row = "treatment",
                                    treatment.scores = c(0,0.2,1.6), outcome.scores = c(0,0.67,0.7)))
n <- tm(env_new$exact.score.sen.IxJ(obs.table = girls, gamma = log(2), delta = c(0,1,1), row = "treatment",
                                    treatment.scores = c(0,0.2,1.6), outcome.scores = c(0,0.67,0.7)))
check(sprintf("exact.score.sen.IxJ girls/profile/Gamma = 2: old %.6f (%.1f s), new %.6f (%.3f s), maximizer old (%s) new (%s)",
              o$v$max.prob, o$s, n$v$max.prob, n$s,
              paste(o$v$maximizer, collapse = ","), paste(n$v$maximizer, collapse = ",")),
      near(o$v$max.prob, n$v$max.prob) && identical(as.numeric(o$v$maximizer), as.numeric(n$v$maximizer)))


# =============================================================================
# possible.table: the same set of tables, not just the same count.
# =============================================================================

# Purpose:  the old (R) and new (C++) enumerations return identical sets
# Input:    the N = 86 table of Appendix F.2, ordinal statistic
#           (0, 1, 2.5) x (0, 1, 2)
# Expected: the sorted cell-vector keys of the two lists are identical
cs <- cases[[7]]; sc <- outer(cs$w, cs$v); stat <- function(tb) sum(sc * tb)
po <- env_old$possible.table(threshold = stat(cs$tb), table = cs$tb, direction = "greater than", transform.fun = stat)
pn <- env_new$possible.table(threshold = stat(cs$tb), table = cs$tb, direction = "greater than", transform.fun = stat)
key <- function(L) sort(sapply(L, function(t) paste(as.vector(t), collapse = ",")))
check(sprintf("possible.table on the N = 86 table: old %d tables, new %d tables, identical sets", length(po), length(pn)),
      identical(key(po), key(pn)))


# =============================================================================
# norm.score.sen.IxJ: the moment code is unchanged between engines.
# =============================================================================

# Purpose:  the normal-approximation worst-case p-value is identical
# Input:    girls' table, gamma = log 2, delta = (0,1,1), profile scores
# Expected: max.prob identical
o <- env_old$norm.score.sen.IxJ(obs.table = girls, gamma = log(2), delta = c(0,1,1),
                                treatment.scores = c(0,0.2,1.6), outcome.scores = c(0,0.67,0.7))
n <- env_new$norm.score.sen.IxJ(obs.table = girls, gamma = log(2), delta = c(0,1,1),
                                treatment.scores = c(0,0.2,1.6), outcome.scores = c(0,0.67,0.7))
check(sprintf("norm.score.sen.IxJ girls/profile/Gamma = 2: old %.8f, new %.8f", o$max.prob, n$max.prob),
      near(o$max.prob, n$max.prob, tol = 1e-12))


# =============================================================================
# sampling.score.sen.IxJ: the importance weights now come from the new
# kernel; with the same seed the same tables are drawn, so the estimates
# must agree.
# =============================================================================

# Purpose:  the Monte Carlo worst-case p-value agrees between engines
# Input:    F.1 block 1 table, gamma = 1, delta = (0,1,1), scores
#           (0,1,2) x (0,1,2), 2000 draws, seed 1
# Expected: max.prob equal within 1e-6
tb_small <- matrix(c(2,3,0, 0,1,4, 0,1,4), 3, byrow = TRUE)
set.seed(1); o <- env_old$sampling.score.sen.IxJ(obs.table = tb_small, gamma = 1, delta = c(0,1,1), mc.iteration = 2000,
                                                 treatment.scores = c(0,1,2), outcome.scores = c(0,1,2))
set.seed(1); n <- env_new$sampling.score.sen.IxJ(obs.table = tb_small, gamma = 1, delta = c(0,1,1), mc.iteration = 2000,
                                                 treatment.scores = c(0,1,2), outcome.scores = c(0,1,2))
check(sprintf("sampling.score.sen.IxJ F.1 block 1, 2000 draws: old %.6f, new %.6f", o$max.prob, n$max.prob),
      near(o$max.prob, n$max.prob, tol = 1e-6))


# =============================================================================
# Summary
# =============================================================================
cat(sprintf("\n%d passed, %d failed\n", sum(results), sum(!results)))
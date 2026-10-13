## edge_case_checks.R
## Exercises the code paths that verify_same_names.R did not: orientation, two-valued
## delta, "less than", threshold ties with fractional scores, gamma = 0 invariance,
## nonmonotone delta, a 2x2 against fisher.test, and degenerate margins.
## Files: sensitivityIxJ_orig.cpp/.R, sensitivityIxJ.cpp/.R in the working directory.

library(Rcpp)
env_old <- new.env(); sourceCpp("sensitivityIxJ_orig.cpp", env = env_old); sys.source("sensitivityIxJ_orig.R", envir = env_old)
env_new <- new.env(); sourceCpp("sensitivityIxJ.cpp", env = env_new); sys.source("sensitivityIxJ.R", envir = env_new)
ok <- function(cond, msg) cat(sprintf("%-6s %s\n", if (isTRUE(cond)) "PASS" else "FAIL", msg))
same <- function(a, b, tol = 1e-8) isTRUE(all.equal(a, b, tolerance = tol))

tb <- matrix(c(12,18,5, 6,12,6, 6,6,15), 3, byrow = TRUE)     # F.2 table, N = 86
w <- c(0,1,2.5); v <- c(0,1,2); sc <- outer(w, v); stat <- function(t) sum(sc * t)

## 1. row = "outcome": pass the transposed table and say so; must equal row = "treatment"
a <- env_new$exact.score.sen.IxJ(obs.table = tb,    gamma = 1, delta = c(0,1,1), row = "treatment", treatment.scores = w, outcome.scores = v)
b <- env_new$exact.score.sen.IxJ(obs.table = t(tb), gamma = 1, delta = c(0,1,1), row = "outcome",   treatment.scores = w, outcome.scores = v)
ok(same(a$max.prob, b$max.prob) && identical(a$maximizer, b$maximizer), "row = 'outcome' gives the same result as row = 'treatment'")

## 2. delta must be binary (model (6)): a two-valued non-binary delta is refused
e <- tryCatch({ env_new$exact.score.sen.IxJ(obs.table = tb, gamma = 0.5, delta = c(1,3,3), treatment.scores = w, outcome.scores = v); "no error" },
              error = function(e) "error")
ok(e == "error", "delta = (1,3,3) is rejected (delta must be binary)")
e <- tryCatch({ env_new$exact.score.sen.IxJ(obs.table = tb, gamma = -1, delta = c(0,1,1), treatment.scores = w, outcome.scores = v); "no error" },
              error = function(e) "error")
ok(e == "error", "gamma = -1 is rejected (gamma must be nonnegative)")

## 3. direction = "less than" in possible.table: old and new give the same set
po <- env_old$possible.table(threshold = stat(tb), table = tb, direction = "less than", transform.fun = stat)
pn <- env_new$possible.table(threshold = stat(tb), table = tb, direction = "less than", transform.fun = stat)
key <- function(L) sort(sapply(L, function(x) paste(as.vector(x), collapse = ",")))
ok(length(po) == length(pn) && identical(key(po), key(pn)), sprintf("possible.table 'less than': identical sets (%d tables)", length(pn)))

## 4. threshold ties with fractional scores: counts must agree table-for-table.
##    Scores chosen so that many distinct tables share the observed statistic.
w2 <- c(0, 0.3, 0.6); v2 <- c(0, 0.7, 1.4); sc2 <- outer(w2, v2); stat2 <- function(t) sum(sc2 * t)
tb2 <- matrix(c(6,4,2, 4,4,4, 2,4,6), 3, byrow = TRUE)
po <- env_old$possible.table(threshold = stat2(tb2), table = tb2, direction = "greater than", transform.fun = stat2)
pn <- env_new$possible.table(threshold = stat2(tb2), table = tb2, direction = "greater than", transform.fun = stat2)
n_tie <- sum(sapply(po, function(x) abs(stat2(x) - stat2(tb2)) < 1e-12))
ok(length(po) == length(pn) && identical(key(po), key(pn)),
   sprintf("ties with fractional scores: identical sets (%d tables, %d at the threshold)", length(pn), n_tie))
a <- env_old$exact.score.sen.IxJ(obs.table = tb2, gamma = 1, delta = c(0,1,1), treatment.scores = w2, outcome.scores = v2)
b <- env_new$exact.score.sen.IxJ(obs.table = tb2, gamma = 1, delta = c(0,1,1), treatment.scores = w2, outcome.scores = v2)
ok(same(a$max.prob, b$max.prob) && b$n.tables == length(po),
   sprintf("exact.score.sen.IxJ with tied fractional scores: old %.8f new %.8f", a$max.prob, b$max.prob))

## 5. gamma = 0: alpha must not depend on u at all
r <- env_new$exact.general.sen.IxJ(u_space = "full", obs.table = tb, gamma = 0, delta = c(0,1,1), transform.fun = stat)
ok(diff(range(r$all.prob)) < 1e-12 && same(r$max.prob, r$rct.prob), "gamma = 0: alpha identical over all of U_C and equal to rct.prob")

## 6. nonmonotone delta (0,1,0) through exact.general.sen.IxJ: old vs new on an explicit u_space
U <- env_new$full.candidate.set(colSums(tb2))
ts <- env_old$possible.table(threshold = stat(tb2), table = tb2, direction = "greater than", transform.fun = stat)
a <- env_old$exact.general.sen.IxJ(u_space = U, obs.table = tb2, table_space = ts, gamma = 1, delta = c(0,1,0))
b <- env_new$exact.general.sen.IxJ(u_space = U, obs.table = tb2, table_space = ts, gamma = 1, delta = c(0,1,0))
ok(same(a$max.prob, b$max.prob) && identical(as.numeric(a$maximizer), as.numeric(b$maximizer)),
   sprintf("nonmonotone delta (0,1,0), full U_C: old %.8f new %.8f", a$max.prob, b$max.prob))
##    and delta with the reference level last, (1,1,0)
a <- env_old$exact.general.sen.IxJ(u_space = U, obs.table = tb2, table_space = ts, gamma = 1, delta = c(1,1,0))
b <- env_new$exact.general.sen.IxJ(u_space = U, obs.table = tb2, table_space = ts, gamma = 1, delta = c(1,1,0))
ok(same(a$max.prob, b$max.prob), sprintf("delta = (1,1,0), full U_C: old %.8f new %.8f", a$max.prob, b$max.prob))

## 7. 2x2 at gamma = 0 against fisher.test (one-sided, greater)
tb3 <- matrix(c(9,3, 4,8), 2, byrow = TRUE)
p_new <- env_new$exact.score.sen.IxJ(obs.table = tb3, gamma = 0, delta = c(0,1), treatment.scores = c(0,1), outcome.scores = c(0,1))$max.prob
## T = N_22, "greater than" on N_22 is Fisher's one-sided test with alternative = "greater" on the odds ratio
p_fisher <- fisher.test(tb3, alternative = "greater")$p.value
ok(same(p_new, p_fisher, 1e-10), sprintf("2x2, gamma = 0: new %.10f vs fisher.test %.10f", p_new, p_fisher))

## 8. table_space given AND transform.fun given: table_space wins and equals the internal enumeration
a <- env_new$exact.general.sen.IxJ(u_space = U, obs.table = tb2, table_space = ts, gamma = 1, delta = c(0,1,1))
b <- env_new$exact.general.sen.IxJ(u_space = U, obs.table = tb2, table_space = NULL, gamma = 1, delta = c(0,1,1), transform.fun = stat)
ok(same(a$max.prob, b$max.prob) && a$n.tables == b$n.tables, "table_space path == internal enumeration path")

## 8b. default u_space of exact.general.sen.IxJ is now the full U_C of Theorem 1:
##     equals the explicit u_space = U run, and is >= the previous release's heuristic default
a <- env_new$exact.general.sen.IxJ(obs.table = tb2, table_space = ts, gamma = 1, delta = c(0,1,1))
b <- env_new$exact.general.sen.IxJ(u_space = U, obs.table = tb2, table_space = ts, gamma = 1, delta = c(0,1,1))
o <- env_old$exact.general.sen.IxJ(obs.table = tb2, table_space = ts, gamma = 1, delta = c(0,1,1))
ok(same(a$max.prob, b$max.prob) && a$n.candidates == prod(colSums(tb2) + 1) && a$max.prob >= o$max.prob - 1e-12,
   sprintf("default u_space = U_C (%d candidates): %.8f; old heuristic default: %.8f", a$n.candidates, a$max.prob, o$max.prob))

## 8c. native chi-square path (statistic = "chisq") == R-callback path (transform.fun = chi-square)
chisq_stat <- function(t) { E <- outer(rowSums(t), colSums(t)) / sum(t); sum((t - E)^2 / E) }
a <- env_new$exact.general.sen.IxJ(u_space = U, obs.table = tb2, gamma = 1, delta = c(0,1,1), transform.fun = chisq_stat)
b <- env_new$exact.general.sen.IxJ(u_space = U, obs.table = tb2, gamma = 1, delta = c(0,1,1), statistic = "chisq")
ok(same(a$max.prob, b$max.prob) && a$n.tables == b$n.tables && identical(as.numeric(a$maximizer), as.numeric(b$maximizer)),
   sprintf("chi-square in C++ vs via R callback: %d vs %d tables, p %.8f vs %.8f", b$n.tables, a$n.tables, b$max.prob, a$max.prob))
g2_stat <- function(t) { E <- outer(rowSums(t), colSums(t)) / sum(t); 2 * sum(ifelse(t > 0, t * log(pmax(t, 1) / E), 0)) }
a <- env_new$exact.general.sen.IxJ(u_space = U, obs.table = tb2, gamma = 1, delta = c(0,1,1), transform.fun = g2_stat)
b <- env_new$exact.general.sen.IxJ(u_space = U, obs.table = tb2, gamma = 1, delta = c(0,1,1), statistic = "G2")
ok(same(a$max.prob, b$max.prob) && a$n.tables == b$n.tables,
   sprintf("G^2 in C++ vs via R callback: %d vs %d tables, p %.8f vs %.8f", b$n.tables, a$n.tables, b$max.prob, a$max.prob))
fisher_stat <- function(t) sum(lgamma(t + 1))
a <- env_new$exact.general.sen.IxJ(u_space = matrix(0L, 1, 3), obs.table = tb2, gamma = 0, delta = c(0,1,1), statistic = "freeman_halton")
ok(same(a$max.prob, fisher.test(tb2)$p.value, 1e-8),
   sprintf("Freeman-Halton in C++ at gamma = 0: %.8f vs fisher.test() (r x c) %.8f", a$max.prob, fisher.test(tb2)$p.value))

## 9. degenerate inputs must error, not mis-answer
e1 <- tryCatch({ env_new$exact.score.sen.IxJ(obs.table = matrix(c(5,0, 3,0), 2, byrow = TRUE), gamma = 1, delta = c(0,1),
                                             treatment.scores = c(0,1), outcome.scores = c(0,1)); "no error" }, error = function(e) "error")
e2 <- tryCatch({ env_new$exact.score.sen.IxJ(obs.table = tb, gamma = 1, delta = c(0,1,2), treatment.scores = w, outcome.scores = v); "no error" },
               error = function(e) "error")
cat(sprintf("INFO   zero outcome margin: %s | three-valued delta: %s  (both should be 'error')\n", e1, e2))

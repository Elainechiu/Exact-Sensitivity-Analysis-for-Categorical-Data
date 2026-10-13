suppressMessages(library(Rcpp)); sourceCpp("sensitivityIxJ.cpp"); sourceCpp("SIS_methods.cpp"); source("sensitivityIxJ.R")
ok <- function(c, m) cat(sprintf("%-6s %s\n", if (isTRUE(c)) "PASS" else "FAIL", m))
all_tables <- function(tb) possible.table(0, tb, "greater than", function(t) 1)

## (A) norm.score.sen.IxJ default u_space == U_O (J >= 3) / u+ (J = 2), and its maximizer
##     is the candidate where the exact worst case is attained (same set), moments equal.
cases <- list(list(tb = matrix(c(12,18,5, 6,12,6, 6,6,15), 3, byrow = TRUE), d = c(0,1,1), w = c(0,1,2.5), v = c(0,1,2)),
              list(tb = matrix(c(40,20, 5,5, 10,10), 3, byrow = TRUE),        d = c(0,0,1), w = c(0,1,2),   v = c(0,1)),
              list(tb = matrix(c(5,2,4,3, 2,6,3,5), 2, byrow = TRUE),         d = c(0,1),   w = c(0,1),     v = c(1,2,3,4)))
for (cs in cases) {
  m <- colSums(cs$tb); J <- ncol(cs$tb); N <- sum(cs$tb)
  U <- if (J == 2) .u_space_sign_score(m) else .u_space_ordinal(m)
  nm <- norm.score.sen.IxJ(obs.table = cs$tb, gamma = 1, delta = cs$d, treatment.scores = cs$w, outcome.scores = cs$v)
  ex <- exact.score.sen.IxJ(obs.table = cs$tb, gamma = 1, delta = cs$d, treatment.scores = cs$w, outcome.scores = cs$v)
  ok(ex$n.candidates == nrow(U), sprintf("(A) %dx%d: exact candidate set has %d rows = default normal set (%s)", nrow(cs$tb), J, ex$n.candidates,
                                         if (J == 2) "u+ of Thm 3" else sprintf("U_O of Thm 2, N+1 = %d", N + 1)))
  ## exact moments of T at EVERY candidate vs Theorem D.1 moments
  tbs <- all_tables(cs$tb); S <- outer(cs$w, cs$v); Tv <- sapply(tbs, function(t) sum(S * as.matrix(t)))
  worst <- 0
  for (k in seq_len(nrow(U))) for (gam in c(0, 1)) {
    p <- sapply(tbs, function(t) point_prob_bd(as.matrix(t) + 0L, as.integer(cs$d), as.integer(U[k, ]), gam))
    me <- sum(Tv * p); ve <- sum(Tv^2 * p) - me^2
    nm1 <- norm_single_u_allocation_p_value(obs.table = cs$tb, gamma = gam, delta = cs$d, u_allocation = U[k, ],
                                            treatment.scores = cs$w, outcome.scores = cs$v)
    worst <- max(worst, abs(me - nm1$mu_T), abs(ve - nm1$var_T))
  }
  ok(worst < 1e-6, sprintf("(A) %dx%d: exact vs D.1 mean/var of T at all %d candidates x gamma in {0,1}: worst diff %.1e",
                           nrow(cs$tb), J, nrow(U), worst))
  cat(sprintf("       exact worst p %.5f at u-bar (%s); normal worst p %.5f at u-bar (%s)\n",
              ex$max.prob, paste(ex$maximizer, collapse = ","), nm$max.prob, paste(nm$maximizer, collapse = ",")))
}

## (B) gamma = 0 through the general path: normal p-value independent of u-bar
tb <- cases[[1]]$tb
p0 <- sapply(1:5, function(k) norm_single_u_allocation_p_value(tb, 0, c(0,1,1), .u_space_ordinal(colSums(tb))[k * 10, ], treatment.scores = c(0,1,2.5), outcome.scores = c(0,1,2))$p_value)
ok(diff(range(p0)) < 1e-12, sprintf("(B) gamma = 0: normal-approx p identical over 5 different u-bar (%.6f)", p0[1]))

## (C) norm.score.sen.IxJ rejects nonmonotone delta / scores (assumption of Theorems 2-3)
e1 <- tryCatch({ norm.score.sen.IxJ(tb, 1, c(0,1,0), treatment.scores = c(0,1,2), outcome.scores = c(0,1,2)); "no error" }, error = function(e) "error")
e2 <- tryCatch({ norm.score.sen.IxJ(tb, 1, c(0,1,1), treatment.scores = c(0,2,1), outcome.scores = c(0,1,2)); "no error" }, error = function(e) "error")
ok(e1 == "error" && e2 == "error", "(C) norm.score.sen.IxJ stops for delta = (0,1,0) and for decreasing treatment scores")

## (D) sampling.general.sen.IxJ default u_space = U_C: small table, chi-square, compare to exact
tb2 <- matrix(c(4,1,1, 1,3,2, 1,1,4), 3, byrow = TRUE)                      # N = 18, |U_C| = 7*5*7 = 245
chisq <- function(t) { E <- outer(rowSums(t), colSums(t)) / sum(t); sum((t - E)^2 / E) }
set.seed(3)
sg <- sampling.general.sen.IxJ(obs.table = tb2, gamma = 1, delta = c(0,1,1), transform.fun = chisq, mc.iteration = 400)
eg <- exact.general.sen.IxJ(obs.table = tb2, gamma = 1, delta = c(0,1,1), transform.fun = chisq)
ok(eg$n.candidates == prod(colSums(tb2) + 1), sprintf("(D) exact U_C has %d candidates = prod(N_.j + 1)", eg$n.candidates))
cat(sprintf("       sampling max p %.4f at (%s) [default u_space = U_C, %d candidates x 400 draws]; exact max p %.4f at (%s)\n",
            sg$max.prob, paste(sg$maximizer, collapse = ","), prod(colSums(tb2) + 1), eg$max.prob, paste(eg$maximizer, collapse = ",")))
ok(abs(sg$max.prob - eg$max.prob) < 0.1, "(D) sampling worst-case p within 0.1 of exact (Monte Carlo, 400 draws)")

## (E) sampling.score.sen.IxJ uses U_O: compare with exact on the 3x3 table
set.seed(4)
ss <- sampling.score.sen.IxJ(obs.table = tb2, gamma = 1, delta = c(0,1,1), treatment.scores = c(0,1,2), outcome.scores = c(0,1,2), mc.iteration = 1000)
es <- exact.score.sen.IxJ(obs.table = tb2, gamma = 1, delta = c(0,1,1), treatment.scores = c(0,1,2), outcome.scores = c(0,1,2))
cat(sprintf("       sampling.score max p %.4f at (%s); exact %.4f at (%s)\n", ss$max.prob, paste(ss$maximizer, collapse = ","), es$max.prob, paste(es$maximizer, collapse = ",")))
ok(abs(ss$max.prob - es$max.prob) < 0.05, "(E) sampling.score worst-case p within 0.05 of exact")

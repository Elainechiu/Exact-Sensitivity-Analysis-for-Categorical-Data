## =============================================================================
## edge_case_check.R
## Check the statements of Theorems 1--4 via checking whether they collapse into the
## usual cases.
##
## Depends only on sensitivityIxJ.cpp and sensitivityIxJ.R (same directory).
##
## Logic.  Theorem 4 (with Corollary C.3) says that, for a binary u, the exact null
## distribution of the table under the generic bias model (6),
##     P_u(N = t)  =  kappa(t) S(t-bar, u) / (kappa_0 S_0(u)),
## and the p-value alpha(T, r, u) = P_u(T(N) >= T_obs).  Theorems 1--3 say which u to
## evaluate (U_C, U_O, or u+).  If the machinery is right it must reproduce three things
## that are known independently of the paper:
##
##   (1) gamma = 0.  The model reduces to random assignment, so P_u(N = t) must be the
##       multivariate hypergeometric probability
##           P(t) = prod_i N_i.! prod_j N_.j! / ( N! prod_ij t_ij! ),
##       for EVERY u (u is irrelevant when gamma = 0), and every worst-case p-value from
##       Theorems 1--3 must equal the plain permutation p-value P(T >= T_obs), which
##       is also what fisher.test() computes for the probability-ordering statistic.
##
##   (2) Total mass.  For every u and every gamma, P_u(. ) is a probability
##       distribution on the tables with the observed margins, so summing the point
##       probabilities over ALL such tables must give exactly 1.  (This exercises the
##       denominator kappa_0 S_0(u) against the numerator over the whole space.)
##
##   (3) Moments.  Theorem D.1 (Supplementary Material, Appendix D) gives closed-form
##       first and second moments of the cell counts under the model at a given u, from
##       which norm.score.sen.IxJ builds the normal approximation.  The mean and
##       variance of an ordinal statistic T computed by brute force from the exact
##       distribution of Theorem 4 -- sum_t T(t) P_u(t) and sum_t T(t)^2 P_u(t) - mean^2
##       over all tables -- must equal the mean and variance that
##       norm_single_u_allocation_p_value() reports.  Two different derivations
##       (kernel counting vs. moment algebra) must agree.
## =============================================================================

library(Rcpp)
sourceCpp("sensitivityIxJ.cpp"); source("sensitivityIxJ.R")
set.seed(20260922)

ok <- function(cond, msg) cat(sprintf("%-6s %s\n", if (isTRUE(cond)) "PASS" else "FAIL", msg))
rand_table <- function(I, J, N) {
  z <- sample.int(I, N, TRUE); r <- sample.int(J, N, TRUE)
  tb <- matrix(0L, I, J); for (s in seq_len(N)) tb[z[s], r[s]] <- tb[z[s], r[s]] + 1L
  if (any(rowSums(tb) == 0) || any(colSums(tb) == 0)) return(rand_table(I, J, N)); tb
}
## multivariate hypergeometric log-pmf (the fixed-margins permutation distribution)
log_mvhyper <- function(t) sum(lgamma(rowSums(t) + 1)) + sum(lgamma(colSums(t) + 1)) -
                           lgamma(sum(t) + 1) - sum(lgamma(t + 1))
## all tables with the margins of tb (threshold 0 on a constant statistic keeps everything)
all_tables <- function(tb) possible.table(threshold = 0, table = tb, direction = "greater than",
                                          transform.fun = function(t) 1)
point_prob <- function(t, u, gamma, delta)
  generic.I.by.J.sensitivity.point.probability(table = as.matrix(t), u_allocation = u, gamma = gamma, delta = delta)

## =============================================================================
## (1) gamma = 0  ==>  multivariate hypergeometric
## =============================================================================
cat("\n==== (1) gamma = 0: the model collapses to random assignment (hypergeometric) ====\n")

## (1a) point probabilities, random tables in seven dimensions, arbitrary u
worst <- 0; n_cases <- 0
for (d in list(c(2,2), c(2,3), c(3,3), c(4,2), c(3,4), c(5,2), c(2,5))) for (rep in 1:20) {
  tb <- rand_table(d[1], d[2], sample(8:30, 1)); delta <- c(0, rep(1, d[1] - 1))
  u  <- sapply(colSums(tb), function(x) sample(0:x, 1))     # any u-bar: irrelevant at gamma = 0
  p_engine <- point_prob(tb, u, 0, delta); p_hyper <- exp(log_mvhyper(tb))
  worst <- max(worst, abs(p_engine - p_hyper) / p_hyper); n_cases <- n_cases + 1
}
ok(worst < 1e-10, sprintf("(1a) P_u(N = t) at gamma = 0 = multivariate hypergeometric pmf on %d random tables (worst rel. diff %.1e)", n_cases, worst))

## (1b) Theorem 1 (any count-based T, max over U_C): at gamma = 0 alpha is the same at
##      every candidate and equals the permutation p-value; illustrate with chi-square
tb <- matrix(c(5,3,2, 4,6,3, 2,3,6), 3, byrow = TRUE); delta <- c(0,1,1)
chisq_stat <- function(t) { E <- outer(rowSums(t), colSums(t)) / sum(t); sum((t - E)^2 / E) }
r1 <- exact.general.sen.IxJ(obs.table = tb, gamma = 0, delta = delta, transform.fun = chisq_stat)   # default u_space = U_C
tabs <- all_tables(tb)
p_perm <- sum(sapply(tabs, function(t) if (chisq_stat(as.matrix(t)) >= chisq_stat(tb) - 1e-9) exp(log_mvhyper(as.matrix(t))) else 0))
ok(diff(range(r1$all.prob)) < 1e-12 && abs(r1$max.prob - p_perm) < 1e-10,
   sprintf("(1b) Theorem 1, chi-square, gamma = 0: alpha constant over all %d candidates of U_C and = permutation p-value %.8f", r1$n.candidates, p_perm))

## (1c) Theorem 2 (ordinal T, max over U_O) and Theorem 3 (sign-score, J = 2) at gamma = 0
tb3 <- matrix(c(6,3,1, 3,4,3, 1,3,6), 3, byrow = TRUE); w <- c(0,1,2); v <- c(0,1,2)
r2 <- exact.score.sen.IxJ(obs.table = tb3, gamma = 0, delta = c(0,1,1), treatment.scores = w, outcome.scores = v)
tabs3 <- all_tables(tb3); T3 <- function(t) sum(outer(w, v) * t)
p_perm3 <- sum(sapply(tabs3, function(t) if (T3(as.matrix(t)) >= T3(tb3) - 1e-9) exp(log_mvhyper(as.matrix(t))) else 0))
ok(abs(r2$max.prob - p_perm3) < 1e-10 && abs(r2$rct.prob - p_perm3) < 1e-10,
   sprintf("(1c) Theorem 2, ordinal, gamma = 0: max over U_O = rct.prob = permutation p-value %.8f", p_perm3))
tb2 <- matrix(c(9,4, 5,10), 2, byrow = TRUE)
r3 <- exact.score.sen.IxJ(obs.table = tb2, gamma = 0, delta = c(0,1), treatment.scores = c(0,1), outcome.scores = c(0,1))
p_fisher <- fisher.test(tb2, alternative = "greater")$p.value     # one-sided Fisher = P(N_22 >= obs)
ok(abs(r3$max.prob - p_fisher) < 1e-10,
   sprintf("(1c) Theorem 3, 2x2 sign-score, gamma = 0: %.10f = fisher.test(alternative = 'greater') %.10f", r3$max.prob, p_fisher))

## (1d) Fisher's r x c exact test (probability ordering) vs fisher.test(), which uses the
##      Mehta-Patel network algorithm (descendant of Pagano & Halvorsen 1981)
fisher_stat <- function(t) sum(lgamma(t + 1))     # P(t) <= P(obs)  <=>  fisher_stat(t) >= fisher_stat(obs)
worst <- 0; n_cases <- 0
for (d in list(c(2,3), c(3,3), c(2,4), c(4,2), c(3,4))) for (rep in 1:6) {
  tbx <- rand_table(d[1], d[2], sample(10:22, 1)); dl <- c(0, rep(1, d[1] - 1))
  pe <- exact.general.sen.IxJ(u_space = matrix(0L, 1, d[2]), obs.table = tbx, gamma = 0, delta = dl, transform.fun = fisher_stat)$max.prob
  pf <- fisher.test(tbx, workspace = 2e7)$p.value
  worst <- max(worst, abs(pe - pf)); n_cases <- n_cases + 1
}
ok(worst < 1e-8, sprintf("(1d) Fisher's r x c exact test at gamma = 0 = fisher.test() on %d random tables (worst abs diff %.1e)", n_cases, worst))

## =============================================================================
## (2) Total probability mass = 1 for every (u, gamma)
## =============================================================================
cat("\n==== (2) P_u(.) is a probability distribution on the tables with the margins ====\n")
cases <- list(list(tb = matrix(c(4,3,1, 2,5,3, 1,2,6), 3, byrow = TRUE), delta = c(0,1,1)),
              list(tb = matrix(c(6,4, 5,5, 3,7, 2,8), 4, byrow = TRUE),   delta = c(0,0,1,1)),
              list(tb = matrix(c(5,2,4,3, 2,6,3,5), 2, byrow = TRUE),     delta = c(0,1)),
              list(tb = matrix(c(3,2,2, 1,4,2, 2,1,3, 2,2,3), 4, byrow = TRUE), delta = c(0,1,0,1)))   # nonmonotone delta
for (cs in cases) {
  tbs <- all_tables(cs$tb); m <- colSums(cs$tb)
  us <- rbind(rep(0, length(m)), m, sapply(m, function(x) sample(0:x, 1)), sapply(m, function(x) sample(0:x, 1)))
  for (gam in c(0, 0.7, 2.5)) for (k in seq_len(nrow(us))) {
    mass <- sum(sapply(tbs, function(t) point_prob(t, us[k, ], gam, cs$delta)))
    ok(abs(mass - 1) < 1e-10, sprintf("(2) %dx%d, |all tables| = %d, gamma = %.1f, u-bar = (%s), delta = (%s): total mass = %.12f",
                                       nrow(cs$tb), ncol(cs$tb), length(tbs), gam, paste(us[k, ], collapse = ","), paste(cs$delta, collapse = ","), mass))
  }
}

## =============================================================================
## (3) Exact moments (Theorem 4, brute force over all tables) = Theorem D.1 moments
##     used by the normal approximation
## =============================================================================
cat("\n==== (3) mean and variance of an ordinal T: exact distribution vs Theorem D.1 (normal approximation) ====\n")
mom_cases <- list(
  list(tb = matrix(c(6,3,1, 3,4,3, 1,3,6), 3, byrow = TRUE), delta = c(0,1,1), w = c(0,1,2),   v = c(0,1,2)),
  list(tb = matrix(c(12,18,5, 6,12,6, 6,6,15), 3, byrow = TRUE), delta = c(0,1,1), w = c(0,1,2.5), v = c(0,1,2)),   # F.2 table
  list(tb = matrix(c(6,4, 5,5, 3,7, 2,8), 4, byrow = TRUE), delta = c(0,1,1,1), w = c(1,2,3,4), v = c(0,1)),
  list(tb = matrix(c(5,2,4,3, 2,6,3,5), 2, byrow = TRUE), delta = c(0,1), w = c(0,1), v = c(1,2,3,4)))
for (cs in mom_cases) {
  tbs <- all_tables(cs$tb); m <- colSums(cs$tb); S <- outer(cs$w, cs$v); Tfun <- function(t) sum(S * t)
  Tvals <- sapply(tbs, function(t) Tfun(as.matrix(t)))
  us <- rbind(rep(0, length(m)), m, sapply(m, function(x) sample(0:x, 1)))
  for (gam in c(0, 1)) for (k in seq_len(nrow(us))) {
    p <- sapply(tbs, function(t) point_prob(t, us[k, ], gam, cs$delta))
    mean_exact <- sum(Tvals * p); var_exact <- sum(Tvals^2 * p) - mean_exact^2
    nm <- norm_single_u_allocation_p_value(obs.table = cs$tb, gamma = gam, delta = cs$delta, u_allocation = us[k, ],
                                           treatment.scores = cs$w, outcome.scores = cs$v)
    ok(abs(mean_exact - nm$mu_T) < 1e-8 && abs(var_exact - nm$var_T) < 1e-6,
       sprintf("(3) %dx%d gamma = %g u-bar = (%s): mean exact %.6f vs D.1 %.6f | var exact %.6f vs D.1 %.6f",
               nrow(cs$tb), ncol(cs$tb), gam, paste(us[k, ], collapse = ","), mean_exact, nm$mu_T, var_exact, nm$var_T))
  }
}
cat("\nDone.\n")

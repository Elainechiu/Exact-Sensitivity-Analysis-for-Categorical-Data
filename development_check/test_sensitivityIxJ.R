## =====================================================================================
## test_sensitivityIxJ.R -- unit tests for sensitivityIxJ.R / sensitivityIxJ.cpp (R2 engine)
##
## ONE file, meant to be run section by section (or block by block) in RStudio.
## Put it in the folder that holds sensitivityIxJ.R, sensitivityIxJ.cpp and (optionally)
## SIS_methods.cpp.  Requires testthat (3rd edition) and Rcpp.
##
##   Interactively : run Section 0 and Section 1 first, then any test_that() block.
##                   A passing block prints "Test passed"; a failing one stops with the
##                   failed expectation.
##   All at once   : testthat::test_file("test_sensitivityIxJ.R")   (report of every block)
##
## Switches (Section 0): RUN_SLOW (validity simulations, Section 13), RUN_KNOWN_ISSUES (Section 12,
## tests that state INTENDED behaviour and are expected to fail until fixed).
## Golden values (Section 11) are stored in sensitivityIxJ_reference_values.rds next to this file.
##
## Requirements and where they are tested
##   R1  Omega_t = tables with the observed margins and T >= c (ties included) ........ Section 2
##   R2  U_C, U_O, u+ are the sets of Theorems 1-3; chunks reproduce U_C ............... Section 3, 5
##   R3  P_u is a probability distribution for every (u-bar, gamma, delta) ............. Section 4
##   R4  gamma = 0 reduces to randomization inference (hypergeometric, Fisher) ......... Sections 4-6
##   R5  engine alpha = definition (brute force over assignments) ...................... Sections 4, 5, 7
##   R6  worst case over u in [0,1]^N attained in the candidate set; alpha depends on u-bar only . Sections 6, 7
##   R7  all input routes for one statistic agree ...................................... Section 5
##   R8  p in [0,1], max >= rct, nondecreasing in Gamma, invariances ................... Sections 5, 6
##   R9  known answers: 2x2 = noncentral hypergeometric; paper Table 3 = 0.021 ......... Section 6
##   R10 normal-approximation moments = exact moments ................................. Section 8
##   R11 invalid input -> informative error; empty Omega_t -> warning + NA ............. Section 10
##   R12 deterministic / reproducible; orientation-invariant .......................... Sections 5, 6, 9
##   R13 valid under hidden bias Gamma* when tested at Gamma >= Gamma* ................. Section 13
##   --  output unchanged from the verified R2 engine (characterization) .............. Section 11
## Section 14 has small debugging tools (not tests).
##
## Open flags found in review (see Section 12 for the regression tests)
##   F1 sampling.*.sen.IxJ, row = "outcome": table transposed twice; statistic applied in the wrong orientation
##   F2 default candidate sets differ by family: U_C (exact.general), U_O / u+ (exact.score), U_O path
##      (sampling), corners {0, N_.j}^J (normal); the last two are heuristics outside Theorem 2's hypotheses
##   F3 norm.score.sen.IxJ returns RCT.prob (capital) and no gamma / delta
##   F4 stale roxygen (u_space default, possible.u.allocation, gamma_delta =, mc.iteration defaults)
##   F5 sampling / normal functions do not use .check_common; single.u.sampling.general.IxJ validates nothing
##   F6 SIS uses seed = draw index, so set.seed() has no effect; no Monte Carlo SE returned
##   F7 exact.score.sen.IxJ accepts constant scores without a warning
##   F9 exact.score.sen.IxJ uses {(0,0), u+} whenever J = 2, for any I and treatment scores; Theorem 3 is
##      stated for the sign-score statistic.  Section 6, "J = 2, I > 2", checks it against the full U_C.
## =====================================================================================


## ==================================================================================
## 0. Settings and loading (run first)
## ==================================================================================

library(testthat)
library(Rcpp)
testthat::local_edition(3)

RUN_SLOW         <- FALSE    # TRUE: run Section 13 (simulation; a few minutes)
RUN_KNOWN_ISSUES <- FALSE    # TRUE: run Section 12 (expected to FAIL until the issues are fixed)

root <- getwd()              # the folder holding the package files (test_file() also runs from there)
stopifnot(file.exists(file.path(root, "sensitivityIxJ.R")), file.exists(file.path(root, "sensitivityIxJ.cpp")))
sourceCpp(file.path(root, "sensitivityIxJ.cpp"))
source(file.path(root, "sensitivityIxJ.R"))
HAS_SIS <- file.exists(file.path(root, "SIS_methods.cpp"))
if (HAS_SIS) sourceCpp(file.path(root, "SIS_methods.cpp"))


## ==================================================================================
## 1. Fixtures and independent reference implementations (run second)
## These do not call the exact engine (except all_tables / pp, which are tested in Sections 2 and 4).
## ==================================================================================

## ---- fixtures ------------------------------------------------------------------
tb_F1   <- matrix(c(2,3,0, 0,1,4, 0,1,4), 3, byrow = TRUE)    # paper F.1 block 1; N = 15; cols (2,5,8)
w_F1    <- c(0, 1, 2); v_F1 <- c(0, 1, 2)
T_F1    <- function(t) sum(outer(w_F1, v_F1) * t)
tb_2x2  <- matrix(c(9,3, 4,8), 2, byrow = TRUE)               # margins (12,12) x (13,11)
tb_tiny <- matrix(c(2,0, 1,1, 1,2), 3, byrow = TRUE)          # N = 7, 210 assignments
T_tiny  <- function(t) sum(outer(c(0, 1, 2), c(0, 1)) * t)
girls   <- matrix(c(12,3,0, 18,12,3, 17,25,4), 3, byrow = TRUE)   # Table 3 data, N = 94
boys    <- matrix(c(10,8,1, 29,11,3, 20,24,6), 3, byrow = TRUE)   # N = 112
tb_F2   <- matrix(c(12,18,5, 6,12,6, 6,6,15), 3, byrow = TRUE)    # F.2, N = 86
chisq_stat <- function(t) { E <- outer(rowSums(t), colSums(t)) / sum(t); sum((t - E)^2 / E) }
g2_stat    <- function(t) { E <- outer(rowSums(t), colSums(t)) / sum(t); 2 * sum(ifelse(t > 0, t * log(pmax(t, 1) / E), 0)) }
ord_stat   <- function(w, v) { S <- outer(w, v); function(t) sum(S * t) }

## ---- reference implementations --------------------------------------------------
## multivariate hypergeometric log-pmf: the fixed-margins permutation distribution
log_mvhyper <- function(t) sum(lgamma(rowSums(t) + 1)) + sum(lgamma(colSums(t) + 1)) -
                           lgamma(sum(t) + 1) - sum(lgamma(t + 1))
## every table with the margins of tb (uses possible.table with a constant statistic)
all_tables <- function(tb) possible.table(threshold = 0, table = tb, direction = "greater than",
                                          transform.fun = function(t) 1)
## shorthand for the point probability
pp <- function(t, u, g, d)
  generic.I.by.J.sensitivity.point.probability(table = as.matrix(t), u_allocation = u, gamma = g, delta = d)
## permutation p-value at gamma = 0, from the hypergeometric pmf (no package kernel involved)
perm_p <- function(tb, stat) {
  tabs <- all_tables(tb); Tobs <- stat(tb)
  sum(sapply(tabs, function(t) { t <- as.matrix(t); if (stat(t) >= Tobs - 1e-9) exp(log_mvhyper(t)) else 0 }))
}
## all distinct orderings of a multiset (rows)
uniq_perms <- function(v) {
  if (length(v) <= 1) return(matrix(v, nrow = 1))
  unname(do.call(rbind, lapply(unique(v), function(x) { i <- match(x, v); cbind(x, uniq_perms(v[-i])) })))
}
## u-bar -> a binary u for subjects sorted by outcome level (first u-bar_j of level j get u = 1)
ubar_to_u <- function(ub, m) unlist(lapply(seq_along(m), function(j) rep(c(1, 0), c(ub[j], m[j] - ub[j]))))
## alpha(T, r, u) straight from the definition, by enumerating every treatment assignment z
##   alpha = sum_z 1{T(N(z, r)) >= T_obs} w(z, u) / sum_z w(z, u),  w = exp(gamma * sum_s delta_{z_s} u_s).
## Feasible for N <= ~8.  Returns alpha at every row of U (default: the full U_C).
brute_alpha <- function(tb, delta, gamma, stat, U = NULL) {
  tb <- as.matrix(tb); I <- nrow(tb); J <- ncol(tb); n <- rowSums(tb); m <- colSums(tb)
  Z  <- uniq_perms(rep(seq_len(I), times = n))
  r  <- rep(seq_len(J), times = m)
  Tz <- apply(Z, 1, function(z) stat(unclass(table(factor(z, levels = seq_len(I)), factor(r, levels = seq_len(J))))))
  hit <- Tz >= stat(tb) - 1e-9
  D  <- matrix(delta[Z], nrow(Z))
  if (is.null(U)) U <- full.candidate.set(m)
  a <- apply(U, 1, function(ub) { w <- exp(gamma * drop(D %*% ubar_to_u(ub, m))); sum(w[hit]) / sum(w) })
  list(U = U, alpha = a, n_assign = nrow(Z), Z = Z, D = D, hit = hit, m = m)
}
## random I x J table with all margins positive
random_table <- function(I, J, N) {
  repeat {
    z <- sample.int(I, N, TRUE); r <- sample.int(J, N, TRUE)
    tb <- matrix(0L, I, J); for (s in seq_len(N)) tb[z[s], r[s]] <- tb[z[s], r[s]] + 1L
    if (all(rowSums(tb) > 0) && all(colSums(tb) > 0)) return(tb)
  }
}
random_delta <- function(I) repeat { d <- sample(0:1, I, TRUE); if (length(unique(d)) == 2) return(d) }

## ---- golden-master (characterization) values -------------------------------------------
## All reference values live in ONE file next to this script: sensitivityIxJ_reference_values.rds.
## First time a name is seen: the value is stored and the test skips.  Afterwards: compared.
## Create the references ONLY with the engine you trust (the verified R2 engine).
## To accept an intentional change: golden_forget("name"), rerun, and record why.
REF_FILE <- file.path(root, "sensitivityIxJ_reference_values.rds")
golden <- function(name, value, tolerance = 1e-10) {
  ref <- if (file.exists(REF_FILE)) readRDS(REF_FILE) else list()
  if (is.null(ref[[name]])) { ref[[name]] <- value; saveRDS(ref, REF_FILE); testthat::skip(paste("reference created:", name)) }
  testthat::expect_equal(value, ref[[name]], tolerance = tolerance)
}
golden_forget <- function(name) { ref <- readRDS(REF_FILE); ref[[name]] <- NULL; saveRDS(ref, REF_FILE) }


## ==================================================================================
## 2. possible.table: the reference set Omega_t  [R1]
## ==================================================================================

test_that("every table in Omega_t has the observed margins and T >= T_obs", {
  Om <- possible.table(T_F1(tb_F1), tb_F1, "greater than", T_F1)
  expect_gt(length(Om), 0)
  for (x in Om) {
    expect_equal(as.numeric(rowSums(x)), as.numeric(rowSums(tb_F1)))
    expect_equal(as.numeric(colSums(x)), as.numeric(colSums(tb_F1)))
    expect_gte(T_F1(x), T_F1(tb_F1) - 1e-9)
  }
})

test_that("the observed table (a tie at the threshold) is in Omega_t", {
  Om <- possible.table(T_F1(tb_F1), tb_F1, "greater than", T_F1)
  keys <- sapply(Om, function(x) paste(as.vector(x), collapse = ","))
  expect_true(paste(as.vector(tb_F1), collapse = ",") %in% keys)
})

test_that("Omega_t has no duplicate tables", {
  Om <- all_tables(tb_F1)
  keys <- sapply(Om, function(x) paste(as.vector(x), collapse = ","))
  expect_false(anyDuplicated(keys) > 0)
})

## ---- R1: known answers ----------------------------------------------------------
test_that("2 x 2 with margins (12,12) x (13,11) has exactly 12 tables (N11 = 1..12)", {
  expect_length(all_tables(tb_2x2), 12)
})

test_that("2 x J: count equals the number of compositions within the margins (brute count)", {
  tb <- matrix(c(2,1,3, 1,2,1), 2, byrow = TRUE)       # rows (6,4), cols (3,3,4)
  m <- colSums(tb); n1 <- rowSums(tb)[1]
  g <- expand.grid(lapply(m, function(x) 0:x))
  expect_length(all_tables(tb), sum(rowSums(g) == n1))
})

## ---- R1: direction ----------------------------------------------------------------
test_that("|>=| + |<=| - |ties| = |all tables|", {
  ge <- possible.table(T_F1(tb_F1), tb_F1, "greater than", T_F1)
  le <- possible.table(T_F1(tb_F1), tb_F1, "less than", T_F1)
  ties <- sum(abs(sapply(ge, T_F1) - T_F1(tb_F1)) < 1e-9)
  expect_equal(length(ge) + length(le) - ties, length(all_tables(tb_F1)))
  expect_true(all(sapply(le, T_F1) <= T_F1(tb_F1) + 1e-9))
})

test_that("ties with fractional scores are kept (tolerance 1e-9)", {
  S <- ord_stat(c(0, 0.3, 0.6), c(0, 0.7, 1.4))
  tb <- matrix(c(6,4,2, 4,4,4, 2,4,6), 3, byrow = TRUE)
  Om <- possible.table(S(tb), tb, "greater than", S)
  brute <- Filter(function(x) S(as.matrix(x)) >= S(tb) - 1e-9, all_tables(tb))
  expect_equal(length(Om), length(brute))
})

test_that("empty Omega_t returns list()", {
  expect_identical(possible.table(1e6, tb_F1, "greater than", T_F1), list())
})


## ==================================================================================
## 3. Candidate sets U_C, U_O, u+  [R2]
## ==================================================================================

test_that("full.candidate.set: |U_C| = prod(N_.j + 1), distinct rows, within margins, integer", {
  for (m in list(c(2, 5, 8), c(13, 11), c(1, 1, 1, 1), c(3, 0, 2))) {
    U <- full.candidate.set(m)
    expect_equal(nrow(U), prod(m + 1))
    expect_equal(ncol(U), length(m))
    expect_false(anyDuplicated(U) > 0)
    expect_true(all(sweep(U, 2, m, "<=")) && all(U >= 0))
    expect_type(U, "integer")
  }
})

test_that("full.candidate.set.rows reproduces full.candidate.set row by row, in the same order", {
  m <- c(2, 5, 8); U <- full.candidate.set(m)
  expect_equal(unname(full.candidate.set.rows(m, seq_len(nrow(U)))), unname(U), ignore_attr = TRUE)
  expect_equal(unname(full.candidate.set.rows(m, c(162, 1))), unname(U[c(162, 1), ]), ignore_attr = TRUE)
})

test_that("full.candidate.set refuses to materialize a U_C larger than max_rows", {
  expect_error(full.candidate.set(c(100, 100, 100, 100), max_rows = 1e6), "exceeds max_rows", fixed = TRUE)
})

test_that(".u_space_ordinal is U_O of Theorem 2: N + 1 rows, nested, fills level J first", {
  m <- c(2, 5, 8); U <- .u_space_ordinal(m)
  expect_equal(nrow(U), sum(m) + 1)
  expect_equal(as.numeric(rowSums(U)), 0:sum(m))
  expect_true(all(U[1, ] == 0)); expect_true(all(U[nrow(U), ] == m))
  expect_true(all(diff(U) >= 0))                                  # each row adds one subject
  expect_true(all(U[U[, 2] > 0, 3] == m[3]))                      # u2 moves only after u3 is full
  expect_true(all(U[U[, 1] > 0, 2] == m[2]))
})

test_that(".u_space_sign_score is {(0,0), u+ = (0, N_.2)} (Theorem 3)", {
  U <- .u_space_sign_score(c(13, 11))
  expect_equal(unname(U), matrix(c(0L, 0L, 0L, 11L), 2, byrow = TRUE), ignore_attr = TRUE)
})


## ==================================================================================
## 4. Point probability  [R3, R4, R5]
## ==================================================================================

test_that("R4: gamma = 0 gives the multivariate hypergeometric pmf for any u (random tables, 7 shapes)", {
  set.seed(101)
  for (d in list(c(2,2), c(2,3), c(3,3), c(4,2), c(3,4), c(5,2), c(2,5))) for (k in 1:5) {
    tb <- random_table(d[1], d[2], sample(8:30, 1)); delta <- random_delta(d[1])
    u  <- sapply(colSums(tb), function(x) sample(0:x, 1))
    expect_equal(pp(tb, u, 0, delta), exp(log_mvhyper(tb)), tolerance = 1e-10)
  }
})

test_that("R3: total mass over all tables = 1 for every (u, gamma), monotone and nonmonotone delta", {
  set.seed(102)
  cases <- list(list(tb = tb_F1, d = c(0,1,1)),
                list(tb = matrix(c(6,4, 5,5, 3,7, 2,8), 4, byrow = TRUE), d = c(0,0,1,1)),
                list(tb = matrix(c(5,2,4,3, 2,6,3,5), 2, byrow = TRUE), d = c(0,1)),
                list(tb = matrix(c(3,2,2, 1,4,2, 2,1,3, 2,2,3), 4, byrow = TRUE), d = c(0,1,0,1)))
  for (cs in cases) {
    tabs <- all_tables(cs$tb); m <- colSums(cs$tb)
    us <- rbind(0 * m, m, sapply(m, function(x) sample(0:x, 1)))
    for (g in c(0, 0.7, 2.5)) for (k in seq_len(nrow(us)))
      expect_equal(sum(sapply(tabs, pp, u = us[k, ], g = g, d = cs$d)), 1, tolerance = 1e-10)
  }
})

test_that("point probability is in [0, 1] and orientation-invariant", {
  p1 <- pp(tb_F1, c(1, 3, 4), 1, c(0,1,1))
  p2 <- generic.I.by.J.sensitivity.point.probability(table = t(tb_F1), row = "outcome",
                                                      u_allocation = c(1, 3, 4), gamma = 1, delta = c(0,1,1))
  expect_true(p1 >= 0 && p1 <= 1)
  expect_equal(p1, p2, tolerance = 1e-12)
})

test_that("Theorem 4 kernel matches brute force over assignments, per table (tiny table)", {
  b <- brute_alpha(tb_tiny, c(0,1,1), 1.2, function(t) 0, U = matrix(c(2, 1), 1))  # stat irrelevant here
  ## P_u(N = t) by brute force = weight of assignments producing t / total weight
  r <- rep(1:2, times = colSums(tb_tiny)); u <- ubar_to_u(c(2, 1), colSums(tb_tiny))
  w <- exp(1.2 * drop(b$D %*% u))
  key <- apply(b$Z, 1, function(z) paste(as.vector(unclass(table(factor(z, 1:3), factor(r, 1:2)))), collapse = ","))
  for (t in all_tables(tb_tiny)) {
    k <- paste(as.vector(as.matrix(t)), collapse = ",")
    expect_equal(pp(t, c(2, 1), 1.2, c(0,1,1)), sum(w[key == k]) / sum(w), tolerance = 1e-10)
  }
})


## ==================================================================================
## 5. exact.general.sen.IxJ  [R4, R5, R7, R8, R12]
## ==================================================================================

eg <- function(...) exact.general.sen.IxJ(...)

## ---- output contract ------------------------------------------------------------
test_that("output has the documented elements, types and ranges", {
  r <- eg(obs.table = tb_F1, gamma = 1, delta = c(0,1,1), statistic = "chisq")
  expect_named(r, c("rct.prob", "max.prob", "maximizer", "gamma", "delta", "n.tables",
                    "n.distinct.tbar", "n.candidates", "all.prob", "time"), ignore.order = TRUE)
  expect_equal(r$n.candidates, prod(colSums(tb_F1) + 1))
  expect_length(r$all.prob, r$n.candidates)
  expect_true(all(r$all.prob >= 0 & r$all.prob <= 1 + 1e-12))
  expect_equal(r$max.prob, max(r$all.prob))
  expect_length(r$maximizer, ncol(tb_F1))
  expect_gte(r$max.prob, r$rct.prob - 1e-12)                       # u = 0 is a candidate
  expect_named(r$time, c("enumerate.W", "sweep"))
})

## ---- R7: every route to the same statistic agrees -------------------------------
test_that("chi-square: C++ built-in = R callback = explicit table_space", {
  a <- eg(obs.table = tb_F1, gamma = 1, delta = c(0,1,1), statistic = "chisq")
  b <- eg(obs.table = tb_F1, gamma = 1, delta = c(0,1,1), transform.fun = chisq_stat)
  cc <- eg(obs.table = tb_F1, gamma = 1, delta = c(0,1,1),
          table_space = possible.table(chisq_stat(tb_F1), tb_F1, "greater than", chisq_stat))
  expect_equal(a$all.prob, b$all.prob, tolerance = 1e-10); expect_equal(a$all.prob, cc$all.prob, tolerance = 1e-10)
  expect_equal(a$n.tables, b$n.tables); expect_equal(a$n.tables, cc$n.tables)
})

test_that("G^2: C++ built-in = R callback", {
  tb <- matrix(c(6,4,2, 4,4,4, 2,4,6), 3, byrow = TRUE)
  a <- eg(obs.table = tb, gamma = 1, delta = c(0,1,1), statistic = "G2")
  b <- eg(obs.table = tb, gamma = 1, delta = c(0,1,1), transform.fun = g2_stat)
  expect_equal(a$all.prob, b$all.prob, tolerance = 1e-10); expect_equal(a$n.tables, b$n.tables)
})

test_that("ordinal: scores path (C++) = transform.fun path (R)", {
  a <- eg(obs.table = tb_F1, gamma = 1, delta = c(0,1,1), treatment.scores = w_F1, outcome.scores = v_F1)
  b <- eg(obs.table = tb_F1, gamma = 1, delta = c(0,1,1), transform.fun = T_F1)
  expect_equal(a$all.prob, b$all.prob, tolerance = 1e-10)
})

test_that("precedence: table_space wins over transform.fun", {
  ts <- possible.table(T_F1(tb_F1), tb_F1, "greater than", T_F1)
  a <- eg(obs.table = tb_F1, gamma = 1, delta = c(0,1,1), table_space = ts, transform.fun = chisq_stat)
  b <- eg(obs.table = tb_F1, gamma = 1, delta = c(0,1,1), transform.fun = T_F1)
  expect_equal(a$max.prob, b$max.prob, tolerance = 1e-10)
})

## ---- u_space options --------------------------------------------------------------
test_that("u_space = NULL and 'full' are the same; a single row matches all.prob; chunks give the same max", {
  m <- colSums(tb_F1); U <- full.candidate.set(m)
  a <- eg(obs.table = tb_F1, gamma = 1, delta = c(0,1,1), treatment.scores = w_F1, outcome.scores = v_F1)
  b <- eg(u_space = "full", obs.table = tb_F1, gamma = 1, delta = c(0,1,1), treatment.scores = w_F1, outcome.scores = v_F1)
  expect_identical(a$all.prob, b$all.prob)
  k <- 77
  one <- eg(u_space = U[k, , drop = FALSE], obs.table = tb_F1, gamma = 1, delta = c(0,1,1),
            treatment.scores = w_F1, outcome.scores = v_F1)
  expect_equal(one$max.prob, a$all.prob[k], tolerance = 1e-12)
  ch <- split(seq_len(nrow(U)), ceiling(seq_len(nrow(U)) / 40))
  mx <- sapply(ch, function(i) eg(u_space = full.candidate.set.rows(m, i), obs.table = tb_F1, gamma = 1,
                                  delta = c(0,1,1), treatment.scores = w_F1, outcome.scores = v_F1)$max.prob)
  expect_equal(max(mx), a$max.prob, tolerance = 1e-12)
})

## ---- R4: reductions at gamma = 0 --------------------------------------------------
test_that("gamma = 0: alpha constant over U_C and = permutation p-value (chi-square, ordinal)", {
  for (s in list(chisq_stat, T_F1)) {
    r <- eg(obs.table = tb_F1, gamma = 0, delta = c(0,1,1), transform.fun = s)
    expect_lt(diff(range(r$all.prob)), 1e-12)
    expect_equal(r$max.prob, perm_p(tb_F1, s), tolerance = 1e-10)
    expect_equal(r$max.prob, r$rct.prob, tolerance = 1e-12)
  }
})

test_that("Freeman-Halton at gamma = 0 = fisher.test() r x c (random tables, 5 shapes)", {
  set.seed(103)
  for (d in list(c(2,3), c(3,3), c(2,4), c(4,2), c(3,4))) for (k in 1:3) {
    tb <- random_table(d[1], d[2], sample(10:20, 1))
    p <- eg(u_space = matrix(0L, 1, d[2]), obs.table = tb, gamma = 0, delta = random_delta(d[1]),
            statistic = "freeman_halton")$max.prob
    expect_equal(p, fisher.test(tb, workspace = 2e7)$p.value, tolerance = 1e-8)
  }
})

## ---- R5: alpha = sum of point probabilities over Omega_t ----------------------------
test_that("alpha(u) = sum over Omega_t of P_u(N = t) (links stage 1/2 to the point probability)", {
  Om <- possible.table(T_F1(tb_F1), tb_F1, "greater than", T_F1)
  for (u in list(c(0,0,0), c(1,3,4), c(2,5,8), c(2,0,1))) {
    a <- eg(u_space = matrix(u, 1), obs.table = tb_F1, gamma = 0.8, delta = c(0,1,1), transform.fun = T_F1)$max.prob
    expect_equal(a, sum(sapply(Om, pp, u = u, g = 0.8, d = c(0,1,1))), tolerance = 1e-10)
  }
})

## ---- R8: structural properties ------------------------------------------------------
test_that("worst-case p nondecreasing in gamma; rct.prob does not depend on gamma (chi-square, nonmonotone delta)", {
  gs <- c(0, 0.25, 0.5, 1, 2, 4)
  r  <- lapply(gs, function(g) eg(obs.table = tb_F1, gamma = g, delta = c(0,1,0), statistic = "chisq"))
  mx <- sapply(r, `[[`, "max.prob"); rc <- sapply(r, `[[`, "rct.prob")
  expect_true(all(diff(mx) >= -1e-12)); expect_lt(diff(range(rc)), 1e-12)
})

test_that("invariance: permuting outcome levels (with u-bar permuted) leaves chi-square results unchanged", {
  perm <- c(3, 1, 2); m <- colSums(tb_F1)
  a <- eg(obs.table = tb_F1, gamma = 1, delta = c(0,1,1), statistic = "chisq")
  b <- eg(obs.table = tb_F1[, perm], gamma = 1, delta = c(0,1,1), statistic = "chisq")
  key <- function(U) apply(U, 1, paste, collapse = ",")
  Ua <- full.candidate.set(m); Ub <- full.candidate.set(m[perm])
  expect_equal(b$all.prob[match(key(Ua[, perm, drop = FALSE]), key(Ub))], a$all.prob, tolerance = 1e-10)
})

test_that("invariance: permuting treatment levels together with delta leaves chi-square results unchanged", {
  perm <- c(2, 3, 1); d <- c(0, 1, 1)
  a <- eg(obs.table = tb_F1, gamma = 1, delta = d, statistic = "chisq")
  b <- eg(obs.table = tb_F1[perm, ], gamma = 1, delta = d[perm], statistic = "chisq")
  expect_equal(b$all.prob, a$all.prob, tolerance = 1e-10)
})

test_that("relabelling: alpha under delta at u-bar = alpha under 1 - delta at N_.j - u-bar", {
  ## Model (6): exp(xi_i + gamma (1 - delta_i) u) is proportional (over i) to
  ## exp(xi_i - gamma delta_i + gamma delta_i (1 - u)); the constant -gamma delta_i cancels with fixed
  ## treatment margins.  So delta -> 1 - delta is the same as u -> 1 - u, i.e. u-bar -> N_.j - u-bar.
  m <- colSums(tb_F1); U <- full.candidate.set(m)
  a <- eg(obs.table = tb_F1, gamma = 1, delta = c(0,1,1), statistic = "chisq")
  b <- eg(obs.table = tb_F1, gamma = 1, delta = c(1,0,0), statistic = "chisq")
  key <- function(V) apply(V, 1, paste, collapse = ",")
  flip <- sweep(-U, 2, m, "+")                                     # N_.j - u-bar
  expect_equal(b$all.prob[match(key(flip), key(U))], a$all.prob, tolerance = 1e-10)
  expect_equal(b$max.prob, a$max.prob, tolerance = 1e-10)
})

## ---- R12 ---------------------------------------------------------------------------
test_that("deterministic; matrix, table object and transposed input with row = 'outcome' agree", {
  a <- eg(obs.table = tb_F1, gamma = 1, delta = c(0,1,1), statistic = "chisq")
  expect_identical(a$all.prob, eg(obs.table = tb_F1, gamma = 1, delta = c(0,1,1), statistic = "chisq")$all.prob)
  expect_equal(eg(obs.table = as.table(tb_F1), gamma = 1, delta = c(0,1,1), statistic = "chisq")$max.prob, a$max.prob)
  b <- eg(obs.table = t(tb_F1), row = "outcome", gamma = 1, delta = c(0,1,1), transform.fun = function(x) chisq_stat(t(x)))
  expect_equal(b$max.prob, a$max.prob, tolerance = 1e-10)
})

test_that("numerical stability: very large gamma gives a finite p in [0, 1]", {
  r <- eg(obs.table = tb_F1, gamma = 50, delta = c(0,1,0), statistic = "chisq")
  expect_true(is.finite(r$max.prob) && r$max.prob >= 0 && r$max.prob <= 1 + 1e-12)
})


## ==================================================================================
## 6. exact.score.sen.IxJ  [R4, R6, R8, R9, R12]
## ==================================================================================

es <- function(tb, g, d, w, v, ...) exact.score.sen.IxJ(obs.table = tb, gamma = g, delta = d,
                                                         treatment.scores = w, outcome.scores = v, ...)

test_that("output contract: obs.stat, obs.table, |U_O| = N + 1 candidates for J >= 3", {
  r <- es(tb_F1, 1, c(0,1,1), w_F1, v_F1)
  expect_equal(r$obs.stat, T_F1(tb_F1))
  expect_equal(r$n.candidates, sum(tb_F1) + 1)
  expect_equal(r$n.tables, length(possible.table(T_F1(tb_F1), tb_F1, "greater than", T_F1)))
  expect_gte(r$max.prob, r$rct.prob - 1e-12)
})

## ---- R6: Theorems 2 and 3 --------------------------------------------------------
test_that("Theorem 2: max over U_O = max over the full U_C (several tables, scores, gammas)", {
  cases <- list(list(tb = tb_F1, d = c(0,1,1), w = w_F1, v = v_F1),
                list(tb = tb_F1, d = c(0,0,1), w = c(0, 1, 2.5), v = c(0, 1, 2)),
                list(tb = matrix(c(5,2,4,3, 2,6,3,5), 2, byrow = TRUE), d = c(0,1), w = c(0,1), v = c(1,2,3,4)))
  for (cs in cases) for (g in c(0.4, 1.5)) {
    a <- es(cs$tb, g, cs$d, cs$w, cs$v)
    b <- exact.general.sen.IxJ(obs.table = cs$tb, gamma = g, delta = cs$d, treatment.scores = cs$w, outcome.scores = cs$v)
    expect_equal(a$max.prob, b$max.prob, tolerance = 1e-10)
  }
})

test_that("Theorem 2 on random 3 x 3 tables with random monotone scores", {
  set.seed(104)
  for (k in 1:8) {
    tb <- random_table(3, 3, sample(9:16, 1))
    w <- sort(round(runif(3, 0, 2), 2)); v <- sort(round(runif(3, 0, 2), 2)); g <- runif(1, 0, 2)
    a <- es(tb, g, c(0,1,1), w, v)
    b <- exact.general.sen.IxJ(obs.table = tb, gamma = g, delta = c(0,1,1), treatment.scores = w, outcome.scores = v)
    expect_equal(a$max.prob, b$max.prob, tolerance = 1e-10)
  }
})

test_that("J = 2, I > 2: the code searches only {(0,0), u+}; check against the full U_C", {
  ## exact.score.sen.IxJ uses .u_space_sign_score whenever J = 2, whatever I and the treatment scores.
  ## Theorem 3 is stated for the sign-score statistic; if this test fails, the J = 2 branch is using
  ## Theorem 3 outside its hypotheses (flag F9, see the file header).
  tb <- matrix(c(6,4, 5,5, 3,7, 2,8), 4, byrow = TRUE)
  for (g in c(0.4, 1.5)) {
    a <- es(tb, g, c(0,1,1,1), 1:4, c(0,1))
    b <- exact.general.sen.IxJ(obs.table = tb, gamma = g, delta = c(0,1,1,1), treatment.scores = 1:4, outcome.scores = c(0,1))
    expect_equal(a$max.prob, b$max.prob, tolerance = 1e-10)
  }
})

## ---- R9: known answer, 2 x 2 -----------------------------------------------------
test_that("2 x 2 worst case = Fisher's noncentral hypergeometric: fisher.test(or = Gamma, 'greater')", {
  ## Under u+ = (0, N_.2), P(N_22 = k) is proportional to C(N_.2,k) C(N_.1, N_2. - k) e^{gamma k},
  ## Fisher's noncentral hypergeometric with odds ratio Gamma = e^gamma.
  for (tb in list(tb_2x2, matrix(c(12,18, 17,3), 2, byrow = TRUE), matrix(c(3,7, 8,2), 2, byrow = TRUE)))
    for (G in c(1, 1.5, 2, 4)) {
      r <- es(tb, log(G), c(0,1), c(0,1), c(0,1))
      expect_equal(r$max.prob, fisher.test(tb, or = G, alternative = "greater")$p.value, tolerance = 1e-10)
      if (r$max.prob > r$rct.prob + 1e-12) expect_equal(r$maximizer, c(0, colSums(tb)[[2]]))
    }
})

test_that("Theorem 3: u+ attains the max over the full U_C in 2 x 2", {
  r <- exact.general.sen.IxJ(obs.table = tb_2x2, gamma = log(2), delta = c(0,1), treatment.scores = c(0,1), outcome.scores = c(0,1))
  expect_equal(r$max.prob, fisher.test(tb_2x2, or = 2, alternative = "greater")$p.value, tolerance = 1e-10)
})

## ---- R9: published number ---------------------------------------------------------
test_that("paper Table 3: girls, profile scores, Gamma = 2 gives 0.021 (to the printed 3 decimals)", {
  r <- es(girls, log(2), c(0,1,1), c(0, 0.2, 1.6), c(0, 0.67, 0.7))
  ## scores here are the 2-decimal rounding of the recorded (0, 0.197, 1.598), (0, 0.668, 0.698)
  expect_equal(round(r$max.prob, 3), 0.021)
})

## ---- R4, R8 ------------------------------------------------------------------------
test_that("gamma = 0: max.prob = rct.prob = permutation p-value", {
  r <- es(tb_F1, 0, c(0,1,1), w_F1, v_F1)
  expect_equal(r$max.prob, r$rct.prob, tolerance = 1e-12)
  expect_equal(r$max.prob, perm_p(tb_F1, T_F1), tolerance = 1e-10)
})

test_that("sensitivity curve is nondecreasing in Gamma (no claim about the limit: it need not be 1 here)", {
  G <- c(1, 1.25, 1.5, 2, 3, 5, 20, 1e4)
  p <- sapply(log(G), function(g) es(tb_F1, g, c(0,1,1), w_F1, v_F1)$max.prob)
  expect_true(all(diff(p) >= -1e-12))
  expect_true(all(is.finite(p)) && all(p <= 1 + 1e-12))
})

test_that("orientation: row = 'outcome' with t(table) gives the same answer", {
  a <- es(tb_F1, 1, c(0,1,1), w_F1, v_F1)
  b <- es(t(tb_F1), 1, c(0,1,1), w_F1, v_F1, row = "outcome")
  expect_equal(b$max.prob, a$max.prob, tolerance = 1e-12)
  expect_equal(b$maximizer, a$maximizer)
})


## ==================================================================================
## 7. Brute force over treatment assignments  [R5, R6]
## ==================================================================================

test_that("R5: engine alpha = brute force at every u-bar in U_C (tiny 3 x 2, ordinal)", {
  b <- brute_alpha(tb_tiny, c(0,1,1), 1.2, T_tiny)
  e <- exact.general.sen.IxJ(obs.table = tb_tiny, gamma = 1.2, delta = c(0,1,1), transform.fun = T_tiny)
  expect_equal(b$n_assign, 210)
  expect_equal(e$all.prob, b$alpha, tolerance = 1e-10)
})

test_that("R5: engine = brute force on random tiny tables (I, J in {2,3}, N <= 7, random delta, gamma, statistic)", {
  set.seed(105)
  for (k in 1:16) {
    I <- sample(2:3, 1); J <- sample(2:3, 1); N <- sample(5:7, 1)
    tb <- random_table(I, J, N); d <- random_delta(I); g <- runif(1, 0, 2)
    stat <- if (k %% 2 == 0) chisq_stat else ord_stat(sort(runif(I)), sort(runif(J)))
    b <- brute_alpha(tb, d, g, stat)
    e <- exact.general.sen.IxJ(obs.table = tb, gamma = g, delta = d, transform.fun = stat)
    expect_equal(e$all.prob, b$alpha, tolerance = 1e-10,
                 label = sprintf("case %d (%dx%d, N=%d, delta=%s)", k, I, J, N, paste(d, collapse = "")))
  }
})

test_that("R6 (Theorem 1): max over ALL binary u (2^N, brute force) = engine max over U_C", {
  b0 <- brute_alpha(tb_tiny, c(0,1,1), 1.2, T_tiny)
  allu <- as.matrix(expand.grid(rep(list(0:1), sum(tb_tiny))))
  a <- apply(allu, 1, function(u) { w <- exp(1.2 * drop(b0$D %*% u)); sum(w[b0$hit]) / sum(w) })
  e <- exact.general.sen.IxJ(obs.table = tb_tiny, gamma = 1.2, delta = c(0,1,1), transform.fun = T_tiny)
  expect_equal(max(a), e$max.prob, tolerance = 1e-10)
})

test_that("R6 (Theorem 1): alpha depends on binary u only through u-bar", {
  b0 <- brute_alpha(tb_tiny, c(0,1,1), 1.2, T_tiny)
  r  <- rep(seq_len(ncol(tb_tiny)), times = colSums(tb_tiny))
  allu <- as.matrix(expand.grid(rep(list(0:1), sum(tb_tiny))))
  a <- apply(allu, 1, function(u) { w <- exp(1.2 * drop(b0$D %*% u)); sum(w[b0$hit]) / sum(w) })
  grp <- apply(allu, 1, function(u) paste(tapply(u, r, sum), collapse = ","))
  expect_true(all(tapply(a, grp, function(x) diff(range(x))) < 1e-12))
})

test_that("R6 (Theorem 1, continuous u): random u in [0,1]^N never beats the max over U_C", {
  set.seed(106)
  b0 <- brute_alpha(tb_tiny, c(0,1,0), 1.5, chisq_stat)
  e  <- exact.general.sen.IxJ(obs.table = tb_tiny, gamma = 1.5, delta = c(0,1,0), transform.fun = chisq_stat)
  a  <- replicate(300, { u <- runif(sum(tb_tiny)); w <- exp(1.5 * drop(b0$D %*% u)); sum(w[b0$hit]) / sum(w) })
  expect_lte(max(a), e$max.prob + 1e-12)
})


## ==================================================================================
## 8. Normal approximation  [R10]
## ==================================================================================

test_that("Theorem D.1 mean/variance = exact moments from Theorem 4 (all supported shapes except 2x5)", {
  cases <- list(
    list(tb = tb_F1, d = c(0,1,1), w = w_F1, v = v_F1),                                                   # 3x3
    list(tb = matrix(c(6,4, 5,5, 3,7, 2,8), 4, byrow = TRUE), d = c(0,1,1,1), w = c(1,2,3,4), v = c(0,1)),  # 4x2
    list(tb = matrix(c(5,2,4,3, 2,6,3,5), 2, byrow = TRUE), d = c(0,1), w = c(0,1), v = c(1,2,3,4)),         # 2x4
    list(tb = matrix(c(4,6,0, 1,3,6), 2, byrow = TRUE), d = c(0,1), w = c(0,1), v = c(0,1,2)),               # 2x3
    list(tb = matrix(c(4,3, 2,5, 3,3), 3, byrow = TRUE), d = c(0,0,1), w = c(0,1,2), v = c(0,1)),            # 3x2
    list(tb = tb_2x2, d = c(0,1), w = c(0,1), v = c(0,1)))                                                    # 2x2
  set.seed(107)
  for (cs in cases) {
    tabs <- all_tables(cs$tb); S <- outer(cs$w, cs$v); Tv <- sapply(tabs, function(t) sum(S * as.matrix(t)))
    m <- colSums(cs$tb)
    us <- rbind(0 * m, m, sapply(m, function(x) sample(0:x, 1)))
    for (g in c(0, 1)) for (k in seq_len(nrow(us))) {
      p <- sapply(tabs, pp, u = us[k, ], g = g, d = cs$d)
      mu <- sum(Tv * p); vr <- sum(Tv^2 * p) - mu^2
      nm <- norm_single_u_allocation_p_value(obs.table = cs$tb, gamma = g, delta = cs$d, u_allocation = us[k, ],
                                             treatment.scores = cs$w, outcome.scores = cs$v)
      lab <- sprintf("%dx%d g=%g u=(%s)", nrow(cs$tb), ncol(cs$tb), g, paste(us[k, ], collapse = ","))
      expect_equal(nm$mu_T, mu, tolerance = 1e-8, label = paste(lab, "mean"))
      expect_equal(nm$var_T, vr, tolerance = 1e-6, label = paste(lab, "var"))
    }
  }
})

test_that("normal p-value is 1 - pnorm(z) with z = (T_obs - mu) / sd", {
  nm <- norm_single_u_allocation_p_value(obs.table = tb_F1, gamma = 1, delta = c(0,1,1), u_allocation = c(1,3,4),
                                         treatment.scores = w_F1, outcome.scores = v_F1)
  expect_equal(nm$T_obs, T_F1(tb_F1))
  expect_equal(nm$z_score, (nm$T_obs - nm$mu_T) / sqrt(nm$var_T))
  expect_equal(nm$p_value, pnorm(nm$z_score, lower.tail = FALSE), tolerance = 1e-12)
})

test_that("norm.score.sen.IxJ: max over its candidate set, >= its u = 0 value, close to exact for N = 94", {
  ns <- norm.score.sen.IxJ(obs.table = girls, gamma = log(2), delta = c(0,1,1),
                           treatment.scores = c(0, 0.2, 1.6), outcome.scores = c(0, 0.67, 0.7))
  ex <- exact.score.sen.IxJ(obs.table = girls, gamma = log(2), delta = c(0,1,1),
                            treatment.scores = c(0, 0.2, 1.6), outcome.scores = c(0, 0.67, 0.7))
  expect_gte(ns$max.prob, ns$RCT.prob - 1e-12)
  expect_lt(abs(ns$max.prob - ex$max.prob), 0.05)   # approximation quality, not a theorem; a failure is a finding, not a bug
})

test_that("unsupported dimension (3 x 4) gives an error", {
  tb <- matrix(c(2,1,1,2, 1,2,2,1, 2,2,1,1), 3, byrow = TRUE)
  expect_error(norm_single_u_allocation_p_value(obs.table = tb, gamma = 1, delta = c(0,1,1), u_allocation = c(0,0,0,0),
                                                treatment.scores = c(0,1,2), outcome.scores = c(0,1,2,3)),
               "not supported", fixed = TRUE)
})


## ==================================================================================
## 9. Sampling (SIS)  [R12]
## ==================================================================================

test_that("SIS estimate at one u is close to the exact alpha (Monte Carlo, 4000 draws)", {
  skip_if_not(HAS_SIS, "SIS_methods.cpp not found")
  ex <- exact.general.sen.IxJ(u_space = matrix(c(1, 3, 4), 1), obs.table = tb_F1, gamma = 1, delta = c(0,1,1),
                              transform.fun = T_F1)$max.prob
  sm <- single.u.sampling.general.IxJ(table = tb_F1, mc.iteration = 4000, u_allocation = c(1, 3, 4),
                                      gamma = 1, delta = c(0,1,1), transform.fun = T_F1)
  expect_length(sm$probability.progression, 4000)
  expect_lt(abs(sm$final.prob.estimate - ex), 0.03)
})

test_that("SIS at gamma = 0 estimates the permutation p-value", {
  skip_if_not(HAS_SIS, "SIS_methods.cpp not found")
  sm <- single.u.sampling.general.IxJ(table = tb_F1, mc.iteration = 4000, u_allocation = c(0, 0, 0),
                                      gamma = 0, delta = c(0,1,1), transform.fun = T_F1)
  expect_lt(abs(sm$final.prob.estimate - perm_p(tb_F1, T_F1)), 0.03)
})

test_that("sampling is reproducible: identical calls give identical output (seed = draw index)", {
  skip_if_not(HAS_SIS, "SIS_methods.cpp not found")
  f <- function() sampling.score.sen.IxJ(obs.table = tb_F1, gamma = 1, delta = c(0,1,1), mc.iteration = 300,
                                         treatment.scores = w_F1, outcome.scores = v_F1)
  set.seed(1); a <- f(); set.seed(999); b <- f()
  expect_identical(a, b)          # documents current behaviour: set.seed() has no effect (flag F6, see the file header)
})

test_that("sampling.general.sen.IxJ with an explicit u_space returns a maximizer from that set", {
  skip_if_not(HAS_SIS, "SIS_methods.cpp not found")
  U <- rbind(c(0,0,8), c(0,5,8), c(2,5,8))
  r <- sampling.general.sen.IxJ(obs.table = tb_F1, gamma = 1, delta = c(0,1,1), mc.iteration = 300,
                                transform.fun = T_F1, u_space = U)
  expect_true(any(apply(U, 1, function(x) all(x == r$maximizer))))
  expect_true(r$max.prob >= 0 && r$max.prob <= 1)
})


## ==================================================================================
## 10. Invalid input and edge cases  [R11]
## Error messages are matched (fixed = TRUE) against the current wording in sensitivityIxJ.R.
## ==================================================================================

es_bad <- function(tb = tb_F1, g = 1, d = c(0,1,1), w = w_F1, v = v_F1, ...)
  exact.score.sen.IxJ(obs.table = tb, gamma = g, delta = d, treatment.scores = w, outcome.scores = v, ...)
eg_bad <- function(tb = tb_F1, g = 1, d = c(0,1,1), ...)
  exact.general.sen.IxJ(obs.table = tb, gamma = g, delta = d, ...)

test_that("gamma: negative or non-scalar", {
  expect_error(es_bad(g = -1), "gamma must be nonnegative", fixed = TRUE)
  expect_error(es_bad(g = c(1, 2)), "gamma must be a scalar", fixed = TRUE)
})

test_that("delta: non-binary, constant, wrong length", {
  expect_error(es_bad(d = c(0,1,2)), "delta must be a binary vector", fixed = TRUE)
  expect_error(es_bad(d = c(1,3,3)), "delta must be a binary vector", fixed = TRUE)
  expect_error(es_bad(d = c(1,1,1)), "delta must not be constant", fixed = TRUE)
  expect_error(es_bad(d = c(0,1)), "'delta' must be a numeric vector of length 3", fixed = TRUE)
})

test_that("table: zero margin, data.frame, 3-d array, bad 'row'", {
  expect_error(es_bad(tb = matrix(c(5,0, 3,0), 2, byrow = TRUE), d = c(0,1), w = c(0,1), v = c(0,1)),
               "zero row or column margin", fixed = TRUE)
  expect_error(eg_bad(tb = as.data.frame(tb_F1), statistic = "chisq"), "must be either a matrix or a table object", fixed = TRUE)
  expect_error(eg_bad(tb = as.table(array(1, c(2,2,2))), d = c(0,1), statistic = "chisq"), "must be two-dimensional", fixed = TRUE)
  expect_error(es_bad(row = "rows"), "'row' must be either 'treatment' or 'outcome'", fixed = TRUE)
})

test_that("exact.score: hypotheses of Theorems 2-3 enforced (monotone delta and scores, lengths)", {
  expect_error(es_bad(d = c(0,1,0)), "'delta' must have a monotone trend", fixed = TRUE)
  expect_error(es_bad(w = c(0,2,1)), "'treatment.scores' must have a monotone trend", fixed = TRUE)
  expect_error(es_bad(v = c(2,1,0)), "'outcome.scores' must have a monotone trend", fixed = TRUE)
  expect_error(es_bad(w = c(0,1)), "'treatment.scores' must be a numeric vector of length 3", fixed = TRUE)
})

test_that("exact.general: u_space and table_space checks; no statistic supplied", {
  expect_error(eg_bad(u_space = matrix(c(3,0,0), 1), statistic = "chisq"), "must not exceed the corresponding outcome margin", fixed = TRUE)
  expect_error(eg_bad(u_space = matrix(0, 1, 2), statistic = "chisq"), "Number of columns in 'u_space'", fixed = TRUE)
  expect_error(eg_bad(u_space = matrix(-1, 1, 3), statistic = "chisq"), "must be non-negative", fixed = TRUE)
  expect_error(eg_bad(u_space = c(0, 0, 0), statistic = "chisq"), "'u_space' must be a numeric matrix", fixed = TRUE)
  wrong <- matrix(c(1,4,0, 0,1,4, 0,1,4), 3, byrow = TRUE)              # column margins (1,6,8), not (2,5,8)
  expect_error(eg_bad(table_space = list(tb_F1, wrong)), "does not have the same margins", fixed = TRUE)
  expect_error(eg_bad(table_space = list()), "'table_space' is empty", fixed = TRUE)
  expect_error(eg_bad(), "supply 'table_space'", fixed = TRUE)
  expect_error(eg_bad(statistic = "kendall"))                                 # match.arg
})

test_that("point probability: u_allocation length / range", {
  expect_error(pp(tb_F1, c(3,0,0), 1, c(0,1,1)), "must not exceed", fixed = TRUE)
  expect_error(pp(tb_F1, c(0,0), 1, c(0,1,1)), "'u_allocation' must be a numeric vector of length 3", fixed = TRUE)
  expect_error(pp(tb_F1, c(-1,0,0), 1, c(0,1,1)), "must be non-negative", fixed = TRUE)
})

test_that("empty Omega_t: warning and NA, not a silent number", {
  expect_warning(r <- eg_bad(transform.fun = T_F1, threshold = 1e6), "No tables found", fixed = TRUE)
  expect_true(is.na(r$max.prob))
})


## ==================================================================================
## 11. Characterization (golden master)
## Freezes current, verified output; does not prove correctness.  First run stores, later runs compare.
## ==================================================================================

paper_cases <- list(
  F1_b1   = list(tb = matrix(c(2,3,0, 0,1,4, 0,1,4),3,byrow=TRUE), u = c(0,0,3),  g = 1,   d = c(0,1,1), w = c(0,1,2),   v = c(0,1,2)),
  F1_b2   = list(tb = matrix(c(2,2,1, 1,2,2, 1,2,2),3,byrow=TRUE), u = c(0,3,5),  g = 1,   d = c(0,1,1), w = c(0,1,2),   v = c(0,1,2)),
  F1_b3   = list(tb = matrix(c(3,2,1, 0,2,4, 0,1,5),3,byrow=TRUE), u = c(0,4,10), g = 1,   d = c(0,1,1), w = c(0,1,2),   v = c(0,1,2)),
  F1_b4   = list(tb = matrix(c(3,3,0, 1,2,3, 2,3,1),3,byrow=TRUE), u = c(0,8,4),  g = 1,   d = c(0,1,1), w = c(0,1,2),   v = c(0,1,2)),
  F1_2x3a = list(tb = matrix(c(4,6,0, 1,3,6),2,byrow=TRUE),        u = c(0,1,6),  g = 1,   d = c(0,1),   w = c(0,1),     v = c(0,1,2)),
  F1_2x3b = list(tb = matrix(c(2,4,4, 2,2,6),2,byrow=TRUE),        u = c(3,6,10), g = 1,   d = c(0,1),   w = c(0,1),     v = c(0,1,2)),
  F2_86   = list(tb = tb_F2, u = c(0,10,20), g = 1,   d = c(0,1,1), w = c(0,1,2.5), v = c(0,1,2)),
  F2_girl = list(tb = girls, u = c(0,40,7),  g = 0.5, d = c(0,1,1), w = c(0,1,2.5), v = c(0,1,2)),
  F2_boy  = list(tb = boys,  u = c(0,20,10), g = 0.5, d = c(0,1,1), w = c(0,1,2.5), v = c(0,1,2)),
  D_girl_prior   = list(tb = girls, u = c(0,15,7),  g = log(2), d = c(0,1,1), w = c(0,0.25,1.5), v = c(0,1,1.5)),
  D_boy_prior    = list(tb = boys,  u = c(0,20,10), g = log(2), d = c(0,1,1), w = c(0,0.25,1.5), v = c(0,1,1.5)),
  D_girl_profile = list(tb = girls, u = c(0,15,7),  g = log(2), d = c(0,1,1), w = c(0,0.2,1.6),  v = c(0,0.67,0.7)),
  D_boy_profile  = list(tb = boys,  u = c(0,20,10), g = log(2), d = c(0,1,1), w = c(0,0.2,1.6),  v = c(0,0.67,0.7)))

test_that("paper tables: p at the verify_same_names u, and |Omega_t|", {
  val <- t(sapply(paper_cases, function(cs) {
    r <- exact.general.sen.IxJ(u_space = matrix(cs$u, 1), obs.table = cs$tb, gamma = cs$g, delta = cs$d,
                               treatment.scores = cs$w, outcome.scores = cs$v)
    c(p = r$max.prob, n.tables = r$n.tables, n.distinct.tbar = r$n.distinct.tbar)
  }))
  golden("paper_tables_single_u", val)
})

test_that("paper tables: worst case over U_O (exact.score) for the education data", {
  val <- t(sapply(paper_cases[grep("^D_", names(paper_cases))], function(cs) {
    r <- exact.score.sen.IxJ(obs.table = cs$tb, gamma = cs$g, delta = cs$d, treatment.scores = cs$w, outcome.scores = cs$v)
    c(max.prob = r$max.prob, rct.prob = r$rct.prob, r$maximizer)
  }))
  golden("education_worst_case", val)
})

test_that("F.2 (N = 86) chi-square at u = (0,10,20), C++ statistic", {
  r <- exact.general.sen.IxJ(u_space = matrix(c(0,10,20), 1), obs.table = tb_F2, gamma = 1, delta = c(0,1,1), statistic = "chisq")
  golden("F2_chisq_single_u", c(p = r$max.prob, n.tables = r$n.tables))
})

test_that("tb_F1: full all.prob vectors (ordinal and chi-square, gamma = 1; nonmonotone delta)", {
  golden("F1_ordinal_allprob", exact.general.sen.IxJ(obs.table = tb_F1, gamma = 1, delta = c(0,1,1),
                                                     treatment.scores = w_F1, outcome.scores = v_F1)$all.prob)
  golden("F1_chisq_allprob", exact.general.sen.IxJ(obs.table = tb_F1, gamma = 1, delta = c(0,1,1), statistic = "chisq")$all.prob)
  golden("F1_chisq_nonmonotone_allprob", exact.general.sen.IxJ(obs.table = tb_F1, gamma = 1, delta = c(0,1,0), statistic = "chisq")$all.prob)
})

test_that("normal approximation on the education data", {
  r <- norm.score.sen.IxJ(obs.table = girls, gamma = log(2), delta = c(0,1,1),
                          treatment.scores = c(0, 0.2, 1.6), outcome.scores = c(0, 0.67, 0.7))
  golden("girls_normal", unlist(r[c("T_obs", "RCT.mean", "RCT.var", "RCT.prob", "max.prob", "max.mean", "max.var")]))
})

test_that("SIS estimate is frozen (deterministic: seed = draw index)", {
  skip_if_not(HAS_SIS, "SIS_methods.cpp not found")
  sm <- single.u.sampling.general.IxJ(table = tb_F1, mc.iteration = 500, u_allocation = c(1, 3, 4),
                                      gamma = 1, delta = c(0,1,1), transform.fun = T_F1)
  golden("F1_sis_500", sm$final.prob.estimate)
})


## ==================================================================================
## 12. Known issues (regression tests; skipped unless RUN_KNOWN_ISSUES = TRUE)
## Each block states the INTENDED behaviour, so it fails until the issue is fixed.  After a fix:
## set RUN_KNOWN_ISSUES = TRUE, confirm the block passes, delete its skip_if_not() line.
##   F1 sampling.*.sen.IxJ with row = 'outcome': table transposed twice, statistic in the wrong orientation
##   F3 norm.score.sen.IxJ output names (RCT.prob) differ from the exact functions (rct.prob)
##   F5 sampling / normal functions skip .check_common; single.u.sampling.general.IxJ validates nothing
##   F7 exact.score.sen.IxJ accepts constant scores without a warning
## ==================================================================================

test_that("F1: sampling.score.sen.IxJ, row = 'outcome' with t(table) = row = 'treatment' (I != J)", {
  skip_if_not(RUN_KNOWN_ISSUES, "known issue F1 (set RUN_KNOWN_ISSUES <- TRUE in Section 0)")
  skip_if_not(HAS_SIS, "SIS_methods.cpp not found")
  tb <- matrix(c(4,6,0, 1,3,6), 2, byrow = TRUE)
  f <- function(x, row) sampling.score.sen.IxJ(obs.table = x, gamma = 1, delta = c(0,1), row = row, mc.iteration = 300,
                                               treatment.scores = c(0,1), outcome.scores = c(0,1,2))
  a <- f(tb, "treatment"); b <- f(t(tb), "outcome")
  expect_equal(b$max.prob, a$max.prob)
})

test_that("F1: sampling.general.sen.IxJ, row = 'outcome' with t(table) = row = 'treatment' (I = J, wrong answer risk)", {
  skip_if_not(RUN_KNOWN_ISSUES, "known issue F1 (set RUN_KNOWN_ISSUES <- TRUE in Section 0)")
  skip_if_not(HAS_SIS, "SIS_methods.cpp not found")
  f_tr <- function(x) T_F1(x)                 # statistic on a treatment-rows table
  f_out <- function(x) T_F1(t(x))             # the same statistic written for an outcome-rows table
  a <- sampling.general.sen.IxJ(obs.table = tb_F1, gamma = 1, delta = c(0,1,1), mc.iteration = 300, transform.fun = f_tr)
  b <- sampling.general.sen.IxJ(obs.table = t(tb_F1), row = "outcome", gamma = 1, delta = c(0,1,1),
                                mc.iteration = 300, transform.fun = f_out)
  expect_equal(b$max.prob, a$max.prob)
})

test_that("F3: norm.score.sen.IxJ uses the same output names as the exact functions (rct.prob, gamma, delta)", {
  skip_if_not(RUN_KNOWN_ISSUES, "known issue F3 (set RUN_KNOWN_ISSUES <- TRUE in Section 0)")
  r <- norm.score.sen.IxJ(obs.table = girls, gamma = log(2), delta = c(0,1,1),
                          treatment.scores = c(0, 0.2, 1.6), outcome.scores = c(0, 0.67, 0.7))
  expect_true(all(c("rct.prob", "max.prob", "maximizer", "gamma", "delta") %in% names(r)))
})

test_that("F5: sampling and normal functions reject a constant delta, like the exact functions", {
  skip_if_not(RUN_KNOWN_ISSUES, "known issue F5 (set RUN_KNOWN_ISSUES <- TRUE in Section 0)")
  expect_error(norm.score.sen.IxJ(obs.table = tb_F1, gamma = 1, delta = c(1,1,1),
                                  treatment.scores = w_F1, outcome.scores = v_F1), "constant")
  skip_if_not(HAS_SIS, "SIS_methods.cpp not found")
  expect_error(sampling.score.sen.IxJ(obs.table = tb_F1, gamma = 1, delta = c(1,1,1), mc.iteration = 10,
                                      treatment.scores = w_F1, outcome.scores = v_F1), "constant")
})

test_that("F5: single.u.sampling.general.IxJ validates u_allocation", {
  skip_if_not(RUN_KNOWN_ISSUES, "known issue F5 (set RUN_KNOWN_ISSUES <- TRUE in Section 0)")
  skip_if_not(HAS_SIS, "SIS_methods.cpp not found")
  expect_error(single.u.sampling.general.IxJ(table = tb_F1, mc.iteration = 10, u_allocation = c(9, 9, 9),
                                             gamma = 1, delta = c(0,1,1), transform.fun = T_F1))
})

test_that("F7: exact.score.sen.IxJ warns on constant scores (as sampling.score does)", {
  skip_if_not(RUN_KNOWN_ISSUES, "known issue F7 (set RUN_KNOWN_ISSUES <- TRUE in Section 0)")
  expect_warning(exact.score.sen.IxJ(obs.table = tb_F1, gamma = 1, delta = c(0,1,1),
                                     treatment.scores = c(1,1,1), outcome.scores = v_F1), "identical")
})


## ==================================================================================
## 13. Validity under hidden bias, by simulation  [R13]  (skipped unless RUN_SLOW = TRUE)
## Data generated EXACTLY from model (6) with hidden bias gamma* and no effect, conditional on the margins:
##   P(z) proportional to exp(gamma* sum_s delta_{z_s} u_s).  With A = subjects in the delta = 1 groups,
##   k = |A cap {u = 1}| is Fisher noncentral hypergeometric; labels uniform within A and its complement.
##   Requirement: at gamma = gamma*, P(worst-case p <= alpha) <= alpha (+ 3 MCSE); also for gamma > gamma*.
## ==================================================================================

sim_table <- function(n, m, delta, u, gstar) {
  N <- sum(n); r <- rep(seq_along(m), times = m); Nbar <- sum(n[delta == 1]); U1 <- sum(u)
  k_rng <- max(0, Nbar - (N - U1)):min(U1, Nbar)
  lw <- lchoose(U1, k_rng) + lchoose(N - U1, Nbar - k_rng) + gstar * k_rng
  k <- if (length(k_rng) == 1) k_rng else sample(k_rng, 1, prob = exp(lw - max(lw)))
  one <- which(u == 1); zero <- which(u == 0)
  A <- c(one[sample.int(length(one), k)], zero[sample.int(length(zero), Nbar - k)])
  z <- integer(N)
  lab1 <- rep(which(delta == 1), times = n[delta == 1]); lab0 <- rep(which(delta == 0), times = n[delta == 0])
  z[A] <- lab1[sample.int(length(lab1))]
  z[setdiff(seq_len(N), A)] <- lab0[sample.int(length(lab0))]
  unclass(table(factor(z, levels = seq_along(n)), factor(r, levels = seq_along(m))))
}

run_size <- function(gstar, gtest, n_sim, n = c(10,10,10), m = c(10,10,10), delta = c(0,1,1), alpha = 0.05) {
  r <- rep(seq_along(m), times = m)
  u <- as.integer(rank(-r, ties.method = "first") <= sum(n[delta == 1]))   # adversarial: u = 1 for the highest outcomes
  p <- replicate(n_sim, {
    tb <- sim_table(n, m, delta, u, gstar)
    exact.score.sen.IxJ(obs.table = tb, gamma = gtest, delta = delta,
                        treatment.scores = c(0,1,2), outcome.scores = c(0,1,2))$max.prob
  })
  c(rate = mean(p <= alpha), mcse = sqrt(alpha * (1 - alpha) / n_sim))
}

test_that("the simulator reproduces model (6): empirical table frequencies ~ point probabilities", {
  skip_if_not(RUN_SLOW, "slow (set RUN_SLOW <- TRUE in Section 0)")
  set.seed(201)
  n <- c(2,2,3); m <- c(4,3); delta <- c(0,1,1); u <- c(0,0,0,0,1,1,1); g <- 1.2
  draws <- replicate(20000, paste(as.vector(sim_table(n, m, delta, u, g)), collapse = ","))
  freq <- table(draws) / length(draws)
  tb0 <- matrix(c(2,0, 1,1, 1,2), 3, byrow = TRUE)
  for (t in all_tables(tb0)) {
    k <- paste(as.vector(as.matrix(t)), collapse = ",")
    p <- pp(t, c(0, 3), g, delta)                                   # u-bar of u above: (0 on level 1, 3 on level 2)
    f <- if (k %in% names(freq)) as.numeric(freq[k]) else 0
    expect_lt(abs(f - p), 4 * sqrt(p * (1 - p) / 20000) + 1e-3)
  }
})

test_that("valid at Gamma = Gamma* (hidden bias present), and no less conservative at Gamma > Gamma*", {
  skip_if_not(RUN_SLOW, "slow (set RUN_SLOW <- TRUE in Section 0)")
  set.seed(202)
  for (gstar in c(0, log(2))) {
    s  <- run_size(gstar, gtest = gstar,      n_sim = 600)
    s2 <- run_size(gstar, gtest = gstar + 0.5, n_sim = 600)
    expect_lte(s[["rate"]], 0.05 + 3 * s[["mcse"]], label = sprintf("size at Gamma = Gamma* = %.2f", exp(gstar)))
    expect_lte(s2[["rate"]], 0.05 + 3 * s2[["mcse"]])
  }
})

test_that("Freeman-Halton at gamma = 0 is valid under random assignment", {
  skip_if_not(RUN_SLOW, "slow (set RUN_SLOW <- TRUE in Section 0)")
  set.seed(203)
  n <- c(8,8,8); m <- c(8,8,8)
  p <- replicate(400, {
    tb <- sim_table(n, m, c(0,1,1), rep(0L, sum(n)), 0)
    exact.general.sen.IxJ(u_space = matrix(0L, 1, 3), obs.table = tb, gamma = 0, delta = c(0,1,1),
                          statistic = "freeman_halton")$max.prob
  })
  expect_lte(mean(p <= 0.05), 0.05 + 3 * sqrt(0.05 * 0.95 / 400))
})


## ==================================================================================
## 14. Debugging tools (not tests)
## Use when a block fails: compare three independent routes to alpha(u-bar), check mass, compare input routes.
## ==================================================================================

## alpha per candidate from (i) the engine, (ii) sum of point probabilities over Omega_t, (iii) brute
## force over assignments (N <= max_brute).  engine != sum_pp -> stage 1/2 (build_W_cpp, alpha_bd);
## sum_pp != brute -> point_prob_bd or Omega_t.
dbg_compare <- function(tb, delta, gamma, stat, U = NULL, max_brute = 8, tol = 1e-9) {
  tb <- as.matrix(tb); m <- colSums(tb)
  if (is.null(U)) U <- full.candidate.set(m)
  e  <- exact.general.sen.IxJ(u_space = U, obs.table = tb, gamma = gamma, delta = delta, transform.fun = stat)
  Om <- possible.table(stat(tb), tb, "greater than", stat)
  sum_pp <- apply(U, 1, function(u) sum(sapply(Om, pp, u = u, g = gamma, d = delta)))
  brute  <- if (sum(tb) <= max_brute) brute_alpha(tb, delta, gamma, stat, U)$alpha else NA
  out <- data.frame(U, engine = e$all.prob, sum_pp = sum_pp, brute = brute)
  out$bad <- abs(out$engine - out$sum_pp) > tol | (!is.na(out$brute) & abs(out$engine - out$brute) > tol)
  cat(sprintf("%d candidates, %d disagree (tol %g); |Omega_t| = %d\n", nrow(out), sum(out$bad), tol, length(Om)))
  invisible(out)
}
## total probability over all tables with the margins; must be 1
dbg_mass <- function(tb, u, gamma, delta) sum(sapply(all_tables(tb), pp, u = u, g = gamma, d = delta))
## stage 1 summary: |Omega_t|, tables visited, distinct t-bar
dbg_stage1 <- function(tb, delta, w, v) {
  n <- as.integer(rowSums(tb)); m <- as.integer(colSums(tb)); S <- outer(w, v)
  W <- build_W_cpp(n, m, as.integer(delta), sum(S * tb), greater = TRUE, score = S, stat_fun = NULL)
  cat(sprintf("|Omega_t| = %d of %d tables visited; %d distinct t-bar\n", W$n_tables, W$n_visited, nrow(W$TDmat)))
  invisible(data.frame(W$TDmat, logW = W$logW))
}
## Examples (uncomment):
# d <- dbg_compare(tb_tiny, c(0,1,1), 1.2, T_tiny); d[d$bad, ]
# dbg_mass(tb_F1, c(1,3,4), 1, c(0,1,1))
# dbg_stage1(tb_F1, c(0,1,1), w_F1, v_F1)
## Base R: traceback(); debugonce(exact.general.sen.IxJ); options(error = recover);
##         trace(.sweep, quote(print(dim(U))), at = 2); untrace(.sweep)

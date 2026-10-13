## =============================================================================
## general_test_feasibility_5x2.R
## General (count-based) test feasibility: chi-square statistic, 5 x 2 tables,
## N in {20, 50, 100, 120, 150, 180, 200}.
##
## Depends only on the sensitivityIxJ package files in the same directory:
##   sensitivityIxJ.cpp, sensitivityIxJ.R
##
## Experiment design
##   * Test statistic: chi-square, T = sum_ij (N_ij - E_ij)^2 / E_ij, E_ij = N_i. N_.j / N
##     (a count-based statistic that is NOT arrangement-increasing, so Theorems 2-3 do not
##     apply and Theorem 1 is the only exact route).  exact.general.sen.IxJ also offers
##     statistic = "G2" (2 sum_ij N_ij log(N_ij / E_ij)) and "freeman_halton" (the r x c
##     extension of Fisher's exact test); this script uses the chi-square only.
##   * Sensitivity parameters: delta = (0, 1, ..., 1); gamma in {0, 1}.
##   * Worst case: Theorem 1, maximum over the FULL candidate set U_C of all u-bar vectors,
##     |U_C| = prod_j (N_.j + 1).  The timed call is exact.general.sen.IxJ with its default
##     u_space (= U_C) and statistic = "chisq".  For |U_C| > CHUNK the sweep is split into
##     checkpointed chunks (see run_sensitivity); the block can be re-run to resume.
##   * Observed tables: fixed below, one per N, with a positive association between
##     treatment level and outcome level; all have chi-square p-value < 0.05 at gamma = 0
##     (values in the comments; block (2) re-computes them).  They differ from the tables
##     of the ordinal scripts because the chi-square test has less power against a
##     monotone association.
##   * For each N, three blocks: (1) the observed table; (2) exact p-value at gamma = 0
##     from exact.general.sen.IxJ evaluated at u = 0; (3) the timed worst case at gamma = 1.
##
## Expected running time.  Stage 2 costs about |U_C| x (distinct t-bar) x N^2 x c; the
## comment on each N = ... block gives |U_C|.  Blocks with |U_C| in the millions take hours
## to days on 4 CPUs: launch the script as an RStudio Background Job and let the
## checkpoints carry it across sessions.
##
## CSV columns (one row per run, appended to general_test_feasibility.csv)
##   test                      "chisq"
##   dimension, I, J, N        table dimension and sample size
##   delta                     sensitivity parameter vector (0, 1, ..., 1)
##   table_rowwise             the observed table, row by row (rows = treatment)
##   row_margins, col_margins  N_i. and N_.j
##   T_obs                     observed chi-square statistic
##   n_tables                  |Omega_t|: tables with the observed margins and T >= T_obs
##   n_distinct_tbar           number of distinct t-bar vectors among them (size of W)
##   p_gamma0                  exact p-value at gamma = 0 (block 2)
##   gamma                     the sensitivity parameter used for the worst case (= gamma1 = 1)
##   candidate_set             "U_C (Theorem 1)": the full candidate set of Theorem 1
##   n_candidates              |U_C| = prod_j (N_.j + 1), all of which are evaluated
##   n_chunks                  1 if the sweep was one call; > 1 if it was split into chunks
##   p_gamma1                  exact worst-case p-value at gamma = 1 (max over U_C)
##   ubar_maximizer            the worst-case u in u-bar coordinates
##   sec_stage1_enumerate_W    seconds for stage 1 of Corollary C.3 (one call)
##   sec_stage2_sweep          seconds for stage 2, summed over chunks
##   sec_total                 wall time of the whole block (all calls, incl. repeated stage 1)
##   mem_MB_before, mem_MB_after, mem_MB_left_behind, mem_MB_peak_increase, peak_is_per_run
##                             as in the ordinal scripts (peak = high-water mark reset per run)
##   R_version, timestamp
## =============================================================================

library(Rcpp)
sourceCpp("sensitivityIxJ.cpp"); source("sensitivityIxJ.R")

CSV_FILE <- "general_test_feasibility.csv"
CKPT_DIR <- "ckpt_general"; dir.create(CKPT_DIR, showWarnings = FALSE)
CHUNK    <- 200000L                   # candidates per exact.general.sen.IxJ call when |U_C| is large
I <- 5; J <- 2
delta <- c(0, rep(1, I - 1))          # level 1 is the reference level
gamma1 <- 1                           # the sensitivity parameter at which the worst case is computed

## ---- helpers ------------------------------------------------------------------------
## Print the setup of one case.
show_setup <- function(tb, delta, label) {
  cat("\n", strrep("=", 78), "\n", label, "\n", strrep("=", 78), "\n", sep = "")
  cat("Test statistic : chi-square, T = sum_ij (N_ij - E_ij)^2 / E_ij, E_ij = N_i. N_.j / N\n")
  cat("delta          :", paste(delta, collapse = ","), "\n")
  cat("Observed table (rows = treatment, columns = outcome):\n"); print(tb)
  cat("Row margins    :", paste(rowSums(tb), collapse = ","),
      "   Column margins :", paste(colSums(tb), collapse = ","), "\n")
  E <- outer(rowSums(tb), colSums(tb)) / sum(tb)
  cat("Observed T     :", round(sum((tb - E)^2 / E), 4), "\n")
  cat("Candidate set  : U_C, |U_C| = prod(N_.j + 1) =", format(prod(colSums(tb) + 1), big.mark = ","), "(Theorem 1)\n")
}

## Memory (Linux / Posit Cloud), read from /proc/self/status:
##   VmRSS = current resident memory of this R process;
##   VmHWM = its peak resident memory ("high-water mark").  VmHWM only ever increases, so
##   before each run we reset it by writing 5 to /proc/self/clear_refs; after the run,
##   VmHWM - (VmRSS before) is the peak memory THIS run held above its starting point,
##   including transients freed before the run ended.  If the reset is refused, the
##   fallback is the session-wide peak, flagged in the CSV (peak_is_per_run = FALSE).
mem_mb <- function() {
  st <- tryCatch(readLines("/proc/self/status"), error = function(e) NULL)
  if (is.null(st)) return(c(rss = NA_real_, peak = NA_real_))
  get <- function(key) as.numeric(sub("^\\S+:\\s+(\\d+).*$", "\\1", grep(key, st, value = TRUE))) / 1024
  c(rss = get("^VmRSS"), peak = get("^VmHWM"))
}
reset_peak <- function() {
  tryCatch({ writeLines("5", "/proc/self/clear_refs"); TRUE }, error = function(e) FALSE, warning = function(w) FALSE)
}

## Timed sensitivity analysis at gamma = gamma1 with the package function
## exact.general.sen.IxJ (Theorem 1: maximum over the FULL candidate set U_C, which is the
## function's default u_space; Corollary C.3 for each evaluation).  The chi-square
## statistic is evaluated inside the C++ enumerator (statistic = "chisq").
##
## Long runs.  |U_C| = prod_j (N_.j + 1) and each candidate costs O(N^{J+1}); for a 3 x 4
## table at N = 200 that is about 7 million candidates.  When |U_C| exceeds CHUNK the sweep
## is done in blocks of CHUNK candidates, each block by its own call to
## exact.general.sen.IxJ with u_space = full.candidate.set.rows(...), and each block's
## maximum is written to a .rds file in CKPT_DIR.  Re-running the block after a
## suspended session resumes from the finished chunks.  Stage 1 (enumeration) is repeated
## by every call, so it is reported for one call; Stage 2 is summed over calls.
run_sensitivity <- function(tb, delta, p0, label) {
  I <- nrow(tb); J <- ncol(tb); N <- sum(tb); m <- colSums(tb)
  n_UC <- prod(m + 1); tag <- sprintf("chisq_%dx%d_N%d", I, J, N)
  gc(); peak_reset_ok <- reset_peak(); m0 <- mem_mb()
  t0 <- proc.time()[["elapsed"]]
  if (n_UC <= CHUNK) {
    res <- exact.general.sen.IxJ(obs.table = tb, gamma = gamma1, delta = delta, row = "treatment", statistic = "chisq")
    best <- res$max.prob; best_u <- res$maximizer; n_tables <- res$n.tables; n_tbar <- res$n.distinct.tbar
    s1 <- res$time[["enumerate.W"]]; s2 <- res$time[["sweep"]]; rct <- res$rct.prob; n_chunks <- 1
  } else {
    starts <- seq(1, n_UC, by = CHUNK); n_chunks <- length(starts)
    best <- -Inf; best_u <- NULL; s1 <- NA; s2 <- 0; n_tables <- NA; n_tbar <- NA; rct <- NA
    for (ci in seq_along(starts)) {
      f <- file.path(CKPT_DIR, sprintf("%s_chunk%04d.rds", tag, ci))
      if (file.exists(f)) { r <- readRDS(f); cat(sprintf("  chunk %d/%d: loaded from checkpoint\n", ci, n_chunks)) }
      else {
        idx <- starts[ci]:min(starts[ci] + CHUNK - 1, n_UC)
        res <- exact.general.sen.IxJ(u_space = full.candidate.set.rows(m, idx), obs.table = tb, gamma = gamma1,
                                     delta = delta, row = "treatment", statistic = "chisq")
        r <- list(max = res$max.prob, u = res$maximizer, s1 = res$time[["enumerate.W"]], s2 = res$time[["sweep"]],
                  n_tables = res$n.tables, n_tbar = res$n.distinct.tbar, rct = res$rct.prob)
        saveRDS(r, f)
        cat(sprintf("  chunk %d/%d: %d candidates, max %.6f, %.1f s\n", ci, n_chunks, length(idx), r$max, r$s1 + r$s2))
      }
      if (r$max > best) { best <- r$max; best_u <- r$u }
      if (is.na(s1)) { s1 <- r$s1; n_tables <- r$n_tables; n_tbar <- r$n_tbar; rct <- r$rct }
      s2 <- s2 + r$s2
    }
  }
  t_total <- proc.time()[["elapsed"]] - t0
  m1 <- mem_mb()
  cat(sprintf("Stage 1 (enumerate Omega_t, build W): |Omega_t| = %s tables, %d distinct t-bar, %.3f s\n",
              format(n_tables, big.mark = ","), n_tbar, s1))
  cat(sprintf("p-value, gamma = 0 : %.6f\n", p0))
  cat(sprintf("Stage 2 (U_C, Theorem 1): p-value, gamma = %g : %.6f  (%s candidates, worst case at u-bar = (%s), %.3f s%s)\n",
              gamma1, best, format(n_UC, big.mark = ","), paste(best_u, collapse = ","), s2,
              if (n_chunks > 1) sprintf(", %d chunks", n_chunks) else ""))
  cat(sprintf("Total time %.3f s | memory: RSS %.0f MB -> %.0f MB (left behind %.0f MB); peak during run +%.0f MB%s\n",
              t_total, m0["rss"], m1["rss"], m1["rss"] - m0["rss"], m1["peak"] - m0["rss"],
              if (peak_reset_ok) "" else "  [peak reset refused: session-wide peak]"))
  E <- outer(rowSums(tb), colSums(tb)) / sum(tb)
  row <- data.frame(
    test = "chisq", dimension = sprintf("%dx%d", I, J), I = I, J = J, N = N,
    delta = paste(delta, collapse = " "),
    table_rowwise = paste(as.vector(t(tb)), collapse = " "),
    row_margins = paste(rowSums(tb), collapse = " "), col_margins = paste(colSums(tb), collapse = " "),
    T_obs = round(sum((tb - E)^2 / E), 6), n_tables = n_tables, n_distinct_tbar = n_tbar,
    p_gamma0 = p0, gamma = gamma1, candidate_set = "U_C (Theorem 1)", n_candidates = n_UC, n_chunks = n_chunks,
    p_gamma1 = best, ubar_maximizer = paste(best_u, collapse = " "),
    sec_stage1_enumerate_W = round(s1, 4), sec_stage2_sweep = round(s2, 4), sec_total = round(t_total, 4),
    mem_MB_before = round(m0["rss"]), mem_MB_after = round(m1["rss"]),
    mem_MB_left_behind = round(m1["rss"] - m0["rss"]),
    mem_MB_peak_increase = round(m1["peak"] - m0["rss"]), peak_is_per_run = peak_reset_ok,
    R_version = R.version.string, timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    stringsAsFactors = FALSE)
  write.table(row, CSV_FILE, sep = ",", row.names = FALSE, col.names = !file.exists(CSV_FILE),
              append = file.exists(CSV_FILE))
  cat("Row appended to", CSV_FILE, "\n")
  invisible(row)
}


## =============================================================================
## Runs.  Three blocks per N; run them in order.
## =============================================================================

## ---- N = 20: (1) observed table   (chi-square p at gamma = 0: 0.01067; |U_C| = 121) ------
tb <- matrix(c(5, 0, 3, 1, 1, 1, 1, 3, 0, 5), nrow = 5, byrow = TRUE)
dimnames(tb) <- list(paste0("Trt", 1:I), paste0("Out", 1:J)); storage.mode(tb) <- "integer"
show_setup(tb, delta, sprintf("chi-square test, %dx%d, N = %d", I, J, sum(tb)))

## ---- N = 20: (2) p-value at gamma = 0 (exact.general.sen.IxJ at u = 0) --------------
p0 <- exact.general.sen.IxJ(u_space = matrix(0L, 1, J), obs.table = tb, gamma = 0, delta = delta,
                            row = "treatment", statistic = "chisq")$max.prob
cat(sprintf("p-value at gamma = 0: %.6f  (%s)\n", p0, if (p0 < 0.05) "below 0.05" else "NOT below 0.05"))

## ---- N = 20: (3) exact.general.sen.IxJ at gamma = 1 over U_C (Theorem 1), timed, to CSV ---
run_sensitivity(tb, delta, p0, sprintf("chi-square test, %dx%d, N = %d", I, J, sum(tb)))

## ---- N = 50: (1) observed table   (chi-square p at gamma = 0: 0.0102; |U_C| = 676) ------
tb <- matrix(c(10, 2, 6, 3, 4, 4, 3, 6, 2, 10), nrow = 5, byrow = TRUE)
dimnames(tb) <- list(paste0("Trt", 1:I), paste0("Out", 1:J)); storage.mode(tb) <- "integer"
show_setup(tb, delta, sprintf("chi-square test, %dx%d, N = %d", I, J, sum(tb)))

## ---- N = 50: (2) p-value at gamma = 0 (exact.general.sen.IxJ at u = 0) --------------
p0 <- exact.general.sen.IxJ(u_space = matrix(0L, 1, J), obs.table = tb, gamma = 0, delta = delta,
                            row = "treatment", statistic = "chisq")$max.prob
cat(sprintf("p-value at gamma = 0: %.6f  (%s)\n", p0, if (p0 < 0.05) "below 0.05" else "NOT below 0.05"))

## ---- N = 50: (3) exact.general.sen.IxJ at gamma = 1 over U_C (Theorem 1), timed, to CSV ---
run_sensitivity(tb, delta, p0, sprintf("chi-square test, %dx%d, N = %d", I, J, sum(tb)))

## ---- N = 100: (1) observed table   (chi-square p at gamma = 0: 0.002876; |U_C| = 2,601) ------
tb <- matrix(c(17, 5, 12, 7, 9, 9, 7, 12, 5, 17), nrow = 5, byrow = TRUE)
dimnames(tb) <- list(paste0("Trt", 1:I), paste0("Out", 1:J)); storage.mode(tb) <- "integer"
show_setup(tb, delta, sprintf("chi-square test, %dx%d, N = %d", I, J, sum(tb)))

## ---- N = 100: (2) p-value at gamma = 0 (exact.general.sen.IxJ at u = 0) --------------
p0 <- exact.general.sen.IxJ(u_space = matrix(0L, 1, J), obs.table = tb, gamma = 0, delta = delta,
                            row = "treatment", statistic = "chisq")$max.prob
cat(sprintf("p-value at gamma = 0: %.6f  (%s)\n", p0, if (p0 < 0.05) "below 0.05" else "NOT below 0.05"))

## ---- N = 100: (3) exact.general.sen.IxJ at gamma = 1 over U_C (Theorem 1), timed, to CSV ---
run_sensitivity(tb, delta, p0, sprintf("chi-square test, %dx%d, N = %d", I, J, sum(tb)))

## ---- N = 120: (1) observed table   (chi-square p at gamma = 0: 0.000505; |U_C| = 3,721) ------
tb <- matrix(c(20, 6, 15, 8, 11, 11, 8, 15, 6, 20), nrow = 5, byrow = TRUE)
dimnames(tb) <- list(paste0("Trt", 1:I), paste0("Out", 1:J)); storage.mode(tb) <- "integer"
show_setup(tb, delta, sprintf("chi-square test, %dx%d, N = %d", I, J, sum(tb)))

## ---- N = 120: (2) p-value at gamma = 0 (exact.general.sen.IxJ at u = 0) --------------
p0 <- exact.general.sen.IxJ(u_space = matrix(0L, 1, J), obs.table = tb, gamma = 0, delta = delta,
                            row = "treatment", statistic = "chisq")$max.prob
cat(sprintf("p-value at gamma = 0: %.6f  (%s)\n", p0, if (p0 < 0.05) "below 0.05" else "NOT below 0.05"))

## ---- N = 120: (3) exact.general.sen.IxJ at gamma = 1 over U_C (Theorem 1), timed, to CSV ---
run_sensitivity(tb, delta, p0, sprintf("chi-square test, %dx%d, N = %d", I, J, sum(tb)))

## ---- N = 150: (1) observed table   (chi-square p at gamma = 0: 0.000625(MC); |U_C| = 5,776) ------
tb <- matrix(c(24, 8, 18, 11, 14, 14, 11, 18, 8, 24), nrow = 5, byrow = TRUE)
dimnames(tb) <- list(paste0("Trt", 1:I), paste0("Out", 1:J)); storage.mode(tb) <- "integer"
show_setup(tb, delta, sprintf("chi-square test, %dx%d, N = %d", I, J, sum(tb)))

## ---- N = 150: (2) p-value at gamma = 0 (exact.general.sen.IxJ at u = 0) --------------
p0 <- exact.general.sen.IxJ(u_space = matrix(0L, 1, J), obs.table = tb, gamma = 0, delta = delta,
                            row = "treatment", statistic = "chisq")$max.prob
cat(sprintf("p-value at gamma = 0: %.6f  (%s)\n", p0, if (p0 < 0.05) "below 0.05" else "NOT below 0.05"))

## ---- N = 150: (3) exact.general.sen.IxJ at gamma = 1 over U_C (Theorem 1), timed, to CSV ---
run_sensitivity(tb, delta, p0, sprintf("chi-square test, %dx%d, N = %d", I, J, sum(tb)))

## ---- N = 180: (1) observed table   (chi-square p at gamma = 0: 0.0002(MC); |U_C| = 8,281) ------
tb <- matrix(c(28, 10, 22, 13, 17, 17, 13, 22, 10, 28), nrow = 5, byrow = TRUE)
dimnames(tb) <- list(paste0("Trt", 1:I), paste0("Out", 1:J)); storage.mode(tb) <- "integer"
show_setup(tb, delta, sprintf("chi-square test, %dx%d, N = %d", I, J, sum(tb)))

## ---- N = 180: (2) p-value at gamma = 0 (exact.general.sen.IxJ at u = 0) --------------
p0 <- exact.general.sen.IxJ(u_space = matrix(0L, 1, J), obs.table = tb, gamma = 0, delta = delta,
                            row = "treatment", statistic = "chisq")$max.prob
cat(sprintf("p-value at gamma = 0: %.6f  (%s)\n", p0, if (p0 < 0.05) "below 0.05" else "NOT below 0.05"))

## ---- N = 180: (3) exact.general.sen.IxJ at gamma = 1 over U_C (Theorem 1), timed, to CSV ---
run_sensitivity(tb, delta, p0, sprintf("chi-square test, %dx%d, N = %d", I, J, sum(tb)))

## ---- N = 200: (1) observed table   (chi-square p at gamma = 0: 0.00035(MC); |U_C| = 10,201) ------
tb <- matrix(c(30, 12, 24, 15, 19, 19, 15, 24, 12, 30), nrow = 5, byrow = TRUE)
dimnames(tb) <- list(paste0("Trt", 1:I), paste0("Out", 1:J)); storage.mode(tb) <- "integer"
show_setup(tb, delta, sprintf("chi-square test, %dx%d, N = %d", I, J, sum(tb)))

## ---- N = 200: (2) p-value at gamma = 0 (exact.general.sen.IxJ at u = 0) --------------
p0 <- exact.general.sen.IxJ(u_space = matrix(0L, 1, J), obs.table = tb, gamma = 0, delta = delta,
                            row = "treatment", statistic = "chisq")$max.prob
cat(sprintf("p-value at gamma = 0: %.6f  (%s)\n", p0, if (p0 < 0.05) "below 0.05" else "NOT below 0.05"))

## ---- N = 200: (3) exact.general.sen.IxJ at gamma = 1 over U_C (Theorem 1), timed, to CSV ---
run_sensitivity(tb, delta, p0, sprintf("chi-square test, %dx%d, N = %d", I, J, sum(tb)))

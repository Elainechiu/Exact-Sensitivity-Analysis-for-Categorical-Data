## =============================================================================
## ordinal_feasibility_3x3.R
## Ordinal test feasibility, 3 x 3 tables, N in {20, 50, 100, 120, 150, 180, 200}.
##
## Depends only on the sensitivityIxJ package files in the same directory:
##   sensitivityIxJ.cpp, sensitivityIxJ.R
##
## Experiment design
##   * Test statistic: ordinal, T = sum_{i,j} w_i v_j N_ij with scores equal to the
##     levels, w = (1, ..., 3), v = (1, ..., 3).  Rows = treatment, columns = outcome.
##   * Sensitivity parameters: delta = (0, 1, ..., 1); gamma in {0, 1}.
##   * Worst case: Theorem 2, maximum over the N + 1 candidates in U_O; the timed call
##     is exact.score.sen.IxJ.
##   * Observed tables: fixed below, one per N, with a positive association between
##     treatment level and outcome level.  All have exact p-value < 0.05 at gamma = 0
##     (values in the comments; block (2) re-computes them).
##   * For each N, three blocks: (1) the observed table; (2) exact p-value at gamma = 0
##     from exact.score.sen.IxJ; (3) the timed worst-case computation at gamma = 1, with
##     memory, appended as one row to ordinal_test_feasibility.csv.
##
## CSV columns (one row per run; the same columns in all four ordinal_feasibility_*.R)
##   test                      "ordinal"
##   dimension, I, J, N        table dimension and sample size
##   treatment_scores          w_1 ... w_I  (= 1 ... I)
##   outcome_scores            v_1 ... v_J  (= 1 ... J)
##   delta                     sensitivity parameter vector (0, 1, ..., 1)
##   table_rowwise             the observed table, row by row (rows = treatment)
##   row_margins, col_margins  N_i. and N_.j
##   T_obs                     observed value of the ordinal statistic
##   n_tables                  |Omega_t|: number of tables with the observed margins and
##                             T >= T_obs (the tables that enter the p-value)
##   n_distinct_tbar           number of distinct t-bar vectors among those tables, i.e. the
##                             size of W(t-bar) in Corollary C.3 (what stage 2 loops over)
##   p_gamma0                  exact p-value at gamma = 0, from exact.score.sen.IxJ (block 2)
##   gamma                     the sensitivity parameter used for the worst case (= gamma1 = 1)
##   candidate_set             "U_O (Theorem 2)" for J >= 3: the N + 1 candidates of Theorem 2;
##                             "u+ (Theorem 3)" for J = 2: the closed-form worst case
##   n_candidates              number of u evaluated by exact.score.sen.IxJ: N + 1 under
##                             Theorem 2; 2 under Theorem 3 (u+ and the reference point u = 0)
##   p_gamma1                  exact worst-case p-value at gamma = 1
##   ubar_maximizer            the worst-case u in u-bar coordinates (u-bar_j = number of
##                             subjects with outcome level j and u = 1)
##   sec_stage1_enumerate_W    seconds for stage 1 of Corollary C.3 inside the package call:
##                             one pass over Omega_t that builds W(t-bar); independent of u, gamma
##   sec_stage2_sweep          seconds for stage 2: evaluating alpha at the candidate(s)
##   sec_total                 wall time of the whole package call (>= sum of the two stages)
##   mem_MB_before, mem_MB_after   resident memory (RSS) of the R process before / after
##   mem_MB_left_behind        mem_MB_after - mem_MB_before: memory still resident after the run
##   mem_MB_peak_increase      (peak RSS during the run) - mem_MB_before: the most memory the
##                             run held above its starting point, transients included.  This is
##                             the number to report as the run's memory requirement.
##   peak_is_per_run           TRUE if the peak counter was reset before the run (per-run peak);
##                             FALSE if the reset was refused and the peak is session-wide
##   R_version, timestamp
##
## Example row (illustrative numbers), shown column = value so that it is readable:
##   test                    = ordinal
##   dimension, I, J, N      = 3x3, 3, 3, 100
##   treatment_scores        = 1 2 3
##   outcome_scores          = 1 2 3
##   delta                   = 0 1 1
##   table_rowwise           = 17 10 7 11 11 10 7 10 17
##   row_margins             = 34 32 34
##   col_margins             = 35 31 34
##   T_obs                   = 601
##   n_tables                = 181016
##   n_distinct_tbar         = 900
##   p_gamma0                = 0.002132
##   gamma                   = 1
##   candidate_set           = U_O (Theorem 2)
##   n_candidates            = 101
##   p_gamma1                = 0.0312
##   ubar_maximizer          = 0 0 34
##   sec_stage1_enumerate_W  = 0.0080
##   sec_stage2_sweep        = 0.0020
##   sec_total               = 0.0110
##   mem_MB_before           = 310
##   mem_MB_after            = 312
##   mem_MB_left_behind      = 2
##   mem_MB_peak_increase    = 9
##   peak_is_per_run         = TRUE
##   R_version, timestamp    = R version 4.6.1 (2026-06-12), 2026-09-22 14:03:11
##
## Reading the example: a 3 x 3 table with N = 100 (rows 17 10 7 / 11 11 10 / 7 10 17),
## ordinal statistic with scores 1..3 on both axes, delta = (0,1,1).  The observed
## statistic is 601.  At gamma = 0 the exact p-value is 0.002 (significant, so the
## sensitivity analysis is meaningful).  There are 181,016 tables with the observed
## margins and T >= 601; they collapse to 900 distinct t-bar vectors, which is what
## stage 2 loops over.  Under Theorem 2 the worst case is searched over the 101
## candidates of U_O; at gamma = 1 it is 0.031, attained when the 34 subjects with
## outcome level 3 have u = 1 and nobody else does (u-bar = (0,0,34)).  Stage 1 took
## 8 ms, stage 2 took 2 ms; the whole exact.score.sen.IxJ call took 11 ms.  The run's
## peak memory was 9 MB above where it started (the run's memory requirement); 2 MB of
## that was still resident when it finished.
## =============================================================================

library(Rcpp)
sourceCpp("sensitivityIxJ.cpp"); source("sensitivityIxJ.R")

CSV_FILE <- "ordinal_test_feasibility.csv"
I <- 3; J <- 3
w <- 1:I; v <- 1:J                    # scores = levels
delta <- c(0, rep(1, I - 1))          # level 1 is the reference level
gamma1 <- 1                           # the sensitivity parameter at which the worst case is computed
CANDIDATE_LABEL <- "U_O (Theorem 2)"
CANDIDATE_NOTE  <- function(tb) cat("Candidate set  : U_O, |U_O| = N + 1 =", sum(tb) + 1, "(Theorem 2)\n")

## ---- helpers ------------------------------------------------------------------------
## Print the setup of one case.
show_setup <- function(tb, w, v, delta, label) {
  cat("\n", strrep("=", 78), "\n", label, "\n", strrep("=", 78), "\n", sep = "")
  cat("Test statistic : ordinal, T = sum w_i v_j N_ij with w =", paste(w, collapse = ","),
      " v =", paste(v, collapse = ","), "\n")
  cat("delta          :", paste(delta, collapse = ","), "\n")
  cat("Observed table (rows = treatment, columns = outcome):\n"); print(tb)
  cat("Row margins    :", paste(rowSums(tb), collapse = ","),
      "   Column margins :", paste(colSums(tb), collapse = ","), "\n")
  cat("Observed T     :", sum(outer(w, v) * tb), "\n")
  CANDIDATE_NOTE(tb)
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

## Append one row to CSV_FILE (same columns in all four ordinal_feasibility_*.R scripts).
write_row <- function(tb, w, v, delta, p0, res, candidate_set, t_total, m0, m1, peak_reset_ok) {
  I <- nrow(tb); J <- ncol(tb); N <- sum(tb)
  row <- data.frame(
    test = "ordinal", dimension = sprintf("%dx%d", I, J), I = I, J = J, N = N,
    treatment_scores = paste(w, collapse = " "), outcome_scores = paste(v, collapse = " "),
    delta = paste(delta, collapse = " "),
    table_rowwise = paste(as.vector(t(tb)), collapse = " "),
    row_margins = paste(rowSums(tb), collapse = " "), col_margins = paste(colSums(tb), collapse = " "),
    T_obs = sum(outer(w, v) * tb), n_tables = res$n.tables, n_distinct_tbar = res$n.distinct.tbar,
    p_gamma0 = p0, gamma = gamma1,
    candidate_set = candidate_set, n_candidates = res$n.candidates,
    p_gamma1 = res$max.prob, ubar_maximizer = paste(res$maximizer, collapse = " "),
    sec_stage1_enumerate_W = round(res$time[["enumerate.W"]], 4), sec_stage2_sweep = round(res$time[["sweep"]], 4),
    sec_total = round(t_total, 4),
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

## Timed sensitivity analysis at gamma = gamma1 with the package function
## exact.score.sen.IxJ.  The function chooses the candidate set from the outcome
## dimension: for J = 2 the closed-form worst case of Theorem 3 (u-bar+ = (0, N_.2),
## together with the reference point (0, 0)); for J >= 3 the N + 1 candidates of U_O
## (Theorem 2).  Records time and memory and appends one row to CSV_FILE.  p0 is the
## gamma = 0 p-value from block (2).
run_sensitivity <- function(tb, w, v, delta, p0, label) {
  gc(); peak_reset_ok <- reset_peak(); m0 <- mem_mb()
  t0 <- proc.time()[["elapsed"]]
  res <- exact.score.sen.IxJ(obs.table = tb, gamma = gamma1, delta = delta, row = "treatment",
                             treatment.scores = w, outcome.scores = v, verbose = FALSE)
  t_total <- proc.time()[["elapsed"]] - t0
  m1 <- mem_mb()
  cat(sprintf("Stage 1 (enumerate Omega_t, build W): |Omega_t| = %s tables, %d distinct t-bar, %.3f s\n",
              format(res$n.tables, big.mark = ","), res$n.distinct.tbar, res$time[["enumerate.W"]]))
  cat(sprintf("p-value, gamma = 0 : %.6f\n", p0))
  cat(sprintf("Stage 2 (%s): p-value, gamma = %g : %.6f  (%d candidates evaluated, worst case at u-bar = (%s), %.3f s)\n",
              CANDIDATE_LABEL, gamma1, res$max.prob, res$n.candidates,
              paste(res$maximizer, collapse = ","), res$time[["sweep"]]))
  cat(sprintf("Total time %.3f s (exact.score.sen.IxJ call) | memory: RSS %.0f MB -> %.0f MB (left behind %.0f MB); peak during run +%.0f MB%s\n",
              t_total, m0["rss"], m1["rss"], m1["rss"] - m0["rss"], m1["peak"] - m0["rss"],
              if (peak_reset_ok) "" else "  [peak reset refused: session-wide peak]"))
  write_row(tb, w, v, delta, p0, res, CANDIDATE_LABEL, t_total, m0, m1, peak_reset_ok)
}

## =============================================================================
## Runs.  Three blocks per N; run them in order.
## =============================================================================

## ---- N = 20: (1) observed table   (p at gamma = 0: 0.0002531) --------------------
tb <- matrix(c(6, 2, 0, 2, 1, 1, 0, 2, 6), nrow = 3, byrow = TRUE)
dimnames(tb) <- list(paste0("Trt", 1:I), paste0("Out", 1:J)); storage.mode(tb) <- "integer"
show_setup(tb, w, v, delta, sprintf("ordinal test, %dx%d, N = %d", I, J, sum(tb)))

## ---- N = 20: (2) p-value at gamma = 0 (exact.score.sen.IxJ) ----------------
p0 <- exact.score.sen.IxJ(obs.table = tb, gamma = 0, delta = delta, row = "treatment",
                          treatment.scores = w, outcome.scores = v, verbose = FALSE)$max.prob
cat(sprintf("p-value at gamma = 0: %.6f  (%s)\n", p0, if (p0 < 0.05) "below 0.05" else "NOT below 0.05"))

## ---- N = 20: (3) exact.score.sen.IxJ at gamma = 1 (Theorem 2, sweep U_O), timed, to CSV
run_sensitivity(tb, w, v, delta, p0, sprintf("ordinal test, %dx%d, N = %d", I, J, sum(tb)))

## ---- N = 50: (1) observed table   (p at gamma = 0: 0.001521) --------------------
tb <- matrix(c(10, 5, 3, 5, 5, 5, 2, 5, 10), nrow = 3, byrow = TRUE)
dimnames(tb) <- list(paste0("Trt", 1:I), paste0("Out", 1:J)); storage.mode(tb) <- "integer"
show_setup(tb, w, v, delta, sprintf("ordinal test, %dx%d, N = %d", I, J, sum(tb)))

## ---- N = 50: (2) p-value at gamma = 0 (exact.score.sen.IxJ) ----------------
p0 <- exact.score.sen.IxJ(obs.table = tb, gamma = 0, delta = delta, row = "treatment",
                          treatment.scores = w, outcome.scores = v, verbose = FALSE)$max.prob
cat(sprintf("p-value at gamma = 0: %.6f  (%s)\n", p0, if (p0 < 0.05) "below 0.05" else "NOT below 0.05"))

## ---- N = 50: (3) exact.score.sen.IxJ at gamma = 1 (Theorem 2, sweep U_O), timed, to CSV
run_sensitivity(tb, w, v, delta, p0, sprintf("ordinal test, %dx%d, N = %d", I, J, sum(tb)))

## ---- N = 100: (1) observed table   (p at gamma = 0: 0.002132) --------------------
tb <- matrix(c(17, 10, 7, 11, 11, 10, 7, 10, 17), nrow = 3, byrow = TRUE)
dimnames(tb) <- list(paste0("Trt", 1:I), paste0("Out", 1:J)); storage.mode(tb) <- "integer"
show_setup(tb, w, v, delta, sprintf("ordinal test, %dx%d, N = %d", I, J, sum(tb)))

## ---- N = 100: (2) p-value at gamma = 0 (exact.score.sen.IxJ) ----------------
p0 <- exact.score.sen.IxJ(obs.table = tb, gamma = 0, delta = delta, row = "treatment",
                          treatment.scores = w, outcome.scores = v, verbose = FALSE)$max.prob
cat(sprintf("p-value at gamma = 0: %.6f  (%s)\n", p0, if (p0 < 0.05) "below 0.05" else "NOT below 0.05"))

## ---- N = 100: (3) exact.score.sen.IxJ at gamma = 1 (Theorem 2, sweep U_O), timed, to CSV
run_sensitivity(tb, w, v, delta, p0, sprintf("ordinal test, %dx%d, N = %d", I, J, sum(tb)))

## ---- N = 120: (1) observed table   (p at gamma = 0: 0.001118) --------------------
tb <- matrix(c(20, 13, 8, 13, 13, 13, 8, 13, 19), nrow = 3, byrow = TRUE)
dimnames(tb) <- list(paste0("Trt", 1:I), paste0("Out", 1:J)); storage.mode(tb) <- "integer"
show_setup(tb, w, v, delta, sprintf("ordinal test, %dx%d, N = %d", I, J, sum(tb)))

## ---- N = 120: (2) p-value at gamma = 0 (exact.score.sen.IxJ) ----------------
p0 <- exact.score.sen.IxJ(obs.table = tb, gamma = 0, delta = delta, row = "treatment",
                          treatment.scores = w, outcome.scores = v, verbose = FALSE)$max.prob
cat(sprintf("p-value at gamma = 0: %.6f  (%s)\n", p0, if (p0 < 0.05) "below 0.05" else "NOT below 0.05"))

## ---- N = 120: (3) exact.score.sen.IxJ at gamma = 1 (Theorem 2, sweep U_O), timed, to CSV
run_sensitivity(tb, w, v, delta, p0, sprintf("ordinal test, %dx%d, N = %d", I, J, sum(tb)))

## ---- N = 150: (1) observed table   (p at gamma = 0: 0.001056) --------------------
tb <- matrix(c(24, 16, 11, 16, 16, 16, 11, 16, 24), nrow = 3, byrow = TRUE)
dimnames(tb) <- list(paste0("Trt", 1:I), paste0("Out", 1:J)); storage.mode(tb) <- "integer"
show_setup(tb, w, v, delta, sprintf("ordinal test, %dx%d, N = %d", I, J, sum(tb)))

## ---- N = 150: (2) p-value at gamma = 0 (exact.score.sen.IxJ) ----------------
p0 <- exact.score.sen.IxJ(obs.table = tb, gamma = 0, delta = delta, row = "treatment",
                          treatment.scores = w, outcome.scores = v, verbose = FALSE)$max.prob
cat(sprintf("p-value at gamma = 0: %.6f  (%s)\n", p0, if (p0 < 0.05) "below 0.05" else "NOT below 0.05"))

## ---- N = 150: (3) exact.score.sen.IxJ at gamma = 1 (Theorem 2, sweep U_O), timed, to CSV
run_sensitivity(tb, w, v, delta, p0, sprintf("ordinal test, %dx%d, N = %d", I, J, sum(tb)))

## ---- N = 180: (1) observed table   (p at gamma = 0: 0.00126) --------------------
tb <- matrix(c(28, 19, 14, 20, 19, 19, 14, 19, 28), nrow = 3, byrow = TRUE)
dimnames(tb) <- list(paste0("Trt", 1:I), paste0("Out", 1:J)); storage.mode(tb) <- "integer"
show_setup(tb, w, v, delta, sprintf("ordinal test, %dx%d, N = %d", I, J, sum(tb)))

## ---- N = 180: (2) p-value at gamma = 0 (exact.score.sen.IxJ) ----------------
p0 <- exact.score.sen.IxJ(obs.table = tb, gamma = 0, delta = delta, row = "treatment",
                          treatment.scores = w, outcome.scores = v, verbose = FALSE)$max.prob
cat(sprintf("p-value at gamma = 0: %.6f  (%s)\n", p0, if (p0 < 0.05) "below 0.05" else "NOT below 0.05"))

## ---- N = 180: (3) exact.score.sen.IxJ at gamma = 1 (Theorem 2, sweep U_O), timed, to CSV
run_sensitivity(tb, w, v, delta, p0, sprintf("ordinal test, %dx%d, N = %d", I, J, sum(tb)))

## ---- N = 200: (1) observed table   (p at gamma = 0: 0.0008933) --------------------
tb <- matrix(c(30, 22, 15, 22, 22, 22, 15, 22, 30), nrow = 3, byrow = TRUE)
dimnames(tb) <- list(paste0("Trt", 1:I), paste0("Out", 1:J)); storage.mode(tb) <- "integer"
show_setup(tb, w, v, delta, sprintf("ordinal test, %dx%d, N = %d", I, J, sum(tb)))

## ---- N = 200: (2) p-value at gamma = 0 (exact.score.sen.IxJ) ----------------
p0 <- exact.score.sen.IxJ(obs.table = tb, gamma = 0, delta = delta, row = "treatment",
                          treatment.scores = w, outcome.scores = v, verbose = FALSE)$max.prob
cat(sprintf("p-value at gamma = 0: %.6f  (%s)\n", p0, if (p0 < 0.05) "below 0.05" else "NOT below 0.05"))

## ---- N = 200: (3) exact.score.sen.IxJ at gamma = 1 (Theorem 2, sweep U_O), timed, to CSV
run_sensitivity(tb, w, v, delta, p0, sprintf("ordinal test, %dx%d, N = %d", I, J, sum(tb)))

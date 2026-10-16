# =============================================================================
# general_test_feasibility_3x3.R
#
# Code for "Exact, Nonparametric Sensitivity Analysis for Observational Studies
# of Contingency Tables"
# Authors: Elaine K. Chiu and Hyunseung Kang
#
# Running-time and memory benchmark of the exact worst-case p-value for the
# chi-square test on 3 x 3 tables, N in {20, 50, 100, 120, 150, 180, 200}.
#
# Setting
#   Statistic  chi-square, T = sum_ij (N_ij - E_ij)^2 / E_ij with
#              E_ij = N_i. N_.j / N. The statistic is not arrangement-
#              increasing, so Theorems 2 and 3 do not apply and the worst case
#              is found by Theorem 1. (exact.general.sen.IxJ also offers
#              statistic = "G2" and "freeman_halton"; this script uses the
#              chi-square only.)
#   delta      (0, 1, 1); gamma in {0, 1}.
#   Worst case Theorem 1: the maximum over the full candidate set U_C of all
#              u-bar vectors, |U_C| = prod_j (N_.j + 1), which is the default
#              u_space of exact.general.sen.IxJ. When |U_C| exceeds CHUNK the
#              sweep is split into checkpointed chunks (see run_sensitivity),
#              so a long run can be resumed after an interruption.
#   Tables     one fixed table per N, with a positive association between
#              treatment and outcome and exact p-value < 0.05 at gamma = 0.
#              They differ from the tables of the ordinal scripts because the
#              chi-square test has less power against a monotone association.
#
# For each N the script prints the table, the exact p-value at gamma = 0
# (exact.general.sen.IxJ evaluated at u = 0), and the timed worst-case
# computation at gamma = 1 with its memory use, and appends one row to
# feasibility/general_test_feasibility.csv. The CSV columns are the same in
# all four general_test_feasibility_*.R scripts and are documented in
# README.txt.
#
# Running time. |U_C| reaches 309,120 at N = 200 and each candidate costs
# O(N^{J+1}), so the blocks for N >= 150 take about an hour or more on 4 CPUs.
# Run the script in the background; the checkpoints in CKPT_DIR let an
# interrupted block resume where it stopped.
#
# Run from the top-level folder (the one with README.txt):
#   source("feasibility/general_test_feasibility_3x3.R")
# =============================================================================
library(Rcpp)

# ---- Setup ------------------------------------------------------------------
sourceCpp("library/sensitivityIxJ.cpp")
source("library/sensitivityIxJ.R")

CSV_FILE <- "feasibility/general_test_feasibility.csv"
CKPT_DIR <- "feasibility/ckpt_general"; dir.create(CKPT_DIR, showWarnings = FALSE)
CHUNK    <- 200000L                   # candidates per exact.general.sen.IxJ call when |U_C| is large
I <- 3; J <- 3
delta <- c(0, rep(1, I - 1))          # level 1 is the reference level
gamma1 <- 1                           # the sensitivity parameter at which the worst case is computed

# ---- Helpers ----------------------------------------------------------------
# Print the setup of one case.
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

# Memory use is read from the operating system, because R's own counters do
# not see the C++ allocations inside build_W_cpp. On Linux the file
# /proc/self/status reports, for this R process, the RAM it is using now
# (VmRSS) and the most RAM it has used so far (VmHWM). Before each run the
# "most so far" counter is reset by writing 5 to /proc/self/clear_refs, the
# Linux command for that reset; then
#     VmHWM after the run  -  VmRSS before the run
# is the most RAM the run needed above its starting point, even if part of it
# was freed again before the run finished. If the reset is refused, the
# session-wide maximum is used and flagged (peak_is_per_run = FALSE). On
# systems without /proc (macOS, Windows) the memory columns are NA; the
# timings are unaffected.
mem_mb <- function() {
  st <- tryCatch(readLines("/proc/self/status"), error = function(e) NULL)
  if (is.null(st)) return(c(rss = NA_real_, peak = NA_real_))
  get <- function(key) as.numeric(sub("^\\S+:\\s+(\\d+).*$", "\\1", grep(key, st, value = TRUE))) / 1024
  c(rss = get("^VmRSS"), peak = get("^VmHWM"))
}
reset_peak <- function() {
  tryCatch({ writeLines("5", "/proc/self/clear_refs"); TRUE }, error = function(e) FALSE, warning = function(w) FALSE)
}

# Timed worst case at gamma = gamma1 via exact.general.sen.IxJ; appends one
# row to CSV_FILE. p0 is the gamma = 0 p-value from block (2).
#
# Chunking. When |U_C| exceeds CHUNK, the sweep is split into blocks of CHUNK
# candidates, each a separate call with u_space = full.candidate.set.rows(...),
# and each block's maximum is saved as a .rds file in CKPT_DIR. Re-running
# after an interruption resumes from the finished chunks. Stage 1 is repeated
# by every call and reported once; stage 2 is summed over calls.
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

# =============================================================================
# Runs. For each N, three blocks, run in order: (1) the observed table;
# (2) the exact p-value at gamma = 0; (3) the timed worst case at gamma = 1,
# appended to the CSV.
# =============================================================================

# ---- N = 20: (1) observed table   (|U_C| = 360) -----------------------------
tb <- matrix(c(7, 1, 0, 2, 1, 1, 0, 1, 7), nrow = 3, byrow = TRUE)
dimnames(tb) <- list(paste0("Trt", 1:I), paste0("Out", 1:J)); storage.mode(tb) <- "integer"
show_setup(tb, delta, sprintf("chi-square test, %dx%d, N = %d", I, J, sum(tb)))

# ---- N = 20: (2) p-value at gamma = 0 (exact.general.sen.IxJ at u = 0) ------
p0 <- exact.general.sen.IxJ(u_space = matrix(0L, 1, J), obs.table = tb, gamma = 0, delta = delta,
                            row = "treatment", statistic = "chisq")$max.prob
cat(sprintf("p-value at gamma = 0: %.6f  (%s)\n", p0, if (p0 < 0.05) "below 0.05" else "NOT below 0.05"))

# ---- N = 20: (3) exact.general.sen.IxJ at gamma = 1 over U_C (Theorem 1), timed, to CSV
run_sensitivity(tb, delta, p0, sprintf("chi-square test, %dx%d, N = %d", I, J, sum(tb)))

# ---- N = 50: (1) observed table   (|U_C| = 5,320) ---------------------------
tb <- matrix(c(12, 5, 2, 5, 4, 4, 2, 4, 12), nrow = 3, byrow = TRUE)
dimnames(tb) <- list(paste0("Trt", 1:I), paste0("Out", 1:J)); storage.mode(tb) <- "integer"
show_setup(tb, delta, sprintf("chi-square test, %dx%d, N = %d", I, J, sum(tb)))

# ---- N = 50: (2) p-value at gamma = 0 (exact.general.sen.IxJ at u = 0) ------
p0 <- exact.general.sen.IxJ(u_space = matrix(0L, 1, J), obs.table = tb, gamma = 0, delta = delta,
                            row = "treatment", statistic = "chisq")$max.prob
cat(sprintf("p-value at gamma = 0: %.6f  (%s)\n", p0, if (p0 < 0.05) "below 0.05" else "NOT below 0.05"))

# ---- N = 50: (3) exact.general.sen.IxJ at gamma = 1 over U_C (Theorem 1), timed, to CSV
run_sensitivity(tb, delta, p0, sprintf("chi-square test, %dx%d, N = %d", I, J, sum(tb)))

# ---- N = 100: (1) observed table   (|U_C| = 40,176) -------------------------
tb <- matrix(c(20, 10, 5, 10, 10, 10, 5, 10, 20), nrow = 3, byrow = TRUE)
dimnames(tb) <- list(paste0("Trt", 1:I), paste0("Out", 1:J)); storage.mode(tb) <- "integer"
show_setup(tb, delta, sprintf("chi-square test, %dx%d, N = %d", I, J, sum(tb)))

# ---- N = 100: (2) p-value at gamma = 0 (exact.general.sen.IxJ at u = 0) -----
p0 <- exact.general.sen.IxJ(u_space = matrix(0L, 1, J), obs.table = tb, gamma = 0, delta = delta,
                            row = "treatment", statistic = "chisq")$max.prob
cat(sprintf("p-value at gamma = 0: %.6f  (%s)\n", p0, if (p0 < 0.05) "below 0.05" else "NOT below 0.05"))

# ---- N = 100: (3) exact.general.sen.IxJ at gamma = 1 over U_C (Theorem 1), timed, to CSV
run_sensitivity(tb, delta, p0, sprintf("chi-square test, %dx%d, N = %d", I, J, sum(tb)))

# ---- N = 120: (1) observed table   (|U_C| = 68,413) -------------------------
tb <- matrix(c(23, 12, 7, 12, 12, 12, 7, 12, 23), nrow = 3, byrow = TRUE)
dimnames(tb) <- list(paste0("Trt", 1:I), paste0("Out", 1:J)); storage.mode(tb) <- "integer"
show_setup(tb, delta, sprintf("chi-square test, %dx%d, N = %d", I, J, sum(tb)))

# ---- N = 120: (2) p-value at gamma = 0 (exact.general.sen.IxJ at u = 0) -----
p0 <- exact.general.sen.IxJ(u_space = matrix(0L, 1, J), obs.table = tb, gamma = 0, delta = delta,
                            row = "treatment", statistic = "chisq")$max.prob
cat(sprintf("p-value at gamma = 0: %.6f  (%s)\n", p0, if (p0 < 0.05) "below 0.05" else "NOT below 0.05"))

# ---- N = 120: (3) exact.general.sen.IxJ at gamma = 1 over U_C (Theorem 1), timed, to CSV
run_sensitivity(tb, delta, p0, sprintf("chi-square test, %dx%d, N = %d", I, J, sum(tb)))

# ---- N = 150: (1) observed table   (|U_C| = 131,652) ------------------------
tb <- matrix(c(28, 15, 9, 16, 15, 15, 9, 15, 28), nrow = 3, byrow = TRUE)
dimnames(tb) <- list(paste0("Trt", 1:I), paste0("Out", 1:J)); storage.mode(tb) <- "integer"
show_setup(tb, delta, sprintf("chi-square test, %dx%d, N = %d", I, J, sum(tb)))

# ---- N = 150: (2) p-value at gamma = 0 (exact.general.sen.IxJ at u = 0) -----
p0 <- exact.general.sen.IxJ(u_space = matrix(0L, 1, J), obs.table = tb, gamma = 0, delta = delta,
                            row = "treatment", statistic = "chisq")$max.prob
cat(sprintf("p-value at gamma = 0: %.6f  (%s)\n", p0, if (p0 < 0.05) "below 0.05" else "NOT below 0.05"))

# ---- N = 150: (3) exact.general.sen.IxJ at gamma = 1 over U_C (Theorem 1), timed, to CSV
run_sensitivity(tb, delta, p0, sprintf("chi-square test, %dx%d, N = %d", I, J, sum(tb)))

# ---- N = 180: (1) observed table   (|U_C| = 226,233) ------------------------
tb <- matrix(c(32, 19, 11, 19, 19, 19, 11, 18, 32), nrow = 3, byrow = TRUE)
dimnames(tb) <- list(paste0("Trt", 1:I), paste0("Out", 1:J)); storage.mode(tb) <- "integer"
show_setup(tb, delta, sprintf("chi-square test, %dx%d, N = %d", I, J, sum(tb)))

# ---- N = 180: (2) p-value at gamma = 0 (exact.general.sen.IxJ at u = 0) -----
p0 <- exact.general.sen.IxJ(u_space = matrix(0L, 1, J), obs.table = tb, gamma = 0, delta = delta,
                            row = "treatment", statistic = "chisq")$max.prob
cat(sprintf("p-value at gamma = 0: %.6f  (%s)\n", p0, if (p0 < 0.05) "below 0.05" else "NOT below 0.05"))

# ---- N = 180: (3) exact.general.sen.IxJ at gamma = 1 over U_C (Theorem 1), timed, to CSV
run_sensitivity(tb, delta, p0, sprintf("chi-square test, %dx%d, N = %d", I, J, sum(tb)))

# ---- N = 200: (1) observed table   (|U_C| = 309,120) ------------------------
tb <- matrix(c(35, 21, 13, 21, 21, 21, 13, 21, 34), nrow = 3, byrow = TRUE)
dimnames(tb) <- list(paste0("Trt", 1:I), paste0("Out", 1:J)); storage.mode(tb) <- "integer"
show_setup(tb, delta, sprintf("chi-square test, %dx%d, N = %d", I, J, sum(tb)))

# ---- N = 200: (2) p-value at gamma = 0 (exact.general.sen.IxJ at u = 0) -----
p0 <- exact.general.sen.IxJ(u_space = matrix(0L, 1, J), obs.table = tb, gamma = 0, delta = delta,
                            row = "treatment", statistic = "chisq")$max.prob
cat(sprintf("p-value at gamma = 0: %.6f  (%s)\n", p0, if (p0 < 0.05) "below 0.05" else "NOT below 0.05"))

# ---- N = 200: (3) exact.general.sen.IxJ at gamma = 1 over U_C (Theorem 1), timed, to CSV
run_sensitivity(tb, delta, p0, sprintf("chi-square test, %dx%d, N = %d", I, J, sum(tb)))
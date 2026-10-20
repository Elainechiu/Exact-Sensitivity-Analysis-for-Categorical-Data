# =============================================================================
# ordinal_feasibility_4x2.R
#
# Code for "Exact, Nonparametric Sensitivity Analysis for Observational Studies
# of Contingency Tables"
#
# Running-time and memory benchmark of the exact worst-case p-value for the
# ordinal test on 4 x 2 tables, N in {20, 50, 100, 120, 150, 180, 200}.
#
# Setting
#   Statistic  T = sum_ij w_i v_j N_ij with w = (1, ..., 4), v = (1, 2). With a
#              binary outcome this is a sign-score statistic: the term
#              sum_i w_i N_i1 is fixed by the margins.
#   delta      (0, 1, 1, 1); gamma in {0, 1}.
#   Worst case Theorem 3 gives it in closed form, u-bar+ = (0, N_.2);
#              exact.score.sen.IxJ evaluates u+ and the reference point u = 0.
#   Tables     one fixed table per N, with a positive association between
#              treatment and outcome and exact p-value < 0.05 at gamma = 0.
#
# For each N the script prints the table, the exact p-value at gamma = 0, and
# the timed worst-case computation at gamma = 1 with its memory use, and
# appends one row to feasibility/ordinal_test_feasibility.csv. The CSV columns
# are the same in all four ordinal_feasibility_*.R scripts and are documented
# in README.txt.
#
# Run from the top-level folder (the one with README.txt):
#   source("feasibility/ordinal_feasibility_4x2.R")
# =============================================================================
library(Rcpp)

# ---- Setup ------------------------------------------------------------------
sourceCpp("library/sensitivityIxJ.cpp")
source("library/sensitivityIxJ.R")

CSV_FILE <- "feasibility/ordinal_test_feasibility.csv"
I <- 4; J <- 2
w <- 1:I; v <- 1:J                    # scores = levels
delta <- c(0, rep(1, I - 1))          # level 1 is the reference level
gamma1 <- 1                           # the sensitivity parameter at which the worst case is computed
CANDIDATE_LABEL <- "u+ (Theorem 3)"
CANDIDATE_NOTE  <- function(tb) cat("Candidate set  : u+ in closed form, u-bar+ = (0, N_.2) = (0,", colSums(tb)[2], ") (Theorem 3); exact.score.sen.IxJ also evaluates the reference point (0,0)\n")

# ---- Helpers ----------------------------------------------------------------
# Print the setup of one case.
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

# Append one row to CSV_FILE (same columns in all four ordinal_feasibility_*.R scripts).
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

# Timed sensitivity analysis at gamma = gamma1 with the package function
# exact.score.sen.IxJ. The function chooses the candidate set from the outcome
# dimension: for J = 2 the closed-form worst case of Theorem 3 (u-bar+ = (0, N_.2),
# together with the reference point (0, 0)); for J >= 3 the N + 1 candidates of U_O
# (Theorem 2). Records time and memory and appends one row to CSV_FILE. p0 is the
# gamma = 0 p-value from block (2).
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

# =============================================================================
# Runs. For each N, three blocks, run in order: (1) the observed table;
# (2) the exact p-value at gamma = 0; (3) the timed worst case at gamma = 1,
# appended to the CSV.
# =============================================================================

# ---- N = 20: (1) observed table   (p at gamma = 0: 0.008043) ----------------
tb <- matrix(c(5, 1, 3, 1, 1, 3, 1, 5), nrow = 4, byrow = TRUE)
dimnames(tb) <- list(paste0("Trt", 1:I), paste0("Out", 1:J)); storage.mode(tb) <- "integer"
show_setup(tb, w, v, delta, sprintf("ordinal test, %dx%d, N = %d", I, J, sum(tb)))

# ---- N = 20: (2) p-value at gamma = 0 (exact.score.sen.IxJ) -----------------
p0 <- exact.score.sen.IxJ(obs.table = tb, gamma = 0, delta = delta, row = "treatment",
                          treatment.scores = w, outcome.scores = v, verbose = FALSE)$max.prob
cat(sprintf("p-value at gamma = 0: %.6f  (%s)\n", p0, if (p0 < 0.05) "below 0.05" else "NOT below 0.05"))

# ---- N = 20: (3) exact.score.sen.IxJ at gamma = 1 (Theorem 3, closed-form u+), timed, to CSV
run_sensitivity(tb, w, v, delta, p0, sprintf("ordinal test, %dx%d, N = %d", I, J, sum(tb)))

# ---- N = 50: (1) observed table   (p at gamma = 0: 0.002947) ----------------
tb <- matrix(c(10, 3, 7, 5, 5, 7, 3, 10), nrow = 4, byrow = TRUE)
dimnames(tb) <- list(paste0("Trt", 1:I), paste0("Out", 1:J)); storage.mode(tb) <- "integer"
show_setup(tb, w, v, delta, sprintf("ordinal test, %dx%d, N = %d", I, J, sum(tb)))

# ---- N = 50: (2) p-value at gamma = 0 (exact.score.sen.IxJ) -----------------
p0 <- exact.score.sen.IxJ(obs.table = tb, gamma = 0, delta = delta, row = "treatment",
                          treatment.scores = w, outcome.scores = v, verbose = FALSE)$max.prob
cat(sprintf("p-value at gamma = 0: %.6f  (%s)\n", p0, if (p0 < 0.05) "below 0.05" else "NOT below 0.05"))

# ---- N = 50: (3) exact.score.sen.IxJ at gamma = 1 (Theorem 3, closed-form u+), timed, to CSV
run_sensitivity(tb, w, v, delta, p0, sprintf("ordinal test, %dx%d, N = %d", I, J, sum(tb)))

# ---- N = 100: (1) observed table   (p at gamma = 0: 0.00177) ----------------
tb <- matrix(c(18, 8, 14, 10, 10, 14, 8, 18), nrow = 4, byrow = TRUE)
dimnames(tb) <- list(paste0("Trt", 1:I), paste0("Out", 1:J)); storage.mode(tb) <- "integer"
show_setup(tb, w, v, delta, sprintf("ordinal test, %dx%d, N = %d", I, J, sum(tb)))

# ---- N = 100: (2) p-value at gamma = 0 (exact.score.sen.IxJ) ----------------
p0 <- exact.score.sen.IxJ(obs.table = tb, gamma = 0, delta = delta, row = "treatment",
                          treatment.scores = w, outcome.scores = v, verbose = FALSE)$max.prob
cat(sprintf("p-value at gamma = 0: %.6f  (%s)\n", p0, if (p0 < 0.05) "below 0.05" else "NOT below 0.05"))

# ---- N = 100: (3) exact.score.sen.IxJ at gamma = 1 (Theorem 3, closed-form u+), timed, to CSV
run_sensitivity(tb, w, v, delta, p0, sprintf("ordinal test, %dx%d, N = %d", I, J, sum(tb)))

# ---- N = 120: (1) observed table   (p at gamma = 0: 0.002345) ---------------
tb <- matrix(c(21, 10, 16, 13, 13, 16, 10, 21), nrow = 4, byrow = TRUE)
dimnames(tb) <- list(paste0("Trt", 1:I), paste0("Out", 1:J)); storage.mode(tb) <- "integer"
show_setup(tb, w, v, delta, sprintf("ordinal test, %dx%d, N = %d", I, J, sum(tb)))

# ---- N = 120: (2) p-value at gamma = 0 (exact.score.sen.IxJ) ----------------
p0 <- exact.score.sen.IxJ(obs.table = tb, gamma = 0, delta = delta, row = "treatment",
                          treatment.scores = w, outcome.scores = v, verbose = FALSE)$max.prob
cat(sprintf("p-value at gamma = 0: %.6f  (%s)\n", p0, if (p0 < 0.05) "below 0.05" else "NOT below 0.05"))

# ---- N = 120: (3) exact.score.sen.IxJ at gamma = 1 (Theorem 3, closed-form u+), timed, to CSV
run_sensitivity(tb, w, v, delta, p0, sprintf("ordinal test, %dx%d, N = %d", I, J, sum(tb)))

# ---- N = 150: (1) observed table   (p at gamma = 0: 0.001228) ---------------
tb <- matrix(c(26, 13, 20, 16, 16, 20, 13, 26), nrow = 4, byrow = TRUE)
dimnames(tb) <- list(paste0("Trt", 1:I), paste0("Out", 1:J)); storage.mode(tb) <- "integer"
show_setup(tb, w, v, delta, sprintf("ordinal test, %dx%d, N = %d", I, J, sum(tb)))

# ---- N = 150: (2) p-value at gamma = 0 (exact.score.sen.IxJ) ----------------
p0 <- exact.score.sen.IxJ(obs.table = tb, gamma = 0, delta = delta, row = "treatment",
                          treatment.scores = w, outcome.scores = v, verbose = FALSE)$max.prob
cat(sprintf("p-value at gamma = 0: %.6f  (%s)\n", p0, if (p0 < 0.05) "below 0.05" else "NOT below 0.05"))

# ---- N = 150: (3) exact.score.sen.IxJ at gamma = 1 (Theorem 3, closed-form u+), timed, to CSV
run_sensitivity(tb, w, v, delta, p0, sprintf("ordinal test, %dx%d, N = %d", I, J, sum(tb)))

# ---- N = 180: (1) observed table   (p at gamma = 0: 0.001437) ---------------
tb <- matrix(c(30, 16, 24, 20, 20, 24, 16, 30), nrow = 4, byrow = TRUE)
dimnames(tb) <- list(paste0("Trt", 1:I), paste0("Out", 1:J)); storage.mode(tb) <- "integer"
show_setup(tb, w, v, delta, sprintf("ordinal test, %dx%d, N = %d", I, J, sum(tb)))

# ---- N = 180: (2) p-value at gamma = 0 (exact.score.sen.IxJ) ----------------
p0 <- exact.score.sen.IxJ(obs.table = tb, gamma = 0, delta = delta, row = "treatment",
                          treatment.scores = w, outcome.scores = v, verbose = FALSE)$max.prob
cat(sprintf("p-value at gamma = 0: %.6f  (%s)\n", p0, if (p0 < 0.05) "below 0.05" else "NOT below 0.05"))

# ---- N = 180: (3) exact.score.sen.IxJ at gamma = 1 (Theorem 3, closed-form u+), timed, to CSV
run_sensitivity(tb, w, v, delta, p0, sprintf("ordinal test, %dx%d, N = %d", I, J, sum(tb)))

# ---- N = 200: (1) observed table   (p at gamma = 0: 0.001026) ---------------
tb <- matrix(c(33, 18, 27, 22, 22, 27, 18, 33), nrow = 4, byrow = TRUE)
dimnames(tb) <- list(paste0("Trt", 1:I), paste0("Out", 1:J)); storage.mode(tb) <- "integer"
show_setup(tb, w, v, delta, sprintf("ordinal test, %dx%d, N = %d", I, J, sum(tb)))

# ---- N = 200: (2) p-value at gamma = 0 (exact.score.sen.IxJ) ----------------
p0 <- exact.score.sen.IxJ(obs.table = tb, gamma = 0, delta = delta, row = "treatment",
                          treatment.scores = w, outcome.scores = v, verbose = FALSE)$max.prob
cat(sprintf("p-value at gamma = 0: %.6f  (%s)\n", p0, if (p0 < 0.05) "below 0.05" else "NOT below 0.05"))

# ---- N = 200: (3) exact.score.sen.IxJ at gamma = 1 (Theorem 3, closed-form u+), timed, to CSV
run_sensitivity(tb, w, v, delta, p0, sprintf("ordinal test, %dx%d, N = %d", I, J, sum(tb)))
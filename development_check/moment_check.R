# =============================================================================
# moment_check.R
#
# Code for "Exact, Nonparametric Sensitivity Analysis for Observational Studies
# of Contingency Tables"
# Authors: Elaine K. Chiu and Hyunseung Kang
#
# Checks the exact means, variances and covariances of the cell counts returned
# by zr_moments_<I>_by_<J> in library/sensitivityIxJ.cpp against the same
# moments computed by brute force from the exact null distribution
# (point_prob_bd over every table with the observed margins). Each of the seven
# supported table dimensions is checked in its own block: the block prints both
# sets of means and variances at one u-bar, then compares all moments at every
# u-bar in the full candidate set for gamma = 0 and gamma = 0.7.
#
# Run from the top-level folder (the one with README.txt):
#   source("development_check/moment_check.R")
# Each check prints PASS or FAIL; a summary is printed at the end.
# =============================================================================
library(Rcpp)
# ---- Setup ------------------------------------------------------------------
sourceCpp("library/sensitivityIxJ.cpp")
source("library/sensitivityIxJ.R")
# Record and print the result of one check.
results <- logical(0)
check <- function(label, expr) {
  ok <- tryCatch(isTRUE(expr), error = function(e) FALSE)
  results[label] <<- ok
  cat(if (ok) "PASS  " else "FAIL  ", label, "\n", sep = "")
}
# Free cells of a table (rows 1 to I - 1, columns 1 to J - 1), read row by row.
free_cells <- function(tab) as.vector(t(tab[-nrow(tab), -ncol(tab), drop = FALSE]))
# Names of the free cells, in the same order: "N_11", "N_12", ...
free_cell_names <- function(tab) {
  I <- nrow(tab); J <- ncol(tab)
  as.vector(t(outer(1:(I - 1), 1:(J - 1), function(i, j) paste0("N_", i, j))))
}
# Every table with the same margins as tab.
all_tables <- function(tab) {
  lapply(possible.table(0, tab, "greater than", function(t) 1), as.matrix)
}
# Moments of the free cells by brute force: each table is weighted by its exact
# probability. Covariances are listed for the pairs (1,2), (1,3), ..., (2,3), ...
brute_force_moments <- function(tables, delta, u, gamma) {
  prob <- sapply(tables, function(t) point_prob_bd(t, delta, u, gamma))
  X    <- do.call(rbind, lapply(tables, free_cells))      # one row per table
  mu   <- colSums(prob * X)
  S    <- t(X) %*% (prob * X) - outer(mu, mu)             # covariance matrix
  list(mean = mu, var = diag(S), cov = S[lower.tri(S)])
}
# Moments of the free cells from zr_moments_<I>_by_<J>, which returns the
# means, then the variances, then the covariances (same pair order).
function_moments <- function(tab, delta, u, gamma) {
  n <- rowSums(tab); m <- colSums(tab); N <- sum(tab); g <- gamma * delta
  out <- switch(paste0(nrow(tab), "x", ncol(tab)),
                "2x2" = zr_moments_two_by_two(n[1], n[2], m[1], m[2], N, g, u[1], u[2], 1),
                "2x3" = zr_moments_two_by_three(n[1], n[2], m[1], m[2], m[3], N, g, u[1], u[2], u[3], 1),
                "2x4" = zr_moments_two_by_four(n[1], n[2], m[1], m[2], m[3], m[4], N, g, u[1], u[2], u[3], u[4], 1),
                "2x5" = zr_moments_two_by_five(n[1], n[2], m[1], m[2], m[3], m[4], m[5], N, g, u[1], u[2], u[3], u[4], u[5], 1),
                "3x2" = zr_moments_three_by_two(n[1], n[2], n[3], m[1], m[2], N, g, u[1], u[2], 1),
                "3x3" = zr_moments_three_by_three(n[1], n[2], n[3], m[1], m[2], m[3], N, g, u[1], u[2], u[3], 1),
                "4x2" = zr_moments_four_by_two(n[1], n[2], n[3], n[4], m[1], m[2], N, g, u[1], u[2], 1))
  out <- as.numeric(out)
  K <- length(free_cells(tab))
  if (length(out) != 2 * K + K * (K - 1) / 2) stop("unexpected number of moments")
  list(mean = out[1:K], var = out[K + 1:K], cov = out[-(1:(2 * K))])
}
# Print the brute-force and function means and variances side by side.
compare_moments <- function(tab, delta, u, gamma) {
  bf <- brute_force_moments(all_tables(tab), delta, u, gamma)
  fm <- function_moments(tab, delta, u, gamma)
  cat(sprintf("\n%d x %d table, gamma = %g, u-bar = (%s)\n",
              nrow(tab), ncol(tab), gamma, toString(u)))
  print(data.frame(cell       = free_cell_names(tab),
                   mean_brute = bf$mean, mean_zr = fm$mean,
                   var_brute  = bf$var,  var_zr  = fm$var),
        digits = 10, row.names = FALSE)
}
# Largest absolute difference between the two sets of moments over every
# u-bar in the full candidate set.
max_error <- function(tab, delta, gamma) {
  tables <- all_tables(tab)
  U   <- full.candidate.set(colSums(tab))
  err <- 0
  for (k in seq_len(nrow(U))) {
    bf  <- brute_force_moments(tables, delta, U[k, ], gamma)
    fm  <- function_moments(tab, delta, U[k, ], gamma)
    err <- max(err, abs(unlist(bf) - unlist(fm)))
  }
  err
}
# =============================================================================
# 2 x 2: zr_moments_two_by_two
# =============================================================================
tab_22   <- matrix(c(3, 2,
                     2, 3), nrow = 2, byrow = TRUE)
delta_22 <- c(0, 1)
compare_moments(tab_22, delta_22, u = c(2, 3), gamma = 0.7)
# Purpose:  means, variances and covariances of the free cells match brute force
# Input:    tab_22, delta = (0, 1), gamma = 0 and 0.7, every candidate u-bar
# Expected: largest absolute difference below 1e-8
check("zr_moments_two_by_two, gamma = 0",   max_error(tab_22, delta_22, 0)   < 1e-8)
check("zr_moments_two_by_two, gamma = 0.7", max_error(tab_22, delta_22, 0.7) < 1e-8)
# =============================================================================
# 2 x 3: zr_moments_two_by_three
# =============================================================================
tab_23   <- matrix(c(3, 2, 1,
                     1, 2, 3), nrow = 2, byrow = TRUE)
delta_23 <- c(0, 1)
compare_moments(tab_23, delta_23, u = c(1, 2, 3), gamma = 0.7)
# Purpose:  means, variances and covariances of the free cells match brute force
# Input:    tab_23, delta = (0, 1), gamma = 0 and 0.7, every candidate u-bar
# Expected: largest absolute difference below 1e-8
check("zr_moments_two_by_three, gamma = 0",   max_error(tab_23, delta_23, 0)   < 1e-8)
check("zr_moments_two_by_three, gamma = 0.7", max_error(tab_23, delta_23, 0.7) < 1e-8)
# =============================================================================
# 2 x 4: zr_moments_two_by_four
# =============================================================================
tab_24   <- matrix(c(2, 1, 2, 1,
                     1, 2, 1, 2), nrow = 2, byrow = TRUE)
delta_24 <- c(0, 1)
compare_moments(tab_24, delta_24, u = c(1, 2, 1, 2), gamma = 0.7)
# Purpose:  means, variances and covariances of the free cells match brute force
# Input:    tab_24, delta = (0, 1), gamma = 0 and 0.7, every candidate u-bar
# Expected: largest absolute difference below 1e-8
check("zr_moments_two_by_four, gamma = 0",   max_error(tab_24, delta_24, 0)   < 1e-8)
check("zr_moments_two_by_four, gamma = 0.7", max_error(tab_24, delta_24, 0.7) < 1e-8)
# =============================================================================
# 2 x 5: zr_moments_two_by_five
# =============================================================================
tab_25   <- matrix(c(2, 1, 1, 1, 2,
                     1, 2, 1, 2, 1), nrow = 2, byrow = TRUE)
delta_25 <- c(0, 1)
compare_moments(tab_25, delta_25, u = c(1, 2, 1, 2, 1), gamma = 0.7)
# Purpose:  means, variances and covariances of the free cells match brute force
# Input:    tab_25, delta = (0, 1), gamma = 0 and 0.7, every candidate u-bar
# Expected: largest absolute difference below 1e-8
check("zr_moments_two_by_five, gamma = 0",   max_error(tab_25, delta_25, 0)   < 1e-8)
check("zr_moments_two_by_five, gamma = 0.7", max_error(tab_25, delta_25, 0.7) < 1e-8)
# =============================================================================
# 3 x 2: zr_moments_three_by_two
# =============================================================================
tab_32   <- matrix(c(3, 1,
                     2, 2,
                     1, 3), nrow = 3, byrow = TRUE)
delta_32 <- c(0, 0, 1)
compare_moments(tab_32, delta_32, u = c(2, 4), gamma = 0.7)
# Purpose:  means, variances and covariances of the free cells match brute force
# Input:    tab_32, delta = (0, 0, 1), gamma = 0 and 0.7, every candidate u-bar
# Expected: largest absolute difference below 1e-8
check("zr_moments_three_by_two, gamma = 0",   max_error(tab_32, delta_32, 0)   < 1e-8)
check("zr_moments_three_by_two, gamma = 0.7", max_error(tab_32, delta_32, 0.7) < 1e-8)
# =============================================================================
# 3 x 3: zr_moments_three_by_three
# =============================================================================
tab_33   <- matrix(c(2, 1, 1,
                     1, 2, 1,
                     1, 1, 2), nrow = 3, byrow = TRUE)
delta_33 <- c(0, 0, 1)
compare_moments(tab_33, delta_33, u = c(1, 2, 3), gamma = 0.7)
# Purpose:  means, variances and covariances of the free cells match brute force
# Input:    tab_33, delta = (0, 0, 1), gamma = 0 and 0.7, every candidate u-bar
# Expected: largest absolute difference below 1e-8
check("zr_moments_three_by_three, gamma = 0",   max_error(tab_33, delta_33, 0)   < 1e-8)
check("zr_moments_three_by_three, gamma = 0.7", max_error(tab_33, delta_33, 0.7) < 1e-8)
# =============================================================================
# 4 x 2: zr_moments_four_by_two
# =============================================================================
tab_42   <- matrix(c(2, 1,
                     1, 2,
                     2, 1,
                     1, 2), nrow = 4, byrow = TRUE)
delta_42 <- c(0, 0, 1, 1)
compare_moments(tab_42, delta_42, u = c(2, 4), gamma = 0.7)
# Purpose:  means, variances and covariances of the free cells match brute force
# Input:    tab_42, delta = (0, 0, 1, 1), gamma = 0 and 0.7, every candidate u-bar
# Expected: largest absolute difference below 1e-8
check("zr_moments_four_by_two, gamma = 0",   max_error(tab_42, delta_42, 0)   < 1e-8)
check("zr_moments_four_by_two, gamma = 0.7", max_error(tab_42, delta_42, 0.7) < 1e-8)
# ---- Summary ----------------------------------------------------------------
cat("\n", sum(results), " / ", length(results), " checks passed\n", sep = "")
if (any(!results)) cat("Failed:\n", paste0("  ", names(results)[!results], collapse = "\n"), "\n", sep = "")
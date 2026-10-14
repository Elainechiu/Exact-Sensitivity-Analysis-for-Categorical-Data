# =============================================================================
# moment_check.R
#
# Code for "Exact, Nonparametric Sensitivity Analysis for Observational Studies
# of Contingency Tables"
# Authors: Elaine K. Chiu and Hyunseung Kang
#
# Checks the exact first and second moments of the cell counts (Theorem D.1,
# implemented in zr_moments_<I>_by_<J> and zu_moments_*_log of
# library/sensitivityIxJ.cpp) against moments computed by brute force from the
# exact null distribution of Theorem 4 (point_prob_bd over every table with
# the observed margins), for every u-bar total 0, 1, ..., N and gamma in
# {0, 0.7}, in each of the seven supported dimensions.
#
# Run from the top-level folder (the one with README.txt):
#   source("development_check/moment_check.R")
# Each check prints PASS or FAIL; a summary is printed at the end.
# =============================================================================
suppressMessages(library(Rcpp))

# ---- Setup ------------------------------------------------------------------
sourceCpp("library/sensitivityIxJ.cpp")
source("library/sensitivityIxJ.R")
set.seed(1)

# Record and print the result of one check.
results <- logical(0)
check <- function(label, expr) {
  ok <- tryCatch(isTRUE(expr), error = function(e) FALSE)
  results[label] <<- ok
  cat(if (ok) "PASS  " else "FAIL  ", label, "\n", sep = "")
}

# Every table with the margins of tb (threshold 0 on a constant statistic).
all_tables <- function(tb) possible.table(0, tb, "greater than", function(t) 1)

# The free cells (i, j) with i < I, j < J, in row-major order.
free_cells <- function(I, J) {
  out <- list()
  for (i in 1:(I - 1)) for (j in 1:(J - 1)) out[[length(out) + 1]] <- c(i, j)
  out
}

# Call the zr_moments function for an I x J table. Every zr function returns
# the means of the free cells (row-major), then their variances, then their
# covariances with the pairs in lexicographic order.
zr_call <- function(I, J, n, m, N, gd, u) {
  switch(paste(I, J),
         "2 2" = zr_moments_two_by_two(n[1], n[2], m[1], m[2], N, gd, u[1], u[2], 1),
         "2 3" = zr_moments_two_by_three(n[1], n[2], m[1], m[2], m[3], N, gd, u[1], u[2], u[3], 1),
         "2 4" = zr_moments_two_by_four(n[1], n[2], m[1], m[2], m[3], m[4], N, gd, u[1], u[2], u[3], u[4], 1),
         "2 5" = zr_moments_two_by_five(n[1], n[2], m[1], m[2], m[3], m[4], m[5], N, gd, u[1], u[2], u[3], u[4], u[5], 1),
         "3 2" = zr_moments_three_by_two(n[1], n[2], n[3], m[1], m[2], N, gd, u[1], u[2], 1),
         "3 3" = zr_moments_three_by_three(n[1], n[2], n[3], m[1], m[2], m[3], N, gd, u[1], u[2], u[3], 1),
         "4 2" = zr_moments_four_by_two(n[1], n[2], n[3], n[4], m[1], m[2], N, gd, u[1], u[2], 1))
}

# Call the zu_moments_*_log function for I treatment levels. The first I - 1
# entries of the result are E[Q_1], ..., E[Q_{I-1}].
zu_call <- function(I, n, N, gd, Us) {
  switch(I, NULL,
         zu_moments_two_treatment_log(n[1], n[2], N, gd, Us, 1),
         zu_moments_three_treatment_log(n[1], n[2], n[3], N, gd, Us, 1),
         zu_moments_four_treatment_log(n[1], n[2], n[3], n[4], N, gd, Us, 1))
}

# Spread a u-bar total over the outcome levels, respecting the margins, so
# that the u-bar vector is interior rather than at a corner.
split_u <- function(ubar, m) {
  u <- integer(length(m)); left <- ubar
  for (j in seq_along(m)) {
    u[j] <- min(m[j], ceiling(left / (length(m) - j + 1))); left <- left - u[j]
  }
  if (left > 0) for (j in seq_along(m)) {
    add <- min(left, m[j] - u[j]); u[j] <- u[j] + add; left <- left - add
  }
  u
}

# One small table per supported dimension: 2x2, 2x3, 2x4, 2x5, 3x2, 3x3, 4x2.
tables <- list(
  matrix(c(3,2, 2,3), 2, byrow = TRUE),
  matrix(c(3,2,1, 1,2,3), 2, byrow = TRUE),
  matrix(c(2,1,2,1, 1,2,1,2), 2, byrow = TRUE),
  matrix(c(2,1,1,1,2, 1,2,1,2,1), 2, byrow = TRUE),
  matrix(c(3,1, 2,2, 1,3), 3, byrow = TRUE),
  matrix(c(2,1,1, 1,2,1, 1,1,2), 3, byrow = TRUE),
  matrix(c(2,1, 1,2, 2,1, 1,2), 4, byrow = TRUE))


# =============================================================================
# zr_moments_<I>_by_<J>: means, variances and covariances of the free cells
# vs brute force over the exact distribution; zu_moments_*_log: E[Q_i] vs the
# brute-force cell means through the identity of Theorem D.1.
# =============================================================================

for (tb in tables) {
  I <- nrow(tb); J <- ncol(tb); n <- rowSums(tb); m <- colSums(tb); N <- sum(tb)
  delta <- c(rep(0, ceiling(I / 2)), rep(1, I - ceiling(I / 2)))     # monotone binary
  tbs <- all_tables(tb)
  cells <- lapply(tbs, function(t) as.vector(t(as.matrix(t))))        # row-major
  fc <- free_cells(I, J); idx <- sapply(fc, function(ij) (ij[1] - 1) * J + ij[2])
  Xf <- matrix(unlist(lapply(cells, function(v) v[idx])), ncol = length(idx), byrow = TRUE)
  
  for (gam in c(0, 0.7)) {
    worst_zr <- 0; worst_zu <- 0
    for (ubar in 0:N) {
      u <- split_u(ubar, m); gd <- gam * delta
      
      # brute-force mean vector and covariance matrix of the free cells
      p <- sapply(tbs, function(t) point_prob_bd(as.matrix(t) + 0L, as.integer(delta), as.integer(u), gam))
      mu <- colSums(Xf * p); Sig <- t(Xf) %*% (Xf * p) - outer(mu, mu)
      truth <- c(mu, diag(Sig)); K <- length(mu)
      if (K > 1) for (a in 1:(K - 1)) for (b in (a + 1):K) truth <- c(truth, Sig[a, b])
      
      got <- as.numeric(zr_call(I, J, n, m, N, gd, u))
      worst_zr <- max(worst_zr, abs(got - truth))
      
      # E[Q_i] from zu_moments: the moments of Q_i are not a function of the
      # table alone, so they are checked through the cell-mean identity
      #   E[N_ij] = rho_{j,1} E[Q_i] + rho_{j,0} (N_i. - E[Q_i]),
      # rho_{j,1} = u-bar_j / u-bar, rho_{j,0} = (N_.j - u-bar_j) / (N - u-bar),
      # which needs 0 < u-bar < N
      if (ubar > 0 && ubar < N) {
        EQ <- as.numeric(zu_call(I, n, N, gd, ubar))[1:(I - 1)]
        rho1 <- u / ubar; rho0 <- (m - u) / (N - ubar)
        EN <- outer(EQ, rho1) + outer(n[1:(I - 1)] - EQ, rho0)        # (I-1) x J
        worst_zu <- max(worst_zu, abs(as.vector(t(EN))[idx] - mu))
      }
    }
    
    # Purpose:  every element returned by zr_moments_<I>_by_<J> equals the
    #           brute-force value from the exact distribution
    # Input:    the I x J table above, delta monotone binary, gamma = gam,
    #           u-bar total 0, 1, ..., N spread over the outcome levels
    # Expected: largest absolute error below 1e-8 over all N + 1 u-bar
    check(sprintf("zr_moments %dx%d, gamma = %.1f: means, variances, covariances at u-bar = 0..%d (worst abs err %.1e)",
                  I, J, gam, N, worst_zr),
          worst_zr < 1e-8)
    
    # Purpose:  E[Q_i] from zu_moments_*_log reproduces the brute-force cell
    #           means through the Theorem D.1 identity
    # Input:    the same table and gamma, u-bar total 1, ..., N - 1
    # Expected: largest absolute error below 1e-8
    check(sprintf("zu_moments I = %d, gamma = %.1f: E[Q_i] via cell-mean identity at u-bar = 1..%d (worst abs err %.1e)",
                  I, gam, N - 1, worst_zu),
          worst_zu < 1e-8)
  }
}


# =============================================================================
# Summary
# =============================================================================
cat(sprintf("\n%d passed, %d failed\n", sum(results), sum(!results)))
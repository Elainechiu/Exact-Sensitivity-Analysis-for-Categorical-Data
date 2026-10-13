## moment_check.R -- every element returned by zr_moments_<I>_by_<J> vs brute force over
## the exact distribution of Theorem 4 (all tables with the margins, point_prob_bd), for
## u-bar totals 0, 1, 2, ..., N-1, N and gamma in {0, 0.7}; plus zu_moments_*_log vs
## brute-force moments of Q_i.  Prints the worst absolute error per (dimension, gamma, ubar).
suppressMessages(library(Rcpp))
args <- commandArgs(TRUE); cpp <- if (length(args)) args[1] else "sensitivityIxJ.cpp"
sourceCpp(cpp); source("sensitivityIxJ.R")

all_tables <- function(tb) possible.table(0, tb, "greater than", function(t) 1)
free_cells <- function(I, J) { out <- list(); for (i in 1:(I-1)) for (j in 1:(J-1)) out[[length(out)+1]] <- c(i, j); out }
## order of the returned vector of each zr function: means of free cells (row-major), variances, covariances (pairs in lexicographic order)
zr_call <- function(I, J, n, m, N, gd, u) {
  switch(paste(I, J),
    "2 2" = zr_moments_two_by_two(n[1],n[2],m[1],m[2],N,gd,u[1],u[2],1),
    "2 3" = zr_moments_two_by_three(n[1],n[2],m[1],m[2],m[3],N,gd,u[1],u[2],u[3],1),
    "2 4" = zr_moments_two_by_four(n[1],n[2],m[1],m[2],m[3],m[4],N,gd,u[1],u[2],u[3],u[4],1),
    "2 5" = zr_moments_two_by_five(n[1],n[2],m[1],m[2],m[3],m[4],m[5],N,gd,u[1],u[2],u[3],u[4],u[5],1),
    "3 2" = zr_moments_three_by_two(n[1],n[2],n[3],m[1],m[2],N,gd,u[1],u[2],1),
    "3 3" = zr_moments_three_by_three(n[1],n[2],n[3],m[1],m[2],m[3],N,gd,u[1],u[2],u[3],1),
    "4 2" = zr_moments_four_by_two(n[1],n[2],n[3],n[4],m[1],m[2],N,gd,u[1],u[2],1))
}
zu_call <- function(I, n, N, gd, Us) {
  switch(I, NULL,
    zu_moments_two_treatment_log(n[1],n[2],N,gd,Us,1),
    zu_moments_three_treatment_log(n[1],n[2],n[3],N,gd,Us,1),
    zu_moments_four_treatment_log(n[1],n[2],n[3],n[4],N,gd,Us,1))
}
## split a total ubar into a u-bar vector (spread over outcome levels, respecting margins)
split_u <- function(ubar, m) { u <- integer(length(m)); left <- ubar
  for (j in seq_along(m)) { u[j] <- min(m[j], ceiling(left / (length(m) - j + 1))); left <- left - u[j] }
  if (left > 0) for (j in seq_along(m)) { add <- min(left, m[j] - u[j]); u[j] <- u[j] + add; left <- left - add }
  u }

set.seed(1)
tables <- list(
  matrix(c(3,2, 2,3), 2, byrow = TRUE),
  matrix(c(3,2,1, 1,2,3), 2, byrow = TRUE),
  matrix(c(2,1,2,1, 1,2,1,2), 2, byrow = TRUE),
  matrix(c(2,1,1,1,2, 1,2,1,2,1), 2, byrow = TRUE),
  matrix(c(3,1, 2,2, 1,3), 3, byrow = TRUE),
  matrix(c(2,1,1, 1,2,1, 1,1,2), 3, byrow = TRUE),
  matrix(c(2,1, 1,2, 2,1, 1,2), 4, byrow = TRUE))
worst_all <- 0
for (tb in tables) {
  I <- nrow(tb); J <- ncol(tb); n <- rowSums(tb); m <- colSums(tb); N <- sum(tb)
  delta <- c(rep(0, ceiling(I/2)), rep(1, I - ceiling(I/2)))          # monotone binary
  tbs <- all_tables(tb); cells <- lapply(tbs, function(t) as.vector(t(as.matrix(t))))   # row-major
  fc <- free_cells(I, J); idx <- sapply(fc, function(ij) (ij[1]-1)*J + ij[2])
  Xf <- matrix(unlist(lapply(cells, function(v) v[idx])), ncol = length(idx), byrow = TRUE)
  for (gam in c(0, 0.7)) for (ubar in 0:N) {
    u <- split_u(ubar, m); gd <- gam * delta
    p <- sapply(tbs, function(t) point_prob_bd(as.matrix(t) + 0L, as.integer(delta), as.integer(u), gam))
    mu <- colSums(Xf * p); Sig <- t(Xf) %*% (Xf * p) - outer(mu, mu)
    truth <- c(mu, diag(Sig)); K <- length(mu)
    if (K > 1) for (a in 1:(K-1)) for (b in (a+1):K) truth <- c(truth, Sig[a, b])
    got <- as.numeric(zr_call(I, J, n, m, N, gd, u))
    err <- max(abs(got - truth))
    ## Q moments
    ## exact Q_i moments need the joint of (t, q); use the kernel formula: E[Q_i] via point probabilities is not
    ## available from tables alone, so check zu via the identity E[N_ij] = rho_j1 E[Q_i] + rho_j0 (N_i. - E[Q_i]).
    errq <- NA
    if (ubar > 0 && ubar < N) {
      zu <- as.numeric(zu_call(I, n, N, gd, ubar)); EQ <- zu[1:(I-1)]
      rho1 <- u / ubar; rho0 <- (m - u) / (N - ubar)
      EN <- outer(EQ, rho1) + outer(n[1:(I-1)] - EQ, rho0)       # (I-1) x J
      errq <- max(abs(as.vector(t(EN))[sapply(fc, function(ij) (ij[1]-1)*J + ij[2])] - mu))
    }
    worst_all <- max(worst_all, err, errq, na.rm = TRUE)
    flag <- if (err < 1e-8 && (is.na(errq) || errq < 1e-8)) "PASS" else "FAIL"
    cat(sprintf("%s %dx%d gamma=%.1f ubar=%2d u=(%s)  zr err %.1e  zu err %s\n", flag, I, J, gam, ubar,
                paste(u, collapse = ","), err, if (is.na(errq)) "  n/a " else sprintf("%.1e", errq)))
  }
}
cat(sprintf("\nWORST ABS ERROR OVERALL: %.2e\n", worst_all))

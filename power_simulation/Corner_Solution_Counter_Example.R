# =============================================================================
# Corner_Solution_Counter_Example.R
#
# Code for "Exact, Nonparametric Sensitivity Analysis for Observational Studies
# of Contingency Tables"
#
# Numerical counterexample of Appendix E.1 (see also Lemma 1 in Appendix B.2
# and its discussion in Appendix E.2). Under the generic bias model, Lemma 1
# shows that the worst-case unmeasured confounder of any count-based test lies
# at a vertex of [0,1]^N. This script considers the more general exposure-dose
# model
#     P(Z_s = i | F) proportional to exp(xi_i(x_s) + gamma * phi(i) * u_s),
# with phi(i) = i, and the ordinal test T = sum_s Z_s r_s, on six-subject
# examples where every treatment assignment can be enumerated:
#   Examples 1-2  binary outcome: the maximizer is still a vertex (consistent
#                 with Theorem 3, which also holds under phi(i) = i);
#   Examples 3-5  non-binary outcome: the maximizer found by a bounded
#                 quasi-Newton search is NOT a vertex, and its p-value exceeds
#                 the p-value at every one of the 2^6 = 64 vertices.
# Example 5 is the one reported in Appendix E.1: treatments (1,1,2,2,3,4),
# r = (1.4, 1.4, 2.1, 2.1, 3.5, 4.7), gamma = 2, c = 40, maximizer
# u+ = (0, 0, 0.4674, 0.4674, 0.9073, 1).
#
# The script is self-contained (it does not use library/); it needs the
# package `arrangements` for the enumeration of treatment assignments.
#
# Run from the top-level folder (the one with README.txt):
#   source("power_simulation/Corner_Solution_Counter_Example.R")
# Results are printed to the console.
# =============================================================================
library(arrangements)

# ---- Objective -------------------------------------------------------------
# All distinct treatment assignments with the given multiset of levels, one per row.
assignments <- function(levels) unique(permutations(levels))

# p-value (10) under the exposure-dose model with phi(i) = i, at confounder u:
#   sum_z 1{z'R >= c} exp(gamma z'u) / sum_z exp(gamma z'u).
sig.fun <- function(R, permute, u, gamma, critical.value) {
  stat   <- as.vector(permute %*% R)
  weight <- exp(gamma * as.vector(permute %*% u))
  sum(weight[stat >= critical.value]) / sum(weight)
}

# optim() minimizes, so return the negative p-value
obj_fun <- function(u, R, permute, gamma, critical.value) -sig.fun(R, permute, u, gamma, critical.value)

# ---- One example -----------------------------------------------------------
# Maximizes the p-value over [0,1]^N by L-BFGS-B from u = (0.5, ..., 0.5), then
# evaluates every vertex of [0,1]^N and reports whether the search ended at a
# vertex and whether any vertex attains a larger p-value.
run_example <- function(label, levels, R, gamma, critical.value) {
  permute <- assignments(levels)
  N <- length(levels)
  stat.range <- range(permute %*% R)              # attainable values of the statistic
  
  opt <- optim(par = rep(0.5, N), fn = obj_fun, R = R, permute = permute, gamma = gamma,
               critical.value = critical.value, method = "L-BFGS-B",
               lower = rep(0, N), upper = rep(1, N))
  u.hat <- opt$par
  p.hat <- -opt$value
  is.vertex <- all(abs(u.hat - round(u.hat)) < 1e-6)
  
  corners  <- as.matrix(expand.grid(rep(list(c(0, 1)), N)))
  corner.p <- apply(corners, 1, function(u) sig.fun(R, permute, u, gamma, critical.value))
  best.corner <- corners[which.max(corner.p), ]
  
  cat("\n", strrep("=", 72), "\n", label, "\n", strrep("=", 72), "\n", sep = "")
  cat("treatments:", levels, "   r:", R, "   gamma =", gamma, "   c =", critical.value, "\n")
  cat("attainable statistic values: [", stat.range[1], ",", stat.range[2], "]\n")
  cat("maximizer found by optim:", round(u.hat, 4), "\n")
  cat("p-value at maximizer:    ", signif(p.hat, 6), "\n")
  cat("best vertex:             ", best.corner, "  p =", signif(max(corner.p), 6), "\n")
  cat(if (is.vertex) "=> maximizer is a vertex of [0,1]^N"
      else if (p.hat > max(corner.p) + 1e-10)
        "=> maximizer is NOT a vertex and beats every vertex: corner solution fails"
      else "=> maximizer is not a vertex but a vertex attains the same p-value (not a counterexample)", "\n")
  invisible(list(u = u.hat, p = p.hat, corner.p = corner.p, is.vertex = is.vertex))
}

# ---- Examples --------------------------------------------------------------
## Binary outcome: maximizer at a vertex (Theorem 3 under phi(i) = i)
ex1 <- run_example("Example 1: binary outcome", levels = c(1,1,2,2,2,3), R = c(1,1,2,2,2,2), gamma = 2, critical.value = 19)
ex2 <- run_example("Example 2: binary outcome, other c and gamma", levels = c(1,1,2,2,2,3), R = c(1,1,2,2,2,2), gamma = 3, critical.value = 20)

## Non-binary integer outcome: maximizer not at a vertex
ex3 <- run_example("Example 3: four outcome levels, integer scores", levels = c(1,1,2,2,2,3), R = c(1,2,2,3,3,4), gamma = 2,   critical.value = 31)
ex4 <- run_example("Example 4: same, other c and gamma",            levels = c(1,1,2,2,2,3), R = c(1,2,2,3,3,4), gamma = 2.5, critical.value = 30)

## Non-binary, non-integer outcome: the example of Appendix E.1
ex5 <- run_example("Example 5 (Appendix E.1): four treatments, four outcome levels",
                   levels = c(1,1,2,2,3,4), R = c(1.4,1.4,2.1,2.1,3.5,4.7), gamma = 2, critical.value = 40)
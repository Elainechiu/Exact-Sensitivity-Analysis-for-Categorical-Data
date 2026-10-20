# =============================================================================
# Kernel_Advantage_Simulation.R
#
# Code for "Exact, Nonparametric Sensitivity Analysis for Observational Studies
# of Contingency Tables"
#
# Computational savings of the kernel representation of the exact null
# distribution (Theorem 4). Two parts:
#
#   Part 1  Exact p-value at a given u: running time of the kernel approach
#           (exact.general.sen.IxJ, which sums over the tables in Omega_t)
#           versus the original formulation (a sum over every treatment
#           assignment z in Z(N_I.), generated one at a time with multicool).
#           Reproduces Table 1 of the main text and Table F.1 of the
#           Supplementary Material; results are appended to
#           power_simulation/output.txt.
#
#   Part 2  Monte Carlo estimation of the same p-value: convergence of
#           importance sampling over tables (sequential importance sampling
#           of Chen, Diaconis, Holmes and Liu (2005) and Eisinger and Chen
#           (2017), proposal implemented in library/SIS_methods.cpp), with and
#           without the known normalizing constant, versus direct sampling of
#           treatment assignments. Reproduces the figures of Section F.2 of
#           the Supplementary Material, saved to power_simulation/.
#
# The script is sequential: the functions defined at the top of each part are
# used by every experiment below them. The permutation approach in Part 1
# enumerates up to 17 million assignments per run and takes several minutes
# per table; run the whole file in the background if needed.
#
# Run from the top-level folder (the one with README.txt):
#   source("power_simulation/Kernel_Advantage_Simulation.R")
# =============================================================================
library(Rcpp)
library(multicool)   # loopless multiset permutations (Williams, 2009)
library(ggplot2)

# ---- Setup ------------------------------------------------------------------
sourceCpp("library/SIS_methods.cpp")
sourceCpp("library/sensitivityIxJ.cpp")
source("library/sensitivityIxJ.R")

OUT_DIR <- "power_simulation"
OUT_TXT <- file.path(OUT_DIR, "output.txt")



###############################################################################
#
#   Part 1: exact p-value at one u, kernel approach vs permutation approach
#
#   Both functions return the same p-value (up to numerical precision); only
#   the running time differs. The permutation approach evaluates (7) of the
#   main text directly, one treatment assignment at a time; the kernel
#   approach evaluates (9) by summing over the tables in Omega_t.
#
###############################################################################

# -----------------------------------------------------------------------------
# time_permute_pvalue: the original formulation (7). Enumerates every
# treatment assignment with the observed treatment margins, computes the
# weight exp(gamma * delta' q) and the table it induces, and accumulates the
# numerator and denominator of (7).
#   table              observed table, rows = treatment levels
#   u.vec              binary u_s for the N subjects, ordered by outcome level
#   transform.fun      count-based statistic, function(table) -> numeric
#   delta, gamma       sensitivity parameters
# -----------------------------------------------------------------------------
time_permute_pvalue <- function(table, u.vec, transform.fun, delta, gamma) {
  obs.stat <- transform.fun(table)
  treatment.margins <- rowSums(table)
  outcome.margins   <- colSums(table)
  if (length(u.vec) != sum(table)) stop("the length of u.vec is incorrect")
  n_treat   <- length(treatment.margins)
  n_outcome <- length(outcome.margins)
  
  treatment.init <- rep(seq_len(n_treat), treatment.margins)      # 1,...,1, 2,...,2, ...
  outcome        <- rep(seq_len(n_outcome), outcome.margins)      # subjects sorted by outcome
  u_int <- if (is.factor(u.vec)) as.integer(u.vec) - 1L else as.integer(u.vec)
  
  mc.object   <- multicool::initMC(treatment.init)
  total_perms <- factorial(length(treatment.init)) / prod(factorial(treatment.margins))
  
  numerator.summation   <- 0
  denominator.summation <- 0
  i <- 0L
  start.time <- Sys.time()
  while (i < total_perms) {
    perm_treat <- nextPerm(mc.object)                    # one assignment z in Z(N_I.)
    
    ## q_i = number of subjects with u_s = 1 assigned to treatment i
    u1.counts <- tabulate(perm_treat[u_int == 1L], nbins = n_treat)
    denominator.term <- exp(gamma * sum(delta * u1.counts))
    
    ## the table induced by z and r
    perm_table_mat <- matrix(0L, nrow = n_treat, ncol = n_outcome)
    for (idx in seq_along(perm_treat))
      perm_table_mat[perm_treat[idx], outcome[idx]] <- perm_table_mat[perm_treat[idx], outcome[idx]] + 1L
    
    numerator.term <- if (transform.fun(perm_table_mat) >= obs.stat) denominator.term else 0
    numerator.summation   <- numerator.summation + numerator.term
    denominator.summation <- denominator.summation + denominator.term
    i <- i + 1L
  }
  run.time <- difftime(Sys.time(), start.time, units = "secs")
  list(permute.pvalue = numerator.summation / denominator.summation,
       run.time = run.time, total_perms = total_perms)
}


# -----------------------------------------------------------------------------
# time_kernel_pvalue: the kernel formulation (9), as implemented in
# exact.general.sen.IxJ. Omega_t is enumerated in C++ (transform.fun is
# evaluated once per table) and the p-value at the single candidate
# u_allocation (u-bar coordinates) is returned.
# -----------------------------------------------------------------------------
time_kernel_pvalue <- function(transform.fun, obs.table, u_allocation, gamma, delta) {
  start.time <- Sys.time()
  sen.result <- exact.general.sen.IxJ(u_space = matrix(u_allocation, nrow = 1),
                                      obs.table = obs.table, gamma = gamma, delta = delta,
                                      transform.fun = transform.fun)
  run.time <- difftime(Sys.time(), start.time, units = "secs")
  list(run.time = run.time, p.value = sen.result$max.prob)
}


# -----------------------------------------------------------------------------
# exact_p_experiment: runs both approaches `run` times on one table and one u
# and returns the p-values and timings. u.vec is given per subject (subjects
# ordered by outcome level); its u-bar coordinates are derived here.
# -----------------------------------------------------------------------------
exact_p_experiment <- function(transform.fun, gamma, delta, obs.table, u.vec, run = 10) {
  outcome.margins <- colSums(obs.table)
  outcome <- rep(seq_along(outcome.margins), outcome.margins)
  u_int   <- if (is.factor(u.vec)) as.integer(u.vec) - 1L else as.integer(u.vec)
  u_allocation <- tabulate(outcome[u_int == 1L], nbins = length(outcome.margins))   # u-bar
  
  permute.time.all.run <- kernel.time.all.run <- permute.p.all.run <- kernel.p.all.run <- rep(NA, run)
  for (i in seq_len(run)) {
    permute.result <- time_permute_pvalue(table = obs.table, u.vec = u.vec,
                                          transform.fun = transform.fun, gamma = gamma, delta = delta)
    kernel.result  <- time_kernel_pvalue(transform.fun = transform.fun, obs.table = obs.table,
                                         u_allocation = u_allocation, gamma = gamma, delta = delta)
    permute.time.all.run[i] <- permute.result$run.time
    kernel.time.all.run[i]  <- kernel.result$run.time
    permute.p.all.run[i]    <- permute.result$permute.pvalue
    kernel.p.all.run[i]     <- kernel.result$p.value
  }
  list(u_allocation = u_allocation,
       permute.p.all.run = permute.p.all.run, kernel.p.all.run = kernel.p.all.run,
       permute.time.all.run = permute.time.all.run, kernel.time.all.run = kernel.time.all.run,
       mean.permute.time = mean(permute.time.all.run), sd.permute.time = sd(permute.time.all.run),
       mean.kernel.time  = mean(kernel.time.all.run),  sd.kernel.time  = sd(kernel.time.all.run))
}


# -----------------------------------------------------------------------------
# run_exact_block: one row of Table F.1. Prints the setup and the result to the
# console and appends them to power_simulation/output.txt.
# -----------------------------------------------------------------------------
run_exact_block <- function(label, obs.table, u.vec, transform.fun, gamma, delta, run = 10) {
  gc()
  res <- exact_p_experiment(transform.fun = transform.fun, gamma = gamma, delta = delta,
                            obs.table = obs.table, u.vec = u.vec, run = run)
  sink(file = OUT_TXT, append = TRUE)
  cat("\n", strrep("=", 70), "\n", label, "\n", strrep("=", 70), "\n", sep = "")
  print(obs.table)
  cat("u-bar =", res$u_allocation, "  gamma =", gamma, "  delta =", delta, "\n")
  print(res)
  sink()
  agree <- all(abs(res$kernel.p.all.run - res$permute.p.all.run) < 1e-8)
  cat(label, ": p =", signif(res$kernel.p.all.run[1], 3),
      "  kernel", signif(res$mean.kernel.time, 3), "s   permutation", signif(res$mean.permute.time, 4), "s",
      if (agree) "  [kernel = permutation: PASS]" else "  [kernel != permutation: FAIL]", "\n")
  invisible(res)
}


# ---- Part 1 experiments -----------------------------------------------------
# Ordinal statistic with scores (0,1,2) x (0,1,2) for the 3 x 3 tables and
# (0,1) x (0,1,2) for the 2 x 3 tables; gamma = 1; delta = (0,1,1) or (0,1).
# u.vec lists u_s for the subjects sorted by outcome level, so, e.g.,
# c(rep(0,12), rep(1,3)) on a table with outcome margins (2,5,8) gives
# u-bar = (0,0,3).

ordinal_stat <- function(w, v) { score <- outer(w, v); function(tb) sum(score * tb) }

## Experiment 1: 3 x 3, treatment margins (5,5,5), Table F.1 rows 1-3
tb1 <- matrix(c(2,3,0, 0,1,4, 0,1,4), nrow = 3, byrow = TRUE)
f3  <- ordinal_stat(c(0,1,2), c(0,1,2)); d3 <- c(0,1,1)
first.result  <- run_exact_block("First result",  tb1, factor(c(rep(0,12), rep(1,3)),  levels = 0:1), f3, 1, d3)
second.result <- run_exact_block("Second result", tb1, factor(c(rep(0,8),  rep(1,7)),  levels = 0:1), f3, 1, d3)
third.result  <- run_exact_block("Third result",  tb1, factor(c(rep(0,2),  rep(1,13)), levels = 0:1), f3, 1, d3)

## Experiment 2: another 3 x 3 with margins (5,5,5), larger rejection region, rows 4-6
tb2 <- matrix(c(2,2,1, 1,2,2, 1,2,2), nrow = 3, byrow = TRUE)
fourth.result <- run_exact_block("Fourth result", tb2, factor(c(rep(0,13), rep(1,2)),  levels = 0:1), f3, 1, d3)
fifth.result  <- run_exact_block("Fifth result",  tb2, factor(c(rep(0,7),  rep(1,8)),  levels = 0:1), f3, 1, d3)
sixth.result  <- run_exact_block("Sixth result",  tb2, factor(c(rep(0,4),  rep(1,11)), levels = 0:1), f3, 1, d3)

## Experiment 3: 3 x 3, treatment margins (6,6,6), rows 7-9 (= Table 1 of the main text)
tb3 <- matrix(c(3,2,1, 0,2,4, 0,1,5), nrow = 3, byrow = TRUE)
seventh.result <- run_exact_block("Seventh result", tb3, factor(c(rep(0,13), rep(1,5)),  levels = 0:1), f3, 1, d3)
eighth.result  <- run_exact_block("Eighth result",  tb3, factor(c(rep(0,8),  rep(1,10)), levels = 0:1), f3, 1, d3)
ninth.result   <- run_exact_block("Ninth result",   tb3, factor(c(rep(0,4),  rep(1,14)), levels = 0:1), f3, 1, d3)

## Experiment 4: another 3 x 3 with margins (6,6,6), rows 10-12
tb4 <- matrix(c(3,3,0, 1,2,3, 2,3,1), nrow = 3, byrow = TRUE)
tenth.result    <- run_exact_block("Tenth result",    tb4, factor(c(rep(0,14), rep(1,4)),  levels = 0:1), f3, 1, d3)
eleventh.result <- run_exact_block("Eleventh result", tb4, factor(c(rep(0,11), rep(1,7)),  levels = 0:1), f3, 1, d3)
twelfth.result  <- run_exact_block("Twelfth result",  tb4, factor(c(rep(0,6),  rep(1,12)), levels = 0:1), f3, 1, d3)

## Experiment 5: 2 x 3, treatment margins (10,10), rows 13-14
tb5 <- matrix(c(4,6,0, 1,3,6), nrow = 2, byrow = TRUE)
f2  <- ordinal_stat(c(0,1), c(0,1,2)); d2 <- c(0,1)
thirteenth.result <- run_exact_block("Thirteenth result", tb5, factor(c(rep(0,1),  rep(1,19)), levels = 0:1), f2, 1, d2)
fourteenth.result <- run_exact_block("Fourteenth result", tb5, factor(c(rep(0,13), rep(1,7)),  levels = 0:1), f2, 1, d2)

## Experiment 6: another 2 x 3, rows 15-16
tb6 <- matrix(c(2,4,4, 2,2,6), nrow = 2, byrow = TRUE)
fifteenth.result <- run_exact_block("Fifteenth result", tb6, factor(c(rep(0,1),  rep(1,19)), levels = 0:1), f2, 1, d2)
sixteenth.result <- run_exact_block("Sixteenth result", tb6, factor(c(rep(0,13), rep(1,7)),  levels = 0:1), f2, 1, d2)



###############################################################################
#
#   Part 2: Monte Carlo estimation of alpha(T, r, u)
#
#   Three estimators of the same p-value at one u:
#     SIS        importance sampling over tables. Tables with the observed
#                margins are drawn from the SIS-G proposal h(t) of Chen et al.
#                (2005) / Eisinger and Chen (2017) (sample_sis_one_table in
#                library/SIS_methods.cpp). Each draw t receives weight
#                P_u(t)/h(t), where P_u(t) is the exact null probability of
#                the table under the sensitivity model (Theorem 4, computed
#                by generic.I.by.J.sensitivity.point.probability), so the
#                average of 1{T(t) >= T_obs} P_u(t)/h(t) is unbiased for
#                alpha(T, r, u).
#     snSIS      the self-normalized version, sum 1{.} w / sum w, which does
#                not use the normalizing constant of P_u(.). This is the
#                estimator implemented in sampling.general.sen.IxJ.
#     PermTreat  importance sampling over treatment assignments: z is drawn
#                uniformly from Z(N_I.) and weighted by exp(gamma delta' q),
#                i.e. a Monte Carlo version of (7).
#   The exact p-value from exact.general.sen.IxJ is the reference line in
#   the convergence plots. All three samplers draw their randomness from R's
#   RNG, so set.seed() before each experiment makes the figures reproducible.
#
###############################################################################

# -----------------------------------------------------------------------------
# importance.sampling.SIS: SIS estimator with the normalizing constant.
#   table         observed table, rows = treatment levels
#   u_allocation  u-bar (length J)
# Returns the running estimate after each draw, the final estimate and the
# running time.
# -----------------------------------------------------------------------------
importance.sampling.SIS <- function(table, u_allocation, gamma, delta, transform.fun,
                                    mc.iteration = 5000) {
  start.time <- Sys.time()
  obs.stat <- transform.fun(table)
  treatment_margins <- rowSums(table)
  outcome_margins   <- colSums(table)
  all.weight.summation <- 0
  mc.probability <- numeric(mc.iteration)
  for (mc in seq_len(mc.iteration)) {
    draw <- sample_sis_one_table(row_s = treatment_margins, col_s = outcome_margins)
    mc.table  <- draw$table
    inverse.h <- draw$weight                       # 1 / h(t)
    if (transform.fun(mc.table) >= obs.stat) {
      p.u <- generic.I.by.J.sensitivity.point.probability(mc.table, u_allocation = u_allocation,
                                                          gamma = gamma, delta = delta)
      all.weight.summation <- all.weight.summation + inverse.h * p.u
    }
    mc.probability[mc] <- all.weight.summation / mc
  }
  list(probability.progression = mc.probability,
       final.prob.estimate = all.weight.summation / mc.iteration,
       run.time = difftime(Sys.time(), start.time, units = "secs"))
}


# -----------------------------------------------------------------------------
# importance.sampling.snis: self-normalized SIS estimator (same estimator as
# single.u.sampling.general.IxJ in the library, repeated here so that the
# running time is measured in the same way as the other two samplers).
# -----------------------------------------------------------------------------
importance.sampling.snis <- function(table, u_allocation, gamma, delta, transform.fun,
                                     mc.iteration = 5000) {
  start.time <- Sys.time()
  obs.stat <- transform.fun(table)
  treatment_margins <- rowSums(table)
  outcome_margins   <- colSums(table)
  numerator.summation <- denominator.summation <- 0
  mc.probability <- numeric(mc.iteration)
  for (mc in seq_len(mc.iteration)) {
    draw <- sample_sis_one_table(row_s = treatment_margins, col_s = outcome_margins)
    mc.table  <- draw$table
    inverse.h <- draw$weight
    p.u <- generic.I.by.J.sensitivity.point.probability(mc.table, u_allocation = u_allocation,
                                                        gamma = gamma, delta = delta)
    weight    <- inverse.h * p.u
    indicator <- as.integer(transform.fun(mc.table) >= obs.stat)
    numerator.summation   <- numerator.summation + indicator * weight
    denominator.summation <- denominator.summation + weight
    mc.probability[mc]    <- numerator.summation / denominator.summation
  }
  list(probability.progression = mc.probability,
       final.prob.estimate = numerator.summation / denominator.summation,
       run.time = difftime(Sys.time(), start.time, units = "secs"),
       total.iteration = mc.iteration)
}


# -----------------------------------------------------------------------------
# permutation.sampling: Monte Carlo version of (7). u.vec is per subject, with
# subjects ordered by outcome level, as in Part 1.
# -----------------------------------------------------------------------------
permutation.sampling <- function(table, u.vec, transform.fun, gamma, delta, mc.iteration = 5000) {
  obs.stat <- transform.fun(table)
  treatment.margins <- rowSums(table)
  outcome.margins   <- colSums(table)
  n_treat   <- length(treatment.margins)
  n_outcome <- length(outcome.margins)
  treatment.init <- rep(seq_len(n_treat), treatment.margins)
  outcome        <- rep(seq_len(n_outcome), outcome.margins)
  u_int <- if (is.factor(u.vec)) as.integer(u.vec) - 1L else as.integer(u.vec)
  N <- length(treatment.init)
  
  numerator.summation <- denominator.summation <- 0
  mc.probability <- numeric(mc.iteration)
  start.time <- Sys.time()
  for (mc in seq_len(mc.iteration)) {
    sample.treatment <- sample(treatment.init, size = N, replace = FALSE)     # uniform on Z(N_I.)
    u1.counts <- tabulate(sample.treatment[u_int == 1L], nbins = n_treat)     # q
    denominator.term <- exp(gamma * sum(delta * u1.counts))
    perm.table <- matrix(0L, n_treat, n_outcome)
    for (idx in seq_len(N))
      perm.table[sample.treatment[idx], outcome[idx]] <- perm.table[sample.treatment[idx], outcome[idx]] + 1L
    numerator.term <- as.integer(transform.fun(perm.table) >= obs.stat) * denominator.term
    numerator.summation   <- numerator.summation + numerator.term
    denominator.summation <- denominator.summation + denominator.term
    mc.probability[mc]    <- numerator.summation / denominator.summation
  }
  list(probability.progression = mc.probability,
       final.prob.estimate = numerator.summation / denominator.summation,
       run.time = difftime(Sys.time(), start.time, units = "secs"))
}


# -----------------------------------------------------------------------------
# run_experiment: exact p-value plus the three Monte Carlo estimators on one
# table at one u. u_allocation is u-bar; the per-subject u.vec used by the
# permutation sampler is derived from it (subjects sorted by outcome level,
# u_s = 1 for the last u-bar_j subjects at each level j).
# -----------------------------------------------------------------------------
run_experiment <- function(data.table, gamma, delta, u_allocation, transform.fun, mc.iteration = 10000) {
  outcome.margins <- colSums(data.table)
  u.vec <- unlist(lapply(seq_along(outcome.margins),
                         function(j) c(rep(0, outcome.margins[j] - u_allocation[j]), rep(1, u_allocation[j]))))
  
  start.time <- Sys.time()
  exact.p <- exact.general.sen.IxJ(u_space = matrix(u_allocation, nrow = 1), obs.table = data.table,
                                   gamma = gamma, delta = delta, transform.fun = transform.fun)$max.prob
  total.enumeration.time <- difftime(Sys.time(), start.time, units = "secs")
  
  sis  <- importance.sampling.SIS(table = data.table, u_allocation = u_allocation, gamma = gamma, delta = delta,
                                  transform.fun = transform.fun, mc.iteration = mc.iteration)
  snis <- importance.sampling.snis(table = data.table, u_allocation = u_allocation, gamma = gamma, delta = delta,
                                   transform.fun = transform.fun, mc.iteration = mc.iteration)
  perm <- permutation.sampling(table = data.table, u.vec = u.vec, transform.fun = transform.fun,
                               gamma = gamma, delta = delta, mc.iteration = mc.iteration)
  
  list(data.table = data.table, u_allocation = u_allocation, gamma = gamma, delta = delta,
       total.enumeration.time = total.enumeration.time,
       table.sample.SIS.time = sis$run.time, table.sample.snSIS.time = snis$run.time,
       permutation.sample.time = perm$run.time,
       exact.p = exact.p, table.sample.SIS.p = sis$final.prob.estimate,
       table.sample.snSIS.p = snis$final.prob.estimate, permutation.p = perm$final.prob.estimate,
       table.SIS.p.progression = sis$probability.progression,
       table.snSIS.p.progression = snis$probability.progression,
       permute.p.progression = perm$probability.progression)
}


# -----------------------------------------------------------------------------
# plot_convergence: running estimates of the three samplers against the exact
# p-value (dashed line).
# -----------------------------------------------------------------------------
plot_convergence <- function(result, skip.initial = 0, ylb = 0, yub = 1) {
  df <- rbind(
    data.frame(Step = seq_along(result$table.SIS.p.progression),   Estimate = result$table.SIS.p.progression,   Method = "SIS"),
    data.frame(Step = seq_along(result$table.snSIS.p.progression), Estimate = result$table.snSIS.p.progression, Method = "snSIS"),
    data.frame(Step = seq_along(result$permute.p.progression),     Estimate = result$permute.p.progression,     Method = "PermTreat"))
  df <- df[df$Step > skip.initial, ]
  df$Method <- factor(df$Method, levels = c("SIS", "snSIS", "PermTreat"))
  
  ggplot(df, aes(x = Step, y = Estimate, color = Method)) +
    geom_line(linewidth = 0.7) +
    geom_hline(yintercept = result$exact.p, linetype = "dashed", color = "black") +
    labs(y = "Estimated p-value", x = "Monte Carlo Iteration") +
    ylim(ylb, yub) +
    scale_color_manual(values = c("SIS" = "red", "snSIS" = "blue", "PermTreat" = "darkgreen"),
                       labels = c(expression(hat(alpha)[SIS]), expression(hat(alpha)[snSIS]),
                                  expression(hat(alpha)[PermTreat]))) +
    theme_minimal() +
    theme(legend.title = element_blank(), legend.position = "bottom",
          legend.spacing.x = unit(1, "cm"),
          axis.title = element_text(size = 15), axis.text = element_text(size = 15),
          legend.text = element_text(size = 15))
}


# -----------------------------------------------------------------------------
# run_sampling_block: one experiment of Section F.2. Saves the convergence
# plot to power_simulation/plot_convergence_experiment_<tag>.png, prints the
# summary (p-values and times, without the progressions) and appends it to
# power_simulation/output.txt.
# -----------------------------------------------------------------------------
run_sampling_block <- function(tag, data.table, gamma, delta, u_allocation, transform.fun,
                               ylb = 0, yub = 1, mc.iteration = 10000) {
  gc(); set.seed(100)
  res <- run_experiment(data.table = data.table, gamma = gamma, delta = delta,
                        u_allocation = u_allocation, transform.fun = transform.fun,
                        mc.iteration = mc.iteration)
  p <- plot_convergence(res, ylb = ylb, yub = yub)
  print(p)
  ggsave(file.path(OUT_DIR, paste0("plot_convergence_experiment_", tag, ".png")),
         plot = p, width = 8, height = 6, dpi = 1200)
  summary.res <- res[!grepl("progression", names(res))]
  sink(file = OUT_TXT, append = TRUE)
  cat("\n", strrep("=", 70), "\nSampling experiment ", tag, "\n", strrep("=", 70), "\n", sep = "")
  print(summary.res)
  sink()
  print(summary.res)
  invisible(res)
}


# ---- Part 2 experiments -----------------------------------------------------
# Ordinal statistic with scores (0,1,2.5) x (0,1,2); delta = (0,1,1).
# Naming: Experiment <table> - <u>.
f_F2 <- ordinal_stat(c(0, 1, 2.5), c(0, 1, 2))
delta <- c(0, 1, 1)

## Table 1 of F.2: outcome margins (24,36,26), N = 86 (matrix() fills column-major)
tbA <- matrix(c(12,6,6, 18,12,6, 5,6,15), ncol = 3)
Experiment.1.1.result <- run_sampling_block("1_1", tbA, gamma = 1, delta, u_allocation = c(0,10,20), f_F2, yub = 0.10)
Experiment.1.2.result <- run_sampling_block("1_2", tbA, gamma = 1, delta, u_allocation = c(0,36,26), f_F2, yub = 0.20)

## Table 2 of F.2: the Black girls table of Section 6, N = 94
tbB <- matrix(c(12,18,17, 3,12,25, 0,3,4), ncol = 3)
Experiment.2.1.result <- run_sampling_block("2_1", tbB, gamma = 0.5, delta, u_allocation = c(0,30,7), f_F2)
Experiment.2.2.result <- run_sampling_block("2_2", tbB, gamma = 0.5, delta, u_allocation = c(0,40,7), f_F2, yub = 0.10)
Experiment.2.3.result <- run_sampling_block("2_3", tbB, gamma = 1,   delta, u_allocation = c(0,30,5), f_F2, yub = 0.15)
Experiment.2.4.result <- run_sampling_block("2_4", tbB, gamma = 1,   delta, u_allocation = c(0,40,7), f_F2, yub = 0.25)

## Table 3 of F.2: the Black boys table of Section 6, N = 112
tbC <- matrix(c(10,29,20, 8,11,24, 1,3,6), ncol = 3)
Experiment.3.1.result <- run_sampling_block("3_1", tbC, gamma = 0.5, delta, u_allocation = c(0,20,10), f_F2, yub = 0.25)
Experiment.3.2.result <- run_sampling_block("3_2", tbC, gamma = 0.5, delta, u_allocation = c(0,30,10), f_F2, yub = 0.25)
Experiment.3.3.result <- run_sampling_block("3_3", tbC, gamma = 1,   delta, u_allocation = c(0,20,10), f_F2)
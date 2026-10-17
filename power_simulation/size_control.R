# =============================================================================
# size_control.R
#
# Code for "Exact, Nonparametric Sensitivity Analysis for Observational Studies
# of Contingency Tables"
# Authors: Elaine K. Chiu and Hyunseung Kang
#
# Size (type I error rate) of the sign-score test on I x 2 tables under the
# worst-case unmeasured confounder, comparing the exact worst-case p-value
# (exact.score.sen.IxJ, Theorem 3) with the normal approximation based on the
# exact first two moments of the cell counts (norm.score.sen.IxJ, Appendix D).
# Reproduces Figures 9-12 of Section F.8 of the Supplementary Material; the
# four PNG files are written to power_simulation/.
#
# Data generation. Under the sharp null and the generic bias model at the
# worst-case confounder u+ of Theorem 3, Corollary 1 states that the counts at
# outcome level 2, (N_12, ..., N_I2), follow a multivariate extended
# hypergeometric distribution with margins N_I., N_.2 and odds exp(gamma delta).
# Tables are therefore drawn with BiasedUrn::rMFNCHypergeo. The first block
# below checks this against direct sampling of treatment assignments from the
# model (both samplers give the same mean and sd of the statistic).
#
# For each setting, 1,000 tables are drawn (set.seed(100) inside the
# simulation function), the two worst-case p-values are computed for each, and
# the empirical rejection rate 1/1000 sum 1{p <= alpha} is plotted against
# alpha in [0, 1]; the 45-degree line is the uniform reference.
#
# Run from the top-level folder (the one with README.txt):
#   source("power_simulation/size_control.R")
# Each exact-vs-normal setting takes a few minutes; the whole file about an
# hour.
# =============================================================================
library(Rcpp)
library(BiasedUrn)    # rMFNCHypergeo: multivariate extended hypergeometric sampler
library(multicool)    # allPerm: all treatment assignments, for the equivalence check only
library(ggplot2)
library(patchwork)

# ---- Setup ------------------------------------------------------------------
sourceCpp("library/sensitivityIxJ.cpp")
source("library/sensitivityIxJ.R")

OUT_DIR <- "power_simulation"



###############################################################################
#
#   Two equivalent ways of sampling the sign-score statistic under the null
#   and the generic bias model at the worst-case confounder u+.
#   Form 1 follows the model: enumerate Z(N_I.) and sample z with probability
#   proportional to exp(gamma delta' q). Form 2 uses Corollary 1 and samples
#   (N_12, ..., N_I2) from the multivariate extended hypergeometric
#   distribution. Both return the sampled statistics and their mean and sd.
#
###############################################################################

## Form 1: direct sampling of treatment assignments (small N only)
biased.ordinal.test.multiple.treatment.binary.outcome.moment = function(treatment.margins, outcome.margins, gamma, delta, treatment.scores, outcome.scores, mc.iteration = 10000){
  if(length(treatment.margins) < 2) stop("this function requires at least two treatments")
  if(length(outcome.margins) != 2) stop("this function only discusses binary outcome")
  if(length(treatment.scores) != length(treatment.margins)) stop("treatment.scores must have the same length as treatment.margins")
  if(length(delta) != length(treatment.margins)) stop("delta must have the same length as treatment.margins")
  
  z.seed = rep(seq_along(treatment.margins), treatment.margins)
  all.permutation = allPerm(initMC(z.seed))                       # every z in Z(N_I.)
  
  worst.case.u = c(rep(0, outcome.margins[1]), rep(1, outcome.margins[2]))   # u+ of Theorem 3
  each.z.probability.weight = apply(all.permutation, 1, function(z) exp(gamma * sum(delta[z] * worst.case.u)))
  
  outcome.seed = c(rep(1, outcome.margins[1]), rep(2, outcome.margins[2]))
  transform.fun = function(tb) sum(tb * outer(treatment.scores, outcome.scores))
  
  obs.stat.vec = rep(NA, mc.iteration)
  for(mc in 1:mc.iteration){
    sampled.z = all.permutation[sample(nrow(all.permutation), size = 1, prob = each.z.probability.weight), ]
    obs.stat.vec[mc] = transform.fun(table(factor(sampled.z, levels = seq_along(treatment.margins)),
                                           factor(outcome.seed, levels = 1:2)))
  }
  list(obs.stat.vec = obs.stat.vec, obs.mean = mean(obs.stat.vec), obs.sd = sd(obs.stat.vec))
}

## Form 2: Corollary 1, multivariate extended hypergeometric (any N)
biased.hyper.ordinal.test.multiple.treatment.binary.outcome.moment = function(treatment.margins, outcome.margins, delta, gamma, treatment.scores, outcome.scores, mc.iteration = 10000){
  if(length(treatment.margins) < 2) stop("this function requires at least two treatments")
  if(length(outcome.margins) != 2) stop("this function only discusses binary outcome")
  if(length(delta) != length(treatment.margins)) stop("delta should have the same length as treatment.margins")
  if(length(treatment.scores) != length(treatment.margins)) stop("treatment.scores should have the same length as treatment.margins")
  if(any(diff(delta) < 0)) stop("delta should be nondecreasing")
  if(any(diff(treatment.scores) < 0)) stop("treatment.scores should be nondecreasing")
  
  transform.fun = function(tb) sum(tb * outer(treatment.scores, outcome.scores))
  weights = exp(gamma * delta)
  hyper.obs.stat.vec = rep(NA, mc.iteration)
  for(mc in 1:mc.iteration){
    trt.outcome.2 <- as.vector(rMFNCHypergeo(nran = 1, m = treatment.margins, n = outcome.margins[2], odds = weights))
    obs.table <- cbind(treatment.margins - trt.outcome.2, trt.outcome.2)
    hyper.obs.stat.vec[mc] = transform.fun(obs.table)
  }
  list(obs.stat.vec = hyper.obs.stat.vec, obs.mean = mean(hyper.obs.stat.vec), obs.sd = sd(hyper.obs.stat.vec))
}

## Equivalence check: the two samplers agree in mean and sd (Monte Carlo error only)
set.seed(100)
print(biased.ordinal.test.multiple.treatment.binary.outcome.moment(treatment.margins = c(2,2,3), outcome.margins = c(3,4), delta = c(0,0,1), gamma = 0.5, treatment.scores = c(0,1,2), outcome.scores = c(0,1), mc.iteration = 100000)[c("obs.mean", "obs.sd")])
print(biased.hyper.ordinal.test.multiple.treatment.binary.outcome.moment(treatment.margins = c(2,2,3), outcome.margins = c(3,4), delta = c(0,0,1), gamma = 0.5, treatment.scores = c(0,1,2), outcome.scores = c(0,1), mc.iteration = 100000)[c("obs.mean", "obs.sd")])
print(biased.ordinal.test.multiple.treatment.binary.outcome.moment(treatment.margins = c(2,3,3), outcome.margins = c(4,4), delta = c(0,1,1), gamma = 1, treatment.scores = c(0,2,3), outcome.scores = c(0,1), mc.iteration = 100000)[c("obs.mean", "obs.sd")])
print(biased.hyper.ordinal.test.multiple.treatment.binary.outcome.moment(treatment.margins = c(2,3,3), outcome.margins = c(4,4), delta = c(0,1,1), gamma = 1, treatment.scores = c(0,2,3), outcome.scores = c(0,1), mc.iteration = 100000)[c("obs.mean", "obs.sd")])



###############################################################################
#
#   Size simulation: draw tables at u+ (Form 2), compute the exact and the
#   normal-approximated worst-case p-value for each, and plot the empirical
#   rejection rate against the nominal level.
#
###############################################################################

## Checks shared by the two simulation functions
.check_size_inputs = function(treatment.margins, outcome.margins, delta, treatment.scores){
  if(length(treatment.margins) < 2) stop("this function requires at least two treatments")
  if(length(outcome.margins) != 2) stop("this function only discusses binary outcome")
  if(length(delta) != length(treatment.margins)) stop("delta should have the same length as treatment.margins")
  if(length(treatment.scores) != length(treatment.margins)) stop("treatment.scores should have the same length as treatment.margins")
  if(any(diff(delta) < 0)) stop("delta should be nondecreasing")
  if(any(diff(treatment.scores) < 0)) stop("treatment.scores should be nondecreasing")
  if(sum(treatment.margins) != sum(outcome.margins)) stop("treatment and outcome margins must sum to the same N")
}

## Rejection-rate plot. Each column of pmat is one method's p-values.
.size_plot = function(pmat, methods, colors, linetypes, base_size){
  alpha_seq <- seq(0, 1, length.out = nrow(pmat))
  plot_df <- do.call(rbind, lapply(seq_along(methods), function(k)
    data.frame(alpha = alpha_seq, type1err = sapply(alpha_seq, function(a) mean(pmat[, k] <= a)), Method = methods[k])))
  plot_df <- rbind(plot_df, data.frame(alpha = alpha_seq, type1err = alpha_seq, Method = "Uniform"))
  plot_df$Method <- factor(plot_df$Method, levels = c(methods, "Uniform"))
  ggplot(plot_df, aes(x = alpha, y = type1err, color = Method, linetype = Method)) +
    geom_line(linewidth = 1.4) +
    scale_color_manual(values = c(colors, "black")) +
    scale_linetype_manual(values = c(linetypes, "dashed")) +
    labs(x = "Nominal Significance Level", y = "Empirical Rejection Rate") +
    theme_minimal(base_size = base_size) +
    theme(legend.title = element_blank(),
          legend.text = element_text(size = base_size), axis.title = element_text(size = base_size),
          axis.text = element_text(size = base_size),
          legend.key.width = unit(2.5, "cm"), legend.position = "bottom")
}

## Exact vs normal approximation (sample sizes where the exact p-value is computed quickly)
normal.exact.multiple.treatment.binary.outcome.size.simu.fun = function(treatment.margins, outcome.margins, delta, gamma,
                                                                        treatment.scores, outcome.scores, mc.iteration = 10000,
                                                                        seed = 100, base_size = 20){
  set.seed(seed)
  .check_size_inputs(treatment.margins, outcome.margins, delta, treatment.scores)
  transform.fun = function(tb) sum(tb * outer(treatment.scores, outcome.scores))
  weights = exp(gamma * delta)
  u.plus  = matrix(c(0, outcome.margins[2]), nrow = 1)        # u+ in u-bar coordinates, for the normal approximation
  
  hyper.obs.stat.vec = exact.p.vec = normal.p.vec = rep(NA, mc.iteration)
  for(mc in 1:mc.iteration){
    trt.outcome.2 <- as.vector(rMFNCHypergeo(nran = 1, m = treatment.margins, n = outcome.margins[2], odds = weights))
    obs.table <- cbind(treatment.margins - trt.outcome.2, trt.outcome.2)
    hyper.obs.stat.vec[mc] = transform.fun(obs.table)
    exact.p.vec[mc]  = exact.score.sen.IxJ(obs.table = obs.table, gamma = gamma, delta = delta,
                                           treatment.scores = treatment.scores, outcome.scores = outcome.scores)$max.prob
    normal.p.vec[mc] = norm.score.sen.IxJ(obs.table = obs.table, gamma = gamma, delta = delta,
                                          treatment.scores = treatment.scores, outcome.scores = outcome.scores,
                                          u_space = u.plus)$max.prob
    if(mc %% 100 == 0) cat("  mc", mc, "of", mc.iteration, "\n")
  }
  p <- .size_plot(cbind(exact.p.vec, normal.p.vec), methods = c("Exact", "Normal"),
                  colors = c("#F8766D", "#00BFC4"), linetypes = c("solid", "longdash"), base_size = base_size)
  list(obs.stat.vec = hyper.obs.stat.vec, type.1.error.graph = p,
       worst.case.mean = mean(hyper.obs.stat.vec), worst.case.sd = sd(hyper.obs.stat.vec),
       exact.p.vec = exact.p.vec, normal.p.vec = normal.p.vec)
}

## Normal approximation only (large sample sizes)
normal.multiple.treatment.binary.outcome.size.simu.fun = function(treatment.margins, outcome.margins, delta, gamma,
                                                                  treatment.scores, outcome.scores, mc.iteration = 10000,
                                                                  seed = 100, base_size = 20){
  set.seed(seed)
  .check_size_inputs(treatment.margins, outcome.margins, delta, treatment.scores)
  transform.fun = function(tb) sum(tb * outer(treatment.scores, outcome.scores))
  weights = exp(gamma * delta)
  u.plus  = matrix(c(0, outcome.margins[2]), nrow = 1)
  
  hyper.obs.stat.vec = normal.p.vec = rep(NA, mc.iteration)
  for(mc in 1:mc.iteration){
    trt.outcome.2 <- as.vector(rMFNCHypergeo(nran = 1, m = treatment.margins, n = outcome.margins[2], odds = weights))
    obs.table <- cbind(treatment.margins - trt.outcome.2, trt.outcome.2)
    hyper.obs.stat.vec[mc] = transform.fun(obs.table)
    normal.p.vec[mc] = norm.score.sen.IxJ(obs.table = obs.table, gamma = gamma, delta = delta,
                                          treatment.scores = treatment.scores, outcome.scores = outcome.scores,
                                          u_space = u.plus)$max.prob
    if(mc %% 100 == 0) cat("  mc", mc, "of", mc.iteration, "\n")
  }
  p <- .size_plot(cbind(normal.p.vec), methods = "Normal",
                  colors = "#00BFC4", linetypes = "longdash", base_size = base_size)
  list(obs.stat.vec = hyper.obs.stat.vec, type.1.error.graph = p,
       worst.case.mean = mean(hyper.obs.stat.vec), worst.case.sd = sd(hyper.obs.stat.vec),
       normal.p.vec = normal.p.vec)
}

## Three panels side by side with a shared legend and x-axis label
combine_plots_patchwork <- function(plot1, plot2, plot3, titles = c("Scenario 1", "Scenario 2", "Scenario 3"), base_size = 18) {
  strip_y <- theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())
  p1 <- plot1 + theme(axis.title.x = element_blank()) + ggtitle(titles[1]) +
    theme(plot.title = element_text(hjust = 0.5, size = base_size)) + coord_fixed(ratio = 1)
  p2 <- plot2 + strip_y + theme(axis.title.x = element_blank()) + ggtitle(titles[2]) +
    theme(plot.title = element_text(hjust = 0.5, size = base_size)) + coord_fixed(ratio = 1)
  p3 <- plot3 + strip_y + theme(axis.title.x = element_blank()) + ggtitle(titles[3]) +
    theme(plot.title = element_text(hjust = 0.5, size = base_size)) + coord_fixed(ratio = 1)
  (p1 | p2 | p3) +
    plot_layout(guides = "collect") +
    plot_annotation(caption = "Nominal Significance Level",
                    theme = theme(plot.caption = element_text(hjust = 0.5, size = base_size))) &
    theme(legend.position = "right",
          legend.text = element_text(size = base_size), legend.title = element_text(size = base_size),
          legend.key.size = unit(0.3, "cm"), legend.key.width = unit(1.4, "cm"), legend.margin = margin(l = -10))
}

## One setting: run, print the summary, show the single-panel plot
run_size <- function(label, fun, treatment.margins, outcome.margins, gamma, delta = c(0,0,1),
                     treatment.scores = c(0,1,2), outcome.scores = c(0,1), mc.iteration = 1000) {
  cat("\n", label, ": treatment margins (", paste(treatment.margins, collapse = ","),
      "), outcome margins (", paste(outcome.margins, collapse = ","), "), gamma = ", gamma, "\n", sep = "")
  res <- fun(treatment.margins = treatment.margins, outcome.margins = outcome.margins, delta = delta, gamma = gamma,
             treatment.scores = treatment.scores, outcome.scores = outcome.scores, mc.iteration = mc.iteration, base_size = 18)
  cat("  worst-case mean", signif(res$worst.case.mean, 4), " sd", signif(res$worst.case.sd, 4), "\n")
  print(res$type.1.error.graph)
  invisible(res)
}



###############################################################################
#
#   Experiments. Settings 1-14 compare the exact and the normal-approximated
#   worst-case p-value; settings 15-20 (large N) report the normal
#   approximation only. delta = (0,0,1), w = (0,1,2), v = (0,1) throughout.
#
###############################################################################
EXACT  <- normal.exact.multiple.treatment.binary.outcome.size.simu.fun
NORMAL <- normal.multiple.treatment.binary.outcome.size.simu.fun

## Figure 9: margins (60,10,20); (15,75), gamma in {0, 1, 1.5}
size_control_res_1 <- run_size("Setting 1", EXACT, c(60,10,20), c(15,75), gamma = 0)
size_control_res_2 <- run_size("Setting 2", EXACT, c(60,10,20), c(15,75), gamma = 1)
size_control_res_3 <- run_size("Setting 3", EXACT, c(60,10,20), c(15,75), gamma = 1.5)

## Figure 11: margins (30,30,40); (15,85), gamma in {0, 1, 1.5}
size_control_res_4 <- run_size("Setting 4", EXACT, c(30,30,40), c(15,85), gamma = 0)
size_control_res_5 <- run_size("Setting 5", EXACT, c(30,30,40), c(15,85), gamma = 1)
size_control_res_6 <- run_size("Setting 6", EXACT, c(30,30,40), c(15,85), gamma = 1.5)

## Figure 12: margins growing proportionally from (60,20,2); (10,72), gamma = 1
size_control_res_7  <- run_size("Setting 7",  EXACT, c(60,20,2),    c(10,72),  gamma = 1)
size_control_res_8  <- run_size("Setting 8",  EXACT, c(120,40,4),   c(20,144), gamma = 1)   # not plotted
size_control_res_9  <- run_size("Setting 9",  EXACT, c(240,80,8),   c(40,288), gamma = 1)
size_control_res_10 <- run_size("Setting 10", EXACT, c(480,160,16), c(80,576), gamma = 1)

## Figure 10: margins growing proportionally from (20,10,10); (10,30), gamma = 1
size_control_res_11 <- run_size("Setting 11", EXACT, c(20,10,10),    c(10,30),   gamma = 1)
size_control_res_12 <- run_size("Setting 12", EXACT, c(40,20,20),    c(20,60),   gamma = 1)   # not plotted
size_control_res_13 <- run_size("Setting 13", EXACT, c(80,40,40),    c(40,120),  gamma = 1)
size_control_res_14 <- run_size("Setting 14", EXACT, c(200,100,100), c(100,300), gamma = 1)

## ---- Combined figures of Section F.8 ----------------------------------------
fig9 <- combine_plots_patchwork(size_control_res_1$type.1.error.graph, size_control_res_2$type.1.error.graph,
                                size_control_res_3$type.1.error.graph,
                                titles = c("\u03B3 = 0, (\u0393 = 1)", "\u03B3 = 1, (\u0393 = 2.718)", "\u03B3 = 1.5, (\u0393 = 4.482)"))
ggsave(file.path(OUT_DIR, "size_control_figure9_margins_60_10_20.png"), fig9, width = 13, height = 5, dpi = 300)

fig10 <- combine_plots_patchwork(size_control_res_11$type.1.error.graph, size_control_res_13$type.1.error.graph,
                                 size_control_res_14$type.1.error.graph,
                                 titles = c("(20,10,10); (10,30)", "(80,40,40); (40,120)", "(200,100,100); (100,300)"))
ggsave(file.path(OUT_DIR, "size_control_figure10_growing_margins_20_10_10.png"), fig10, width = 13, height = 5, dpi = 300)

fig11 <- combine_plots_patchwork(size_control_res_4$type.1.error.graph, size_control_res_5$type.1.error.graph,
                                 size_control_res_6$type.1.error.graph,
                                 titles = c("\u03B3 = 0, (\u0393 = 1)", "\u03B3 = 1, (\u0393 = 2.718)", "\u03B3 = 1.5, (\u0393 = 4.482)"))
ggsave(file.path(OUT_DIR, "size_control_figure11_margins_30_30_40.png"), fig11, width = 13, height = 5, dpi = 300)

fig12 <- combine_plots_patchwork(size_control_res_7$type.1.error.graph, size_control_res_9$type.1.error.graph,
                                 size_control_res_10$type.1.error.graph,
                                 titles = c("(60,20,2); (10,72)", "(240,80,8); (40,288)", "(480,160,16); (80,576)"))
ggsave(file.path(OUT_DIR, "size_control_figure12_growing_margins_60_20_2.png"), fig12, width = 13, height = 5, dpi = 300)



###############################################################################
#
#   Large sample sizes: normal approximation only. Not reported in the paper;
#   kept to show that the anti-conservativeness of the normal approximation
#   vanishes as N grows.
#
###############################################################################
size_control_res_15 <- run_size("Setting 15", NORMAL, c(35,225,240),     c(25,475),    gamma = 0)
size_control_res_16 <- run_size("Setting 16", NORMAL, c(35,225,240)*10,  c(25,475)*10, gamma = 1)
size_control_res_17 <- run_size("Setting 17", NORMAL, c(35,225,240),     c(25,475),    gamma = 1.5)
size_control_res_18 <- run_size("Setting 18", NORMAL, c(35,225,240)*10,  c(25,475)*10, gamma = 1.5)
size_control_res_19 <- run_size("Setting 19", NORMAL, c(110,195,195),    c(230,270),   gamma = 1)
size_control_res_20 <- run_size("Setting 20", NORMAL, c(110,195,195)*10, c(230,270)*10, gamma = 1)
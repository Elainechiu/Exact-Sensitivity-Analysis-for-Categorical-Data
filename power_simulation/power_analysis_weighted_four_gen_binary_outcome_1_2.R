#!/usr/bin/env Rscript
# =============================================================================
# power_analysis_weighted_four_gen_binary_outcome_1_2.R
#
# Code for "Exact, Nonparametric Sensitivity Analysis for Observational Studies
# of Contingency Tables"
# Authors: Elaine K. Chiu and Hyunseung Kang
#
# Power of the 3 x 2 (V1) test: the ordinal test with the data-generating
# treatment scores w* after collapsing outcome levels 1 and 2, i.e. outcome
# scores (0, 1) on the collapsed 3 x 2 table (equivalently v = (0, 0, 1) on the
# 3 x 3 table). All three treatment levels are kept, with delta = (0, 1, 1).
# Tables of size N = 60 (treatment margins 20, 20, 20) are drawn from the
# linear-by-linear association model; the exact worst-case p-value at the
# given gamma is computed with exact.score.sen.IxJ and the number of
# rejections at level 0.05 is printed. The data-generating parameters are set
# below and are identical to those in the other power_analysis_*.R scripts.
#
# Usage, from the top-level folder (the one with README.txt):
#   Rscript power_simulation/power_analysis_weighted_four_gen_binary_outcome_1_2.R <seed> <gamma> <replications>
# Replication k uses set.seed(seed + k - 1), so the same tables are generated
# by every power_analysis_*.R script given the same (seed, gamma). The last
# line of output, "how many rejections? <n>", is read by summarize_power.sh.
# =============================================================================
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) stop("usage: <seed> <gamma> <replications>")
arg1 <- args[1]   # starting seed
arg2 <- args[2]   # gamma
arg3 <- args[3]   # number of Monte Carlo replications in this job

library(Rcpp)
sourceCpp("library/sensitivityIxJ.cpp")
source("library/sensitivityIxJ.R")

#################################################################################
#                                                                               #
#             This is where to define the data generating process               #
#              by specifying the parameters in the log-linear models.           #
#                                                                               #
#################################################################################

## four weight model
# Define the row and column scores
w <- c(0, 1.7, 2.45)  # Row scores from treatment
v <- c(0, 1.25, 1.4)  # Column scores from outcome

# Specify model parameters
lambda <- 0  # Overall effect
lambda_Z <- c(1, 0, 0)  # Row (Z) effects
lambda_R <- c(1, 0.2, 0)  # Column (R) effects
beta <- 1  # Association parameter (positive association)

# Compute expected cell counts (mu_ij)
mu <- matrix(0, nrow = length(w), ncol = length(v))  # Matrix to store expected counts
for (i in 1:length(w)) {
  for (j in 1:length(v)) {
    mu[i, j] <- exp(lambda + lambda_Z[i] + lambda_R[j] + beta * w[i] * v[j])
  }
}

# Convert expected counts to conditional probabilities for each row (Z = i)
conditional_prob <- matrix(0, nrow = length(w), ncol = length(v))  # Matrix to store conditional probabilities
for (i in 1:length(w)) {
  row_sum <- sum(mu[i, ])  # Sum of expected counts for row i
  for (j in 1:length(v)) {
    conditional_prob[i, j] <- mu[i, j] / row_sum  # Normalize to get conditional probabilities
  }
}

## print model setup
cat("w",w,"\n")
cat("v",v,"\n")
cat("lambda",lambda,"\n")
cat("row effect", lambda_Z,"\n")
cat("column effect",lambda_R,"\n")
cat("effect size beta", beta,"\n")
out = outer(w,v)
print("outer product of scores")
print(out)

# Print conditional probabilities
cat("Conditional Probabilities P(R = j | Z = i):\n")
print(conditional_prob)


trt1.prob = conditional_prob[1,]
trt2.prob = conditional_prob[2,]
trt3.prob = conditional_prob[3,]




###############################################################################
#                                                                             #
#       Start the Monte Carlo Simulations                                     #
#                                                                             #
###############################################################################

gamma.rejection = 0

mc = as.numeric(arg1)
gamma = as.numeric(arg2)
end.count = as.numeric(arg3)
cat("starting seed number is:", mc,"\n")
cat("gamma in this file is:", gamma, "\n")
count = 0
while(count < end.count) {
  set.seed(mc)
  ## generate a table with treatment effect
  ## assume 20 people getting trt1, 20 getting 2, 20 getting trt3
  outcome.trt1 = sample(x=c("o1","o2","o3"),size=20,replace=T,prob=trt1.prob)
  outcome.trt2 = sample(x=c("o1","o2","o3"),size=20,replace=T,prob=trt2.prob)
  outcome.trt3 = sample(x=c("o1","o2","o3"),size=20,replace=T,prob=trt3.prob)
  n11=sum(outcome.trt1=="o1")
  n12=sum(outcome.trt1=="o2")
  n13=sum(outcome.trt1=="o3")
  n21=sum(outcome.trt2=="o1")
  n22=sum(outcome.trt2=="o2")
  n23=sum(outcome.trt2=="o3")
  n31=sum(outcome.trt3=="o1")
  n32=sum(outcome.trt3=="o2")
  n33=sum(outcome.trt3=="o3")
  ## collapse to 3 x 2: column 1 = outcome levels 1 and 2, column 2 = outcome level 3;
  ## all three treatment levels are kept
  observed.table = matrix(data=c(n11+n12,n13,n21+n22,n23,n31+n32,n33),nrow = 3,byrow=T)
  
  p.gamma = exact.score.sen.IxJ(obs.table = observed.table,treatment.scores = w, outcome.scores = c(0,1),
                                gamma = gamma,
                                delta=c(0,1,1),
                                row="treatment")$max.prob
  
  gamma.rejection = gamma.rejection + ifelse(p.gamma < 0.05, 1, 0)
  count = count + 1
  mc = mc + 1
  
}

cat("number of rejections for the 3 x 2 (V1) ordinal test (outcome levels 1-2 collapsed) with starting seed:",
    as.numeric(arg1), "total iterations", end.count, "gamma", gamma, "\n")
cat("how many rejections?",gamma.rejection,"\n")
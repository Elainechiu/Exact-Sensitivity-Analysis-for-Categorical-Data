## verify_same_names.R
## Previous release (sensitivityIxJ_orig.cpp + sensitivityIxJ_orig.R) vs the new engine under the
## SAME public names (sensitivityIxJ.cpp + sensitivityIxJ.R), on the paper's own tables:
##   (1) |Omega_t| from build_W_cpp == length(possible.table(...)),
##   (2) p-value at the given u equal, (3) timings.
## Files in the working directory: sensitivityIxJ_orig.cpp, sensitivityIxJ_orig.R,
##   sensitivityIxJ.cpp, sensitivityIxJ.R, SIS_methods_orig.cpp.  Then: source("verify_same_names.R")

library(Rcpp)
env_old <- new.env(); sourceCpp("sensitivityIxJ_orig.cpp", env = env_old); sys.source("sensitivityIxJ_orig.R", envir = env_old)
env_new <- new.env(); sourceCpp("sensitivityIxJ.cpp", env = env_new); sys.source("sensitivityIxJ.R", envir = env_new)
sourceCpp("SIS_methods_orig.cpp")   # table sampler used by sampling.*.sen.IxJ in both stacks (unchanged file)
tm <- function(expr) { t0 <- proc.time()[["elapsed"]]; v <- expr; list(v = v, s = proc.time()[["elapsed"]] - t0) }

cases <- list(
  list(src="F.1 block 1", tb=matrix(c(2,3,0, 0,1,4, 0,1,4),3,byrow=TRUE), u=c(0,0,3),  gamma=1,  delta=c(0,1,1), w=c(0,1,2),   v=c(0,1,2)),
  list(src="F.1 block 2", tb=matrix(c(2,2,1, 1,2,2, 1,2,2),3,byrow=TRUE), u=c(0,3,5),  gamma=1,  delta=c(0,1,1), w=c(0,1,2),   v=c(0,1,2)),
  list(src="F.1 block 3", tb=matrix(c(3,2,1, 0,2,4, 0,1,5),3,byrow=TRUE), u=c(0,4,10), gamma=1,  delta=c(0,1,1), w=c(0,1,2),   v=c(0,1,2)),
  list(src="F.1 block 4", tb=matrix(c(3,3,0, 1,2,3, 2,3,1),3,byrow=TRUE), u=c(0,8,4),  gamma=1,  delta=c(0,1,1), w=c(0,1,2),   v=c(0,1,2)),
  list(src="F.1 2x3 b1",  tb=matrix(c(4,6,0, 1,3,6),2,byrow=TRUE),        u=c(0,1,6),  gamma=1,  delta=c(0,1),   w=c(0,1),     v=c(0,1,2)),
  list(src="F.1 2x3 b2",  tb=matrix(c(2,4,4, 2,2,6),2,byrow=TRUE),        u=c(3,6,10), gamma=1,  delta=c(0,1),   w=c(0,1),     v=c(0,1,2)),
  list(src="F.2 exp 1.1 (N=86)",  tb=matrix(c(12,18,5, 6,12,6, 6,6,15),3,byrow=TRUE),  u=c(0,10,20), gamma=1,   delta=c(0,1,1), w=c(0,1,2.5), v=c(0,1,2)),
  list(src="F.2 girls (N=94)",    tb=matrix(c(12,3,0, 18,12,3, 17,25,4),3,byrow=TRUE), u=c(0,40,7),  gamma=0.5, delta=c(0,1,1), w=c(0,1,2.5), v=c(0,1,2)),
  list(src="F.2 boys (N=112)",    tb=matrix(c(10,8,1, 29,11,3, 20,24,6),3,byrow=TRUE), u=c(0,20,10), gamma=0.5, delta=c(0,1,1), w=c(0,1,2.5), v=c(0,1,2)),
  list(src="Data girls, prior",   tb=matrix(c(12,3,0, 18,12,3, 17,25,4),3,byrow=TRUE), u=c(0,15,7),  gamma=log(2), delta=c(0,1,1), w=c(0,0.25,1.5), v=c(0,1,1.5)),
  list(src="Data boys, prior",    tb=matrix(c(10,8,1, 29,11,3, 20,24,6),3,byrow=TRUE), u=c(0,20,10), gamma=log(2), delta=c(0,1,1), w=c(0,0.25,1.5), v=c(0,1,1.5)),
  list(src="Data girls, profile", tb=matrix(c(12,3,0, 18,12,3, 17,25,4),3,byrow=TRUE), u=c(0,15,7),  gamma=log(2), delta=c(0,1,1), w=c(0,0.2,1.6),  v=c(0,0.67,0.7)),
  list(src="Data boys, profile",  tb=matrix(c(10,8,1, 29,11,3, 20,24,6),3,byrow=TRUE), u=c(0,20,10), gamma=log(2), delta=c(0,1,1), w=c(0,0.2,1.6),  v=c(0,0.67,0.7)),
  ## a NON-linear (chi-square) statistic through the R-callback path
  list(src="F.2 (N=86) chi-sq",   tb=matrix(c(12,18,5, 6,12,6, 6,6,15),3,byrow=TRUE),  u=c(0,10,20), gamma=1,   delta=c(0,1,1), chisq=TRUE)
)
chisq_stat <- function(tb) { E <- outer(rowSums(tb), colSums(tb)) / sum(tb); sum((tb - E)^2 / E) }

rows <- list()
for (cs in cases) {
  if (isTRUE(cs$chisq)) { stat <- chisq_stat; sc <- NULL } else { sc <- outer(cs$w, cs$v); stat <- function(tb) sum(sc * tb) }
  ## OLD: enumerate + single-u p-value
  o_enum <- tm(env_old$possible.table(threshold = stat(cs$tb), table = cs$tb, direction = "greater than", transform.fun = stat))
  o_p    <- tm(env_old$exact.general.sen.IxJ(u_space = t(as.matrix(cs$u)), obs.table = cs$tb, table_space = o_enum$v,
                                             gamma = cs$gamma, delta = cs$delta))
  ## NEW (same function name): internal enumeration + single-u p-value
  n_p <- tm(env_new$exact.general.sen.IxJ(u_space = t(as.matrix(cs$u)), obs.table = cs$tb, table_space = NULL,
                                          gamma = cs$gamma, delta = cs$delta, transform.fun = stat))
  rows[[length(rows)+1]] <- data.frame(
    case = cs$src, N = sum(cs$tb),
    tables_old = length(o_enum$v), tables_new = n_p$v$n.tables, distinct_TD = n_p$v$n.distinct.tbar,
    p_old = signif(o_p$v$max.prob, 7), p_new = signif(n_p$v$max.prob, 7),
    same = isTRUE(all.equal(o_p$v$max.prob, n_p$v$max.prob, tolerance = 1e-8)) && length(o_enum$v) == n_p$v$n.tables,
    enum_old_s = round(o_enum$s, 2), enum_new_s = round(n_p$v$time[["enumerate.W"]], 3),
    p_old_s = round(o_p$s, 2), p_new_s = round(n_p$v$time[["sweep"]], 3))
  print(rows[[length(rows)]], row.names = FALSE)
}
out <- do.call(rbind, rows); rownames(out) <- NULL
cat("\n"); print(out, row.names = FALSE)
cat("\n", if (all(out$same)) "NEW STACK MATCHES OLD on all cases (table counts and p-values)."
          else "MISMATCH: see 'same' column.", "\n")

## ---- exact.score.sen.IxJ: old vs new on the Table 3 data (girls, profile score, Gamma = 2) --
girls <- matrix(c(12,3,0, 18,12,3, 17,25,4), 3, byrow = TRUE)
o <- tm(env_old$exact.score.sen.IxJ(obs.table = girls, gamma = log(2), delta = c(0,1,1), row = "treatment",
                                    treatment.scores = c(0,0.2,1.6), outcome.scores = c(0,0.67,0.7)))
n <- tm(env_new$exact.score.sen.IxJ(obs.table = girls, gamma = log(2), delta = c(0,1,1), row = "treatment",
                                    treatment.scores = c(0,0.2,1.6), outcome.scores = c(0,0.67,0.7)))
cat(sprintf("\nexact.score.sen.IxJ girls/profile/Gamma=2 (Table 3 says 0.021): old=%.6f (%.1fs)  new=%.6f (%.3fs)  maximizer old=(%s) new=(%s)  %s\n",
            o$v$max.prob, o$s, n$v$max.prob, n$s, paste(o$v$maximizer, collapse=","), paste(n$v$maximizer, collapse=","),
            if (isTRUE(all.equal(o$v$max.prob, n$v$max.prob, tolerance = 1e-8))) "[same]" else "[DIFFER]"))

## ---- possible.table: same set of tables? -----------------------------------------
cs <- cases[[7]]; sc <- outer(cs$w, cs$v); stat <- function(tb) sum(sc * tb)
po <- env_old$possible.table(threshold = stat(cs$tb), table = cs$tb, direction = "greater than", transform.fun = stat)
pn <- env_new$possible.table(threshold = stat(cs$tb), table = cs$tb, direction = "greater than", transform.fun = stat)
key <- function(L) sort(sapply(L, function(t) paste(as.vector(t), collapse = ",")))
cat(sprintf("possible.table on the N=86 table: old %d tables, new %d tables, identical sets: %s\n",
            length(po), length(pn), identical(key(po), key(pn))))

## ---- full U_C on the N=86 table, chi-square, edge (default) vs full ---------------
re <- tm(env_new$exact.general.sen.IxJ(obs.table = cs$tb, gamma = cs$gamma, delta = cs$delta, transform.fun = chisq_stat))
rf <- tm(env_new$exact.general.sen.IxJ(u_space = "full", obs.table = cs$tb, gamma = cs$gamma, delta = cs$delta, transform.fun = chisq_stat))
cat(sprintf("chi-square N=86: default edge set (%d u) %.6f in %.2fs; full U_C (%d u) %.6f in %.1fs  %s\n",
            re$v$n.candidates, re$v$max.prob, re$s, rf$v$n.candidates, rf$v$max.prob, rf$s,
            if (abs(re$v$max.prob - rf$v$max.prob) < 1e-10) "[edge = full]" else "[EDGE SET MISSES THE MAXIMUM]"))

## ---- normal approximation (code unchanged; must be identical) --------------------
o <- env_old$norm.score.sen.IxJ(obs.table = girls, gamma = log(2), delta = c(0,1,1),
                                treatment.scores = c(0,0.2,1.6), outcome.scores = c(0,0.67,0.7))
n <- env_new$norm.score.sen.IxJ(obs.table = girls, gamma = log(2), delta = c(0,1,1),
                                treatment.scores = c(0,0.2,1.6), outcome.scores = c(0,0.67,0.7))
cat(sprintf("norm.score.sen.IxJ girls/profile/Gamma=2: old=%.8f new=%.8f  %s\n", o$max.prob, n$max.prob,
            if (isTRUE(all.equal(o$max.prob, n$max.prob))) "[identical]" else "[DIFFER]"))

## ---- sampling estimator (weights now from the new kernel; same seed -> same tables) -
tb_small <- matrix(c(2,3,0, 0,1,4, 0,1,4), 3, byrow = TRUE)
set.seed(1); o <- env_old$sampling.score.sen.IxJ(obs.table = tb_small, gamma = 1, delta = c(0,1,1), mc.iteration = 2000,
                                                  treatment.scores = c(0,1,2), outcome.scores = c(0,1,2))
set.seed(1); n <- env_new$sampling.score.sen.IxJ(obs.table = tb_small, gamma = 1, delta = c(0,1,1), mc.iteration = 2000,
                                                  treatment.scores = c(0,1,2), outcome.scores = c(0,1,2))
cat(sprintf("sampling.score.sen.IxJ F.1 block 1, 2000 draws, seed 1: old=%.6f new=%.6f  %s\n", o$max.prob, n$max.prob,
            if (isTRUE(all.equal(o$max.prob, n$max.prob, tolerance = 1e-6))) "[same]" else "[DIFFER: check weight scaling]"))

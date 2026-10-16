############################################################################################################
#                                                                                                          #
# sensitivityIxJ.R                                                                                         #
#                                                                                                          #
# Code for "Exact, Nonparametric Sensitivity Analysis for Observational Studies                            #
# of Contingency Tables"                                                                                   #
# Authors: Elaine K. Chiu and Hyunseung Kang                                                               #
#                                                                                                          #
# R interface to the methods in the paper: sensitivity analysis for an I x J                               #
# contingency table under the generic bias model, with an I-level treatment and                            #
# a J-level categorical outcome.  Three approaches are provided.                                           #
#                                                                                                          #
# 1. EXACT METHODS                                                                                         #
#    - Worst-case p-value of any count-based test statistic, and of the ordinal                            #
#      test, with no asymptotic approximation                                                              #
#    - Any I and J; the cost grows with the number of tables sharing the                                   #
#      observed margins                                                                                    #
#    - exact.general.sen.IxJ, exact.score.sen.IxJ                                                          #
#                                                                                                          #
# 2. SAMPLING-BASED METHODS                                                                                #
#    - Monte Carlo estimate of the same worst-case p-value, by sequential                                  #
#      importance sampling of tables with the observed margins                                             #
#    - For tables too large to enumerate                                                                   #
#    - sampling.general.sen.IxJ, sampling.score.sen.IxJ                                                    #
#                                                                                                          #
# 3. NORMAL APPROXIMATION                                                                                  #
#    - Approximate p-value of the ordinal test from the exact mean and variance                            #
#      of the cell counts                                                                                  #
#    - Fastest of the three; available for 2x2, 2x3, 2x4, 2x5, 3x2, 3x3 and                                #
#      4x2 tables, the dimensions for which the moments are implemented                                    #
#    - norm.score.sen.IxJ                                                                                  #
#                                                                                                          #
# Under the generic bias model, an unmeasured confounder u_s in [0, 1] shifts                              #
# the treatment assignment probabilities by exp(gamma * delta_i * u_s), where                              #
# gamma >= 0 measures the strength of the confounding and the binary vector                                #
# delta marks the treatment levels it favours.                                                             #
#                                                                                                          #
# Dependencies, loaded before this file:                                                                   #
# - Rcpp                                                                                                   #
# - library/sensitivityIxJ.cpp   exact computation and moments                                             #
# - library/SIS_methods.cpp      table sampling, for the sampling methods only                             #
#                                                                                                          #
############################################################################################################



require(Rcpp)

## =============================================================================
## FUNCTION  .as_treatment_rows
## PURPOSE   Put the observed table in the orientation used throughout this file:
##           rows are treatment levels i = 1..I, columns are outcome levels
##           j = 1..J, so that the entry (i, j) counts the subjects with
##           treatment i and outcome j.  A table supplied with outcomes in rows
##           (row = "outcome") is transposed.  The storage mode is set to integer
##           because the C++ routines take an integer matrix.
## INPUT     tb   matrix or table object
##           row  "treatment" or "outcome", describing the supplied table
## OUTPUT    integer matrix with treatment levels in rows.
## EXAMPLE
##   tb <- matrix(c(12,18,5, 6,12,6, 6,6,15), 3, byrow = TRUE)   # rows = treatment
##   identical(.as_treatment_rows(tb, "treatment"),
##             .as_treatment_rows(t(tb), "outcome"))             # TRUE
## =============================================================================

.as_treatment_rows <- function(tb, row) {
  tb <- as.matrix(tb); if (row == "outcome") tb <- t(tb)
  storage.mode(tb) <- "integer"; tb
}

## =============================================================================
## FUNCTION  .check_common
## PURPOSE   Check the arguments shared by the exact routines, and stop with a
##           message if any of them is outside what the sensitivity model allows:
##             row     is "treatment" or "outcome"
##             gamma   is a single number and is nonnegative
##             delta   is a vector of 0s and 1s, one entry per treatment level
##             delta   is not constant: if every delta_i is the same, the
##                     confounder shifts every treatment level equally and the
##                     model reduces to random assignment
##             margins are all nonzero: an empty treatment or outcome level
##                     should be dropped before calling
## INPUT     tb     observed table, rows = treatment levels
##           gamma  sensitivity parameter
##           delta  binary vector of length I
##           row    "treatment" or "outcome", as supplied by the caller
## OUTPUT    delta as an integer vector, returned invisibly.
## EXAMPLE
##   tb <- matrix(c(12,18,5, 6,12,6, 6,6,15), 3, byrow = TRUE)
##   d <- .check_common(tb, gamma = 1, delta = c(0,1,1), row = "treatment")
##   d                                                             # 0 1 1
##   try(.check_common(tb, gamma = 1, delta = c(1,3,3), row = "treatment"))
##   ## error: delta must be a binary vector
## =============================================================================


.check_common <- function(tb, gamma, delta, row) {
  if (!(is.character(row) && length(row) == 1 && row %in% c("treatment", "outcome")))
    stop("'row' must be either 'treatment' or 'outcome'.")
  if (!(is.numeric(gamma) && length(gamma) == 1)) stop("gamma must be a scalar.")
  if (gamma < 0) stop("gamma must be nonnegative.")
  if (!is.numeric(delta) || length(delta) != nrow(tb))
    stop(paste0("'delta' must be a numeric vector of length ", nrow(tb), " (one entry per treatment level)."))
  if (!all(delta %in% c(0, 1))) stop("delta must be a binary vector (entries 0 or 1).")
  if (all(delta == delta[1]))
    stop("delta must not be constant (delta_1 = ... = delta_I reduces the model to random assignment).")
  if (any(rowSums(tb) == 0) || any(colSums(tb) == 0))
    stop("'obs.table' has a zero row or column margin; drop the empty level first.")
  invisible(as.integer(delta))
}

## =============================================================================
## FUNCTION  .check_u_space
## PURPOSE   Check a user-supplied matrix of candidate u-bar vectors, and stop
##           with a message if it is not usable:
##             U        is a numeric matrix, one candidate per row
##             columns  one per outcome level, so ncol(U) = J
##             entries  nonnegative, and u-bar_j <= N_.j, since u-bar_j counts
##                      the subjects at outcome level j with u_s = 1 and cannot
##                      exceed the number of subjects at that level
## INPUT     U                numeric matrix of candidates, one per row
##           outcome_margins  N_.1, ..., N_.J
## OUTPUT    U as an integer matrix.
## EXAMPLE
##   m <- c(24, 36, 26)
##   .check_u_space(matrix(c(0, 0, 10,  5, 5, 5), 2, byrow = TRUE), m)
##   try(.check_u_space(matrix(c(0, 0, 30), 1), m))   # error: 30 exceeds N_.3 = 26
## =============================================================================

.check_u_space <- function(U, outcome_margins) {
  if (!is.matrix(U) || !is.numeric(U)) stop("'u_space' must be a numeric matrix with one u_allocation per row.")
  if (ncol(U) != length(outcome_margins))
    stop(paste0("Number of columns in 'u_space' (", ncol(U), ") must equal the number of outcome levels (", length(outcome_margins), ")."))
  if (any(U < 0)) stop("All elements of 'u_space' must be non-negative.")
  if (any(sweep(U, 2, outcome_margins, ">"))) stop("Each u_allocation entry must not exceed the corresponding outcome margin N_.j.")
  storage.mode(U) <- "integer"; U
}

## ---- Candidate sets ---------------------------------------------------------
## Three sets of candidate u-bar vectors are used, depending on the test:
##   full.candidate.set     every u-bar vector; the default for a general
##                          count-based statistic
##   .u_space_ordinal       the N + 1 vectors used by an ordinal test, J >= 3
##   .u_space_sign_score    the two vectors used when the outcome is binary
## =============================================================================
## FUNCTION  full.candidate.set
## PURPOSE   Build the set of all candidate u-bar vectors: every
##           (u-bar_1, ..., u-bar_J) with 0 <= u-bar_j <= N_.j, so the set has
##           prod_j (N_.j + 1) rows.  This is the default candidate set of
##           exact.general.sen.IxJ, and the set in which the worst case of a
##           general count-based statistic lies.
##           The set grows quickly with the margins, so max_rows stops the call
##           rather than attempting to build a matrix too large for memory; for
##           such a table, evaluate the candidates in chunks with
##           full.candidate.set.rows.
## INPUT     outcome_margins  N_.1, ..., N_.J
##           max_rows         largest number of candidates to build at once
## OUTPUT    integer matrix with prod_j (N_.j + 1) rows and columns u1, ..., uJ.
##           Rows run through column 1 fastest, then column 2, and so on.
## EXAMPLE
##   U <- full.candidate.set(c(24, 36, 26))
##   nrow(U)                        # 25 * 37 * 27 = 24,975
##   head(U, 3)                     # (0,0,0), (1,0,0), (2,0,0)
##   U[nrow(U), ]                   # (24, 36, 26)
## =============================================================================


full.candidate.set <- function(outcome_margins, max_rows = 2e7) {
  M <- prod(outcome_margins + 1)
  if (M > max_rows) stop("|U_C| = ", format(M, big.mark = ","), " candidates exceeds max_rows = ", format(max_rows, big.mark = ","),
                         "; evaluate u_space in chunks (see full.candidate.set.rows).")
  U <- as.matrix(expand.grid(lapply(outcome_margins, function(x) 0:x)))
  storage.mode(U) <- "integer"; colnames(U) <- paste0("u", seq_along(outcome_margins)); U
}
## =============================================================================
## FUNCTION  full.candidate.set.rows
## PURPOSE   Return selected rows of the full candidate set without building the
##           whole matrix, so that a set too large to hold in memory can be
##           swept in chunks.
##           The candidates are ordered as in full.candidate.set: column 1 runs
##           through 0, ..., N_.1 fastest, then column 2 advances, and so on.
##           Row k is recovered from the 1-based index k by taking successive
##           remainders, column j using base N_.j + 1.
## INPUT     outcome_margins  N_.1, ..., N_.J
##           k                row indices, between 1 and prod_j (N_.j + 1)
## OUTPUT    integer matrix with length(k) rows and columns u1, ..., uJ.
## EXAMPLE
##   m <- c(24, 36, 26)
##   full.candidate.set.rows(m, 1:5)              # same as the first 5 rows of
##                                                # full.candidate.set(m)
##   full.candidate.set.rows(m, prod(m + 1))      # (24, 36, 26), the last row
##
##   ## sweep a large candidate set in chunks of 10,000
##   total <- prod(m + 1)
##   for (start in seq(1, total, by = 10000)) {
##     U <- full.candidate.set.rows(m, start:min(start + 9999, total))
##     ## evaluate the candidates in U
##   }
## =============================================================================
full.candidate.set.rows <- function(outcome_margins, k) {
  J <- length(outcome_margins); bases <- outcome_margins + 1; idx <- as.numeric(k) - 1
  U <- matrix(0L, length(k), J)
  for (j in seq_len(J)) { U[, j] <- as.integer(idx %% bases[j]); idx <- idx %/% bases[j] }
  colnames(U) <- paste0("u", seq_len(J)); U
}

## =============================================================================
## FUNCTION  .u_space_ordinal
## PURPOSE   Build the candidate set used by an ordinal test when there are
##           three or more outcome levels.  Candidate k, for k = 0, ..., N,
##           gives u_s = 1 to the k subjects with the highest outcome levels, so
##           the set has N + 1 rows.  In u-bar coordinates this fills outcome
##           level J first, then level J - 1, and so on: row k + 1 is the vector
##           whose entries sum to k and whose last entries are at their margins.
##           Far smaller than the full candidate set, which has
##           prod_j (N_.j + 1) rows, and large enough for an ordinal test whose
##           scores and delta are nondecreasing.
## INPUT     outcome_margins  N_.1, ..., N_.J, with the outcome levels in
##                            increasing order
## OUTPUT    integer matrix with N + 1 rows and columns u1, ..., uJ.
## EXAMPLE
##   U <- .u_space_ordinal(c(3, 2, 2))      # N = 7
##   nrow(U)                                # 8
##   U[c(1, 3, 5, 8), ]                     # (0,0,0), (0,0,2), (0,2,2), (3,2,2)
## =============================================================================

.u_space_ordinal <- function(outcome_margins) {
  m <- outcome_margins; J <- length(m); N <- sum(m)
  U <- t(sapply(0:N, function(k) { out <- integer(J); left <- k
    for (j in J:1) { out[j] <- min(m[j], left); left <- left - out[j] }; out }))
  storage.mode(U) <- "integer"; colnames(U) <- paste0("u", seq_len(J)); U
}

## =============================================================================
## FUNCTION  .u_space_sign_score
## PURPOSE   Candidate matrix used by an ordinal test when the outcome is binary
##           (J = 2).  No search is needed here: the worst case is known in
##           closed form and gives u_s = 1 to exactly the subjects at outcome
##           level 2, that is u-bar = (0, N_.2).  The matrix also carries the
##           reference point (0, 0), at which no subject has u_s = 1, so that a
##           sweep over these two rows returns both the worst-case p-value and
##           the p-value under random assignment.
## INPUT     outcome_margins  c(N_.1, N_.2)
## OUTPUT    2 x 2 integer matrix: row 1 is (0, 0), row 2 is (0, N_.2).
## EXAMPLE
##   .u_space_sign_score(c(15, 75))     # rows (0, 0) and (0, 75)
## =============================================================================


.u_space_sign_score <- function(outcome_margins) {
  U <- matrix(c(0L, 0L, 0L, as.integer(outcome_margins[2])), ncol = 2, byrow = TRUE)
  colnames(U) <- c("u1", "u2"); U
}

## ---- The two computational stages -------------------------------------------
## Stage 1 summarizes the reference set Omega_t into W(t-bar), the part of the
##   p-value that depends on the observed table and the test statistic but not
##   on u or gamma.  Done in C++ by build_W_cpp or build_W_from_tables, once.
## Stage 2 evaluates the p-value at each candidate u-bar from that summary.
##   Done in C++ by alpha_bd, once per candidate.
## .prep_from_W packages the stage-1 output with the constants stage 2 needs;
## .sweep runs stage 2 over a set of candidates; .rct runs it at u-bar = 0.
## =============================================================================
## FUNCTION  .prep_from_W
## PURPOSE   Collect everything stage 2 needs into one list, so that .sweep and
##           .rct can be called without assembling the arguments of alpha_bd
##           each time.  Three fields are copied from the stage-1 result, two
##           from the arguments, and three are derived here.
## INPUT     res  list returned by build_W_cpp or build_W_from_tables
##           n    treatment margins N_i.
##           m    outcome margins N_.j
##           d    delta, binary, length I
##           J    number of outcome levels
## OUTPUT    list with
##             TDmat     distinct t-bar vectors, from res
##             logW      log W(t-bar), from res
##             n_tables  number of tables in Omega_t, from res
##             m, J      as supplied
##             N         total number of subjects, sum of n
##             ND        number of subjects in the treatment levels with
##                       delta_i = 1
##             logc0     log of the factor of the denominator that does not
##                       depend on u, from log_c_0
## EXAMPLE
##   tb   <- matrix(c(12,18,5, 6,12,6, 6,6,15), 3, byrow = TRUE)
##   sc   <- outer(0:2, 0:2)
##   res  <- build_W_cpp(rowSums(tb), colSums(tb), c(0,1,1), sum(sc * tb),
##                       TRUE, sc, NULL)
##   prep <- .prep_from_W(res, rowSums(tb), colSums(tb), c(0,1,1), 3)
##   c(prep$N, prep$ND)                     # 86, 51
## =============================================================================


.prep_from_W <- function(res, n, m, d, J)
  list(TDmat = res$TDmat, logW = res$logW, m = m, N = sum(n), ND = sum(n[d == 1L]),
       logc0 = log_c_0(n, d), J = J, n_tables = res$n_tables)


## =============================================================================
## FUNCTION  .sweep
## PURPOSE   Evaluate alpha(T, r, u) at every candidate u-bar in U (rows) for one
##           gamma, by calling the C++ function alpha_bd, and record the elapsed
##           time. Each candidate costs O(N^{J+1}); W is reused across all
##           candidates and all gamma.
## INPUT     prep  from .prep_from_W;  U  integer candidate matrix;  gamma >= 0.
## OUTPUT    list(all.prob = alpha at each row of U, seconds = wall time).
## EXAMPLE
##   sw <- .sweep(prep, .u_space_ordinal(colSums(tb)), gamma = 1)
##   max(sw$all.prob)                       # worst-case p-value over U_O
## =============================================================================
.sweep <- function(prep, U, gamma) {
  t0 <- proc.time()[["elapsed"]]
  a <- alpha_bd(prep$TDmat, prep$logW, U, prep$m, gamma, prep$N, prep$ND, prep$logc0)
  list(all.prob = a, seconds = proc.time()[["elapsed"]] - t0)
}

## =============================================================================
## FUNCTION  .rct
## PURPOSE   alpha(T, r, u) at u = 0 (all u-bar_j = 0), i.e. the p-value under
##           random assignment. At u = 0 the weights exp(gamma delta' q) are all 1,
##           so this equals the multivariate hypergeometric p-value for any gamma.
## INPUT     prep  from .prep_from_W;  gamma (does not affect the value).
## OUTPUT    scalar p-value.
## EXAMPLE
##   .rct(prep, gamma = 1)                  # equals .rct(prep, gamma = 0)
## =============================================================================

.rct <- function(prep, gamma)
  alpha_bd(prep$TDmat, prep$logW, matrix(0L, 1, prep$J), prep$m, gamma, prep$N, prep$ND, prep$logc0)[1]

## =============================================================================
## FUNCTION  possible.table
## PURPOSE   Enumerate Omega_t: every I x J table with the row and column margins
##           of `table` whose statistic transform.fun(t) is >= threshold
##           ("greater than") or <= threshold ("less than"). Enumeration is done
##           in C++ (build_W_cpp with keep_tables = TRUE); transform.fun is an R
##           callback evaluated once per table with the given margins. For a
##           large Omega_t, prefer exact.general.sen.IxJ with transform.fun or
##           scores, which never stores the tables.
## INPUT     threshold  critical value c;  table  observed table (rows = treatment);
##           direction;  transform.fun  function(table) -> numeric statistic.
## OUTPUT    list of tables (class "table"); list() if Omega_t is empty.
## EXAMPLE
##   tb   <- matrix(c(3,2,1, 0,2,4, 0,1,5), 3, byrow = TRUE)
##   stat <- function(t) sum(outer(0:2, 0:2) * t)
##   Omega <- possible.table(threshold = stat(tb), table = tb,
##                           direction = "greater than", transform.fun = stat)
##   length(Omega)
##   all(sapply(Omega, function(t) stat(t) >= stat(tb)))
## =============================================================================

possible.table <- function(threshold, table, direction = c("greater than", "less than"), transform.fun) {
  direction <- match.arg(direction)
  tb <- as.matrix(table); storage.mode(tb) <- "integer"
  I <- nrow(tb); J <- ncol(tb)
  res <- build_W_cpp(as.integer(rowSums(tb)), as.integer(colSums(tb)), rep(1L, I), threshold,
                     greater = (direction == "greater than"), score = NULL,
                     stat_fun = function(t) transform.fun(t), keep_tables = TRUE)
  if (res$n_tables == 0) return(list())
  lapply(seq_len(nrow(res$tables)), function(k) as.table(matrix(res$tables[k, ], nrow = I, byrow = TRUE)))
}



## =============================================================================
## FUNCTION  is_binary_vector                                    [internal helper]
## PURPOSE   TRUE if every entry of x is 0 or 1.  Used to check delta in the
##           sampling and normal-approximation routines (model (6) needs delta binary).
## INPUT     x  numeric vector.
## OUTPUT    logical scalar.
## EXAMPLE
##   is_binary_vector(c(0, 1, 1))   # TRUE
##   is_binary_vector(c(0, 2, 1))   # FALSE
## =============================================================================
is_binary_vector <- function(x) {
  all(x %in% c(0, 1))
}


## =============================================================================
## FUNCTION  is_scalar_numeric                                   [internal helper]
## PURPOSE   TRUE if x is a single numeric value.  Used to check gamma.
## INPUT     x.
## OUTPUT    logical scalar.
## EXAMPLE
##   is_scalar_numeric(0.5)       # TRUE
##   is_scalar_numeric(c(0, 1))   # FALSE
## =============================================================================
is_scalar_numeric <- function(x) {
  is.numeric(x) && length(x) == 1
}



# ------------------------------------------------------------------------------------------------------------------

################################################################################
#                                                                              
#                        EXACT COMPUTATION METHODS                             
#                                                                              
# This section implements exact methods for sensitivity analysis in I×J         
# contingency tables under the generic bias model. The methods compute exact   
# p-values by enumerating and evaluating all possible tables in the reference  
# set.                                          
#                                                                              
# Main Functions:                                                              
# - generic.I.by.J.sensitivity.point.probability()                           
#    - Core function that computes exact probability for a single table       
#    - Routes to dimension-specific C++ implementations for efficiency        
#                                                                                
# - exact.general.sen.IxJ()                                                  
#    - Main interface for general test statistics                             
#    - Computes exact p-values by summing over all relevant tables            
#                                                                                
# -  exact.score.sen.IxJ()                                                    
#    - Specialized for ordinal tests with monotone trends                       
#             
#################################################################################






## =============================================================================
## FUNCTION  generic.I.by.J.sensitivity.point.probability
## PURPOSE   Exact null probability P_u(N = table) of ONE table under Fisher's sharp
##           null and the generic bias model at (gamma, delta), for one u given in
##           u-bar coordinates (u_allocation[j] = number of subjects at outcome level
##           j with u_s = 1). The probability equals kappa(t) S(t-bar, u) divided by
##           kappa_0 S_0(u); at gamma = 0, or when every u_allocation entry is 0 or
##           equal to its outcome margin, this is the multivariate hypergeometric
##           probability of the table. Validates the inputs and calls the C++
##           function point_prob_bd.
## INPUT     table, row, u_allocation (length J), gamma >= 0, delta (binary, length I);
##           `...` absorbs the shared_divisor argument of a previous release, no
##           longer needed since the computation is done in log space.
## OUTPUT    scalar probability.
## EXAMPLE
##   tb <- matrix(c(12,18,5, 6,12,6, 6,6,15), 3, byrow = TRUE)
##   p0 <- generic.I.by.J.sensitivity.point.probability(tb, u_allocation = c(0,0,0), gamma = 0, delta = c(0,1,1))
##   p1 <- generic.I.by.J.sensitivity.point.probability(tb, u_allocation = c(0,10,20), gamma = 1, delta = c(0,1,1))
## =============================================================================

generic.I.by.J.sensitivity.point.probability <- function(table, row = "treatment", u_allocation,
                                                        gamma, delta, ...) {
  if (!is.matrix(table) && !inherits(table, "table")) stop("'table' must be either a matrix or a table object.")
  if (length(dim(table)) != 2) stop("'table' must be two-dimensional.")
  tb <- .as_treatment_rows(table, row)
  d  <- .check_common(tb, gamma, delta, row)
  if (!is.numeric(u_allocation) || length(u_allocation) != ncol(tb))
    stop(paste0("'u_allocation' must be a numeric vector of length ", ncol(tb), " (one entry per outcome level)."))
  if (any(u_allocation < 0)) stop("All elements of 'u_allocation' must be non-negative.")
  if (any(u_allocation > colSums(tb))) stop("Each u_allocation entry must not exceed the corresponding outcome margin N_.j.")
  point_prob_bd(tb, d, as.integer(u_allocation), gamma)
}




#' Exact worst-case p-value for a general test statistic in an I x J table
#'
#' Computes the exact worst-case p-value, over all values of an unmeasured
#' confounder, for any count-based test statistic in an I x J contingency
#' table under the generic bias sensitivity model. The computation has two
#' stages: first, the reference set of tables at least as extreme as the
#' observed table is summarized once, independently of the confounder and of
#' gamma; second, the p-value is evaluated at every candidate confounder
#' allocation, and the maximum is returned.
#'
#' @param obs.table A matrix or table object giving the observed contingency
#'   table.
#' @param gamma A nonnegative scalar, the sensitivity parameter.
#' @param delta A binary vector marking the treatment levels favoured by the
#'   confounder. Its length must equal the number of treatment levels (rows
#'   of \code{obs.table} if \code{row = "treatment"}, or columns if
#'   \code{row = "outcome"}), and it must not be constant.
#' @param row Either \code{"treatment"} or \code{"outcome"}, indicating which
#'   dimension of \code{obs.table} holds the treatment levels. Default is
#'   \code{"treatment"}.
#' @param u_space Candidate confounder allocations to evaluate: either
#'   \code{NULL} or \code{"full"} for the full candidate set of all u-bar
#'   vectors, which is guaranteed to contain the worst case for any
#'   count-based statistic; or a matrix with one allocation per row (one
#'   column per outcome level), to evaluate only a chosen subset. Default is
#'   \code{NULL}.
#' @param table_space A list of tables with the margins of \code{obs.table}
#'   whose test statistic is at least as extreme as observed, as returned by
#'   \code{\link{possible.table}}. If supplied, this defines the reference
#'   set directly.
#' @param transform.fun A function of one table returning a numeric test
#'   statistic, used to build the reference set internally when
#'   \code{table_space}, \code{treatment.scores}/\code{outcome.scores} and
#'   \code{statistic} are all \code{NULL}. The statistic is evaluated once
#'   per candidate table and no table is stored.
#' @param threshold The critical value defining the reference set, on the
#'   scale of \code{transform.fun} or \code{statistic}. Defaults to the
#'   observed statistic.
#' @param treatment.scores,outcome.scores Numeric score vectors for the
#'   treatment and outcome levels. If both are supplied (and
#'   \code{table_space} is \code{NULL}), the reference set is built around
#'   the ordinal statistic \code{sum_ij treatment.scores_i * outcome.scores_j
#'   * obs.table_ij}, evaluated in C++.
#' @param statistic One of \code{"chisq"}, \code{"G2"} or
#'   \code{"freeman_halton"}, selecting a built-in count-based statistic
#'   evaluated in C++ (the chi-square statistic, the likelihood-ratio
#'   statistic, or the r x c extension of Fisher's exact test). Used only
#'   when \code{table_space} and the score arguments are all \code{NULL}.
#' @param verbose Logical; if \code{TRUE}, prints the p-value at every
#'   candidate allocation. Default is \code{FALSE}.
#'
#' @return A list containing:
#' \describe{
#'   \item{rct.prob}{The p-value under random assignment (no confounding).}
#'   \item{max.prob}{The worst-case p-value, maximized over \code{u_space}.}
#'   \item{maximizer}{The allocation achieving \code{max.prob}.}
#'   \item{gamma}{The supplied \code{gamma}.}
#'   \item{delta}{The supplied \code{delta}.}
#'   \item{n.tables}{The size of the reference set.}
#'   \item{n.distinct.tbar}{The number of distinct treatment-margin
#'     summaries in the reference set.}
#'   \item{n.candidates}{The number of allocations evaluated.}
#'   \item{all.prob}{The p-value at every candidate, in the order of
#'     \code{u_space}.}
#'   \item{time}{Elapsed seconds for building the reference set and for
#'     evaluating the candidates, as \code{c(enumerate.W, sweep)}.}
#' }
#'
#' @details
#' The reference set is taken from the first of the following that is
#' supplied: \code{table_space}; \code{treatment.scores} together with
#' \code{outcome.scores}; \code{statistic}; \code{transform.fun}. At least
#' one of these must be given.
#'
#' Unlike \code{\link{exact.score.sen.IxJ}}, this function places no
#' monotonicity requirement on \code{delta} or on any scores, at the cost of
#' evaluating the p-value at every allocation in the full candidate set
#' rather than a small subset.
#'
#' @examples
#' ## chi-square statistic, built in, worst case over the full candidate set
#' obs.table <- matrix(c(12, 18, 5, 6, 12, 6, 6, 6, 15), nrow = 3, byrow = TRUE)
#' res1 <- exact.general.sen.IxJ(obs.table = obs.table, gamma = 1,
#'                                delta = c(0, 1, 1), statistic = "chisq")
#' c(res1$rct.prob, res1$max.prob)
#'
#' ## a user-defined statistic, via possible.table
#' transform.fun <- function(tb) 4 * tb[3, 3] + 3 * tb[2, 3]
#' table.set <- possible.table(threshold = transform.fun(obs.table), table = obs.table,
#'                              direction = "greater than", transform.fun = transform.fun)
#' res2 <- exact.general.sen.IxJ(obs.table = obs.table, table_space = table.set,
#'                                gamma = 0.5, delta = c(0, 1, 1))
#'
#' ## evaluating a single candidate allocation
#' res3 <- exact.general.sen.IxJ(u_space = matrix(c(0, 10, 20), nrow = 1),
#'                                obs.table = obs.table, gamma = 1, delta = c(0, 1, 1),
#'                                treatment.scores = c(0, 1, 2.5), outcome.scores = c(0, 1, 2))
#'
#' @seealso
#' \code{\link{exact.score.sen.IxJ}} for ordinal tests with monotone scores,
#' \code{\link{possible.table}} for building a reference set
#'
#' @export
exact.general.sen.IxJ <- function(u_space = NULL, obs.table, table_space = NULL, gamma, delta,
                                  row = "treatment", verbose = FALSE,
                                  transform.fun = NULL, threshold = NULL,
                                  treatment.scores = NULL, outcome.scores = NULL,
                                  statistic = NULL) {
  ## ---- input checks (model (6): gamma >= 0, delta binary; table orientation) ----
  if (!is.matrix(obs.table) && !inherits(obs.table, "table")) stop("'obs.table' must be either a matrix or a table object.")
  if (length(dim(obs.table)) != 2) stop("'obs.table' must be two-dimensional.")
  tb <- .as_treatment_rows(obs.table, row)
  d  <- .check_common(tb, gamma, delta, row)
  I <- nrow(tb); J <- ncol(tb)
  n <- as.integer(rowSums(tb)); m <- as.integer(colSums(tb))     # N_i. and N_.j

  t0 <- proc.time()[["elapsed"]]
  if (!is.null(table_space)) {
    if (length(table_space) == 0) stop("'table_space' is empty.")
    tabs <- t(sapply(table_space, function(x) {
      x <- .as_treatment_rows(x, row)
      if (!all(rowSums(x) == n) || !all(colSums(x) == m))
        stop("some table in the table_space does not have the same margins as the observed table")
      as.vector(t(x))
    }))
    storage.mode(tabs) <- "integer"
    res <- build_W_from_tables(tabs, I, J, d)
  } else if (!is.null(treatment.scores) && !is.null(outcome.scores)) {
    if (length(treatment.scores) != I || length(outcome.scores) != J)
      stop("'treatment.scores' / 'outcome.scores' must have length I / J.")
    score <- outer(treatment.scores, outcome.scores)
    if (is.null(threshold)) threshold <- sum(score * tb)
    res <- build_W_cpp(n, m, d, threshold, greater = TRUE, score = score, stat_fun = NULL)
  } else if (!is.null(statistic)) {
    statistic <- match.arg(statistic, c("chisq", "G2", "freeman_halton"))
    ## reduced monotone forms (margins fixed); threshold = observed value of the reduced form
    E <- outer(n, m) / sum(n)
    code <- switch(statistic, chisq = 2L, G2 = 3L, freeman_halton = 4L)
    threshold <- switch(statistic,
      chisq  = sum(tb^2 / E),                                   # chi-square + N
      G2     = sum(ifelse(tb > 0, tb * log(pmax(tb, 1)), 0)),   # G^2 / 2 + constant
      freeman_halton = sum(lgamma(tb + 1)))                     # -log P(t) + constant
    res <- build_W_cpp(n, m, d, threshold, greater = TRUE, score = NULL, stat_fun = NULL, statistic = code)
  } else {
    if (is.null(transform.fun)) stop("supply 'table_space', 'transform.fun', treatment/outcome scores, or a built-in statistic.")
    stat <- if (row == "outcome") function(x) transform.fun(t(x)) else transform.fun
    if (is.null(threshold)) threshold <- stat(tb)
    res <- build_W_cpp(n, m, d, threshold, greater = TRUE, score = NULL, stat_fun = stat)
  }
  t_enum <- proc.time()[["elapsed"]] - t0
  if (res$n_tables == 0) {
    warning("No tables found in the reference set. Check your test statistic and table.")
    return(list(rct.prob = NA, max.prob = NA, maximizer = NA, gamma = gamma, delta = delta))
  }
  prep <- .prep_from_W(res, n, m, d, J)

  U <- if (is.null(u_space) || identical(u_space, "full")) full.candidate.set(m)
       else .check_u_space(u_space, m)

  sw <- .sweep(prep, U, gamma)
  k <- which.max(sw$all.prob)
  if (verbose) for (i in seq_len(nrow(U)))
    cat("u_allocation", paste(U[i, ], collapse = ","), " prob", sw$all.prob[i], "\n")
  list(rct.prob = .rct(prep, gamma), max.prob = sw$all.prob[k], maximizer = as.numeric(U[k, ]),
       gamma = gamma, delta = delta,
       n.tables = res$n_tables, n.distinct.tbar = nrow(res$TDmat), n.candidates = nrow(U),
       all.prob = sw$all.prob, time = c(enumerate.W = t_enum, sweep = sw$seconds))
}




#' Exact worst-case p-value for an ordinal test in an I x J table
#'
#' Computes the exact worst-case p-value, over all values of an unmeasured
#' confounder, for the ordinal test statistic T = sum_ij w_i v_j N_ij, where w
#' and v are scores assigned to the treatment and outcome levels. Both
#' treatment and outcome are assumed ordinal, with scores and the confounder
#' indicator increasing together. Under that assumption the worst case is
#' known to lie in a small set of candidate confounder allocations -- just two
#' when the outcome is binary, and N + 1 otherwise -- so this function is far
#' faster than \code{\link{exact.general.sen.IxJ}} for ordinal data, without
#' approximation.
#'
#' @param obs.table A matrix or table object giving the observed contingency
#'   table.
#' @param gamma A nonnegative scalar, the sensitivity parameter.
#' @param delta A binary, nondecreasing vector marking the treatment levels
#'   favoured by the confounder. Its length must equal the number of
#'   treatment levels (rows of \code{obs.table} if \code{row = "treatment"},
#'   or columns if \code{row = "outcome"}).
#' @param row Either \code{"treatment"} or \code{"outcome"}, indicating which
#'   dimension of \code{obs.table} holds the treatment levels. Default is
#'   \code{"treatment"}.
#' @param treatment.scores A nondecreasing numeric vector of scores for the
#'   treatment levels, one per row of the oriented table.
#' @param outcome.scores A nondecreasing numeric vector of scores for the
#'   outcome levels, one per column of the oriented table.
#' @param verbose Logical; if \code{TRUE}, prints the p-value at every
#'   candidate allocation. Default is \code{FALSE}.
#'
#' @return A list containing:
#' \describe{
#'   \item{rct.prob}{The p-value under random assignment (no confounding).}
#'   \item{max.prob}{The worst-case p-value, maximized over the candidate
#'     allocations.}
#'   \item{maximizer}{The allocation achieving \code{max.prob}.}
#'   \item{gamma}{The supplied \code{gamma}.}
#'   \item{delta}{The supplied \code{delta}.}
#'   \item{obs.stat}{The observed value of T.}
#'   \item{obs.table}{The observed table, as supplied.}
#'   \item{n.tables}{The size of the reference set of tables at least as
#'     extreme as observed.}
#'   \item{n.distinct.tbar}{The number of distinct treatment-margin summaries
#'     in the reference set.}
#'   \item{n.candidates}{The number of allocations evaluated (2 when the
#'     outcome is binary, N + 1 otherwise).}
#'   \item{all.prob}{The p-value at every candidate allocation.}
#'   \item{time}{Elapsed seconds for building the reference set and for
#'     evaluating the candidates, as \code{c(enumerate.W, sweep)}.}
#' }
#'
#' @details
#' This function stops with an error if \code{treatment.scores},
#' \code{outcome.scores} or \code{delta} is not nondecreasing, since the
#' restriction to a small candidate set is not valid otherwise; in that case
#' use \code{\link{exact.general.sen.IxJ}}, which places no monotonicity
#' requirement on \code{delta} or on any scores but searches a larger
#' candidate set.
#'
#' @examples
#' ## 2x2 table: the outcome is binary, so only two allocations are evaluated
#' obs.table <- matrix(c(12, 18, 17, 3), nrow = 2, byrow = TRUE)
#' result_2x2 <- exact.score.sen.IxJ(obs.table = obs.table, gamma = 0.5,
#'                                    delta = c(0, 1),
#'                                    treatment.scores = c(0, 1),
#'                                    outcome.scores = c(0, 1))
#' c(result_2x2$rct.prob, result_2x2$max.prob)
#'
#' ## 3x3 table with an ordinal outcome
#' obs.table <- matrix(c(12, 18, 17, 3, 12, 25, 0, 3, 4), nrow = 3, byrow = TRUE)
#' result_3x3 <- exact.score.sen.IxJ(obs.table = obs.table, gamma = 0.5,
#'                                    delta = c(0, 1, 1),
#'                                    treatment.scores = c(0, 1, 2),
#'                                    outcome.scores = c(1, 2, 3))
#' c(result_3x3$rct.prob, result_3x3$max.prob); result_3x3$maximizer
#'
#' @seealso
#' \code{\link{exact.general.sen.IxJ}} for test statistics or treatment
#' effects that are not monotone,
#' \code{\link{possible.table}} for building a reference set
#'
#' @export
exact.score.sen.IxJ <- function(obs.table, gamma, delta, row = "treatment", verbose = FALSE,
                                treatment.scores, outcome.scores) {
  ## ---- input checks (model (6)) ---------------------------------------------------
  if (!is.matrix(obs.table) && !inherits(obs.table, "table")) stop("'obs.table' must be either a matrix or a table object.")
  if (length(dim(obs.table)) != 2) stop("'obs.table' must be two-dimensional.")
  tb <- .as_treatment_rows(obs.table, row)
  d  <- .check_common(tb, gamma, delta, row)
  I <- nrow(tb); J <- ncol(tb)
  if (!is.numeric(treatment.scores) || length(treatment.scores) != I)
    stop(paste0("'treatment.scores' must be a numeric vector of length ", I, "."))
  if (!is.numeric(outcome.scores) || length(outcome.scores) != J)
    stop(paste0("'outcome.scores' must be a numeric vector of length ", J, "."))
  ## ---- assumptions of Theorems 2-3: simultaneous monotonicity ---------------------
  if (!all(diff(treatment.scores) >= 0)) stop("'treatment.scores' must have a monotone trend.")
  if (!all(diff(outcome.scores)   >= 0)) stop("'outcome.scores' must have a monotone trend.")
  if (!all(diff(delta) >= 0))            stop("'delta' must have a monotone trend.")
  n <- as.integer(rowSums(tb)); m <- as.integer(colSums(tb))
  score <- outer(treatment.scores, outcome.scores)
  obs.stat <- sum(score * tb)

  t0 <- proc.time()[["elapsed"]]
  res <- build_W_cpp(n, m, d, obs.stat, greater = TRUE, score = score, stat_fun = NULL)
  t_enum <- proc.time()[["elapsed"]] - t0
  if (res$n_tables == 0) {
    warning("No tables found in the reference set. Check your test statistic and table.")
    return(list(rct.prob = NA, max.prob = NA, maximizer = NA, gamma = gamma, delta = delta,
                obs.stat = obs.stat, obs.table = obs.table))
  }
  prep <- .prep_from_W(res, n, m, d, J)

  ## ---- candidate set: Theorem 3 (J = 2) or Theorem 2 (J >= 3) ----------------------
  U <- if (J == 2) .u_space_sign_score(m) else .u_space_ordinal(m)

  sw <- .sweep(prep, U, gamma)
  k <- which.max(sw$all.prob)
  if (verbose) for (i in seq_len(nrow(U)))
    cat("u_allocation", paste(U[i, ], collapse = ","), " prob", sw$all.prob[i], "\n")
  list(rct.prob = .rct(prep, gamma), max.prob = sw$all.prob[k], maximizer = as.numeric(U[k, ]),
       gamma = gamma, delta = delta, obs.stat = obs.stat, obs.table = obs.table,
       n.tables = res$n_tables, n.distinct.tbar = nrow(res$TDmat), n.candidates = nrow(U),
       all.prob = sw$all.prob, time = c(enumerate.W = t_enum, sweep = sw$seconds))
}




################################################################################
#
#                      SAMPLING-BASED COMPUTATION METHODS
#
# This section implements sampling-based methods for sensitivity analysis in
# I×J contingency tables under the generic bias model, for tables too large
# for exact enumeration. Tables with the observed margins are drawn by
# sequential importance sampling, and the worst-case p-value is estimated by
# self-normalized importance sampling.
#
# Main Functions:
# - compute_sampling_target_v()
#    - The importance-sampling target for one sampled table: its exact null
#      probability under the generic bias model, for any I and J
#
# - single.u.sampling.general.IxJ()
#    - Monte Carlo estimate of the p-value at a single u-allocation, with the
#      running estimate after each draw
#
# - sampling.general.sen.IxJ()
#    - Main interface for a general test statistic: searches over a u-space
#      and returns the largest Monte Carlo estimate
#
# - sampling.score.sen.IxJ()
#    - Specialized for the ordinal test statistic, searching only the small
#      candidate set valid under monotone scores and delta
#
# Supporting Functions (C++ Backend):
# - sample_sis_one_table()   : draws one table with the given margins by
#                              sequential importance sampling (SIS_methods.cpp)
#
################################################################################

## =============================================================================
## FUNCTION  compute_sampling_target_v
## PURPOSE   Importance-sampling target for one sampled table: the exact null
##           probability P_u(N = mc.table) under the generic bias model at
##           (gamma, delta, u), computed by point_prob_bd. Any positive multiple
##           of this quantity would serve equally well, since the SIS estimator
##           is a ratio of weighted sums and any shared constant cancels; this
##           implementation uses the normalized probability directly, and the
##           computation does not depend on I or J.
## INPUT     I, J, mc.table (rows = treatment), treatment_margins, outcome_margins,
##           u_allocation (u-bar), N, gamma_delta = gamma * delta.
## OUTPUT    scalar target value.
## EXAMPLE
##   tb <- matrix(c(12,18,5, 6,12,6, 6,6,15), 3, byrow = TRUE)
##   compute_sampling_target_v(3, 3, tb, rowSums(tb), colSums(tb), u_allocation = c(0,10,20),
##                             N = 86, gamma_delta = 1 * c(0,1,1), shared_divisor = 1)
##   ## equals generic.I.by.J.sensitivity.point.probability(tb, u_allocation = c(0,10,20), gamma = 1, delta = c(0,1,1))
## =============================================================================

compute_sampling_target_v <- function(I, J, mc.table, treatment_margins,
                                      outcome_margins, u_allocation,
                                      N, gamma_delta, shared_divisor) {
  ## The caller passes gamma_delta = gamma * delta with delta binary and gamma >= 0
  ## (validated by the sampling functions), so delta = (gamma_delta != 0) and
  ## gamma = max(gamma_delta); at gamma = 0 both are zero, giving uniform weights.
  d <- as.integer(gamma_delta != 0)
  g <- if (any(gamma_delta != 0)) max(gamma_delta) else 0
  tb <- as.matrix(mc.table); storage.mode(tb) <- "integer"
  point_prob_bd(tb, d, as.integer(u_allocation), g)
}




## =============================================================================
## FUNCTION  single.u.sampling.general.IxJ
## PURPOSE   Self-normalized importance-sampling estimate of alpha(T, r, u) at ONE
##           u-bar. Tables with the observed margins are drawn by sequential
##           importance sampling (sample_sis_one_table from SIS_methods.cpp,
##           proposal h); each draw t gets weight P_u(t) / h(t), and the estimate is
##             sum_t 1{T(t) >= T_obs} P_u(t)/h(t)  /  sum_t P_u(t)/h(t).
## INPUT     table (matrix), row, mc.iteration, u_allocation (u-bar), gamma, delta,
##           transform.fun. Every caller in this file supplies table already
##           oriented with treatment levels in rows and passes row = "treatment";
##           with row = "outcome", the observed statistic is computed before the
##           table is transposed, so transform.fun must be invariant to the
##           table's orientation for the result to be correct.
## OUTPUT    list(probability.progression = running estimate after each draw,
##                final.prob.estimate, total.iteration).
## EXAMPLE
##   tb   <- matrix(c(12,18,5, 6,12,6, 6,6,15), 3, byrow = TRUE)
##   stat <- function(t) sum(outer(c(0,1,2.5), c(0,1,2)) * t)
##   s <- single.u.sampling.general.IxJ(tb, mc.iteration = 2000, u_allocation = c(0,10,20),
##                                      gamma = 1, delta = c(0,1,1), transform.fun = stat)
##   s$final.prob.estimate
## =============================================================================

single.u.sampling.general.IxJ = function(table, row = "treatment", mc.iteration = 20000,
                                         u_allocation, gamma, delta, transform.fun) {
  obs.stat = transform.fun(table)
  
  if (row == "outcome") {
    table = t(table)
  }
  
  N = sum(table)
  I = dim(table)[1]
  J = dim(table)[2]
  
  treatment_margins = rowSums(table)
  outcome_margins   = colSums(table)
  
  # Precompute shared values
  shared_divisor = .Machine$double.xmax^0.9
  
  # Initialize storage
  numerator.summation = 0
  denominator.summation = 0
  mc.probability = numeric(mc.iteration)
  
  for (mc in seq_len(mc.iteration)) {
    mc.table.result = sample_sis_one_table(row_s = treatment_margins,
                                           col_s = outcome_margins)
    inverse.h = mc.table.result$weight
    mc.table  = mc.table.result$table
    
    target.v = compute_sampling_target_v(I = I, J = J,
                                         mc.table = mc.table,
                                         treatment_margins = treatment_margins,
                                         outcome_margins = outcome_margins,
                                         u_allocation = u_allocation,
                                         N = N,
                                         gamma_delta = gamma*delta,
                                         shared_divisor = shared_divisor)
    
    weight = inverse.h * target.v
    indicator = as.integer(transform.fun(mc.table) >= obs.stat)
    
    numerator.summation   = numerator.summation + indicator * weight
    denominator.summation = denominator.summation + weight
    mc.probability[mc]    = numerator.summation / denominator.summation
  }
  
  return(list(
    probability.progression = mc.probability,
    final.prob.estimate = numerator.summation / denominator.summation,
    total.iteration = mc.iteration
  ))
}


#' Monte Carlo worst-case p-value for an ordinal test in an I x J table
#'
#' Monte Carlo counterpart of \code{\link{exact.score.sen.IxJ}} for the ordinal
#' test statistic T = sum_ij w_i v_j N_ij, for use when the reference set of
#' tables at least as extreme as observed is too large to enumerate exactly.
#' Tables with the observed margins are drawn by sequential importance
#' sampling, and the p-value at each candidate confounder allocation is
#' estimated by self-normalized importance sampling.
#'
#' @param obs.table A matrix or table object giving the observed contingency
#'   table.
#' @param gamma A nonnegative scalar, the sensitivity parameter.
#' @param delta A binary, nondecreasing vector marking the treatment levels
#'   favoured by the confounder. Its length must equal the number of
#'   treatment levels (rows of \code{obs.table} if \code{row = "treatment"},
#'   or columns if \code{row = "outcome"}).
#' @param row Either \code{"treatment"} or \code{"outcome"}, indicating which
#'   dimension of \code{obs.table} holds the treatment levels. Default is
#'   \code{"treatment"}.
#' @param mc.iteration The number of Monte Carlo draws used to estimate the
#'   p-value at each candidate allocation. Higher values increase accuracy at
#'   the cost of computation time. Default is 20000.
#' @param treatment.scores A nondecreasing numeric vector of scores for the
#'   treatment levels, one per row of the oriented table.
#' @param outcome.scores A nondecreasing numeric vector of scores for the
#'   outcome levels, one per column of the oriented table.
#' @param verbose Logical; if \code{TRUE}, prints the estimated p-value at
#'   every candidate allocation. Default is \code{FALSE}.
#'
#' @return A list containing:
#' \describe{
#'   \item{rct.prob}{A Monte Carlo estimate of the p-value under random
#'     assignment (no confounding), estimated separately from the candidates
#'     below.}
#'   \item{max.prob}{The largest Monte Carlo estimate across the candidate
#'     allocations.}
#'   \item{maximizer}{The allocation achieving \code{max.prob}.}
#'   \item{obs.stat}{The observed value of T.}
#'   \item{obs.table}{The observed table, as supplied.}
#'   \item{gamma}{The supplied \code{gamma}.}
#'   \item{delta}{The supplied \code{delta}.}
#' }
#'
#' @details
#' This function requires \code{treatment.scores}, \code{outcome.scores} and
#' \code{delta} to be nondecreasing, for the same reason as
#' \code{\link{exact.score.sen.IxJ}}: only then is the candidate set of
#' allocations guaranteed to contain the worst case. The candidate set is the
#' same small set used there -- a single allocation when the outcome is
#' binary, or N allocations otherwise -- excluding the allocation u = 0,
#' which is instead estimated separately and returned as \code{rct.prob}.
#'
#' Being a Monte Carlo method, the returned values are estimates that
#' approach the exact values of \code{\link{exact.score.sen.IxJ}} as
#' \code{mc.iteration} increases, and will vary between calls unless the
#' random seed is fixed with \code{set.seed()} beforehand.
#'
#' @examples
#' \dontrun{
#' obs.table <- matrix(c(12, 3, 0, 18, 12, 3, 17, 25, 4), nrow = 3, byrow = TRUE)
#' result <- sampling.score.sen.IxJ(obs.table = obs.table, gamma = log(2),
#'                                   delta = c(0, 1, 1),
#'                                   treatment.scores = c(0, 0.25, 1.5),
#'                                   outcome.scores = c(0, 1, 1.5),
#'                                   mc.iteration = 2000)
#' c(result$rct.prob, result$max.prob)
#' }
#'
#' @seealso
#' \code{\link{exact.score.sen.IxJ}} for the exact computation,
#' \code{\link{sampling.general.sen.IxJ}} for test statistics or treatment
#' effects that are not monotone
#'
#' @export
sampling.score.sen.IxJ = function(obs.table, gamma,delta, row = "treatment", mc.iteration = 20000,
                                  treatment.scores, outcome.scores, verbose = FALSE) {
  
  # Validate 'obs.table'
  if (!is.matrix(obs.table) && !inherits(obs.table, "table")) {
    stop("'obs.table' must be either a matrix or a table object.")
  }
  if (is.matrix(obs.table)) {
    if (!is.numeric(obs.table) || length(dim(obs.table)) != 2) {
      stop("If 'obs.table' is a matrix, it must be a two-dimensional numeric matrix.")
    }
  }
  if (inherits(obs.table, "table") && length(dim(obs.table)) != 2) {
    stop("'obs.table' must be a two-dimensional table object.")
  }
  obs.table <- as.matrix(obs.table)
  
  # Validate 'row'
  if (!is.character(row) || length(row) != 1 || !(row %in% c("outcome", "treatment"))) {
    stop("'row' must be a single string: either 'outcome' or 'treatment'.")
  }
  
  # Define I, J and margins
  if (row == "treatment") {
    I <- nrow(obs.table)
    J <- ncol(obs.table)
    outcome_margins <- colSums(obs.table)
    treatment_margins <- rowSums(obs.table)
  } else {
    I <- ncol(obs.table)
    J <- nrow(obs.table)
    outcome_margins <- rowSums(obs.table)
    treatment_margins <- colSums(obs.table)
  }
  
  # Validate gamma
  if(!is_scalar_numeric(gamma)){
    stop("gamma must be a scalar")
  }
  if(gamma < 0){
    stop("gamma must be nonnegative")
  }
  if(!is_binary_vector(delta)){
    stop("delta must be binary")
  }
  if (length(delta) != I) {
    stop(paste0("'delta' must be a numeric vector of length ", I, "."))
  }
  
  # Validate scores
  if (!is.numeric(treatment.scores) || length(treatment.scores) != I) {
    stop(paste0("'treatment.scores' must be a numeric vector of length ", I, "."))
  }
  if (!is.numeric(outcome.scores) || length(outcome.scores) != J) {
    stop(paste0("'outcome.scores' must be a numeric vector of length ", J, "."))
  }
  
  # Check for non-trivial scores
  if (length(unique(treatment.scores)) == 1) {
    warning("All treatment scores are identical. The test may not be meaningful.")
  }
  if (length(unique(outcome.scores)) == 1) {
    warning("All outcome scores are identical. The test may not be meaningful.")
  }
  
  # Validate monotonicity of scores
  if (!all(diff(treatment.scores) >= 0)) {
    stop("'treatment.scores' must have a monotone trend.")
  }
  if (!all(diff(outcome.scores) >= 0)) {
    stop("'outcome.scores' must have a monotone trend.")
  }
  
  # Validate monotonicity of delta
  if (!all(diff(delta) >= 0)) {
    stop("'delta' must have a monotone trend.")
  }
  
  # Orient the table so that rows are treatment levels; everything below, including
  # transform.fun and the tables drawn by the sampler, uses this orientation.
  if (row == "outcome") {
    obs.table <- t(obs.table)
  }
  
  transform.fun <- function(tb) sum(tb * outer(treatment.scores, outcome.scores))
  obs.stat <- transform.fun(obs.table)
  
  I <- nrow(obs.table)
  J <- ncol(obs.table)
  outcome_margins <- colSums(obs.table)
  treatment_margins <- rowSums(obs.table)
  
  u_space <- if (J == 2) matrix(c(0, as.numeric(outcome_margins[2])), ncol = 2) else .u_space_ordinal(outcome_margins)[-1, , drop = FALSE]
  
  # Iterate through all candidate u vectors and evaluate via sampling
  max.prob <- -Inf
  maximizer <- NULL
  prob.vector <- numeric(nrow(u_space))
  
  for (idx in seq_len(nrow(u_space))) {
    u_allocation <- as.numeric(u_space[idx, ])
    result <- single.u.sampling.general.IxJ(
      table = obs.table,
      row = "treatment",
      mc.iteration = mc.iteration,
      u_allocation = u_allocation,
      gamma = gamma,
      delta = delta,
      transform.fun = transform.fun
    )
    prob.vector[idx] <- result$final.prob.estimate
    if (prob.vector[idx] > max.prob) {
      max.prob <- prob.vector[idx]
      maximizer <- u_allocation
    }
    
    if (verbose) {
      cat("Current u =", u_allocation, "\n")
      cat("Current estimated p =", prob.vector[idx], "\n")
    }
  }
  
  # RCT (u = 0) probability
  rct.result <- single.u.sampling.general.IxJ(
    table = obs.table,
    row = "treatment",
    mc.iteration = mc.iteration,
    u_allocation = rep(0, J),
    gamma = gamma,
    delta = delta,
    transform.fun = transform.fun
  )
  if (!is.finite(max.prob)) max.prob <- 0
  
  return(list(
    rct.prob = rct.result$final.prob.estimate,
    max.prob = max.prob,
    maximizer = maximizer,
    obs.stat = obs.stat,
    obs.table = obs.table,
    delta = delta,
    gamma = gamma
  ))
}


#----------------------------------------------------------------------------------------------



#' Monte Carlo worst-case p-value for a general test statistic in an I x J table
#'
#' Monte Carlo counterpart of \code{\link{exact.general.sen.IxJ}} for any
#' count-based test statistic, for use when the reference set of tables at
#' least as extreme as observed is too large to enumerate exactly. Tables
#' with the observed margins are drawn by sequential importance sampling, and
#' the p-value at each candidate confounder allocation is estimated by
#' self-normalized importance sampling.
#'
#' @param obs.table A matrix or table object giving the observed contingency
#'   table.
#' @param gamma A nonnegative scalar, the sensitivity parameter.
#' @param delta A binary vector marking the treatment levels favoured by the
#'   confounder. Its length must equal the number of treatment levels (rows
#'   of \code{obs.table} if \code{row = "treatment"}, or columns if
#'   \code{row = "outcome"}).
#' @param row Either \code{"treatment"} or \code{"outcome"}, indicating which
#'   dimension of \code{obs.table} holds the treatment levels. Default is
#'   \code{"treatment"}.
#' @param mc.iteration The number of Monte Carlo draws used to estimate the
#'   p-value at each candidate allocation. Higher values increase accuracy at
#'   the cost of computation time. Default is 20000.
#' @param transform.fun A function of one table returning a numeric test
#'   statistic, invariant to permuting subjects within a treatment level.
#' @param verbose Logical; if \code{TRUE}, prints the estimated p-value at
#'   every candidate allocation. Default is \code{FALSE}.
#' @param u_space Candidate confounder allocations to evaluate: a matrix with
#'   one allocation per row (one column per outcome level). If \code{NULL}
#'   (default), the full candidate set of all u-bar vectors is used, which is
#'   guaranteed to contain the worst case for any count-based statistic, but
#'   can be costly since every row costs \code{mc.iteration} draws.
#'
#' @return A list containing:
#' \describe{
#'   \item{rct.prob}{A Monte Carlo estimate of the p-value under random
#'     assignment (no confounding), estimated separately from the candidates
#'     below.}
#'   \item{max.prob}{The largest Monte Carlo estimate across the candidate
#'     allocations in \code{u_space}.}
#'   \item{maximizer}{The allocation achieving \code{max.prob}.}
#'   \item{obs.stat}{The observed value of \code{transform.fun}.}
#'   \item{obs.table}{The observed table, as supplied.}
#'   \item{gamma}{The supplied \code{gamma}.}
#'   \item{delta}{The supplied \code{delta}.}
#' }
#'
#' @details
#' Unlike \code{\link{sampling.score.sen.IxJ}}, this function places no
#' monotonicity requirement on \code{delta} or on any scores, and accepts an
#' arbitrary \code{transform.fun} rather than only the ordinal statistic, at
#' the cost of a much larger default candidate set: every allocation in the
#' full candidate set costs \code{mc.iteration} draws, so for a large table,
#' either pass a smaller \code{u_space} (e.g. a chunk of
#' \code{\link{full.candidate.set}}), or use
#' \code{\link{sampling.score.sen.IxJ}} if \code{transform.fun} is the
#' ordinal statistic with nondecreasing scores and \code{delta}.
#'
#' Being a Monte Carlo method, the returned values are estimates that
#' approach the exact values of \code{\link{exact.general.sen.IxJ}} as
#' \code{mc.iteration} increases, and will vary between calls unless the
#' random seed is fixed with \code{set.seed()} beforehand.
#'
#' @examples
#' \dontrun{
#' obs.table <- matrix(c(12, 18, 5, 6, 12, 6, 6, 6, 15), nrow = 3, byrow = TRUE)
#' chisq <- function(t) { E <- outer(rowSums(t), colSums(t)) / sum(t); sum((t - E)^2 / E) }
#'
#' ## worst case over a chunk of the full candidate set
#' result <- sampling.general.sen.IxJ(obs.table = obs.table, gamma = 1,
#'                                     delta = c(0, 1, 1), transform.fun = chisq,
#'                                     mc.iteration = 2000,
#'                                     u_space = full.candidate.set(colSums(obs.table))[seq(1, 24975, by = 500), ])
#' c(result$rct.prob, result$max.prob); result$maximizer
#' }
#'
#' @seealso
#' \code{\link{exact.general.sen.IxJ}} for the exact computation,
#' \code{\link{sampling.score.sen.IxJ}} for ordinal tests with monotone scores
#'
#' @export
sampling.general.sen.IxJ = function(obs.table, gamma, delta, row = "treatment",
                                    mc.iteration = 20000,
                                    transform.fun,
                                    verbose = FALSE,
                                    u_space = NULL) {
  # Validate 'obs.table' is a matrix or table object
  if (!is.matrix(obs.table) && !inherits(obs.table, "table")) {
    stop("'obs.table' must be either a matrix or a table object.")
  }
  
  # If 'obs.table' is a matrix, ensure it's numeric and two-dimensional
  if (is.matrix(obs.table)) {
    if (!is.numeric(obs.table)) {
      stop("If 'obs.table' is a matrix, it must be numeric.")
    }
    if (length(dim(obs.table)) != 2) {
      stop("'obs.table' must be a two-dimensional matrix.")
    }
  }
  
  # If 'obs.table' is a table object, ensure it's two-dimensional
  if (inherits(obs.table, "table")) {
    if (length(dim(obs.table)) != 2) {
      stop("'obs.table' must be a two-dimensional table object.")
    }
  }
  
  # Ensure obs.table is a matrix for consistent handling
  obs.table <- as.matrix(obs.table)
  
  # Validate 'row' argument
  if (!is.character(row) || length(row) != 1) {
    stop("'row' must be a single string, either 'outcome' or 'treatment'.")
  }
  
  if (!(row %in% c("outcome", "treatment"))) {
    stop("'row' must be either 'outcome' or 'treatment'.")
  }
  
  # Validate mc.iteration
  if (!is.numeric(mc.iteration) || length(mc.iteration) != 1 || mc.iteration <= 0 || mc.iteration != round(mc.iteration)) {
    stop("'mc.iteration' must be a positive integer.")
  }
  
  # Validate transform.fun
  if (!is.function(transform.fun)) {
    stop("'transform.fun' must be a function.")
  }
  
  # Apply the transformation function to get the observed statistic, in the
  # orientation the user supplied
  obs.stat <- transform.fun(obs.table)
  
  # Orient the table so that rows are treatment levels.  stat.fun evaluates the
  # user's statistic on a table in that orientation, so the sampled tables and the
  # observed table are scored the same way.
  stat.fun <- transform.fun
  if (row == "outcome") {
    obs.table <- t(obs.table)
    stat.fun <- function(tb) transform.fun(t(tb))
  }
  
  # Get dimensions and margins
  I <- nrow(obs.table)
  J <- ncol(obs.table)
  outcome_margins <- colSums(obs.table)
  treatment_margins <- rowSums(obs.table)
  
  # Validate 'gamma'
  if (!is_scalar_numeric(gamma)) {
    stop("'gamma' must be a scalar.")
  }
  if(gamma<0){
    stop("gamma must be nonnegative")
  }
  if(!is_binary_vector(delta)){
    stop("delta must be a binary vector")
  }
  
  if (length(delta) != I) {
    stop(paste0("'delta' must have the same length as the number of treatments (", I, ")."))
  }
  
  if (is.null(u_space)) {
    u_space <- full.candidate.set(outcome_margins)
  } else {
    # Validate provided u_space
    if (!is.matrix(u_space)) {
      stop("'u_space' must be a matrix where each row is a 'u_allocation' vector.")
    }
    
    if (ncol(u_space) != J) {
      stop(paste0("Number of columns in 'u_space' (", ncol(u_space),
                  ") must match the number of outcomes (", J, ")."))
    }
    
    # Check that all elements are numeric and within [0, outcome_margins]
    if (!is.numeric(u_space)) {
      stop("'u_space' must be a numeric matrix.")
    }
    
    if (any(u_space < 0)) {
      stop("All elements of 'u_space' must be non-negative.")
    }
    
    # Validate that each u element is <= the corresponding outcome margin
    margin_mtx <- matrix(
      rep(outcome_margins, times = nrow(u_space)),
      nrow = nrow(u_space),
      byrow = TRUE
    )
    
    check_mat <- (u_space > margin_mtx)
    
    if (any(check_mat)) {
      # Find problem rows
      problem_idx <- which(check_mat, arr.ind = TRUE)
      problem_rows <- unique(problem_idx[, 1])
      
      stop(
        paste0(
          "Each element of 'u_space' must be <= the corresponding outcome margin.\n",
          "Problem row(s) in 'u_space': ",
          paste(problem_rows, collapse = ", ")
        )
      )
    }
  }
  
  # Iterate through all candidate u vectors and evaluate via sampling
  max.prob <- -Inf
  maximizer <- NULL
  prob.vector <- numeric(nrow(u_space))
  
  for (idx in seq_len(nrow(u_space))) {
    u_allocation <- as.numeric(u_space[idx, ])
    
    # Call the sampling function for this u_allocation
    result <- single.u.sampling.general.IxJ(
      table = obs.table,
      row = "treatment",
      mc.iteration = mc.iteration,
      u_allocation = u_allocation,
      gamma = gamma,
      delta = delta,
      transform.fun = stat.fun
    )
    
    prob.vector[idx] <- result$final.prob.estimate
    
    if (prob.vector[idx] > max.prob) {
      max.prob <- prob.vector[idx]
      maximizer <- u_allocation
    }
    
    if (verbose) {
      cat("Current u =", u_allocation, "\n")
      cat("Current estimated p =", prob.vector[idx], "\n")
    }
  }
  
  # RCT (u = 0) probability
  rct.result <- single.u.sampling.general.IxJ(
    table = obs.table,
    row = "treatment",
    mc.iteration = mc.iteration,
    u_allocation = rep(0, J),
    gamma = gamma,
    delta = delta,
    transform.fun = stat.fun
  )

  if (!is.finite(max.prob)) max.prob <- 0
  
  return(list(
    rct.prob = rct.result$final.prob.estimate,
    max.prob = max.prob,
    maximizer = maximizer,
    obs.stat = obs.stat,
    obs.table = obs.table,
    gamma = gamma,
    delta = delta
  ))
}



################################################################################
#
#                    NORMAL APPROXIMATION METHODS
#
# This section implements a normal approximation for sensitivity analysis in
# I×J contingency tables under the generic bias model, for the ordinal test
# statistic T = sum_ij w_i v_j N_ij. The exact mean and variance of T are
# computed from the exact moments of the cell counts, and the p-value is
# approximated from the resulting z-score, without enumerating or sampling
# any tables.
#
# Main Functions:
# - norm_single_u_allocation_p_value()
#    - The approximate p-value at a single u-allocation, from the exact mean
#      and variance of T and a normal approximation to its distribution
#
# - norm.score.sen.IxJ()
#    - Main interface: searches over the candidate set valid for an ordinal
#      test with monotone scores and delta, and returns the largest
#      approximate p-value together with the p-value under random assignment
#
# Supporting Functions (C++ Backend):
# - zr_moments_two_by_two()    : exact moments of the free cells, 2x2 tables
# - zr_moments_two_by_three()  : exact moments of the free cells, 2x3 tables
# - zr_moments_two_by_four()   : exact moments of the free cells, 2x4 tables
# - zr_moments_two_by_five()   : exact moments of the free cells, 2x5 tables
# - zr_moments_three_by_two()  : exact moments of the free cells, 3x2 tables
# - zr_moments_three_by_three(): exact moments of the free cells, 3x3 tables
# - zr_moments_four_by_two()   : exact moments of the free cells, 4x2 tables
#
################################################################################



## =============================================================================
## FUNCTION  construct_A_matrix
## PURPOSE   Linear map from the (I-1)(J-1) "free" cells (i < I, j < J) to all IJ
##           cells of a table with fixed margins, so that
##             vec(N) = A %*% vec(free cells) + (a constant from the margins).
##           Given the covariance matrix Sigma_free of the free cells (from the
##           C++ zr_moments_* functions), the covariance of all IJ cells is
##           A %*% Sigma_free %*% t(A). Cells are flattened row-major: index
##           (i-1)*J + j corresponds to cell (i, j), and rows i = I or columns
##           j = J are recovered from the margins N_i., N_.j.
## INPUT     I, J.
## OUTPUT    IJ x (I-1)(J-1) matrix of 0, 1, -1.
## EXAMPLE
##   A <- construct_A_matrix(2, 2)          # 4 x 1: c(1, -1, -1, 1)
##   ## for a 2 x 2 table, Var(N_ij) = Var(N_11) for every cell, and
##   ## Cov(N_11, N_12) = -Var(N_11): A %*% matrix(v) %*% t(A)
##
##   ## 3 x 3 table: free cells are N_11, N_12, N_21, N_22; given their
##   ## covariance matrix Sigma_free, A %*% Sigma_free %*% t(A) recovers the
##   ## covariance of all nine cells. Entry [9, 3] of the result is
##   ## Cov(N_33, N_13), since flattened index 9 is cell (3,3) and index 3 is
##   ## cell (1,3).
##   A <- construct_A_matrix(3, 3)
##   Sigma_free <- matrix(c(1,   -0.4, -0.3,  0.7,
##                          -0.4,  2,    0.7, -0.5,
##                          -0.3,  0.7,  0.5, -0.15,
##                           0.7, -0.5, -0.15, 0.4), nrow = 4)
##   var_mat <- A %*% Sigma_free %*% t(A)
## =============================================================================
construct_A_matrix <- function(I, J) {
  # A matrix relates the free cells to all cells via marginal constraints
  # The matrix has dimension IJ x (I-1)(J-1)
  df <- (I-1) * (J-1)  # degrees of freedom
  A <- matrix(0, nrow = I*J, ncol = df)
  
  # For each free cell (i,j) with i < I, j < J
  for (ell in 1:df) {
    # Map from linear index to (i,j)
    i <- floor((ell - 1) / (J - 1)) + 1
    j <- ((ell - 1) %% (J - 1)) + 1
    
    # Linear index in the full IJ vector
    idx <- (i - 1) * J + j
    A[idx, ell] <- 1
  }
  
  # Now handle the constrained cells
  for (idx in 1:(I*J)) {
    # Convert linear index to (i,j)
    i <- floor((idx - 1) / J) + 1
    j <- ((idx - 1) %% J) + 1
    
    # Skip free cells (already handled above)
    if (i < I && j < J) next
    
    if (i < I && j == J) {
      # Right column (except corner): z_{i,J} = n_i - sum(z_{i,j} for j=1 to J-1)
      for (jj in 1:(J-1)) {
        ell_jj <- (i - 1) * (J - 1) + jj
        A[idx, ell_jj] <- -1
      }
    } else if (i == I && j < J) {
      # Bottom row (except corner): z_{I,j} = m_j - sum(z_{i,j} for i=1 to I-1)
      for (ii in 1:(I-1)) {
        ell_ii <- (ii - 1) * (J - 1) + j
        A[idx, ell_ii] <- -1
      }
    } else if (i == I && j == J) {
      # Corner cell: z_{I,J} = N - sum(all other cells)
      # This equals N - sum(n_i for i=1 to I-1) - sum(m_j for j=1 to J-1) + sum(all free cells)
      for (ell in 1:df) {
        A[idx, ell] <- 1
      }
    }
  }
  
  return(A)
}



#--------------------------------------------------------------------------------------------------


## =============================================================================
## FUNCTION  create_A_vector
## PURPOSE   Flatten the score matrix (w_i v_j) row-major into a vector a of length
##           IJ, so that the ordinal statistic is T = a' vec(N), its mean is
##           a' E[vec(N)] and its variance a' Cov(vec(N)) a. Uses the same
##           row-major flattening as construct_A_matrix, index (i-1)*J + j
##           corresponding to cell (i, j), so the two can be combined directly.
## INPUT     I, J, treatment.scores (w), outcome.scores (v).
## OUTPUT    numeric vector of length IJ; entry (i-1)*J + j equals w_i * v_j.
## EXAMPLE
##   create_A_vector(2, 2, treatment.scores = c(1, 2), outcome.scores = c(1, 3))
##   ## c(1, 3, 2, 6)
##   identical(create_A_vector(3, 2, c(0, 1, 2), c(0, 1)),
##             as.vector(t(outer(c(0, 1, 2), c(0, 1)))))
##   ## TRUE
## =============================================================================

create_A_vector <- function(I, J, treatment.scores, outcome.scores) {
  # Full vector for all IJ cells
  A_full <- numeric(I*J)
  
  # For each cell in the table
  for (idx in 1:(I*J)) {
    i <- floor((idx - 1) / J) + 1
    j <- ((idx - 1) %% J) + 1
    
    # A_idx = αᵢ * βⱼ
    A_full[idx] <- treatment.scores[i] * outcome.scores[j]
  }
  
  return(A_full)
}




#-----------------------------------------------------------------------------------------------------------


#' Normal-approximation p-value at one confounder allocation, for an ordinal test
#'
#' Computes the exact mean and variance of the ordinal test statistic
#' T = sum_ij w_i v_j N_ij under the generic bias model, at one allocation of
#' an unmeasured confounder, and returns the upper-tail p-value from the
#' normal approximation z = (T_obs - mu_T) / sqrt(var_T). Supported
#' dimensions are 2x2, 2x3, 2x4, 2x5, 3x2, 3x3 and 4x2, the tables for which
#' the exact moments are implemented in C++.
#'
#' @param obs.table A matrix or table object giving the observed contingency
#'   table.
#' @param gamma A nonnegative scalar, the sensitivity parameter.
#' @param delta A binary, nondecreasing vector marking the treatment levels
#'   favored by the confounder. Its length must equal the number of
#'   treatment levels (rows of \code{obs.table} if \code{row = "treatment"},
#'   or columns if \code{row = "outcome"}).
#' @param u_allocation A numeric vector giving the confounder allocation
#'   (u-bar), with one nonnegative entry per outcome level, each no larger
#'   than its outcome margin.
#' @param row Either \code{"treatment"} or \code{"outcome"}, indicating which
#'   dimension of \code{obs.table} holds the treatment levels. Default is
#'   \code{"treatment"}.
#' @param treatment.scores A nondecreasing numeric vector of scores for the
#'   treatment levels, one per row of the oriented table.
#' @param outcome.scores A nondecreasing numeric vector of scores for the
#'   outcome levels, one per column of the oriented table.
#' @param shared_divisor A scaling constant passed to the underlying C++
#'   moment functions for numerical stability; it cancels in every ratio
#'   computed here and does not affect the result. Default is 1e6.
#'
#' @return A list containing:
#' \describe{
#'   \item{T_obs}{The observed value of T.}
#'   \item{mu_T}{The exact mean of T under the generic bias model at
#'     \code{u_allocation}.}
#'   \item{var_T}{The exact variance of T under the generic bias model at
#'     \code{u_allocation}.}
#'   \item{z_score}{\code{(T_obs - mu_T) / sqrt(var_T)}.}
#'   \item{p_value}{The upper-tail p-value \code{1 - pnorm(z_score)}.}
#'   \item{treatment.scores}{The supplied \code{treatment.scores}.}
#'   \item{outcome.scores}{The supplied \code{outcome.scores}.}
#' }
#'
#' @details
#' The mean and variance of the free cells (those not determined by the
#' margins) are computed exactly by a dimension-specific C++ function
#' (\code{zr_moments_<I>_by_<J>}); \code{\link{construct_A_matrix}} recovers
#' the mean vector and covariance matrix of all IJ cells from these, and
#' \code{\link{create_A_vector}} gives the score vector so that
#' \code{mu_T = a' mean} and \code{var_T = a' Sigma a}. At \code{gamma = 0},
#' these are the moments of the multivariate hypergeometric distribution.
#'
#' @examples
#' ## 3x2 table, u_allocation = (0, 35) is the allocation N_.2 = 35 at its margin
#' obs.table <- matrix(c(40, 20, 5, 5, 10, 10), nrow = 3, byrow = TRUE)
#' result <- norm_single_u_allocation_p_value(obs.table = obs.table, gamma = 1,
#'                                             delta = c(0, 0, 1),
#'                                             u_allocation = c(0, 35),
#'                                             treatment.scores = c(0, 1, 2),
#'                                             outcome.scores = c(0, 1))
#' c(result$mu_T, result$var_T, result$p_value)
#'
#' @seealso
#' \code{\link{norm.score.sen.IxJ}}, which calls this function at every
#' candidate allocation and reports the worst case
#'
#' @export
norm_single_u_allocation_p_value <- function(
    obs.table,
    gamma,
    delta,
    u_allocation,
    row = "treatment",
    treatment.scores, outcome.scores,
    shared_divisor = 1e6
) {
  # Validate 'obs.table' is a matrix or table object
  if (!is.matrix(obs.table) && !inherits(obs.table, "table")) {
    stop("'obs.table' must be either a matrix or a table object.")
  }
  
  # If 'obs.table' is a matrix, ensure it's numeric and two-dimensional
  if (is.matrix(obs.table)) {
    if (!is.numeric(obs.table)) {
      stop("If 'obs.table' is a matrix, it must be numeric.")
    }
    if (length(dim(obs.table)) != 2) {
      stop("'obs.table' must be a two-dimensional matrix.")
    }
  }
  
  # If 'obs.table' is a table object, ensure it's two-dimensional
  if (inherits(obs.table, "table")) {
    if (length(dim(obs.table)) != 2) {
      stop("'obs.table' must be a two-dimensional table object.")
    }
  }
  
  # Validate 'row' argument
  if (!is.character(row) || length(row) != 1) {
    stop("'row' must be a single string, either 'outcome' or 'treatment'.")
  }
  if (!(row %in% c("outcome", "treatment"))) {
    stop("'row' must be either 'outcome' or 'treatment'.")
  }
  
  if (row == "outcome") {
    ## our following operation assumes the treatments are in rows, so if the row=="outcome", we
    ## take the transpose, and then the following stays the same 
    obs.table <- t(obs.table)  
  }
  
  # Now get I, J
  I <- nrow(obs.table)
  J <- ncol(obs.table)
  
  # Margins
  N <- sum(obs.table)
  treatment_margins <- rowSums(obs.table)
  outcome_margins   <- colSums(obs.table)
  
  # Check delta, u_allocation lengths
  if (!is_scalar_numeric(gamma)) {
    stop("'gamma' must be a scalar.")
  }
  if(gamma < 0){
    stop("gamma must be nonnegative")
  }
  
  # (gamma = 0 needs no special treatment: the C++ moment functions return the
  #  multivariate hypergeometric moments for every u-bar in that case.)
  if(!is_binary_vector(delta)){
    stop("delta must be binary")
  }
  if (length(delta) != I) {
    stop(paste0("'delta' must have the same length as the number of treatments (", I, ")."))
  }
  
  
  if (!is.numeric(u_allocation) || !is.vector(u_allocation)) {
    stop("'u_allocation' must be a numeric vector.")
  }
  if (length(u_allocation) != J) {
    stop(paste0("'u_allocation' must have the same length as the number of outcomes (", J, ")."))
  }
  if (any(u_allocation < 0)) {
    stop("All elements of 'u_allocation' must be non-negative.")
  }
  if (any(u_allocation > outcome_margins)) {
    stop("Each 'u_allocation' entry must be <= the corresponding outcome margin.")
  }
  
  # Validate scores
  if (!is.numeric(treatment.scores) || length(treatment.scores) != I) {
    stop(paste0("'treatment.scores' must be a numeric vector of length ", I, "."))
  }
  if (!is.numeric(outcome.scores) || length(outcome.scores) != J) {
    stop(paste0("'outcome.scores' must be a numeric vector of length ", J, "."))
  }
  
  # Check for non-trivial scores
  if (length(unique(treatment.scores)) == 1) {
    warning("All treatment scores are identical. The test may not be meaningful.")
  }
  if (length(unique(outcome.scores)) == 1) {
    warning("All outcome scores are identical. The test may not be meaningful.")
  }
  
  # Validate monotonicity of scores
  if (!all(diff(treatment.scores) >= 0)) {
    stop("'treatment.scores' must have a monotone trend.")
  }
  if (!all(diff(outcome.scores) >= 0)) {
    stop("'outcome.scores' must have a monotone trend.")
  }
  
  # Validate monotonicity of gamma_delta
  if (!all(diff(delta) >= 0)) {
    stop("'delta' must have a monotone trend.")
  }
  
  
  
  gamma_delta = gamma*delta
  rcpp_out <- NULL
  
  if (I == 2 && J == 2) {
    # parse n1,n2,m1,m2,u1,u2
    n1 <- treatment_margins[1]
    n2 <- treatment_margins[2]
    m1 <- outcome_margins[1]
    m2 <- outcome_margins[2]
    u1 <- u_allocation[1]
    u2 <- u_allocation[2]
    
    # Call Rcpp
    rcpp_out <- zr_moments_two_by_two(
      n1, n2,
      m1, m2,
      N, gamma_delta,
      u1, u2,
      shared_divisor
    )
    # rcpp_out => c(mean, var)
    
    # Extract values
    mean_z1r1 <- rcpp_out[1]
    var_z1r1 <- rcpp_out[2]
    
    # Recover full mean vector
    mean_vec <- c(mean_z1r1, n1-mean_z1r1, m1-mean_z1r1, N-n1-m1+mean_z1r1)
    
    # Recover full covariance matrix
    A_matrix <- construct_A_matrix(I, J)
    Sigma_free <- matrix(var_z1r1, 1, 1)
    var_mat <- A_matrix %*% Sigma_free %*% t(A_matrix)
    
    # Create full A vector with scores for all cells
    A_vec <- create_A_vector(I, J, treatment.scores, outcome.scores)
    
    # Compute mean and variance using matrix operations
    mu_T <- as.numeric(t(A_vec) %*% mean_vec)
    var_T <- as.numeric(t(A_vec) %*% var_mat %*% A_vec)
    
  } else if (I == 2 && J == 3) {
    # parse n1,n2,m1,m2,m3,u1,u2,u3
    n1 <- treatment_margins[1]
    n2 <- treatment_margins[2]
    m1 <- outcome_margins[1]
    m2 <- outcome_margins[2]
    m3 <- outcome_margins[3]
    u1 <- u_allocation[1]
    u2 <- u_allocation[2]
    u3 <- u_allocation[3]
    
    rcpp_out <- zr_moments_two_by_three(
      n1, n2,
      m1, m2, m3,
      N, gamma_delta,
      u1, u2, u3,
      shared_divisor
    )
    # => c(z1r1_mean, z1r2_mean, z1r1_var, z1r2_var, cov(z1r1,z1r2))
    
    # Extract values
    mean_z1r1 <- rcpp_out[1]
    mean_z1r2 <- rcpp_out[2]
    var_z1r1 <- rcpp_out[3]
    var_z1r2 <- rcpp_out[4]
    cov_z1r1_z1r2 <- rcpp_out[5]
    
    # Recover full mean vector
    mean_vec <- c(
      mean_z1r1,
      mean_z1r2,
      n1 - mean_z1r1 - mean_z1r2,
      m1 - mean_z1r1,
      m2 - mean_z1r2,
      N - n1 - m1 - m2 + mean_z1r1 + mean_z1r2
    )
    
    # Recover full covariance matrix
    A_matrix <- construct_A_matrix(I, J)
    Sigma_free <- matrix(c(var_z1r1, cov_z1r1_z1r2,
                           cov_z1r1_z1r2, var_z1r2), nrow = 2)
    var_mat <- A_matrix %*% Sigma_free %*% t(A_matrix)
    
    # Create full A vector with scores for all cells
    A_vec <- create_A_vector(I, J, treatment.scores, outcome.scores)
    
    # Compute mean and variance using matrix operations
    mu_T <- as.numeric(t(A_vec) %*% mean_vec)
    var_T <- as.numeric(t(A_vec) %*% var_mat %*% A_vec)
    
  } else if (I == 2 && J == 4) {
    # parse n1,n2,m1,m2,m3,m4, u1,u2,u3,u4
    n1 <- treatment_margins[1]
    n2 <- treatment_margins[2]
    m1 <- outcome_margins[1]
    m2 <- outcome_margins[2]
    m3 <- outcome_margins[3]
    m4 <- outcome_margins[4]
    u1 <- u_allocation[1]
    u2 <- u_allocation[2]
    u3 <- u_allocation[3]
    u4 <- u_allocation[4]
    
    rcpp_out <- zr_moments_two_by_four(
      n1, n2,
      m1, m2, m3, m4,
      N, gamma_delta,
      u1, u2, u3, u4,
      shared_divisor
    )
    # => length=9:
    # [1] z1r1_mean, z1r2_mean, z1r3_mean,
    # [4] var(z1r1), var(z1r2), var(z1r3),
    # [7] cov(z1r1,z1r2), cov(z1r1,z1r3), cov(z1r2,z1r3)
    
    # Extract values
    mean_z1r1 <- rcpp_out[1]
    mean_z1r2 <- rcpp_out[2]
    mean_z1r3 <- rcpp_out[3]
    var_z1r1 <- rcpp_out[4]
    var_z1r2 <- rcpp_out[5]
    var_z1r3 <- rcpp_out[6]
    cov_z1r1_z1r2 <- rcpp_out[7]
    cov_z1r1_z1r3 <- rcpp_out[8]
    cov_z1r2_z1r3 <- rcpp_out[9]
    
    # Recover full mean vector
    mean_vec <- c(
      mean_z1r1,
      mean_z1r2,
      mean_z1r3,
      n1 - mean_z1r1 - mean_z1r2 - mean_z1r3,
      m1 - mean_z1r1,
      m2 - mean_z1r2,
      m3 - mean_z1r3,
      N - n1 - m1 - m2 - m3 + mean_z1r1 + mean_z1r2 + mean_z1r3
    )
    
    # Recover full covariance matrix
    A_matrix <- construct_A_matrix(I, J)
    Sigma_free <- matrix(c(var_z1r1, cov_z1r1_z1r2, cov_z1r1_z1r3,
                           cov_z1r1_z1r2, var_z1r2, cov_z1r2_z1r3,
                           cov_z1r1_z1r3, cov_z1r2_z1r3, var_z1r3), nrow = 3)
    var_mat <- A_matrix %*% Sigma_free %*% t(A_matrix)
    
    # Create full A vector with scores for all cells
    A_vec <- create_A_vector(I, J, treatment.scores, outcome.scores)
    
    # Compute mean and variance using matrix operations
    mu_T <- as.numeric(t(A_vec) %*% mean_vec)
    var_T <- as.numeric(t(A_vec) %*% var_mat %*% A_vec)
    
  } else if (I == 2 && J == 5) {
    # parse n1,n2,m1,m2,m3,m4,m5, u1..u5
    n1 <- treatment_margins[1]
    n2 <- treatment_margins[2]
    m1 <- outcome_margins[1]
    m2 <- outcome_margins[2]
    m3 <- outcome_margins[3]
    m4 <- outcome_margins[4]
    m5 <- outcome_margins[5]
    u1 <- u_allocation[1]
    u2 <- u_allocation[2]
    u3 <- u_allocation[3]
    u4 <- u_allocation[4]
    u5 <- u_allocation[5]
    
    rcpp_out <- zr_moments_two_by_five(
      n1, n2,
      m1, m2, m3, m4, m5,
      N, gamma_delta,
      u1, u2, u3, u4, u5,
      shared_divisor
    )
    # => length=14:
    #  means: z1r1, z1r2, z1r3, z1r4
    #  vars:  var1, var2, var3, var4
    #  covs:  c12, c13, c14, c23, c24, c34
    
    # Extract values
    mean_z1r1 <- rcpp_out[1]
    mean_z1r2 <- rcpp_out[2]
    mean_z1r3 <- rcpp_out[3]
    mean_z1r4 <- rcpp_out[4]
    var_z1r1 <- rcpp_out[5]
    var_z1r2 <- rcpp_out[6]
    var_z1r3 <- rcpp_out[7]
    var_z1r4 <- rcpp_out[8]
    cov_z1r1_z1r2 <- rcpp_out[9]
    cov_z1r1_z1r3 <- rcpp_out[10]
    cov_z1r1_z1r4 <- rcpp_out[11]
    cov_z1r2_z1r3 <- rcpp_out[12]
    cov_z1r2_z1r4 <- rcpp_out[13]
    cov_z1r3_z1r4 <- rcpp_out[14]
    
    # Recover full mean vector
    mean_vec <- c(
      mean_z1r1,
      mean_z1r2,
      mean_z1r3,
      mean_z1r4,
      n1 - mean_z1r1 - mean_z1r2 - mean_z1r3 - mean_z1r4,
      m1 - mean_z1r1,
      m2 - mean_z1r2,
      m3 - mean_z1r3,
      m4 - mean_z1r4,
      N - n1 - m1 - m2 - m3 - m4 + mean_z1r1 + mean_z1r2 + mean_z1r3 + mean_z1r4
    )
    
    # Recover full covariance matrix
    A_matrix <- construct_A_matrix(I, J)
    Sigma_free <- matrix(c(var_z1r1, cov_z1r1_z1r2, cov_z1r1_z1r3, cov_z1r1_z1r4,
                           cov_z1r1_z1r2, var_z1r2, cov_z1r2_z1r3, cov_z1r2_z1r4,
                           cov_z1r1_z1r3, cov_z1r2_z1r3, var_z1r3, cov_z1r3_z1r4,
                           cov_z1r1_z1r4, cov_z1r2_z1r4, cov_z1r3_z1r4, var_z1r4), nrow = 4)
    var_mat <- A_matrix %*% Sigma_free %*% t(A_matrix)
    
    # Create full A vector with scores for all cells
    A_vec <- create_A_vector(I, J, treatment.scores, outcome.scores)
    

    # Compute mean and variance using matrix operations
    mu_T <- as.numeric(t(A_vec) %*% mean_vec)
    var_T <- as.numeric(t(A_vec) %*% var_mat %*% A_vec)
    
  } else if (I == 3 && J == 2) {
    # parse n1,n2,n3,m1,m2,u1,u2
    n1 <- treatment_margins[1]
    n2 <- treatment_margins[2]
    n3 <- treatment_margins[3]
    m1 <- outcome_margins[1]
    m2 <- outcome_margins[2]
    u1 <- u_allocation[1]
    u2 <- u_allocation[2]
    
    
    rcpp_out <- zr_moments_three_by_two(
      n1, n2, n3,
      m1, m2,
      N, gamma_delta,
      u1, u2,
      shared_divisor
    )
    # => length=5:
    #  means: z1r1, z2r1
    #  vars:  var(z1r1), var(z2r1)
    #  covs:  cov(z1r1,z2r1)
    
    # Extract values
    mean_z1r1 <- rcpp_out[1]
    mean_z2r1 <- rcpp_out[2]
    var_z1r1 <- rcpp_out[3]
    var_z2r1 <- rcpp_out[4]
    cov_z1r1_z2r1 <- rcpp_out[5]
    
    # Recover full mean vector
    mean_vec <- c(
      mean_z1r1,
      n1 - mean_z1r1,
      mean_z2r1,
      n2 - mean_z2r1,
      m1 - mean_z1r1 - mean_z2r1,
      N - n1 - n2 - m1 + mean_z1r1 + mean_z2r1
    )
    
    # Recover full covariance matrix
    A_matrix <- construct_A_matrix(I, J)
    Sigma_free <- matrix(c(var_z1r1, cov_z1r1_z2r1,
                           cov_z1r1_z2r1, var_z2r1), nrow = 2)
    var_mat <- A_matrix %*% Sigma_free %*% t(A_matrix)
    
    # Create full A vector with scores for all cells
    A_vec <- create_A_vector(I, J, treatment.scores, outcome.scores)
    

    # Compute mean and variance using matrix operations
    mu_T <- as.numeric(t(A_vec) %*% mean_vec)
    var_T <- as.numeric(t(A_vec) %*% var_mat %*% A_vec)
   
    
  } else if (I == 3 && J == 3) {
    # parse n1,n2,n3,m1,m2,m3,u1,u2,u3
    n1 <- treatment_margins[1]
    n2 <- treatment_margins[2]
    n3 <- treatment_margins[3]
    m1 <- outcome_margins[1]
    m2 <- outcome_margins[2]
    m3 <- outcome_margins[3]
    u1 <- u_allocation[1]
    u2 <- u_allocation[2]
    u3 <- u_allocation[3]
    
    rcpp_out <- zr_moments_three_by_three(
      n1, n2, n3,
      m1, m2, m3,
      N, gamma_delta,
      u1, u2, u3,
      shared_divisor
    )
    # => length=14:
    #  means: z1r1, z1r2, z2r1, z2r2
    #  vars:  var(z1r1), var(z1r2), var(z2r1), var(z2r2)
    #  covs:  cov(z1r1,z1r2), cov(z1r1,z2r1), cov(z1r1,z2r2)
    #         cov(z1r2,z2r1), cov(z1r2,z2r2), cov(z2r1,z2r2)
    
    # Extract values
    mean_z1r1 <- rcpp_out[1]
    mean_z1r2 <- rcpp_out[2]
    mean_z2r1 <- rcpp_out[3]
    mean_z2r2 <- rcpp_out[4]
    var_z1r1 <- rcpp_out[5]
    var_z1r2 <- rcpp_out[6]
    var_z2r1 <- rcpp_out[7]
    var_z2r2 <- rcpp_out[8]
    cov_z1r1_z1r2 <- rcpp_out[9]
    cov_z1r1_z2r1 <- rcpp_out[10]
    cov_z1r1_z2r2 <- rcpp_out[11]
    cov_z1r2_z2r1 <- rcpp_out[12]
    cov_z1r2_z2r2 <- rcpp_out[13]
    cov_z2r1_z2r2 <- rcpp_out[14]
    
    # Recover full mean vector
    mean_vec <- c(
      mean_z1r1,
      mean_z1r2,
      n1 - mean_z1r1 - mean_z1r2,
      mean_z2r1,
      mean_z2r2,
      n2 - mean_z2r1 - mean_z2r2,
      m1 - mean_z1r1 - mean_z2r1,
      m2 - mean_z1r2 - mean_z2r2,
      N - n1 - n2 - m1 - m2 + mean_z1r1 + mean_z1r2 + mean_z2r1 + mean_z2r2
    )
    
    # Recover full covariance matrix
    A_matrix <- construct_A_matrix(I, J)
    Sigma_free <- matrix(c(var_z1r1, cov_z1r1_z1r2, cov_z1r1_z2r1, cov_z1r1_z2r2,
                           cov_z1r1_z1r2, var_z1r2, cov_z1r2_z2r1, cov_z1r2_z2r2,
                           cov_z1r1_z2r1, cov_z1r2_z2r1, var_z2r1, cov_z2r1_z2r2,
                           cov_z1r1_z2r2, cov_z1r2_z2r2, cov_z2r1_z2r2, var_z2r2), nrow = 4)
    var_mat <- A_matrix %*% Sigma_free %*% t(A_matrix)
    
    # Create full A vector with scores for all cells
    A_vec <- create_A_vector(I, J, treatment.scores, outcome.scores)
    

    # Compute mean and variance using matrix operations
    mu_T <- as.numeric(t(A_vec) %*% mean_vec)
    var_T <- as.numeric(t(A_vec) %*% var_mat %*% A_vec)
    
  } else if (I == 4 && J == 2) {
    # parse n1,n2,n3,n4,m1,m2,u1,u2
    n1 <- treatment_margins[1]
    n2 <- treatment_margins[2]
    n3 <- treatment_margins[3]
    n4 <- treatment_margins[4]
    m1 <- outcome_margins[1]
    m2 <- outcome_margins[2]
    u1 <- u_allocation[1]
    u2 <- u_allocation[2]
    
    rcpp_out <- zr_moments_four_by_two(
      n1, n2, n3, n4,
      m1, m2,
      N, gamma_delta,
      u1, u2,
      shared_divisor
    )
    
    # Extract values
    mean_z1r1 <- rcpp_out[1]
    mean_z2r1 <- rcpp_out[2]
    mean_z3r1 <- rcpp_out[3]
    var_z1r1 <- rcpp_out[4]
    var_z2r1 <- rcpp_out[5]
    var_z3r1 <- rcpp_out[6]
    cov_z1r1_z2r1 <- rcpp_out[7]
    cov_z1r1_z3r1 <- rcpp_out[8]
    cov_z2r1_z3r1 <- rcpp_out[9]
    
    # Recover full mean vector
    mean_vec <- c(
      mean_z1r1,
      n1 - mean_z1r1,
      mean_z2r1,
      n2 - mean_z2r1,
      mean_z3r1,
      n3 - mean_z3r1,
      m1 - mean_z1r1 - mean_z2r1 - mean_z3r1,
      N - n1 - n2 - n3 - m1 + mean_z1r1 + mean_z2r1 + mean_z3r1
    )
    
    # Recover full covariance matrix
    A_matrix <- construct_A_matrix(I, J)
    Sigma_free <- matrix(c(var_z1r1, cov_z1r1_z2r1, cov_z1r1_z3r1,
                           cov_z1r1_z2r1, var_z2r1, cov_z2r1_z3r1,
                           cov_z1r1_z3r1, cov_z2r1_z3r1, var_z3r1), nrow = 3)
    var_mat <- A_matrix %*% Sigma_free %*% t(A_matrix)
    
    # Create full A vector with scores for all cells
    A_vec <- create_A_vector(I, J, treatment.scores, outcome.scores)
    

    # Compute mean and variance using matrix operations
    mu_T <- as.numeric(t(A_vec) %*% mean_vec)
    var_T <- as.numeric(t(A_vec) %*% var_mat %*% A_vec)
    
  } else {
    stop("Dimension (I x J) not supported by current package.")
  }
  
  ## compute the T_obs 
  transform.fun = function(tb){
    scores  = outer(treatment.scores, outcome.scores, FUN="*")
    return(sum(tb*scores))
  }
  T_obs = transform.fun(obs.table)
  
  z_score <- (T_obs - mu_T) / sqrt(var_T)
  # One-sided upper-tail
  p_value <- 1 - pnorm(z_score)
  
  #### Return result ####
  return(
    list(
      T_obs   = T_obs,
      mu_T    = mu_T,
      var_T   = var_T,
      z_score = z_score,
      p_value = p_value,
      treatment.scores = treatment.scores,
      outcome.scores = outcome.scores
    )
  )
}


#----------------------------------------------------------------------------------

#' Normal-approximation worst-case p-value for an ordinal test in an I x J table
#'
#' Normal-approximation counterpart of \code{\link{exact.score.sen.IxJ}} for
#' the ordinal test statistic T = sum_ij w_i v_j N_ij: evaluates
#' \code{\link{norm_single_u_allocation_p_value}} at every candidate
#' confounder allocation and reports the largest approximate p-value, together
#' with the p-value under random assignment.
#'
#' @param obs.table A matrix or table object giving the observed contingency
#'   table.
#' @param gamma A nonnegative scalar, the sensitivity parameter.
#' @param delta A binary, nondecreasing vector marking the treatment levels
#'   favoured by the confounder. Its length must equal the number of
#'   treatment levels (rows of \code{obs.table} if \code{row = "treatment"},
#'   or columns if \code{row = "outcome"}).
#' @param row Either \code{"treatment"} or \code{"outcome"}, indicating which
#'   dimension of \code{obs.table} holds the treatment levels. Default is
#'   \code{"treatment"}.
#' @param treatment.scores A nondecreasing numeric vector of scores for the
#'   treatment levels, one per row of the oriented table.
#' @param outcome.scores A nondecreasing numeric vector of scores for the
#'   outcome levels, one per column of the oriented table.
#' @param shared_divisor A scaling constant passed to the underlying C++
#'   moment functions; it cancels in every ratio and does not affect the
#'   result. Default is 1e6.
#' @param u_space Candidate confounder allocations to evaluate: a matrix with
#'   one allocation per row. If \code{NULL} (default), the candidate set
#'   valid for an ordinal test with nondecreasing scores and delta is used
#'   (a single allocation when the outcome is binary, or N + 1 allocations
#'   otherwise), built the same way as in \code{\link{exact.score.sen.IxJ}}.
#' @param verbose Logical; if \code{TRUE}, prints the allocation and
#'   approximate p-value at every candidate. Default is \code{FALSE}.
#'
#' @return A list containing:
#' \describe{
#'   \item{T_obs}{The observed value of T.}
#'   \item{RCT.mean}{The exact mean of T under random assignment (u = 0).}
#'   \item{RCT.var}{The exact variance of T under random assignment (u = 0).}
#'   \item{RCT.prob}{The approximate p-value under random assignment.}
#'   \item{maximizer}{The allocation achieving \code{max.prob}.}
#'   \item{max.prob}{The largest approximate p-value across the candidate
#'     allocations.}
#'   \item{max.mean}{The exact mean of T at \code{maximizer}.}
#'   \item{max.var}{The exact variance of T at \code{maximizer}.}
#'   \item{treatment.scores}{The supplied \code{treatment.scores}.}
#'   \item{outcome.scores}{The supplied \code{outcome.scores}.}
#' }
#'
#' @details
#' This function requires \code{treatment.scores}, \code{outcome.scores} and
#' \code{delta} to be nondecreasing, for the same reason as
#' \code{\link{exact.score.sen.IxJ}}: only then is the candidate set of
#' allocations guaranteed to contain the worst case. At each candidate, the
#' exact mean and variance of T are computed and the p-value is approximated
#' from the resulting z-score; see \code{\link{norm_single_u_allocation_p_value}}
#' for that computation.
#'
#' @examples
#' ## 2x3 table with ordinal scores
#' obs.table <- matrix(c(10, 20, 30, 15, 25, 10), nrow = 2, byrow = TRUE)
#' result <- norm.score.sen.IxJ(obs.table = obs.table, gamma = 0.5,
#'                               delta = c(0, 1),
#'                               treatment.scores = c(0, 1),
#'                               outcome.scores = c(0, 1, 2))
#' c(result$RCT.prob, result$max.prob); result$maximizer
#'
#' @seealso
#' \code{\link{exact.score.sen.IxJ}} for the exact computation,
#' \code{\link{norm_single_u_allocation_p_value}} for the computation at a
#' single allocation
#'
#' @export
norm.score.sen.IxJ <- function(
    obs.table,
    gamma,
    delta,
    row = "treatment",
    treatment.scores,
    outcome.scores,
    shared_divisor = 1e6,
    u_space = NULL,
    verbose = FALSE
) {

  ## ----------------------------------------------------------------
  ## 1) Validate obs.table
  ## ----------------------------------------------------------------
  if (!is.matrix(obs.table) && !inherits(obs.table, "table")) {
    stop("'obs.table' must be either a matrix or a table object.")
  }
  if (is.matrix(obs.table)) {
    if (!is.numeric(obs.table)) {
      stop("If 'obs.table' is a matrix, it must be numeric.")
    }
    if (length(dim(obs.table)) != 2) {
      stop("'obs.table' must be two-dimensional.")
    }
  }
  if (inherits(obs.table, "table")) {
    if (length(dim(obs.table)) != 2) {
      stop("'obs.table' must be a two-dimensional table object.")
    }
  }
  
  # Convert to matrix to ensure consistent handling
  obs.table <- as.matrix(obs.table)
  
  ## ----------------------------------------------------------------
  ## 2) Validate row orientation; determine I, J, margins
  ## ----------------------------------------------------------------
  if (!is.character(row) || length(row) != 1) {
    stop("'row' must be a single string, either 'outcome' or 'treatment'.")
  }
  if (!(row %in% c("outcome", "treatment"))) {
    stop("'row' must be either 'outcome' or 'treatment'.")
  }
  
  # Extract dimensions and margins based on orientation
  if (row == "treatment") {
    I <- nrow(obs.table)
    J <- ncol(obs.table)
    outcome_margins    <- colSums(obs.table)
    treatment_margins  <- rowSums(obs.table)
  } else {
    # row == "outcome"
    I <- ncol(obs.table)
    J <- nrow(obs.table)
    outcome_margins    <- rowSums(obs.table)
    treatment_margins  <- colSums(obs.table)
  }
  
  ## ----------------------------------------------------------------
  ## 3) Validate gamma
  ## ----------------------------------------------------------------
  if (!is_scalar_numeric(gamma)) {
    stop("'gamma' must be a scalar.")
  }
  if(gamma < 0){
    stop("gamma must be nonnegative")
  }
  if(!is_binary_vector(delta)){
    stop("delta must be a binary vector")
  }
  if (length(delta) != I) {
    stop(paste0("'delta' must have the same length as the number of treatments (", I, ")."))
  }
  
  
  ## ----------------------------------------------------------------
  ## 4) Validate scores
  ## ----------------------------------------------------------------
  if (!is.numeric(treatment.scores) || length(treatment.scores) != I) {
    stop(paste0("'treatment.scores' must be a numeric vector of length ", I, "."))
  }
  if (!is.numeric(outcome.scores) || length(outcome.scores) != J) {
    stop(paste0("'outcome.scores' must be a numeric vector of length ", J, "."))
  }
  
  # Check for monotonicity
  if (!all(diff(treatment.scores) >= 0)) {
    stop("'treatment.scores' must have a monotone trend.")
  }
  if (!all(diff(outcome.scores) >= 0)) {
    stop("'outcome.scores' must have a monotone trend.")
  }
  if(!all(diff(delta) >= 0)){
    stop("delta must have a monotone trend")
  }
  
  ## ----------------------------------------------------------------
  ## 5) If u_space is NULL, use the candidate set of the theorems for ordinal tests
  ##    (the normal approximation is only provided for the ordinal statistic
  ##    T = sum_ij w_i v_j N_ij, and the monotonicity of the scores and of delta was
  ##    checked in step 4):
  ##      J = 2  -> Theorem 3: u-bar+ = (0, N_.2), plus the reference point (0, 0);
  ##      J >= 3 -> Theorem 2: U_O, the N + 1 vectors filling outcome level J first.
  ## ----------------------------------------------------------------
  if (is.null(u_space)) {
    u_space <- if (J == 2) .u_space_sign_score(outcome_margins) else .u_space_ordinal(outcome_margins)
    
    if (verbose) {
      cat("Using the candidate set of", if (J == 2) "Theorem 3 (u+)" else "Theorem 2 (U_O)", "for J =", J, "\n")
      cat("Total allocations to check:", nrow(u_space), "\n")
    }
  } else {
    ## ----------------------------------------------------------------
    ## 6) If u_space is provided, validate it
    ## ----------------------------------------------------------------
    if (!is.matrix(u_space)) {
      stop("'u_space' must be a matrix where each row is a 'u_allocation' vector.")
    }
    if (ncol(u_space) != J) {
      stop(paste0("Number of columns in 'u_space' (", ncol(u_space),
                  ") must match the number of outcomes (", J, ")."))
    }
    if (!is.numeric(u_space)) {
      stop("'u_space' must be a numeric matrix.")
    }
    if (any(u_space < 0)) {
      stop("All elements of 'u_space' must be non-negative.")
    }
    
    # Compare each row to outcome_margins
    margin_mtx <- matrix(
      rep(outcome_margins, times = nrow(u_space)),
      nrow = nrow(u_space),
      byrow = TRUE
    )
    check_mat <- (u_space > margin_mtx)
    if (any(check_mat)) {
      problem_idx <- which(check_mat, arr.ind = TRUE)
      problem_rows <- unique(problem_idx[, 1])
      stop(
        paste0(
          "Each element of 'u_space' must be <= the corresponding outcome margin.\n",
          "Problem row(s) in 'u_space': ",
          paste(problem_rows, collapse = ", ")
        )
      )
    }
  }
  
  ## ----------------------------------------------------------------
  ## 7) Main computation loop: evaluate normal approximation for each u
  ## ----------------------------------------------------------------
  
  # Initialize tracking variables
  max.mean <- NA
  max.var <- NA
  max.prob <- -Inf
  maximizer <- NULL
  
  # Iterate through each u-allocation
  for (i in seq_len(nrow(u_space))) {
    
    current_u <- as.numeric(u_space[i, ])
    
    # Compute normal approximation p-value for current u-allocation
    out <- norm_single_u_allocation_p_value(
      obs.table               = obs.table,
      gamma               =gamma,
      delta         = delta,
      u_allocation        = current_u,
      row                 = row,
      treatment.scores    = treatment.scores,  
      outcome.scores      = outcome.scores,    
      shared_divisor      = shared_divisor
    )
    
    # Update maximum if current p-value is larger
    if (out$p_value > max.prob) {
      maximizer <- current_u
      max.prob <- out$p_value
      max.mean <- out$mu_T
      max.var  <- out$var_T
    }
    
    # Print progress if verbose mode is enabled
    if (verbose) {
      cat("Processed u_allocation #", i, ": ", current_u, "\n")
      cat("Current maximizer:", maximizer, "\n")
      cat("Current max p-value:", max.prob, "\n\n")
    }
  }
  
  ## ----------------------------------------------------------------
  ## 8) Compute RCT baseline (all u = 0) for comparison
  ## ----------------------------------------------------------------
  RCT.out <- norm_single_u_allocation_p_value(
    obs.table               = obs.table,
    gamma               = gamma,
    delta         = delta,
    u_allocation        = rep(0, J),
    row                 = row,
    treatment.scores    = treatment.scores,  
    outcome.scores      = outcome.scores,    
    shared_divisor      = shared_divisor
  )
  
  # Extract RCT results
  RCT.mean <- RCT.out$mu_T
  RCT.var  <- RCT.out$var_T
  RCT.prob <- RCT.out$p_value
  T_obs    <- RCT.out$T_obs
  
  
  if (!is.finite(max.prob)) max.prob <- 0
  
  # Return comprehensive results
  return(list(
    T_obs             = T_obs,              # Observed test statistic
    RCT.mean          = RCT.mean,           # Mean under RCT
    RCT.var           = RCT.var,            # Variance under RCT
    RCT.prob          = RCT.prob,           # P-value under RCT
    maximizer         = maximizer,          # u-allocation maximizing p-value
    max.prob          = max.prob,           # Maximum p-value (sensitivity bound)
    max.mean          = max.mean,           # Mean at maximizer
    max.var           = max.var,            # Variance at maximizer
    treatment.scores  = treatment.scores,   # Treatment scores used
    outcome.scores    = outcome.scores      # Outcome scores used
  ))
}







########## End of script ######################################################
###############################################################################






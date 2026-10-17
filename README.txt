README - Code and Data for "Exact, Nonparametric Sensitivity Analysis for
         Observational Studies of Contingency Tables"
Authors: Elaine K. Chiu and Hyunseung Kang

------------------------------------------------------------
DESCRIPTION
------------------------------------------------------------
This repository contains the R and C++ code and the derived data set needed to
reproduce the numerical results reported in the paper "Exact, Nonparametric
Sensitivity Analysis for Observational Studies of Contingency Tables",
including the figures and tables in the main text and in the Supplementary
Material. The repository is self-contained: the methods are implemented in
the /library/ folder and are loaded directly by every script, so no package
installation is needed.

All scripts are run with the R working directory set to the top-level folder
of this repository, i.e. the folder that contains this README.txt. All file
paths in this README.txt and in the scripts are relative to that folder.

The material is organized into six folders. The /library/ folder contains
the implementation of the methods in the paper (the sensitivityIxJ code): the
exact worst-case p-value, the sampling-based estimators, and the normal
approximation based on the exact moments of the cell counts. The
/data_analysis/ folder contains the analysis of the Early Childhood
Longitudinal Study, Kindergarten Class of 1998-99 (ECLS-K). The
/power_simulation/ folder contains the power and size simulations, the study
of the computational savings of the kernel representation, and a
counterexample showing that, under a more general sensitivity model, the
maximizer of the p-value need not be a binary vector. The /feasibility/
folder contains the running-time and memory benchmarks of the exact
computation. The /development_check/ folder contains the scripts that verify
the code. The /legacy/ folder contains the previous version of the code,
kept for comparison.

------------------------------------------------------------
REQUIREMENTS
------------------------------------------------------------
- R version 4.4.3 or higher
- A C++ compiler for Rcpp. The .cpp files are compiled at run time by
  Rcpp::sourceCpp; no package installation is needed.
- Required R packages:
    - Rcpp            (all scripts; compiles sensitivityIxJ.cpp and
                       SIS_methods.cpp)
    - multicool       (Kernel_Advantage_Simulation.R, size_control.R:
                       enumerates treatment permutations for the
                       permutation-based comparison)
    - BiasedUrn       (size_control.R: samples the multivariate extended
                       hypergeometric distribution and sensitivityIxJ_cpp_check.R)
    - arrangements    (Corner_Solution_Counter_Example.R)
    - sensitivitymv   (education_data_analysis.R: truncated product for the
                       joint p-value)
    - dplyr           (education_data_analysis.R)
    - ggplot2, gridExtra, patchwork
                      (Kernel_Advantage_Simulation.R, size_control.R: figures)
- Optional, only for re-deriving the profile-likelihood scores in
  education_data_analysis.R:
    - R package reticulate, pointing to a Python 3 installation with the
      Python packages gurobipy and numpy, and a Gurobi license. The optimized
      scores are recorded in the script and this step is off by default; see
      HOW TO REPRODUCE THE RESULTS, item 4.
- Code files in /library/ (paths relative to the top-level folder):
    - library/sensitivityIxJ.R      (sourced by every script)
    - library/sensitivityIxJ.cpp    (sourced by every script)
    - library/SIS_methods.cpp       (sourced by the scripts that use the
                                     sampling-based estimators)

------------------------------------------------------------
FOLDER STRUCTURE
------------------------------------------------------------
/library/
    sensitivityIxJ.R
        The R interface. The main functions are:
          exact.score.sen.IxJ     exact worst-case p-value of an ordinal
                                  test (the sign-score test when the outcome
                                  is binary); returns the worst-case p-value
                                  (max.prob), the maximizing u (maximizer),
                                  and the p-value at gamma = 0 (rct.prob).
          exact.general.sen.IxJ   exact worst-case p-value of a general
                                  count-based statistic (chi-square, G^2,
                                  Freeman-Halton, or a user-supplied
                                  function), maximizing over the full
                                  candidate set of u or over a user-supplied
                                  u_space.
          possible.table          enumerates the tables with the observed
                                  margins whose statistic is at or beyond
                                  the observed value.
          generic.I.by.J.sensitivity.point.probability
                                  the null probability of one table at one u.
          sampling.score.sen.IxJ, sampling.general.sen.IxJ
                                  the sampling-based estimators of the
                                  worst-case p-value (sequential importance
                                  sampling and its self-normalized version).
          norm.score.sen.IxJ      the normal approximation based on the exact
                                  first and second moments of the cell counts.
          full.candidate.set      constructs the full candidate set of u,
                                  of size prod_j (N_.j + 1).

    sensitivityIxJ.cpp
        The computational engine, compiled by Rcpp::sourceCpp. Implements the
        kernel representation of the exact null distribution (point_prob_bd,
        build_W_cpp, alpha_bd) and the exact moments of the cell counts
        (zr_moments_*, zu_moments_*). Exact computation works for any I x J
        table; the moment routines cover 2x2 to 2x5, 3x2, 3x3, and 4x2
        tables.

    SIS_methods.cpp
        The sequential-importance-sampling table sampler adapted from
        Chen, Diaconis, Holmes and Liu (2005) and Eisinger and Chen (2017),
        used as the proposal distribution by sampling.score.sen.IxJ and
        sampling.general.sen.IxJ and by the sampling comparison in
        Kernel_Advantage_Simulation.R.

/data_analysis/
    education_data_analysis.R
        The ECLS-K data analysis. It first documents the ECLS-K variables,
        cleans the data, and constructs the 3x3 tables for Black girls and
        Black boys in low-income, single-parent families. It then computes
        the sensitivity analysis of the ordinal test under the prior scores
        and the profile-likelihood scores, including the joint p-value over
        the two subgroups.

    education.data.csv
        The extract of the ECLS-K public-use child file used by the analysis
        (one row per child). Columns are the ECLS-K variables named in the
        script; codes are as in the ECLS-K codebook.

/power_simulation/
    power_analysis_weighted_four_gen_oracle.R
    power_analysis_weighted_four_gen_binary_outcome_1_2.R
    power_analysis_weighted_four_gen_binary_outcome_2_3.R
    power_analysis_weighted_four_gen_Fisher_outcome_1_2.R
    power_analysis_weighted_four_gen_Fisher_outcome_2_3.R
        Power of the tests in the sensitivity analysis. One script per test:
          oracle               3x3 (Opt), the ordinal test on the full table
                               whose scores equal the data-generating scores
          binary_outcome_1_2   3x2 (V1), outcome levels 1 and 2 collapsed
          binary_outcome_2_3   3x2 (V2), outcome levels 2 and 3 collapsed
          Fisher_outcome_1_2   2x2 (V1), Fisher's exact test after collapsing
                               treatment levels 2-3 and outcome levels 1-2
          Fisher_outcome_2_3   2x2 (V2), Fisher's exact test after collapsing
                               treatment levels 2-3 and outcome levels 2-3
        Each script takes three command-line arguments (starting seed, gamma,
        number of Monte Carlo replications), generates 3x3 tables from the
        linear-by-linear association model, computes the exact worst-case
        p-value, and prints the number of rejections at level 0.05. The
        data-generating parameters (w*, v*, lambda^Z, lambda^r, beta) are set
        near the top of each script; the treatment margins are the size
        arguments of the three sample() calls in the Monte Carlo loop.

    submit.txt, params.txt, command_generating_summary.txt
        The HTCondor submit description, its (seed, gamma) grid (10 seeds x
        11 gammas; 100 replications per job, so 1,000 per gamma), and the
        bash script that sums the rejection counts of one test over seeds
        into power_simulation/summary_<test>.txt. Power at a gamma is the
        total in that file divided by 1,000.

    size_control.R
        For each configuration of margins and gamma, samples 3x2 tables from
        the null distribution under the worst-case u (via the multivariate
        extended hypergeometric distribution, package BiasedUrn) and compares
        the empirical size of the exact worst-case p-value with that of the
        normal approximation.

    Kernel_Advantage_Simulation.R
        The computational savings of the kernel representation of the exact
        null distribution.

    Corner_Solution_Counter_Example.R
        A numerical counterexample for a sensitivity model more general than
        the generic bias model. Enumerates all treatment assignments for a
        six-subject example and shows that, under the exposure-dose model,
        the maximizer of the ordinal test is not a binary vector u (it is not
        at a corner of [0,1]^N).

/feasibility/
    ordinal_feasibility_3x3.R, ordinal_feasibility_3x4.R,
    ordinal_feasibility_4x2.R, ordinal_feasibility_5x2.R
    general_test_feasibility_3x3.R, general_test_feasibility_3x4.R,
    general_test_feasibility_4x2.R, general_test_feasibility_5x2.R
        End-to-end running time and memory of the exact worst-case p-value
        for increasing sample sizes N, for the ordinal test and for the
        chi-square statistic (maximum over the full candidate set of u). Each
        run appends one row to ordinal_test_feasibility.csv or
        general_test_feasibility.csv; the column layout is documented in the
        header of each script. The chi-square sweeps are chunked and
        checkpointed to .rds files so that an interrupted run can be resumed.

    ordinal_test_feasibility.csv, general_test_feasibility.csv
        Example benchmark results. Running the feasibility scripts appends
        new rows to these files, so users can examine the computational
        burden under their own CPU and RAM setting.

/development_check/

    sanity_check.R
        Checks the engine against results known independently of the paper,
        e.g. at gamma = 0 the null distribution is multivariate
        hypergeometric and the worst-case p-value equals fisher.test(), and
        the point probabilities sum to one over all tables with the observed
        margins for every (u, gamma).

    moment_check.R
        Checks the exact first and second moments of the cell counts used by
        the normal approximation against moments computed from the exact
        null distribution.

    verify_same_names.R
        Runs the previous engine (/legacy/) and the current engine on the
        tables that appear in the paper and checks that the reference sets
        and p-values agree; prints timings of both engines.
    
    sensitivityIxJ_cpp_check.R
        Unit tests of the C++ building blocks in library/sensitivityIxJ.cpp.

    sensitivityIxJ_R_check.R
        Unit tests of the R interface in library/sensitivityIxJ.R.

/legacy/
    sensitivityIxJ_orig.R, sensitivityIxJ_orig.cpp
        The previous release of the engine (one numerator per table
        dimension, enumeration in R). Kept so that the verification scripts
        can compare the two engines. Not used by any results script.

------------------------------------------------------------
DATA AVAILABILITY
------------------------------------------------------------
The data analysis uses the Early Childhood Longitudinal Study, Kindergarten
Class of 1998-99 (ECLS-K), Kindergarten-Eighth Grade public-use data file,
released by the National Center for Education Statistics (NCES), U.S.
Department of Education. The public-use file may be downloaded by anyone,
without registration or a data use agreement. No component of the analysis
relies on restricted-use data.

Study landing page:
    https://nces.ed.gov/ecls/kindergarten.asp

Public-use data and online codebook (NCES DataLab, "ECLS-K:1998-99
Kindergarten-Eighth Grade Public-Use File"):
    https://nces.ed.gov/datalab/onlinecodebook

User's manuals (variable definitions and coding):
    https://nces.ed.gov/ecls/data/eclsk_k8_manual_part1.pdf
    https://nces.ed.gov/ecls/data/eclsk_k8_manual_part2.pdf

Variables used:
    Treatment            P1PRIMPK, primary type of nonparental care before
                         kindergarten, recoded to three levels: center-based
                         care (code 6), relative care (codes 1-2), no
                         nonparental care (code 0)
    Outcome              C1R4MPF, highest mathematics proficiency level
                         mastered in fall kindergarten, levels 1-3
    Subgroup definition  GENDER, RACE (Black), P1HFAMIL (single-parent
                         household, codes 3-4), WKINCOME (household income
                         at or below $50,000)
    Study population     Black girls and Black boys in low-income,
                         single-parent families with complete data on the
                         variables above: N = 94 girls and N = 112 boys

Obtaining the file:
    From the DataLab online codebook, select the ECLS-K:1998-99 K-8 public-use
    file and download the child-level CSV (ECLSK_98_99_K8_CHILD_v1_0.csv).
    The beginning of data_analysis/education_data_analysis.R shows the
    extraction of the analysis columns and the write of education.data.csv.

Derived analysis data:
    education.data.csv is the extracted file and is included in
    /data_analysis/, so the data analysis can be reproduced directly, without
    downloading the ECLS-K file.

------------------------------------------------------------
HOW TO REPRODUCE THE RESULTS
------------------------------------------------------------
1. Set the R working directory to the top-level folder of this repository
   (the folder that contains this README.txt), e.g.
       setwd("path/to/this/folder")
   Every script is run from there, and every script that uses the methods
   begins with
       library(Rcpp)
       sourceCpp("library/sensitivityIxJ.cpp")
       source("library/sensitivityIxJ.R")
   (Scripts that use the sampling-based estimators also compile
   library/SIS_methods.cpp; Corner_Solution_Counter_Example.R is
   self-contained.) Scripts are run by giving their path relative to the
   top-level folder, e.g.
       source("power_simulation/size_control.R")

2. To reproduce the computational comparison of the kernel representation
   and the sampling-based estimators, run
       source("power_simulation/Kernel_Advantage_Simulation.R")
   or step through the script from the top. The script is sequential: the
   helper functions defined at the beginning of each part are used by every
   experiment below them.

3. To reproduce the power results. Each power_analysis_*.R script takes
   (seed, gamma, replications) as arguments and prints, as its last line,
   "how many rejections? <n>".
     (a) On an HTCondor cluster, from the top-level folder:
             chmod +x power_simulation/power_analysis_*.R
             mkdir -p power_simulation/logs
         In power_simulation/submit.txt, set container_image to your own
         image (R >= 4.4, Rcpp, and a C++ compiler) and set
             test = oracle
         to the test of interest (oracle, binary_outcome_1_2,
         binary_outcome_2_3, Fisher_outcome_1_2, Fisher_outcome_2_3). Then
             condor_submit power_simulation/submit.txt
         queues one job per row of params.txt, with output in
         power_simulation/logs/<test>_gamma_<gamma>_seed_<seed>.out. When
         all jobs have finished,
             bash power_simulation/command_generating_summary.txt oracle
         writes power_simulation/summary_oracle.txt; power = total / 1,000.
     (b) Without HTCondor, keep the folder structure and run each cell in a
         terminal from the top-level folder (about one to two minutes per
         cell). The command below reproduces one cell of Table 3 of the
         Supplementary Material:
             INPUT   for s in 100 200 300 400 500 600 700 800 900 1000; do
                       Rscript power_simulation/power_analysis_weighted_four_gen_oracle.R $s 1 100 | tail -1
                     done | awk '{s += $4} END {print s/1000}'
                     
             OUTPUT  0.829          (Table 3, gamma = 1.00, column 3x3 (Opt))
         The same loop with another script gives the other columns:
             power_analysis_weighted_four_gen_binary_outcome_1_2.R   0.133   (3x2 (V1))
             power_analysis_weighted_four_gen_binary_outcome_2_3.R   0.795   (3x2 (V2))
             power_analysis_weighted_four_gen_Fisher_outcome_1_2.R   0.074   (2x2 (V1))
             power_analysis_weighted_four_gen_Fisher_outcome_2_3.R   0.694   (2x2 (V2))
         For another row, replace the "1" after "$s" by that row's gamma.
         Seeds are fixed and the p-values exact, so the values are
         reproduced exactly. The Cross-Cut column is not produced by these
         scripts.
     For other data-generating settings, edit the data-generating parameters
     near the top of each script and the treatment margins in the three
     sample() calls; set beta to 0 for the settings with no average
     treatment effect. The ordinal test's scores equal the data-generating
     scores (w*, v*) in every run.

4. To reproduce the data analysis, run
       source("data_analysis/education_data_analysis.R")
   The script reads data_analysis/education.data.csv and requires no ECLS-K
   download. The block that re-derives the profile-likelihood scores with
   Gurobi through reticulate is off by default, because the resulting scores,
   (0, 0.197, 1.598) for treatment and (0, 0.668, 0.698) for outcome, are
   entered directly where they are used. To re-derive them, set
   rederive_scores <- TRUE in the script; this requires the optional
   packages listed under REQUIREMENTS. To start from the raw ECLS-K file
   instead, uncomment the extraction lines at the beginning of the script
   after downloading ECLSK_98_99_K8_CHILD_v1_0.csv as described under DATA
   AVAILABILITY.

5. To reproduce the size-control results, run
       source("power_simulation/size_control.R")
   The experiments draw tables with BiasedUrn::rMFNCHypergeo and set their
   seed internally.

6. To reproduce the counterexample, run
       source("power_simulation/Corner_Solution_Counter_Example.R")
   The results are printed to the console.

7. To run the feasibility benchmarks, run each script in /feasibility/, e.g.
       source("feasibility/ordinal_feasibility_3x3.R")
       source("feasibility/general_test_feasibility_3x3.R")
   The ordinal runs complete in seconds to minutes; the chi-square runs for
   larger N take hours and should be launched as a background job, e.g.
       Rscript feasibility/general_test_feasibility_3x3.R
   Each script appends to its CSV in /feasibility/ and skips chunks already
   checkpointed, so an interrupted run can be restarted.

8. To verify the code, run the scripts in /development_check/, e.g.
      source("development_check/sensitivityIxJ_cpp_check.R")
   All checks should report PASS. verify_same_names.R also loads the
   previous engine from legacy/.

9. Random seeds are fixed within the scripts. The power scripts set the seed
   from the replication index, so the same simulated tables are used across
   the five tests at each gamma; comparisons across tests therefore reflect
   the tests rather than Monte Carlo variation between data sets. Exact
   p-values involve no randomness and are reproduced to the printed
   precision on any platform.
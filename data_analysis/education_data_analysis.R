#################################################################################
#                                                                               #
#        This file include the codes of data cleaning and data analysis for     #
#        the paper "Exact Sensitivity Analysis of Observational Studies of      #
#        Contingency Tables", Section 8, data analysis for the results based on # 
#        three by three ordinal test, and Appendix F for the results based on   #
#        tests shown to be less powerful in Section 7. The scores are chosen in #
#        two ways, see also discussion and the mentioned                        #
#        literature in Section 7 and Section 8 about various                    #
#        ways of choosing scores.                                               #
#                                                                               #
#                                                                               # 
#                                                                               #
#################################################################################




###############################################################################
#
#     Data Cleaning 
# 
##############################################################################

# Assume the user obtains the "ECLSK_98_99_K8_CHILD_v1_0.csv" downloaded from website, National Center 
# for Education Statistics, online codebook.

# data.t= read.csv("ECLSK_98_99_K8_CHILD_v1_0.csv")
# VartoKeep = c("P1PRIMPK","CHILDID","CREGION","GENDER","RACE","C1R4MSCL","WKBLACK",
# "WKHISP","P1HFAMIL","C1CMOTOR","R1_KAGE","P1HIG_1","W1INCCAT","WKINCOME","A2TXSPEN",
# "C1R4MPF","C1R4RPF","S2HVISIT","S2ORIENT","C7R4RPF","C7R4MPF")

# education.data = data.t[,VartoKeep]
# write.csv(x=education.data,file="education.data.csv")
## Here first check the distributions of variables are as the one from website:
## https://nces.ed.gov/datalab/onlinecodebook/session/codebook/f8b2403e-1679-45c0-96ff-5fbec50a799b
## Also adds description from the user manual:
## https://nces.ed.gov/ecls/data/eclsk_k8_manual_part1.pdf
## https://nces.ed.gov/ecls/data/eclsk_k8_manual_part2.pdf

# -----------------------------------
library(Rcpp)
sourceCpp("sensitivityIxJ.cpp")
source("sensitivityIxJ.R")



# -----------------------------------
# Read and clean the data
# -----------------------------------
## if the variable is $-1$ or $-9$ means the variable value is not provided by the subject, so we do >=0 to filter

# Read the CSV file
data.t <- read.csv("education.data.csv")

## Data Description from Manual page 1-1
## The ECLS-K focuses on children’s early school experiences beginning with kindergarten
## and ending with eighth grade. It is a multisource, multimethod study that includes
## interviews with parents, the collection of data from principals and teachers,
## and student records abstracts, as well as direct child assessments.
## In the eighth-grade data collection, a student paper-and-pencil questionnaire was added.
## The ECLS-K was developed under the sponsorship of the U.S. Department of Education,
## Institute of Education Sciences, National Center for Education Statistics (NCES).
## Westat conducted this study with assistance provided by Educational Testing Service (ETS) in Princeton, New Jersey.
## The ECLS-K followed a nationally representative cohort of children from kindergarten into middle school.
## The base-year data were collected in the fall and spring of the 1998–99 school year when the sampled children
## were in kindergarten. A total of 21,260 kindergartners throughout the nation participated.
## Two more waves of data were collected in the fall and spring of the 1999–2000 school year when most,
## but not all, of the base-year children were in first grade.3 The fall-first grade data collection was limited to a 30 percent
## subsample of schools.
## A fifth wave of data was collected in the spring of the 2001–02 school year when most,
## but not all, of the sampled children were in third grade.
## A sixth wave of data was collected in the spring of the 2003–04 school year when most, but not all,
## of the sampled children were in fifth grade.
## A seventh wave of data was collected in the spring of the 2006–07 school year when most,
## but not all, of the sampled children were in eighth grade.

## Exhibit 1-1                                     Date of Collection            Sample
## Fall-kindergarten                               Fall 1998                     Full sample
## Spring-kindergarten                             Spring 1999                   Full sample
## Fall-first grade                                Fall 1999                     30 percent subsample
## Spring-first grade                              Spring 2000                   Full sample plus freshening
## Spring-third grade                              Spring 2002                   Full sample
## Spring-fifth grade                              Spring 2004                   Full sample
## Spring-eighth grade                             Spring 2007                   Full sample


## Explanation of the measure of math ability - Manual 2-4 & Manual 3-11
## The kindergarten through eighth-grade mathematics proficiency levels include
## (1) Number and Shape—identifying some one-digit numerals, recognizing geometric shapes,
## and one-to-one counting up to 10 objects;
## (2) Relative Size—reading all one-digit numerals, counting beyond 10,
## recognizing a sequence of patterns, and using nonstandard units of length to compare the size of objects;
## (3) Ordinality and Sequence—reading two-digit numerals, recognizing the next
## number in a sequence, identifying the ordinal position of an object, and solving a simple word problem;
## (4) Addition and Subtraction—solving simple addition and subtraction problems;
## (5) Multiplication and Division—solving simple multiplication and division problems and recognizing more complex number patterns;
## (6) Place Value—demonstrating understanding of place value in integers to hundreds’ place;
## (7)Rate and Measurement—using knowledge of measurement and rate to solve word problems;
## (8) Fractions—solving problems using fractions; and (9)Area and Volume—solving word problems involving area and volume.
## No new mathematics proficiency level was added at the eighth grade because
## it was not warranted. Previously defined proficiency levels were sufficiently “difficult” to allow for the demonstration of growth in the higher proficiency levels at eighth grade.

## The proficiency levels were assumed to follow a Guttman model,
## that is, a child passing a particular skill level was expected to have mastered all lower levels;
## a failure should be consistent with nonmastery at higher levels. Only a very small percentage of
## children in kindergarten through eighth grade had response patterns that did not follow the Guttman model,
## that is, a failing score at a lower level followed by a pass on a more difficult item cluster.
## For the first six rounds of data collection, less than 7 percent of reading response patterns,
## and about 3 percent of mathematics assessment results, failed to follow the expected hierarchical pattern;
## Mastery of a proficiency level was defined as answering correctly at least three of the four questions in a cluster.
## This definition results in a very low probability of guessing enough right answers by chance,
## generally less than 2 percent. At least two incorrect or “don’t know” responses indicated lack of mastery of a cluster.


## Explanation of related variables - Manual 2-7
## Parent Interview
## The eighth-grade parent interview was conducted using a computer-assisted interview (CAI).
## The parent interview was conducted primarily in English, but provisions were made
## to interview parents who spoke other languages with bilingual English-Spanish interviewers or
## interpreters for other languages. Most of the interviews were conducted by telephone, but a small percentage
## (2.2 percent) were conducted in person.
## Data collection for the eighth-grade parent interview started in fall 2006.
## The parent interview lasted on average 46 minutes and contained approximately 300 questions concerning eighth-
## grade school experiences, parent characteristics, and child health.

## Assessment Format - Manual 3-2
## The format of the eighth-grade assessment was modified from that of prior
## rounds to accommodate administration differences for the older sample.
## In all previous rounds, an assessor presented the questions to the child and
## entered responses into a computer for each individually administered assessment. In spring-eighth grade, groups of ECLS-K sampled children who attended the same school
## were assessed in a single, proctored group administration.

## Variable Naming Convention - Manual 3-3
## The name and description for each variable in the tables begin with a “C,”
## indicating that it is a child variable, and a data collection round number: 1 (fall- kindergarten),
## 2 (spring-kindergarten), 3 (fall-first grade), 4 (spring-first grade), 5 (spring-third grade),
## 6 (spring-fifth grade), or 7 (spring-eighth grade).

## The main variable we use in this data analysis is "C1R4MPF"
## Stands for "C1 RC4 Math Highest Prof Lvl Mastered"
## R4 stands for, there are four question with right or wrong for each level
## a kid showed the mastery of a particular level if more than 3 questions out of
## 4 were answered correctly.
## Table 3-7 from manual shows that most kids are in level 1 to 3
##
## Level             weighted percentage
## Below Level 1     6
## Level 1           32
## Level 2           27
## Level 3           20
## Level 4           4
## Level 5 and above 0

## The distribution
# Category	Label            	                            UnW1	UnW1
# 0        	NON-MASTERY OF THE LOWEST PROFICIENCY LEVEL	 1,385	6.47
# 1	        NUMBER AND SHAPE	                           6,456	30.16
# 2	        RELATIVE SIZE	                               6,613	30.89
# 3	        ORDINALITY, SEQUENCE	                       3,016	14.09
# 4	        ADDITION/SUBTRACTION	                         614	2.87
# 5	        MULTIPLICATION/DIVISION	                        63	0.29
# 6	              PLACE VALUE	                              3	0.01
# 7	              RATE AND MEASUREMENT	                   0	0.00
# 8	              FRACTIONS	                                0	0.00
# 9	             AREA AND VOLUME	                          0	0.00
# -9	             NOT ASCERTAINED	                        561	2.62
# -1	         NOT APPLICABLE	                              415	1.94
# —		                                                     2,283	10.66
#                TOTAL	                               	21,409	100.00


## Family variable P1HFMIL - Manual 7-43
# Category	Label            	        UnW1	UnW1
#  1	     2 PARENTS PLUS SIBLINGS	  11,916	55.66
#  2	     2 PARENTS NO SIBLING	       1,744	8.15
#  3	     1 PARENT PLUS SIBLINGS	     2,837	13.25
#  4	     1 PARENT NO SIBLING	       1,252	5.85
#  5	       OTHER	                  348	  1.63
#  - 	                                3,312	15.47
#                           TOTAL		21,409	100.00



## P1 PRIMARY TYPE NONPARENTAL CARE PRE-K (P1PRIMPK)
##     Category	Label            	UnW1	UnW1
## 0	NO NON-PARENTAL CARE	      3,362	15.70
## 1	RELATIVE CARE, CHILD'S HOME	1,054	4.92
## 2	RELATIVE CARE, OTHER'S HOME	1,386	6.47
## 3	NON-REL CARE, CHILD'S HOME	332	  1.55
## 4	NON-REL CARE, OTHER HOME	  1,423	6.65
## 5	HEAD START PROGRAM	        1,701	7.95
## 6	CENTER-BASED PROGRAM	      7,645	35.71
## 7	2 OR MORE PROGRAMS	        700	  3.27
## 8	LOCATION VARIES	            208	  0.97
##-9	NOT ASCERTAINED	            286	 1.34


## CHILD COMPOSITE GENDER (GENDER)
##     Category	Label            	UnW1	UnW1
## 1	   MALE	                 10,950	51.15
## 2	  FEMALE	               10,446	48.79
## -9	  NOT ASCERTAINED	          13	0.06


## CHILD COMPOSITE RACE (RACE)
## Category	      Label                   	               UnW1	UnW1
## 1	       WHITE, NON-HISPANIC	                        11,788	55.06
## 2	    BLACK OR AFRICAN AMERICAN, NON-HISPANIC	        3,224	  15.06
## 3	    HISPANIC, RACE SPECIFIED	                      1,839	  8.59
## 4	     HISPANIC, RACE NOT SPECIFIED	                  1,987	  9.28
## 5	           ASIAN	                                  1,366	  6.38
## 6	  NATIVE HAWAIIAN, OTHER PACIFIC ISLANDER	           224	  1.05
## 7	  AMERICAN INDIAN OR ALASKA NATIVE	                 381	  1.78
## 8	     MORE THAN ONE RACE, NON HISPANIC	               549	  2.56
## -9	          NOT ASCERTAINED	                           51	    0.24
##                   TOTAL		                            21,409	100.00


## WK INCOME (IMPUTED) (WKINCOME)
## Category	   Min	 Max	          UnW1	   UnW1
## Continuous	0.00	1000000.00	    20,141	 94.08
##   —			                        1,268	   5.92
##TOTAL			                        21,409	100.00


################################################################################

## Verifying the distribution of variables after focusing on the levels we want


# -----------------------------------
# Create binary variables
# -----------------------------------
## FEMALE
# female variable is $1$ if girl, otherwise, boy, missing data, put NA
data.t$FEMALE <- ifelse(data.t$GENDER == 2, 1, ifelse(data.t$GENDER==1,0,NA))



# Create binary variable for two-parent families, two parent families, 1, single parent familes 0, no parent figure, NA
data.t$two.parent <- ifelse(data.t$P1HFAMIL == 1 | data.t$P1HFAMIL == 2, 1, ifelse(data.t$P1HFAMIL==3 | data.t$P1HFAMIL==4,0,NA))
data.t$two.parent <- as.factor(data.t$two.parent)

# -----------------------------------
# Create categorical race variable with consistent data types
# -----------------------------------
### RACE, only consider white, black, and combine two hispanic categories together,
## the other racial groups not the focus so put NA
data.t$race <- ifelse(data.t$RACE == 1, "White",
                      ifelse(data.t$RACE == 2, "Black", ifelse(data.t$RACE==3|data.t$RACE==4,"Hispanic",NA)))
data.t$race <- as.factor(data.t$race)

# ------------------------------------------------------------------------------
## combine some prek category and name it treatment.prek
data.t$treatment.prek = factor(ifelse(data.t$P1PRIMPK==6,"center-based",
                                      ifelse(data.t$P1PRIMPK==1|data.t$P1PRIMPK==2,"relative care",
                                             ifelse(data.t$P1PRIMPK==0,"no care",NA))), levels = c("no care","relative care","center-based"))

## focus on some math level only,
data.t$early.math.prof = ifelse(data.t$C1R4MPF==1,1,ifelse(data.t$C1R4MPF==2,2,ifelse(data.t$C1R4MPF==3,3,NA)))



## only keep the data with variables we want to discuss, also filter out kids without data
## in the last round, we don't talk about the missing mechanism of this longitudinal data set
analysis.data = data.t[,c("P1PRIMPK","GENDER","RACE","FEMALE","P1HFAMIL","WKINCOME","C1R4MPF",
                          "two.parent","race","treatment.prek","early.math.prof",
                          "C7R4MPF")]

### Examine the distribution is right - prek
table(analysis.data$P1PRIMPK)
table(analysis.data$treatment.prek)


## Examine the distribution, early math
table(analysis.data$C1R4MPF)
table(analysis.data$early.math.prof)


## examine the distribution - race
table(analysis.data$RACE)
table(analysis.data$race)


## examine the distribution - family
table(analysis.data$P1HFAMIL)
table(analysis.data$two.parent)




##############################################################################
# Create a new column indicating complete cases
analysis.data$sample <- complete.cases(analysis.data)


analysis.data.complete = subset(analysis.data, sample==TRUE & C7R4MPF>0)

## calculate the median of the income after retaining all the data we are going
## to use
summary(analysis.data.complete$WKINCOME) ## the median is 50000
analysis.data.complete$low.income= ifelse(analysis.data.complete$WKINCOME<=50000,1,0)
head(analysis.data.complete)
# -----------------------------------
# Define Functions
# -----------------------------------


RCT_probability_three_by_three = function(table){
  n11 = table[1,1]
  n12 = table[1,2]
  n21 = table[2,1]
  n22 = table[2,2]
  treatment_margins = rowSums(table)
  outcome_margins = colSums(table)
  n1 = as.numeric(treatment_margins[1])
  n2 = as.numeric(treatment_margins[2])
  n3 = as.numeric(treatment_margins[3])
  m1 = as.numeric(outcome_margins[1])
  m2 = as.numeric(outcome_margins[2])
  m3 = as.numeric(outcome_margins[3])
  N = n1+n2+n3
  ## here computes the probability of an observed table fixing the treatment and outcome margins
  denominator = choose(N,n1)*choose(N-n1,n2)
  numerator = choose(m1,n11)*choose(m1-n11,n21)*choose(m2,n12)*choose(m2-n12,n22)*choose(m3,n1-n11-n12)*choose(m3-(n1-n11-n12),n2-n21-n22)
  probability = numerator/denominator
  return(probability)
  
}




#-----------------------------------
# 
## Black girl low income single parent table
subtable.data = dplyr::filter(analysis.data.complete,FEMALE==1,race=="Black",two.parent==0,low.income==1)
Black.girl.subtable = with(subtable.data, table(treatment.prek, early.math.prof))


## Black boy low income single parent table 


subtable.data = dplyr::filter(analysis.data.complete,FEMALE==0,race=="Black",two.parent==0,low.income==1)
Black.boy.subtable = with(subtable.data, table(treatment.prek, early.math.prof))









#######################################################################################
#
#           Sensitivity analysis with some prior-chosen scores - Black Girl Table 
#
#
#######################################################################################




Gamma = c(1,1.5,2,2.5,3, 3.5,4,4.5)
gamma =log(Gamma)


alpha = c(0,0.25,1.5)
beta  = c(0,1.0,1.5)


## RCT p-value:0.005726717

## Gamma 1.5 p value: 0.01529282
prior.score.girl.Gamma.1.5.result = exact.score.sen.IxJ(obs.table = Black.girl.subtable, gamma=gamma[2], delta = c(0,1,1),
                                                        row = "treatment",treatment.scores=alpha, outcome.scores=beta,
                                                        verbose = FALSE)


## Gamma 2 p value: 0.02755185
prior.score.girl.Gamma.2.result = exact.score.sen.IxJ(obs.table = Black.girl.subtable, gamma=gamma[3], delta = c(0,1,1),
                                                      row = "treatment",treatment.scores=alpha, outcome.scores=beta,
                                                      verbose = FALSE)
## Gamma 2.5 p value: 0.04088577
prior.score.girl.Gamma.2.5.result = exact.score.sen.IxJ(obs.table = Black.girl.subtable, gamma=gamma[4], delta = c(0,1,1),
                                                        row = "treatment",treatment.scores=alpha, outcome.scores=beta,
                                                        verbose = FALSE)

## Gamma 3.0 p value:0.05430924
prior.score.girl.Gamma.3.result = exact.score.sen.IxJ(obs.table = Black.girl.subtable, gamma=gamma[5], delta= c(0,1,1),
                                                      row = "treatment",treatment.scores=alpha, outcome.scores=beta,
                                                      verbose = FALSE)


## Gamma 3.5 p value:0.06729827
prior.score.girl.Gamma.3.5.result = exact.score.sen.IxJ(obs.table = Black.girl.subtable, gamma=gamma[6], delta= c(0,1,1),
                                                        row = "treatment",treatment.scores=alpha, outcome.scores=beta,
                                                        verbose = FALSE)


## Gamma 4.0 p value: 0.07960623
prior.score.girl.Gamma.4.0.result = exact.score.sen.IxJ(obs.table = Black.girl.subtable, gamma=gamma[7], delta= c(0,1,1),
                                                        row = "treatment",treatment.scores=alpha, outcome.scores=beta,
                                                        verbose = FALSE)


## Gamma 4.5 p value: 0.09114105 
prior.score.girl.Gamma.4.5.result = exact.score.sen.IxJ(obs.table = Black.girl.subtable, gamma=gamma[8], delta= c(0,1,1),
                                                        row = "treatment",treatment.scores=alpha, outcome.scores=beta,
                                                        verbose = FALSE)







## second table is the Black boy low income single parent table


subtable.data = dplyr::filter(analysis.data.complete,FEMALE==0,race=="Black",two.parent==0,low.income==1)
Black.boy.subtable = with(subtable.data, table(treatment.prek, early.math.prof))



################################################################################
#  Sensitivity Analysis for Black Boy Table                                    #
################################################################################

Gamma = c(1,1.5,2,2.5,3, 3.5, 4.0, 4.5)
gamma =log(Gamma)


alpha = c(0,0.25,1.5)
beta  = c(0,1.0,1.5)


## RCT p-value:0.01250206

## Gamma 1.5 p value: 0.03220748
prior.score.boy.Gamma.1.5.result = exact.score.sen.IxJ(obs.table = Black.boy.subtable, gamma=gamma[2], delta = c(0,1,1),
                                                       row = "treatment",treatment.scores=alpha, outcome.scores=beta,
                                                       verbose = FALSE)


## Gamma 2 p value: 0.0563288
prior.score.boy.Gamma.2.result = exact.score.sen.IxJ(obs.table = Black.boy.subtable, gamma=gamma[3], delta = c(0,1,1),
                                                     row = "treatment",treatment.scores=alpha, outcome.scores=beta,
                                                     verbose = FALSE)
## Gamma 2.5 p value: 0.08167532
prior.score.boy.Gamma.2.5.result = exact.score.sen.IxJ(obs.table = Black.boy.subtable, gamma=gamma[4], delta = c(0,1,1),
                                                       row = "treatment",treatment.scores=alpha, outcome.scores=beta,
                                                       verbose = FALSE)

## Gamma 3.0 p value: 0.1064913
prior.score.boy.Gamma.3.result = exact.score.sen.IxJ(obs.table = Black.boy.subtable, gamma=gamma[5], delta = c(0,1,1),
                                                     row = "treatment",treatment.scores=alpha, outcome.scores=beta,
                                                     verbose = FALSE)


## Gamma 3.5 p value: 0.1299347
prior.score.boy.Gamma.3.5.result = exact.score.sen.IxJ(obs.table = Black.boy.subtable, gamma=gamma[6], delta = c(0,1,1),
                                                       row = "treatment",treatment.scores=alpha, outcome.scores=beta,
                                                       verbose = FALSE)

## Gamma 4.0 p value: 0.1516706
prior.score.boy.Gamma.4.0.result = exact.score.sen.IxJ(obs.table = Black.boy.subtable, gamma=gamma[7], delta = c(0,1,1),
                                                       row = "treatment",treatment.scores=alpha, outcome.scores=beta,
                                                       verbose = FALSE)
## Gamma 4.5 p value: 0.1716283
prior.score.boy.Gamma.4.5.result = exact.score.sen.IxJ(obs.table = Black.boy.subtable, gamma=gamma[8], delta = c(0,1,1),
                                                       row = "treatment",treatment.scores=alpha, outcome.scores=beta,
                                                       verbose = FALSE)


###############################################################################
#
#     multiple testing based on the prior chosen scores
#         
#
################################################################################


## the multiple testing procedure based on truncated p value method, cut-off 0.2
sensitivitymv::truncatedP(p=c(prior.score.girl.Gamma.1.5.result$rct.prob,prior.score.boy.Gamma.1.5.result$rct.prob),trunc = 0.2)
## the combined p value 0.0006390349

sensitivitymv::truncatedP(p=c(prior.score.girl.Gamma.1.5.result$max.prob,prior.score.boy.Gamma.1.5.result$max.prob),trunc = 0.2)
## the combined p value:0.00344635

sensitivitymv::truncatedP(p=c(prior.score.girl.Gamma.2.result$max.prob,prior.score.boy.Gamma.2.result$max.prob), trunc=0.2)
## the combined p value: 0.009077988

sensitivitymv::truncatedP(p=c(prior.score.girl.Gamma.2.5.result$max.prob,prior.score.boy.Gamma.2.5.result$max.prob), trunc=0.2)
## the combined p value: 0.0169743


sensitivitymv::truncatedP(p=c(prior.score.girl.Gamma.3.result$max.prob,prior.score.boy.Gamma.3.result$max.prob), trunc=0.2)
## the combined p value: 0.02622151



sensitivitymv::truncatedP(p=c(prior.score.girl.Gamma.3.5.result$max.prob,prior.score.boy.Gamma.3.5.result$max.prob), trunc=0.2)
## the combined p value: 0.03603094




sensitivitymv::truncatedP(p=c(prior.score.girl.Gamma.4.0.result$max.prob,prior.score.boy.Gamma.4.0.result$max.prob), trunc=0.2)
## the combined p value: 0.04585473




sensitivitymv::truncatedP(p=c(prior.score.girl.Gamma.4.5.result$max.prob,prior.score.boy.Gamma.4.5.result$max.prob), trunc=0.2)
## the combined p value:0.05535676 













############################################################################################################################
############################################################################################################################
############################################################################################################################
############################################################################################################################
## The following uses the data to construct test scores, with a profile likelihood method, it is required 
## to have Gurobi licence and library installed, also, Gurobi is easier to operate in python, so we use the 
## R library reticulate to run python script. The results in the comment are from gurobi v12.0.2, with reticulate 
## points to python 3.13.3, the user should also install gurobipy and numpy 

library(reticulate)
library(gurobi)





############################################################
#
#    Helper function to estimate the scores
#
###########################################################

fit_linear_by_linear_shared_scores_two_tables <- function(
    n1, n2, 
    alpha_init = NULL, 
    beta_init = NULL,
    max_iter = 100,
    tol = 1e-12,
    verbose = TRUE
) {
  py_run_string('
import numpy as np
import gurobipy as gp
from gurobipy import GRB
import time

def fit_alpha_given_beta_shared_scores(n1, n2, beta_fixed):
    n1 = np.asarray(n1)
    n2 = np.asarray(n2)
    beta_fixed = np.asarray(beta_fixed)

    if n1.shape != n2.shape:
        raise ValueError("The dimensions of n1 and n2 must be the same")

    I, J = n1.shape

    model = gp.Model("fit_alpha_given_beta")
    model.setParam("OutputFlag", 0)

    lambda0_1 = model.addVar(lb=-GRB.INFINITY, name="lambda0_1")
    lambda0_2 = model.addVar(lb=-GRB.INFINITY, name="lambda0_2")

    lambda_z1 = [model.addVar(lb=-GRB.INFINITY, name=f"lambda_z1_{i}") for i in range(1, I)]
    lambda_z2 = [model.addVar(lb=-GRB.INFINITY, name=f"lambda_z2_{i}") for i in range(1, I)]

    lambda_R1 = [model.addVar(lb=-GRB.INFINITY, name=f"lambda_R1_{j}") for j in range(1, J)]
    lambda_R2 = [model.addVar(lb=-GRB.INFINITY, name=f"lambda_R2_{j}") for j in range(1, J)]

    delta_alpha = [model.addVar(lb=0, name=f"delta_alpha_{i}") for i in range(1, I)]

    alpha = [0]
    for i in range(1, I):
        alpha_i = model.addVar(lb=0, name=f"alpha_{i}")
        model.addConstr(alpha_i == sum(delta_alpha[:i]), name=f"alpha_constr_{i}")
        alpha.append(alpha_i)

    eta1 = {}
    eta2 = {}
    mu1 = {}
    mu2 = {}
    ll_terms = []

    for i in range(I):
        for j in range(J):
            eta1[i, j] = model.addVar(lb=-700, ub=700, name=f"eta1_{i}_{j}")
            eta2[i, j] = model.addVar(lb=-700, ub=700, name=f"eta2_{i}_{j}")

            expr1 = lambda0_1
            if i > 0:
                expr1 += lambda_z1[i-1]
            if j > 0:
                expr1 += lambda_R1[j-1]
            expr1 += alpha[i] * beta_fixed[j]
            model.addConstr(eta1[i, j] == expr1, name=f"eta1_def_{i}_{j}")

            expr2 = lambda0_2
            if i > 0:
                expr2 += lambda_z2[i-1]
            if j > 0:
                expr2 += lambda_R2[j-1]
            expr2 += alpha[i] * beta_fixed[j]
            model.addConstr(eta2[i, j] == expr2, name=f"eta2_def_{i}_{j}")

            mu1[i, j] = model.addVar(name=f"mu1_{i}_{j}")
            mu2[i, j] = model.addVar(name=f"mu2_{i}_{j}")

            model.addGenConstrExp(eta1[i, j], mu1[i, j], name=f"exp1_{i}_{j}")
            model.addGenConstrExp(eta2[i, j], mu2[i, j], name=f"exp2_{i}_{j}")

            if n1[i, j] > 0:
                ll_term1 = model.addVar(lb=-GRB.INFINITY, name=f"ll_term1_{i}_{j}")
                model.addConstr(ll_term1 == n1[i, j] * eta1[i, j] - mu1[i, j], name=f"ll_constr1_{i}_{j}")
                ll_terms.append(ll_term1)
            else:
                ll_terms.append(-mu1[i, j])

            if n2[i, j] > 0:
                ll_term2 = model.addVar(lb=-GRB.INFINITY, name=f"ll_term2_{i}_{j}")
                model.addConstr(ll_term2 == n2[i, j] * eta2[i, j] - mu2[i, j], name=f"ll_constr2_{i}_{j}")
                ll_terms.append(ll_term2)
            else:
                ll_terms.append(-mu2[i, j])

    loglik = model.addVar(lb=-GRB.INFINITY, name="loglik")
    model.addConstr(loglik == gp.quicksum(ll_terms), name="loglik_constr")
    model.setObjective(loglik, GRB.MAXIMIZE)

    model.optimize()

    if model.status != GRB.OPTIMAL:
        print(f"Optimization failed. Status code: {model.status}")
        return None

    alpha_values = np.zeros(I)
    for i in range(I):
        if i == 0:
            alpha_values[i] = 0
        else:
            alpha_values[i] = alpha[i].X

    return {
        "alpha": alpha_values,
        "loglik": loglik.X
    }

def fit_beta_given_alpha_shared_scores(n1, n2, alpha_fixed):
    result = fit_alpha_given_beta_shared_scores(n1.T, n2.T, alpha_fixed)
    return {
        "beta": result["alpha"],
        "loglik": result["loglik"]
    }

def fit_linear_by_linear_shared_scores_two_tables_alternating_full(
    n1, n2,
    alpha_init=None,
    beta_init=None,
    max_iter=50,
    tol=1e-6,
    verbose=True
):
    n1 = np.asarray(n1)
    n2 = np.asarray(n2)

    if n1.shape != n2.shape:
        raise ValueError("The dimensions of n1 and n2 must be the same")

    I, J = n1.shape

    if alpha_init is None:
        alpha_init = np.linspace(0, 1, I)
    if beta_init is None:
        beta_init = np.linspace(0, 1, J)

    alpha = np.asarray(alpha_init)
    beta = np.asarray(beta_init)
    loglik_old = float("-inf")

    start_time = time.time()

    for iter in range(1, max_iter + 1):
        if verbose:
            print(f"Iteration {iter}")

        fit_alpha_result = fit_alpha_given_beta_shared_scores(n1, n2, beta)
        if fit_alpha_result is None:
            print("Failed to optimize row scores. Terminating.")
            break

        alpha = fit_alpha_result["alpha"]

        fit_beta_result = fit_beta_given_alpha_shared_scores(n1, n2, alpha)
        if fit_beta_result is None:
            print("Failed to optimize column scores. Terminating.")
            break

        beta = fit_beta_result["beta"]
        loglik_new = fit_beta_result["loglik"]

        if verbose:
            print(f"  Log-likelihood: {loglik_new}")

        if abs(loglik_new - loglik_old) < tol:
            if verbose:
                print("Converged")
            break

        loglik_old = loglik_new

    end_time = time.time()

    return {
        "alpha": alpha.tolist(),
        "beta": beta.tolist(),
        "loglik": float(loglik_new),
        "converged": iter < max_iter,
        "iterations": iter,
        "computation_time": end_time - start_time
    }
  ')
  
  n1_py <- as.matrix(n1)
  n2_py <- as.matrix(n2)
  
  if (is.null(alpha_init)) {
    alpha_init_py <- NULL
  } else {
    alpha_init_py <- as.numeric(alpha_init)
  }
  
  if (is.null(beta_init)) {
    beta_init_py <- NULL
  } else {
    beta_init_py <- as.numeric(beta_init)
  }
  
  result_py <- py$fit_linear_by_linear_shared_scores_two_tables_alternating_full(
    n1_py, n2_py, 
    alpha_init = alpha_init_py, 
    beta_init = beta_init_py,
    max_iter = as.integer(max_iter),
    tol = as.numeric(tol),
    verbose = as.logical(verbose)
  )
  
  result <- list(
    alpha = as.numeric(result_py$alpha),
    beta = as.numeric(result_py$beta),
    loglik = as.numeric(result_py$loglik),
    converged = as.logical(result_py$converged),
    iterations = as.integer(result_py$iterations),
    computation_time = as.numeric(result_py$computation_time)
  )
  
  return(result)
}



result <- fit_linear_by_linear_shared_scores_two_tables(
  Black.girl.subtable, Black.boy.subtable, verbose = TRUE
)

cat("\nFinal results:\n")
cat("Converged:", result$converged, "\n")
cat("Log-likelihood:", result$loglik, "\n")
cat("Row scores (alpha):", round(result$alpha, 3), "\n")
cat("Column scores (beta):", round(result$beta, 3), "\n")
cat(sprintf("Computation time: %.2f seconds\n", result$computation_time))


alpha = result$alpha
beta = result$beta


## alpha =  0.000 0.197 1.598
alpha = round(alpha,3)

## beta = 0.000 0.668 0.698
beta = round(beta,3)

##################################################################
#
#    Sensitivity Analysis based on profile likelihood method 
#
##################################################################


Gamma = c(1,1.5,2,2.5,3,3.5,4.0,4.5)
gamma =log(Gamma)


## RCT p-value: 0.00395886

## Gamma 1.5 p value: 0.01095246
profile.score.girl.Gamma.1.5.result = exact.score.sen.IxJ(obs.table = Black.girl.subtable, gamma=gamma[2], delta = c(0,1,1),
                                                          row = "treatment",treatment.scores=alpha, outcome.scores=beta,
                                                          verbose = FALSE)


## Gamma 2 p value: 0.02110251
profile.score.girl.Gamma.2.result = exact.score.sen.IxJ(obs.table = Black.girl.subtable, gamma=gamma[3], delta = c(0,1,1),
                                                        row = "treatment",treatment.scores=alpha, outcome.scores=beta,
                                                        verbose = FALSE)
## Gamma 2.5 p value: 0.03319327
profile.score.girl.Gamma.2.5.result = exact.score.sen.IxJ(obs.table = Black.girl.subtable, gamma=gamma[4], delta = c(0,1,1),
                                                          row = "treatment",treatment.scores=alpha, outcome.scores=beta,
                                                          verbose = FALSE)

## Gamma 3.0 p value: 0.04608662
profile.score.girl.Gamma.3.result = exact.score.sen.IxJ(obs.table = Black.girl.subtable, gamma=gamma[5], delta = c(0,1,1),
                                                        row = "treatment",treatment.scores=alpha, outcome.scores=beta,
                                                        verbose = FALSE)


## Gamma 3.5 p value: 0.05899391
profile.score.girl.Gamma.3.5.result = exact.score.sen.IxJ(obs.table = Black.girl.subtable, gamma=gamma[6], delta = c(0,1,1),
                                                          row = "treatment",treatment.scores=alpha, outcome.scores=beta,
                                                          verbose = FALSE)


## Gamma 4.0 p value: 0.07144845
profile.score.girl.Gamma.4.0.result = exact.score.sen.IxJ(obs.table = Black.girl.subtable, gamma=gamma[7], delta = c(0,1,1),
                                                          row = "treatment",treatment.scores=alpha, outcome.scores=beta,
                                                          verbose = FALSE)



## Gamma 4.5 p value: 0.08320904
profile.score.girl.Gamma.4.5.result = exact.score.sen.IxJ(obs.table = Black.girl.subtable, gamma=gamma[8], delta = c(0,1,1),
                                                          row = "treatment",treatment.scores=alpha, outcome.scores=beta,
                                                          verbose = FALSE)




########################################################################
#
#  Sensitivity Analysis for Black Boy with profile likelihood score 
#
########################################################################

################################################################################
#  Sensitivity Analysis for Black Boy Table                                    #
################################################################################

Gamma = c(1,1.5,2,2.5,3,3.5,4.0,4.5)
gamma =log(Gamma)



## RCT p-value: 0.01156032

## Gamma 1.5 p value: 0.02753322
profile.score.boy.Gamma.1.5.result = exact.score.sen.IxJ(obs.table = Black.boy.subtable, gamma=gamma[2], 
                                                         delta = c(0,1,1),
                                                         row = "treatment",treatment.scores=alpha, outcome.scores=beta,
                                                         verbose = FALSE)


## Gamma 2 p value: 0.0453901
profile.score.boy.Gamma.2.result = exact.score.sen.IxJ(obs.table = Black.boy.subtable, gamma=gamma[3], delta = c(0,1,1),
                                                       row = "treatment",treatment.scores=alpha, outcome.scores=beta,
                                                       verbose = FALSE)
## Gamma 2.5 p value: 0.06365891
profile.score.boy.Gamma.2.5.result = exact.score.sen.IxJ(obs.table = Black.boy.subtable, gamma=gamma[4], delta = c(0,1,1),
                                                         row = "treatment",treatment.scores=alpha, outcome.scores=beta,
                                                         verbose = FALSE)

## Gamma 3.0 p value: 0.08177947
profile.score.boy.Gamma.3.result = exact.score.sen.IxJ(obs.table = Black.boy.subtable, gamma=gamma[5],delta = c(0,1,1),
                                                       row = "treatment",treatment.scores=alpha, outcome.scores=beta,
                                                       verbose = FALSE)


## Gamma 3.5 p value: 0.09947659
profile.score.boy.Gamma.3.5.result = exact.score.sen.IxJ(obs.table = Black.boy.subtable, gamma=gamma[6],delta = c(0,1,1),
                                                         row = "treatment",treatment.scores=alpha, outcome.scores=beta,
                                                         verbose = FALSE)


## Gamma 4.0 p value: 0.1165848
profile.score.boy.Gamma.4.0.result = exact.score.sen.IxJ(obs.table = Black.boy.subtable, gamma=gamma[7],delta = c(0,1,1),
                                                         row = "treatment",treatment.scores=alpha, outcome.scores=beta,
                                                         verbose = FALSE)


## Gamma 4.5 p value: 0.1329985 
profile.score.boy.Gamma.4.5.result = exact.score.sen.IxJ(obs.table = Black.boy.subtable, gamma=gamma[8],delta = c(0,1,1),
                                                         row = "treatment",treatment.scores=alpha, outcome.scores=beta,
                                                         verbose = FALSE)




###############################################################################
#
#     multiple testing based on the profile likelihood method scores
#         
#
################################################################################


## the multiple testing procedure based on truncated p value method, cut-off 0.2
sensitivitymv::truncatedP(p=c(profile.score.girl.Gamma.1.5.result$rct.prob,profile.score.boy.Gamma.1.5.result$rct.prob),trunc = 0.2)
## the combined p value 0.0004289664

sensitivitymv::truncatedP(p=c(profile.score.girl.Gamma.1.5.result$max.prob,profile.score.boy.Gamma.1.5.result$max.prob),trunc = 0.2)
## the combined p value:0.002257957

sensitivitymv::truncatedP(p=c(profile.score.girl.Gamma.2.result$max.prob,profile.score.boy.Gamma.2.result$max.prob), trunc=0.2)
## the combined p value: 0.006065025

sensitivitymv::truncatedP(p=c(profile.score.girl.Gamma.2.5.result$max.prob,profile.score.boy.Gamma.2.5.result$max.prob), trunc=0.2)
## the combined p value: 0.01170786


sensitivitymv::truncatedP(p=c(profile.score.girl.Gamma.3.result$max.prob,profile.score.boy.Gamma.3.result$max.prob), trunc=0.2)
## the combined p value: 0.0187018


sensitivitymv::truncatedP(p=c(profile.score.girl.Gamma.3.5.result$max.prob,profile.score.boy.Gamma.3.5.result$max.prob), trunc=0.2)
## the combined p value: 0.02652144 



sensitivitymv::truncatedP(p=c(profile.score.girl.Gamma.4.0.result$max.prob,profile.score.boy.Gamma.4.0.result$max.prob), trunc=0.2)
## the combined p value: 0.03472728


sensitivitymv::truncatedP(p=c(profile.score.girl.Gamma.4.5.result$max.prob,profile.score.boy.Gamma.4.5.result$max.prob), trunc=0.2)
## the combined p value: 0.04299339























#####################################################################################################################################
#
#  Here conducts experiment to examine the p-value if we apply either Fisher's exact test, two
#  versions as in the main text, 
#  2x2, V1 combines the second and the third treatment levels into one, and combines the first and second outcome
#  2x2, V2 combines the second and the third treatment levels into one, and combines the second and the third outcome
#  cross-cut, keep only the corner [1,1], [1,3], [3,1], [3,3]
#
#
#
#
#
#
########################################################################################################




############################################
#
#  2x2, V1, Black girl testing
#
#
###########################################

Black.girl.2x2.v1.table = matrix(data=c(Black.girl.subtable[1,1]+Black.girl.subtable[1,2],Black.girl.subtable[1,3],
                                        Black.girl.subtable[2,1]+Black.girl.subtable[2,2]+Black.girl.subtable[3,1]+Black.girl.subtable[3,2],
                                        Black.girl.subtable[2,3]+Black.girl.subtable[3,3]),nrow=2,byrow = TRUE)


Gamma = c(1,1.5,2,2.5,3,3.5,4.0,4.5)
gamma =log(Gamma)


## RCT p-value: 0.2831957
## this RCT p-value can be found by setting gamma = 0, or use fisher.test
## we would obtain the same RCT p-value as 0.2832, and the 95% confidence interval
## covers "1"
fisher.test(Black.girl.2x2.v1.table, alternative = "greater")
## Gamma 1.5 p value: 0.4167711
Fisher.V1.girl.Gamma.1.5.result = exact.score.sen.IxJ(obs.table = Black.girl.2x2.v1.table, gamma=gamma[2],delta = c(0,1),
                                                      row = "treatment",treatment.scores=c(0,1), outcome.scores=c(0,1),
                                                      verbose = FALSE)


## Gamma 2 p value:0.5114623
Fisher.V1.girl.Gamma.2.result = exact.score.sen.IxJ(obs.table = Black.girl.2x2.v1.table, gamma=gamma[3],delta = c(0,1),
                                                    row = "treatment",treatment.scores=c(0,1), outcome.scores=c(0,1),
                                                    verbose = FALSE)
## Gamma 2.5 p value:0.5807063
Fisher.V1.girl.Gamma.2.5.result = exact.score.sen.IxJ(obs.table = Black.girl.2x2.v1.table, gamma=gamma[4], delta = c(0,1),
                                                      row = "treatment",treatment.scores=c(0,1), outcome.scores=c(0,1),
                                                      verbose = FALSE)

## Gamma 3.0 p value:0.6331651
Fisher.V1.girl.Gamma.3.result = exact.score.sen.IxJ(obs.table = Black.girl.2x2.v1.table, gamma=gamma[5],delta = c(0,1),
                                                    row = "treatment",treatment.scores=c(0,1), outcome.scores=c(0,1),
                                                    verbose = FALSE)

## Gamma 3.5 p value: 0.6741471
Fisher.V1.girl.Gamma.3.5.result = exact.score.sen.IxJ(obs.table = Black.girl.2x2.v1.table, gamma=gamma[6],delta = c(0,1),
                                                      row = "treatment",treatment.scores=c(0,1), outcome.scores=c(0,1),
                                                      verbose = FALSE)


## Gamma 4.0 p value: 0.7069909
Fisher.V1.girl.Gamma.4.0.result = exact.score.sen.IxJ(obs.table = Black.girl.2x2.v1.table, gamma=gamma[7],delta = c(0,1),
                                                      row = "treatment",treatment.scores=c(0,1), outcome.scores=c(0,1),
                                                      verbose = FALSE)


## Gamma 4.5 p value: 0.7338751
Fisher.V1.girl.Gamma.4.5.result = exact.score.sen.IxJ(obs.table = Black.girl.2x2.v1.table, gamma=gamma[8],delta = c(0,1),
                                                      row = "treatment",treatment.scores=c(0,1), outcome.scores=c(0,1),
                                                      verbose = FALSE)


############################################
#
#  2x2, V2, Black girl testing
#
#
###########################################


Black.girl.2x2.v2.table = matrix(data=c(Black.girl.subtable[1,1], Black.girl.subtable[1,2]+Black.girl.subtable[1,3],
                                        Black.girl.subtable[2,1]+Black.girl.subtable[3,1],
                                        Black.girl.subtable[2,2]+Black.girl.subtable[3,2]+Black.girl.subtable[2,3]+Black.girl.subtable[3,3]),
                                 nrow=2,byrow = TRUE)




Gamma = c(1,1.5,2,2.5,3,3.5,4.0,4.5)
gamma =log(Gamma)


## RCT p-value:0.01093 
## this RCT p-value can be found by setting gamma = 0, or use fisher.test
## we would obtain the same RCT p-value as 0.2832, and the 95% confidence interval
## covers "1"
fisher.test(Black.girl.2x2.v2.table, alternative = "greater")
## Gamma 1.5 p value: 0.0572564
Fisher.V2.girl.Gamma.1.5.result = exact.score.sen.IxJ(obs.table = Black.girl.2x2.v2.table, gamma=gamma[2],delta = c(0,1),
                                                      row = "treatment",treatment.scores=c(0,1), outcome.scores=c(0,1),
                                                      verbose = FALSE)


## Gamma 2 p value:0.1396436
Fisher.V2.girl.Gamma.2.result = exact.score.sen.IxJ(obs.table = Black.girl.2x2.v2.table, gamma=gamma[3],delta = c(0,1),
                                                    row = "treatment",treatment.scores=c(0,1), outcome.scores=c(0,1),
                                                    verbose = FALSE)
## Gamma 2.5 p value: 0.2403418
Fisher.V2.girl.Gamma.2.5.result = exact.score.sen.IxJ(obs.table = Black.girl.2x2.v2.table, gamma=gamma[4], delta = c(0,1),
                                                      row = "treatment",treatment.scores=c(0,1), outcome.scores=c(0,1),
                                                      verbose = FALSE)

## Gamma 3.0 p value:0.3432144
Fisher.V2.girl.Gamma.3.result = exact.score.sen.IxJ(obs.table = Black.girl.2x2.v2.table, gamma=gamma[5],delta = c(0,1),
                                                    row = "treatment",treatment.scores=c(0,1), outcome.scores=c(0,1),
                                                    verbose = FALSE)


## Gamma 3.5 p value: 0.4389258
Fisher.V2.girl.Gamma.3.5.result = exact.score.sen.IxJ(obs.table = Black.girl.2x2.v2.table, gamma=gamma[6],delta = c(0,1),
                                                      row = "treatment",treatment.scores=c(0,1), outcome.scores=c(0,1),
                                                      verbose = FALSE)


## Gamma 4.0 p value: 0.5235199 
Fisher.V2.girl.Gamma.4.0.result = exact.score.sen.IxJ(obs.table = Black.girl.2x2.v2.table, gamma=gamma[7],delta = c(0,1),
                                                      row = "treatment",treatment.scores=c(0,1), outcome.scores=c(0,1),
                                                      verbose = FALSE)

## Gamma 4.5 p value:0.5961697  
Fisher.V2.girl.Gamma.4.5.result = exact.score.sen.IxJ(obs.table = Black.girl.2x2.v2.table, gamma=gamma[8],delta = c(0,1),
                                                      row = "treatment",treatment.scores=c(0,1), outcome.scores=c(0,1),
                                                      verbose = FALSE)





#################################
#
# 2x2, V1, Black boy testing 
#
##################################



Black.boy.2x2.v1.table = matrix(data=c(Black.boy.subtable[1,1]+Black.boy.subtable[1,2],Black.boy.subtable[1,3],
                                       Black.boy.subtable[2,1]+Black.boy.subtable[2,2]+Black.boy.subtable[3,1]+Black.boy.subtable[3,2],
                                       Black.boy.subtable[2,3]+Black.boy.subtable[3,3]),nrow=2,byrow = TRUE)


Gamma = c(1,1.5,2,2.5,3,3.5,4.0,4.5)
gamma =log(Gamma)


## RCT p-value: 0.465673
## this RCT p-value can be found by setting gamma = 0, or use fisher.test
## we would obtain the same RCT p-value as 0.2832, and the 95% confidence interval
## covers "1"
fisher.test(Black.boy.2x2.v1.table, alternative = "greater")
## Gamma 1.5 p value:0.6454492 
Fisher.V1.boy.Gamma.1.5.result = exact.score.sen.IxJ(obs.table = Black.boy.2x2.v1.table, gamma=gamma[2],delta = c(0,1),
                                                     row = "treatment",treatment.scores=c(0,1), outcome.scores=c(0,1),
                                                     verbose = FALSE)


## Gamma 2 p value: 0.7510768
Fisher.V1.boy.Gamma.2.result = exact.score.sen.IxJ(obs.table = Black.boy.2x2.v1.table, gamma=gamma[3],delta = c(0,1),
                                                   row = "treatment",treatment.scores=c(0,1), outcome.scores=c(0,1),
                                                   verbose = FALSE)
## Gamma 2.5 p value: 0.8166104
Fisher.V1.boy.Gamma.2.5.result = exact.score.sen.IxJ(obs.table = Black.boy.2x2.v1.table, gamma=gamma[4], delta = c(0,1),
                                                     row = "treatment",treatment.scores=c(0,1), outcome.scores=c(0,1),
                                                     verbose = FALSE)

## Gamma 3.0 p value: 0.8596273
Fisher.V1.boy.Gamma.3.result = exact.score.sen.IxJ(obs.table = Black.boy.2x2.v1.table, gamma=gamma[5],delta = c(0,1),
                                                   row = "treatment",treatment.scores=c(0,1), outcome.scores=c(0,1),
                                                   verbose = FALSE)



## Gamma 3.5 p value: 0.8892446
Fisher.V1.boy.Gamma.3.5.result = exact.score.sen.IxJ(obs.table = Black.boy.2x2.v1.table, gamma=gamma[6],delta = c(0,1),
                                                     row = "treatment",treatment.scores=c(0,1), outcome.scores=c(0,1),
                                                     verbose = FALSE)


## Gamma 4.0 p value: 0.9104518 
Fisher.V1.boy.Gamma.4.0.result = exact.score.sen.IxJ(obs.table = Black.boy.2x2.v1.table, gamma=gamma[7],delta = c(0,1),
                                                     row = "treatment",treatment.scores=c(0,1), outcome.scores=c(0,1),
                                                     verbose = FALSE)

## Gamma 4.5 p value: 0.9261353  
Fisher.V1.boy.Gamma.4.5.result = exact.score.sen.IxJ(obs.table = Black.boy.2x2.v1.table, gamma=gamma[8],delta = c(0,1),
                                                     row = "treatment",treatment.scores=c(0,1), outcome.scores=c(0,1),
                                                     verbose = FALSE)








#################################
#
# 2x2, V2, Black boy testing 
#
##################################


## Since some results give pvalues so large to believe, we also use the Fisher's exact test from rbounds which implements Rosenbaum (2002) section 4.4,
## to verify our results, interested readers may find the package "rbounds" and the function FisherSens, we put the running in the comments





Black.boy.2x2.v2.table = matrix(data=c(Black.boy.subtable[1,1], Black.boy.subtable[1,2]+Black.boy.subtable[1,3],
                                       Black.boy.subtable[2,1]+Black.boy.subtable[3,1],
                                       Black.boy.subtable[2,2]+Black.boy.subtable[3,2]+Black.boy.subtable[2,3]+Black.boy.subtable[3,3]),nrow=2,byrow = TRUE)


Gamma = c(1,1.5,2,2.5,3,3.5,4.0,4.5)
gamma =log(Gamma)


#### use Rosenbaum and Keele's function in rbound 
## install.packages(rbounds)
## library(rbounds)
## FisherSens(totalN = sum(Black.boy.2x2.v2.table), treatedN = colSums(Black.boy.2x2.v2.table)[1], totalSuccesses = rowSums(Black.boy.2x2.v2.table)[1], treatedSuccesses = Black.boy.2x2.v2.table[1,1],Gammas = Gamma)

## RCT p-value:0.6021 
## this RCT p-value can be found by setting gamma = 0, or use fisher.test
## we would obtain the same RCT p-value as 0.2832, and the 95% confidence interval
## covers "1"
fisher.test(Black.boy.2x2.v2.table, alternative = "greater")
## Gamma 1.5 p value:0.85878 
Fisher.V2.boy.Gamma.1.5.result = exact.score.sen.IxJ(obs.table = Black.boy.2x2.v2.table, gamma=gamma[2],delta = c(0,1),
                                                     row = "treatment",treatment.scores=c(0,1), outcome.scores=c(0,1),
                                                     verbose = FALSE)


## Gamma 2 p value: 0.9502186
Fisher.V2.boy.Gamma.2.result = exact.score.sen.IxJ(obs.table = Black.boy.2x2.v2.table, gamma=gamma[3],delta = c(0,1),
                                                   row = "treatment",treatment.scores=c(0,1), outcome.scores=c(0,1),
                                                   verbose = FALSE)
## Gamma 2.5 p value: 0.9814398
Fisher.V2.boy.Gamma.2.5.result = exact.score.sen.IxJ(obs.table = Black.boy.2x2.v2.table, gamma=gamma[4], delta = c(0,1),
                                                     row = "treatment",treatment.scores=c(0,1), outcome.scores=c(0,1),
                                                     verbose = FALSE)

## Gamma 3.0 p value: 0.9925898
Fisher.V2.boy.Gamma.3.result = exact.score.sen.IxJ(obs.table = Black.boy.2x2.v2.table, gamma=gamma[5],delta = c(0,1),
                                                   row = "treatment",treatment.scores=c(0,1), outcome.scores=c(0,1),
                                                   verbose = FALSE)


## Gamma 3.5 p value: 0.9968379
Fisher.V2.boy.Gamma.3.5.result = exact.score.sen.IxJ(obs.table = Black.boy.2x2.v2.table, gamma=gamma[6],delta = c(0,1),
                                                     row = "treatment",treatment.scores=c(0,1), outcome.scores=c(0,1),
                                                     verbose = FALSE)



## Gamma 4.0 p value: 0.9985662
Fisher.V2.boy.Gamma.4.0.result = exact.score.sen.IxJ(obs.table = Black.boy.2x2.v2.table, gamma=gamma[7],delta = c(0,1),
                                                     row = "treatment",treatment.scores=c(0,1), outcome.scores=c(0,1),
                                                     verbose = FALSE)


## Gamma 4.5 p value:0.9993136 
Fisher.V2.boy.Gamma.4.5.result = exact.score.sen.IxJ(obs.table = Black.boy.2x2.v2.table, gamma=gamma[8],delta = c(0,1),
                                                     row = "treatment",treatment.scores=c(0,1), outcome.scores=c(0,1),
                                                     verbose = FALSE)










##############################################################################
#                                                                         
#     Examine the p value if we ignore the middle level and only keep 
#     treatment level 1 and 3, outcome level 1 and 3, in other words,
#     we only keep the treatment no care and center-based care, and outcome level
#     Number and shape, Ordinality and sequence
#
#
#############################################################################


Gamma = c(1,1.5,2,2.5,3,3.5,4.0,4.5)
gamma =log(Gamma)

Black.girl.cut.subtble = matrix(data=c(Black.girl.subtable[1,1],Black.girl.subtable[1,3], Black.girl.subtable[3,1],Black.girl.subtable[3,3]),nrow=2, byrow=2)


## RCT p-value: 0.146261 

## Gamma 1.5 p value: 0.2467872
crosscut.girl.Gamma.1.5.result = exact.score.sen.IxJ(obs.table = Black.girl.cut.subtble, gamma=gamma[2],delta = c(0,1),
                                                     row = "treatment",treatment.scores=c(0,1), outcome.scores=c(0,1),
                                                     verbose = FALSE)


## Gamma 2 p value: 0.3317915
crosscut.girl.Gamma.2.result = exact.score.sen.IxJ(obs.table = Black.girl.cut.subtble, gamma=gamma[3],delta = c(0,1),
                                                   row = "treatment",treatment.scores=c(0,1), outcome.scores=c(0,1),
                                                   verbose = FALSE)
## Gamma 2.5 p value: 0.401814
crosscut.girl.Gamma.2.5.result = exact.score.sen.IxJ(obs.table = Black.girl.cut.subtble, gamma=gamma[4], delta = c(0,1),
                                                     row = "treatment",treatment.scores=c(0,1), outcome.scores=c(0,1),
                                                     verbose = FALSE)

## Gamma 3.0 p value: 0.459599
crosscut.girl.Gamma.3.result = exact.score.sen.IxJ(obs.table = Black.girl.cut.subtble, gamma=gamma[5],delta = c(0,1),
                                                   row = "treatment",treatment.scores=c(0,1), outcome.scores=c(0,1),
                                                   verbose = FALSE)


## Gamma 3.5 p value:0.507742
crosscut.girl.Gamma.3.5.result = exact.score.sen.IxJ(obs.table = Black.girl.cut.subtble, gamma=gamma[6],delta = c(0,1),
                                                     row = "treatment",treatment.scores=c(0,1), outcome.scores=c(0,1),
                                                     verbose = FALSE)

## Gamma 4.0 p value: 0.5483093
crosscut.girl.Gamma.4.0.result = exact.score.sen.IxJ(obs.table = Black.girl.cut.subtble, gamma=gamma[7],delta = c(0,1),
                                                     row = "treatment",treatment.scores=c(0,1), outcome.scores=c(0,1),
                                                     verbose = FALSE)


## Gamma 4.5 p value: 0.5828778
crosscut.girl.Gamma.4.5.result = exact.score.sen.IxJ(obs.table = Black.girl.cut.subtble, gamma=gamma[8],delta = c(0,1),
                                                     row = "treatment",treatment.scores=c(0,1), outcome.scores=c(0,1),
                                                     verbose = FALSE)





################################################################################
#
#
#  Cross cut test for Black boy and the sensitivity analysis
#
################################################################################




Gamma = c(1,1.5,2,2.5,3,3.5,4.0,4.5)
gamma =log(Gamma)

Black.boy.cut.subtble = matrix(data=c(Black.boy.subtable[1,1],Black.boy.subtable[1,3], Black.boy.subtable[3,1],Black.boy.subtable[3,3]),nrow=2, byrow=2)


## RCT p-value: 0.309877 

## Gamma 1.5 p value: 0.48123
crosscut.boy.Gamma.1.5.result = exact.score.sen.IxJ(obs.table = Black.boy.cut.subtble, gamma=gamma[2], delta = c(0,1),
                                                    row = "treatment",treatment.scores=c(0,1), outcome.scores=c(0,1),
                                                    verbose = FALSE)


## Gamma 2 p value: 0.6021375
crosscut.boy.Gamma.2.result = exact.score.sen.IxJ(obs.table = Black.boy.cut.subtble, gamma=gamma[3], delta = c(0,1),
                                                  row = "treatment",treatment.scores=c(0,1), outcome.scores=c(0,1),
                                                  verbose = FALSE)
## Gamma 2.5 p value: 0.6872658
crosscut.boy.Gamma.2.5.result = exact.score.sen.IxJ(obs.table = Black.boy.cut.subtble, gamma=gamma[4], delta = c(0,1),
                                                    row = "treatment",treatment.scores=c(0,1), outcome.scores=c(0,1),
                                                    verbose = FALSE)

## Gamma 3.0 p value: 0.7485355
crosscut.boy.Gamma.3.result = exact.score.sen.IxJ(obs.table = Black.boy.cut.subtble, gamma=gamma[5], delta = c(0,1),
                                                  row = "treatment",treatment.scores=c(0,1), outcome.scores=c(0,1),
                                                  verbose = FALSE)



## Gamma 3.5 p value: 0.7937769 
crosscut.boy.Gamma.3.5.result = exact.score.sen.IxJ(obs.table = Black.boy.cut.subtble, gamma=gamma[6], delta = c(0,1),
                                                  row = "treatment",treatment.scores=c(0,1), outcome.scores=c(0,1),
                                                  verbose = FALSE)





## Gamma 4.0 p value: 0.828001 
crosscut.boy.Gamma.4.0.result = exact.score.sen.IxJ(obs.table = Black.boy.cut.subtble, gamma=gamma[7], delta = c(0,1),
                                                    row = "treatment",treatment.scores=c(0,1), outcome.scores=c(0,1),
                                                    verbose = FALSE)




## Gamma 4.5 p value:  0.854457
crosscut.boy.Gamma.4.5.result = exact.score.sen.IxJ(obs.table = Black.boy.cut.subtble, gamma=gamma[8], delta = c(0,1),
                                                    row = "treatment",treatment.scores=c(0,1), outcome.scores=c(0,1),
                                                    verbose = FALSE)












#######################################################################################
#
# Here considers only binarizing the outcome, using either the prior treatment scores 
# or the treatment scores we obtain by profile likelihood method, and consider 
# the 3x2 V1, and 3x2 V2 test
#
#######################################################################################

Black.girl.3x2.V1.table = matrix(data = c(Black.girl.subtable[1,1]+Black.girl.subtable[1,2], Black.girl.subtable[1,3],
                                          Black.girl.subtable[2,1]+Black.girl.subtable[2,2], Black.girl.subtable[2,3],
                                          Black.girl.subtable[3,1]+Black.girl.subtable[3,2], Black.girl.subtable[3,3]),nrow=3,byrow = TRUE)

Black.girl.3x2.V2.table = matrix(data = c(Black.girl.subtable[1,1], Black.girl.subtable[1,2]+Black.girl.subtable[1,3],
                                          Black.girl.subtable[2,1], Black.girl.subtable[2,2]+Black.girl.subtable[2,3],
                                          Black.girl.subtable[3,1], Black.girl.subtable[3,2]+Black.girl.subtable[3,3]),nrow=3,byrow = TRUE)


Black.boy.3x2.V1.table = matrix(data = c(Black.boy.subtable[1,1]+Black.boy.subtable[1,2], Black.boy.subtable[1,3],
                                         Black.boy.subtable[2,1]+Black.boy.subtable[2,2], Black.boy.subtable[2,3],
                                         Black.boy.subtable[3,1]+Black.boy.subtable[3,2], Black.boy.subtable[3,3]),nrow=3,byrow = TRUE)

Black.boy.3x2.V2.table = matrix(data = c(Black.boy.subtable[1,1], Black.boy.subtable[1,2]+Black.boy.subtable[1,3],
                                         Black.boy.subtable[2,1], Black.boy.subtable[2,2]+Black.boy.subtable[2,3],
                                         Black.boy.subtable[3,1], Black.boy.subtable[3,2]+Black.boy.subtable[3,3]),nrow=3,byrow = TRUE)





#########################################
# prior treatment scores, 3x2, V1, girls
########################################


Gamma = c(1,1.5,2,2.5,3, 3.5,4,4.5)
gamma =log(Gamma)


## this is the treatment scores we come up with, non-data-driven 
alpha = c(0,0.25,1.5)


## RCT p-value:0.287194


## Gamma 1.5 pvalue:0.3719068
prior.three.by.two.V1.girl.Gamma.1.5.result = exact.score.sen.IxJ(obs.table = Black.girl.3x2.V1.table, gamma=gamma[2],delta = c(0,1,1),
                                                                  row = "treatment",treatment.scores=alpha, outcome.scores=c(0,1),
                                                                  verbose = FALSE)


## Gamma 2 p value: 0.4273817
prior.three.by.two.V1.girl.Gamma.2.result = exact.score.sen.IxJ(obs.table = Black.girl.3x2.V1.table, gamma=gamma[3],delta = c(0,1,1),
                                                                row = "treatment",treatment.scores=alpha, outcome.scores=c(0,1),
                                                                verbose = FALSE)
## Gamma 2.5 p value:0.46624 
prior.three.by.two.V1.girl.Gamma.2.5.result = exact.score.sen.IxJ(obs.table = Black.girl.3x2.V1.table, gamma=gamma[4],delta = c(0,1,1),
                                                                  row = "treatment",treatment.scores=alpha, outcome.scores=c(0,1),
                                                                  verbose = FALSE)

## Gamma 3.0 p value: 0.4948949
prior.three.by.two.V1.girl.Gamma.3.0.result = exact.score.sen.IxJ(obs.table = Black.girl.3x2.V1.table, gamma=gamma[5],delta = c(0,1,1),
                                                                  row = "treatment",treatment.scores=alpha, outcome.scores=c(0,1),
                                                                  verbose = FALSE)



## Gamma 3.5 p value:0.5168703 
prior.three.by.two.V1.girl.Gamma.3.5.result = exact.score.sen.IxJ(obs.table = Black.girl.3x2.V1.table, gamma=gamma[6],delta = c(0,1,1),
                                                                  row = "treatment",treatment.scores=alpha, outcome.scores=c(0,1),
                                                                  verbose = FALSE)


## Gamma 4.0 p value: 0.5342459  
prior.three.by.two.V1.girl.Gamma.4.0.result = exact.score.sen.IxJ(obs.table = Black.girl.3x2.V1.table, gamma=gamma[7],delta = c(0,1,1),
                                                                  row = "treatment",treatment.scores=alpha, outcome.scores=c(0,1),
                                                                  verbose = FALSE)

## Gamma 4.5 p value: 0.5483236 
prior.three.by.two.V1.girl.Gamma.4.5.result = exact.score.sen.IxJ(obs.table = Black.girl.3x2.V1.table, gamma=gamma[8],delta = c(0,1,1),
                                                                  row = "treatment",treatment.scores=alpha, outcome.scores=c(0,1),
                                                                  verbose = FALSE)



#########################################
# prior treatment scores, 3x2, V2, girls
########################################


Gamma = c(1,1.5,2,2.5,3, 3.5,4,4.5)
gamma =log(Gamma)


## this is the treatment scores we come up with, non-data-driven 
alpha = c(0,0.25,1.5)


## RCT p-value: 0.004103672


## Gamma 1.5 pvalue:0.01185651
prior.three.by.two.V2.girl.Gamma.1.5.result = exact.score.sen.IxJ(obs.table = Black.girl.3x2.V2.table, gamma=gamma[2],delta = c(0,1,1),
                                                                  row = "treatment",treatment.scores=alpha, outcome.scores=c(0,1),
                                                                  verbose = FALSE)


## Gamma 2 p value: 0.02315253
prior.three.by.two.V2.girl.Gamma.2.result = exact.score.sen.IxJ(obs.table = Black.girl.3x2.V2.table, gamma=gamma[3],delta = c(0,1,1),
                                                                row = "treatment",treatment.scores=alpha, outcome.scores=c(0,1),
                                                                verbose = FALSE)
## Gamma 2.5 p value: 0.03642567
prior.three.by.two.V2.girl.Gamma.2.5.result = exact.score.sen.IxJ(obs.table = Black.girl.3x2.V2.table, gamma=gamma[4],delta = c(0,1,1),
                                                                  row = "treatment",treatment.scores=alpha, outcome.scores=c(0,1),
                                                                  verbose = FALSE)

## Gamma 3.0 p value: 0.05030819 
prior.three.by.two.V2.girl.Gamma.3.0.result = exact.score.sen.IxJ(obs.table = Black.girl.3x2.V2.table, gamma=gamma[5],delta = c(0,1,1),
                                                                  row = "treatment",treatment.scores=alpha, outcome.scores=c(0,1),
                                                                  verbose = FALSE)



## Gamma 3.5 p value: 0.06393221
prior.three.by.two.V2.girl.Gamma.3.5.result = exact.score.sen.IxJ(obs.table = Black.girl.3x2.V2.table, gamma=gamma[6],delta = c(0,1,1),
                                                                  row = "treatment",treatment.scores=alpha, outcome.scores=c(0,1),
                                                                  verbose = FALSE)


## Gamma 4.0 p value: 0.07683972
prior.three.by.two.V2.girl.Gamma.4.0.result = exact.score.sen.IxJ(obs.table = Black.girl.3x2.V2.table, gamma=gamma[7],delta = c(0,1,1),
                                                                  row = "treatment",treatment.scores=alpha, outcome.scores=c(0,1),
                                                                  verbose = FALSE)

## Gamma 4.5 p value: 0.08883371 
prior.three.by.two.V2.girl.Gamma.4.5.result = exact.score.sen.IxJ(obs.table = Black.girl.3x2.V2.table, gamma=gamma[8],delta = c(0,1,1),
                                                                  row = "treatment",treatment.scores=alpha, outcome.scores=c(0,1),
                                                                  verbose = FALSE)











#################################################################
#
#     profile treatment score 3x2 v1 for the girls 
#
#################################################################

Gamma = c(1,1.5,2,2.5,3, 3.5,4,4.5)
gamma =log(Gamma)

alpha =  c(0.000, 0.197, 1.598)



## RCT p-value:0.287194


## Gamma 1.5 pvalue:0.3719068
profile.three.by.two.V1.girl.Gamma.1.5.result = exact.score.sen.IxJ(obs.table = Black.girl.3x2.V1.table, gamma=gamma[2],delta = c(0,1,1),
                                                                    row = "treatment",treatment.scores=alpha, outcome.scores=c(0,1),
                                                                    verbose = FALSE)


## Gamma 2 p value: 0.4273817
profile.three.by.two.V1.girl.Gamma.2.result = exact.score.sen.IxJ(obs.table = Black.girl.3x2.V1.table, gamma=gamma[3],delta = c(0,1,1),
                                                                  row = "treatment",treatment.scores=alpha, outcome.scores=c(0,1),
                                                                  verbose = FALSE)
## Gamma 2.5 p value:0.46624 
profile.three.by.two.V1.girl.Gamma.2.5.result = exact.score.sen.IxJ(obs.table = Black.girl.3x2.V1.table, gamma=gamma[4],delta = c(0,1,1),
                                                                    row = "treatment",treatment.scores=alpha, outcome.scores=c(0,1),
                                                                    verbose = FALSE)

## Gamma 3.0 p value: 0.4948949 
profile.three.by.two.V1.girl.Gamma.3.0.result = exact.score.sen.IxJ(obs.table = Black.girl.3x2.V1.table, gamma=gamma[5],delta = c(0,1,1),
                                                                    row = "treatment",treatment.scores=alpha, outcome.scores=c(0,1),
                                                                    verbose = FALSE)



## Gamma 3.5 p value: 0.5168703
profile.three.by.two.V1.girl.Gamma.3.5.result = exact.score.sen.IxJ(obs.table = Black.girl.3x2.V1.table, gamma=gamma[6],delta = c(0,1,1),
                                                                    row = "treatment",treatment.scores=alpha, outcome.scores=c(0,1),
                                                                    verbose = FALSE)


## Gamma 4.0 p value: 0.5342459  
profile.three.by.two.V1.girl.Gamma.4.0.result = exact.score.sen.IxJ(obs.table = Black.girl.3x2.V1.table, gamma=gamma[7],delta = c(0,1,1),
                                                                    row = "treatment",treatment.scores=alpha, outcome.scores=c(0,1),
                                                                    verbose = FALSE)

## Gamma 4.5 p value: 0.5483236  
profile.three.by.two.V1.girl.Gamma.4.5.result = exact.score.sen.IxJ(obs.table = Black.girl.3x2.V1.table, gamma=gamma[8],delta = c(0,1,1),
                                                                    row = "treatment",treatment.scores=alpha, outcome.scores=c(0,1),
                                                                    verbose = FALSE)





#################################################################
#
#     profile treatment score 3x2 v2 for the girls 
#
#################################################################

Gamma = c(1,1.5,2,2.5,3, 3.5,4,4.5)
gamma =log(Gamma)

alpha =  c(0.000, 0.197, 1.598)



## RCT p-value: 0.004153878


## Gamma 1.5 pvalue: 0.01187648
profile.three.by.two.V2.girl.Gamma.1.5.result = exact.score.sen.IxJ(obs.table = Black.girl.3x2.V2.table, gamma=gamma[2],delta = c(0,1,1),
                                                                    row = "treatment",treatment.scores=alpha, outcome.scores=c(0,1),
                                                                    verbose = FALSE)


## Gamma 2 p value: 0.02316021
profile.three.by.two.V2.girl.Gamma.2.result = exact.score.sen.IxJ(obs.table = Black.girl.3x2.V2.table, gamma=gamma[3],delta = c(0,1,1),
                                                                  row = "treatment",treatment.scores=alpha, outcome.scores=c(0,1),
                                                                  verbose = FALSE)
## Gamma 2.5 p value: 0.03642878 
profile.three.by.two.V2.girl.Gamma.2.5.result = exact.score.sen.IxJ(obs.table = Black.girl.3x2.V2.table, gamma=gamma[4],delta = c(0,1,1),
                                                                    row = "treatment",treatment.scores=alpha, outcome.scores=c(0,1),
                                                                    verbose = FALSE)

## Gamma 3.0 p value: 0.05030953 
profile.three.by.two.V2.girl.Gamma.3.0.result = exact.score.sen.IxJ(obs.table = Black.girl.3x2.V2.table, gamma=gamma[5],delta = c(0,1,1),
                                                                    row = "treatment",treatment.scores=alpha, outcome.scores=c(0,1),
                                                                    verbose = FALSE)



## Gamma 3.5 p value: 0.06393283
profile.three.by.two.V2.girl.Gamma.3.5.result = exact.score.sen.IxJ(obs.table = Black.girl.3x2.V2.table, gamma=gamma[6],delta = c(0,1,1),
                                                                    row = "treatment",treatment.scores=alpha, outcome.scores=c(0,1),
                                                                    verbose = FALSE)


## Gamma 4.0 p value: 0.07684002  
profile.three.by.two.V2.girl.Gamma.4.0.result = exact.score.sen.IxJ(obs.table = Black.girl.3x2.V2.table, gamma=gamma[7],delta = c(0,1,1),
                                                                    row = "treatment",treatment.scores=alpha, outcome.scores=c(0,1),
                                                                    verbose = FALSE)

## Gamma 4.5 p value: 0.08883386  
profile.three.by.two.V2.girl.Gamma.4.5.result = exact.score.sen.IxJ(obs.table = Black.girl.3x2.V2.table, gamma=gamma[8],delta = c(0,1,1),
                                                                    row = "treatment",treatment.scores=alpha, outcome.scores=c(0,1),
                                                                    verbose = FALSE)





############################################################################
#
#    prior treatment scores, 3x2, V1, for boys 
#########################################################################


Gamma = c(1,1.5,2,2.5,3, 3.5,4,4.5)
gamma =log(Gamma)


## this is the treatment scores we come up with, non-data-driven 
alpha = c(0,0.25,1.5)


## RCT p-value:0.1881044


## Gamma 1.5 pvalue:0.2593915
prior.three.by.two.V1.boy.Gamma.1.5.result = exact.score.sen.IxJ(obs.table = Black.boy.3x2.V1.table, gamma=gamma[2],delta = c(0,1,1),
                                                                 row = "treatment",treatment.scores=alpha, outcome.scores=c(0,1),
                                                                 verbose = FALSE)


## Gamma 2 p value: 0.3044207
prior.three.by.two.V1.boy.Gamma.2.result = exact.score.sen.IxJ(obs.table = Black.boy.3x2.V1.table, gamma=gamma[3],delta = c(0,1,1),
                                                               row = "treatment",treatment.scores=alpha, outcome.scores=c(0,1),
                                                               verbose = FALSE)
## Gamma 2.5 p value:0.334582 
prior.three.by.two.V1.boy.Gamma.2.5.result = exact.score.sen.IxJ(obs.table = Black.boy.3x2.V1.table, gamma=gamma[4],delta = c(0,1,1),
                                                                 row = "treatment",treatment.scores=alpha, outcome.scores=c(0,1),
                                                                 verbose = FALSE)

## Gamma 3.0 p value: 0.3559135
prior.three.by.two.V1.boy.Gamma.3.0.result = exact.score.sen.IxJ(obs.table = Black.boy.3x2.V1.table, gamma=gamma[5],delta = c(0,1,1),
                                                                 row = "treatment",treatment.scores=alpha, outcome.scores=c(0,1),
                                                                 verbose = FALSE)



## Gamma 3.5 p value: 0.3716813
prior.three.by.two.V1.boy.Gamma.3.5.result = exact.score.sen.IxJ(obs.table = Black.boy.3x2.V1.table, gamma=gamma[6],delta = c(0,1,1),
                                                                 row = "treatment",treatment.scores=alpha, outcome.scores=c(0,1),
                                                                 verbose = FALSE)


## Gamma 4.0 p value: 0.3837562  
prior.three.by.two.V1.boy.Gamma.4.0.result = exact.score.sen.IxJ(obs.table = Black.boy.3x2.V1.table, gamma=gamma[7],delta = c(0,1,1),
                                                                 row = "treatment",treatment.scores=alpha, outcome.scores=c(0,1),
                                                                 verbose = FALSE)

## Gamma 4.5 p value: 0.3932709  
prior.three.by.two.V1.boy.Gamma.4.5.result = exact.score.sen.IxJ(obs.table = Black.boy.3x2.V1.table, gamma=gamma[8],delta = c(0,1,1),
                                                                 row = "treatment",treatment.scores=alpha, outcome.scores=c(0,1),
                                                                 verbose = FALSE)






#########################################
# prior treatment scores, 3x2, V2, boys
########################################


Gamma = c(1,1.5,2,2.5,3, 3.5,4,4.5)
gamma =log(Gamma)


## this is the treatment scores we come up with, non-data-driven 
alpha = c(0,0.25,1.5)


## RCT p-value:0.01223853 


## Gamma 1.5 pvalue:0.02795916
prior.three.by.two.V2.boy.Gamma.1.5.result = exact.score.sen.IxJ(obs.table = Black.boy.3x2.V2.table, gamma=gamma[2],delta = c(0,1,1),
                                                                 row = "treatment",treatment.scores=alpha, outcome.scores=c(0,1),
                                                                 verbose = FALSE)


## Gamma 2 p value: 0.04528917
prior.three.by.two.V2.boy.Gamma.2.result = exact.score.sen.IxJ(obs.table = Black.boy.3x2.V2.table, gamma=gamma[3],delta = c(0,1,1),
                                                               row = "treatment",treatment.scores=alpha, outcome.scores=c(0,1),
                                                               verbose = FALSE)
## Gamma 2.5 p value: 0.06286778 
prior.three.by.two.V2.boy.Gamma.2.5.result = exact.score.sen.IxJ(obs.table = Black.boy.3x2.V2.table, gamma=gamma[4],delta = c(0,1,1),
                                                                 row = "treatment",treatment.scores=alpha, outcome.scores=c(0,1),
                                                                 verbose = FALSE)

## Gamma 3.0 p value: 0.08017721 
prior.three.by.two.V2.boy.Gamma.3.0.result = exact.score.sen.IxJ(obs.table = Black.boy.3x2.V2.table, gamma=gamma[5],delta = c(0,1,1),
                                                                 row = "treatment",treatment.scores=alpha, outcome.scores=c(0,1),
                                                                 verbose = FALSE)



## Gamma 3.5 p value: 0.09701119
prior.three.by.two.V2.boy.Gamma.3.5.result = exact.score.sen.IxJ(obs.table = Black.boy.3x2.V2.table, gamma=gamma[6],delta = c(0,1,1),
                                                                 row = "treatment",treatment.scores=alpha, outcome.scores=c(0,1),
                                                                 verbose = FALSE)


## Gamma 4.0 p value: 0.1132757 
prior.three.by.two.V2.boy.Gamma.4.0.result = exact.score.sen.IxJ(obs.table = Black.boy.3x2.V2.table, gamma=gamma[7],delta = c(0,1,1),
                                                                 row = "treatment",treatment.scores=alpha, outcome.scores=c(0,1),
                                                                 verbose = FALSE)

## Gamma 4.5 p value: 0.1289191 
prior.three.by.two.V2.boy.Gamma.4.5.result = exact.score.sen.IxJ(obs.table = Black.boy.3x2.V2.table, gamma=gamma[8],delta = c(0,1,1),
                                                                 row = "treatment",treatment.scores=alpha, outcome.scores=c(0,1),
                                                                 verbose = FALSE)





#################################################################
#
#     profile treatment score 3x2 v1 for the boys 
#
#################################################################

Gamma = c(1,1.5,2,2.5,3, 3.5,4,4.5)
gamma =log(Gamma)

alpha =  c(0.000, 0.197, 1.598)



## RCT p-value: 0.1881044


## Gamma 1.5 pvalue: 0.2593915
profile.three.by.two.V1.boy.Gamma.1.5.result = exact.score.sen.IxJ(obs.table = Black.boy.3x2.V1.table, gamma=gamma[2],delta = c(0,1,1),
                                                                   row = "treatment",treatment.scores=alpha, outcome.scores=c(0,1),
                                                                   verbose = FALSE)


## Gamma 2 p value: 0.3044207 
profile.three.by.two.V1.boy.Gamma.2.result = exact.score.sen.IxJ(obs.table = Black.boy.3x2.V1.table, gamma=gamma[3],delta = c(0,1,1),
                                                                 row = "treatment",treatment.scores=alpha, outcome.scores=c(0,1),
                                                                 verbose = FALSE)
## Gamma 2.5 p value: 0.334582
profile.three.by.two.V1.boy.Gamma.2.5.result = exact.score.sen.IxJ(obs.table = Black.boy.3x2.V1.table, gamma=gamma[4],delta = c(0,1,1),
                                                                   row = "treatment",treatment.scores=alpha, outcome.scores=c(0,1),
                                                                   verbose = FALSE)

## Gamma 3.0 p value: 0.3559135 
profile.three.by.two.V1.boy.Gamma.3.0.result = exact.score.sen.IxJ(obs.table = Black.boy.3x2.V1.table, gamma=gamma[5],delta = c(0,1,1),
                                                                   row = "treatment",treatment.scores=alpha, outcome.scores=c(0,1),
                                                                   verbose = FALSE)



## Gamma 3.5 p value: 0.3716813
profile.three.by.two.V1.boy.Gamma.3.5.result = exact.score.sen.IxJ(obs.table = Black.boy.3x2.V1.table, gamma=gamma[6],delta = c(0,1,1),
                                                                   row = "treatment",treatment.scores=alpha, outcome.scores=c(0,1),
                                                                   verbose = FALSE)


## Gamma 4.0 p value: 0.3837562  
profile.three.by.two.V1.boy.Gamma.4.0.result = exact.score.sen.IxJ(obs.table = Black.boy.3x2.V1.table, gamma=gamma[7],delta = c(0,1,1),
                                                                   row = "treatment",treatment.scores=alpha, outcome.scores=c(0,1),
                                                                   verbose = FALSE)

## Gamma 4.5 p value: 0.3932709  
profile.three.by.two.V1.boy.Gamma.4.5.result = exact.score.sen.IxJ(obs.table = Black.boy.3x2.V1.table, gamma=gamma[8],delta = c(0,1,1),
                                                                   row = "treatment",treatment.scores=alpha, outcome.scores=c(0,1),
                                                                   verbose = FALSE)




#################################################################
#
#     profile treatment score 3x2 v2 for the boys 
#
#################################################################

Gamma = c(1,1.5,2,2.5,3, 3.5,4,4.5)
gamma =log(Gamma)

alpha =  c(0.000, 0.197, 1.598)



## RCT p-value: 0.01223853


## Gamma 1.5 pvalue: 0.02795916
profile.three.by.two.V2.boy.Gamma.1.5.result = exact.score.sen.IxJ(obs.table = Black.boy.3x2.V2.table, gamma=gamma[2],delta = c(0,1,1),
                                                                   row = "treatment",treatment.scores=alpha, outcome.scores=c(0,1),
                                                                   verbose = FALSE)


## Gamma 2 p value:0.04528917 
profile.three.by.two.V2.boy.Gamma.2.result = exact.score.sen.IxJ(obs.table = Black.boy.3x2.V2.table, gamma=gamma[3],delta = c(0,1,1),
                                                                 row = "treatment",treatment.scores=alpha, outcome.scores=c(0,1),
                                                                 verbose = FALSE)
## Gamma 2.5 p value: 0.06286778  
profile.three.by.two.V2.boy.Gamma.2.5.result = exact.score.sen.IxJ(obs.table = Black.boy.3x2.V2.table, gamma=gamma[4],delta = c(0,1,1),
                                                                   row = "treatment",treatment.scores=alpha, outcome.scores=c(0,1),
                                                                   verbose = FALSE)

## Gamma 3.0 p value: 0.08017721 
profile.three.by.two.V2.boy.Gamma.3.0.result = exact.score.sen.IxJ(obs.table = Black.boy.3x2.V2.table, gamma=gamma[5],delta = c(0,1,1),
                                                                   row = "treatment",treatment.scores=alpha, outcome.scores=c(0,1),
                                                                   verbose = FALSE)



## Gamma 3.5 p value:0.09701119 
profile.three.by.two.V2.boy.Gamma.3.5.result = exact.score.sen.IxJ(obs.table = Black.boy.3x2.V2.table, gamma=gamma[6],delta = c(0,1,1),
                                                                   row = "treatment",treatment.scores=alpha, outcome.scores=c(0,1),
                                                                   verbose = FALSE)

## Gamma 4.0 p value: 0.1132757  
profile.three.by.two.V2.boy.Gamma.4.0.result = exact.score.sen.IxJ(obs.table = Black.boy.3x2.V2.table, gamma=gamma[7],delta = c(0,1,1),
                                                                   row = "treatment",treatment.scores=alpha, outcome.scores=c(0,1),
                                                                   verbose = FALSE)


## Gamma 4.5 p value: 0.1289191   
profile.three.by.two.V2.boy.Gamma.4.5.result = exact.score.sen.IxJ(obs.table = Black.boy.3x2.V2.table, gamma=gamma[8],delta = c(0,1,1),
                                                                   row = "treatment",treatment.scores=alpha, outcome.scores=c(0,1),
                                                                   verbose = FALSE)
















#----------------------------------------------------------------------------

# End of Script
# -----------------------------------



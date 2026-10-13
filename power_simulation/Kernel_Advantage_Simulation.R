
#######################################################################################################
#                      Advantage of using kernels                                                     #
#               This file presents simulation results for the paper "Exact Sensitivity Analysis       #
#               for Observational Studies of Contingency tables", Section 4, Benefits of kernels      #
#               and additional simulations in Appendix F.                                             #
#                                                                                                     #
#                                                                                                     #
#               This file compares the speed of computing the exact p-value                           #
#               based on either permutation of treatment or kernel expression;                        #
#               and develop and compares the convergence behavior of sampling                         #
#               estimation of p-value based on kernels and permutation of treatments.                 #
#                                                                                                     #
#######################################################################################################


#######################################################################################################
# Method Brief: 
# The kernel expression can be combined with modern table sampling method. We adapt the C code from 
# "Sampling for Conditional Inference on Contingency Tables” (2017)", which uses the hypergeometric
# conditional probability as the proposal probability to sample tables. The adapted codes is in the 
# file "SIS_methods.cpp", which will be used in combination of kernels computation in the file 
# "sensitivityIxJ.cpp" and "sensitivityIxJ.R" to develop the sampling method. For the permutation 
# treatment method, we use the multicool library developed by developed by James Curran, Aaron Williams, 
# Jerome Kelleher, Dave Barber, see https://github.com/jmcurran/multicool. This library implements in Cpp 
# and generates the permutation one at a time. The algorithm is described in A. Loopless Generation of 
# Multiset Permutations by Prefix Shifts.
#
########################################################################################################


############### Load the Libraries ###################################

library(Rcpp)
sourceCpp("SIS_methods.cpp")
sourceCpp("sensitivityIxJ.cpp")
source("sensitivityIxJ.R")
library(multicool)




####################################################################################################### 
#   Comparison 1: The exact computation of p-value, using kernel and permutation, the speed test 
#   In the following we develop experiment functions to compute the pvalue based on the permutation 
#   of test label space and the kernel space, and compare the speeds
#######################################################################################################


time_permute_pvalue <- function(
    table,             # the original contingency table
    treatment.margins, # numeric vector of length 2 or 3, e.g. c(5,10,2)
    outcome.margins,   # numeric vector of length 2 or 3
    u.vec,             # 0/1 factor or numeric
    transform.fun,     # function that takes a (treatment x outcome) table and returns a statistic
    delta,             # vector of length = # of treatment categories, used in the exponent
    gamma
) {
  ## compute the observed test statistic on the *original* table
  obs.stat <- transform.fun(table)
  
  ## length check
  if (length(u.vec) != sum(table)) {
    stop("the length of u.vec is incorrect")
  }
  
  ## Build 'treatment.init' (as integers, no factor)
  if (length(treatment.margins) == 2) {
    treatment.init <- c(
      rep(1L, treatment.margins[1]),
      rep(2L, treatment.margins[2])
    )
  } else if (length(treatment.margins) == 3) {
    treatment.init <- c(
      rep(1L, treatment.margins[1]),
      rep(2L, treatment.margins[2]),
      rep(3L, treatment.margins[3])
    )
  } else {
    stop("Currently only supports 2 or 3 treatment margins.")
  }
  n_treat <- length(treatment.margins)  # 2 or 3
  
  ## Build 'outcome' (but we'll store it as an integer for speed, 1-based)
  if (length(outcome.margins) == 2) {
    outcome <- c(
      rep(1L, outcome.margins[1]),
      rep(2L, outcome.margins[2])
    )
  } else if (length(outcome.margins) == 3) {
    outcome <- c(
      rep(1L, outcome.margins[1]),
      rep(2L, outcome.margins[2]),
      rep(3L, outcome.margins[3])
    )
  } else {
    stop("Currently only supports 2 or 3 outcome margins.")
  }
  n_outcome <- length(outcome.margins)  # 2 or 3
  
  ## Convert u.vec to an integer 0/1 vector (avoid factor overhead)
  # If u.vec is a factor with levels=c("0","1"), as.integer() yields (1,2).
  # We'll shift by 1 to make it (0,1).
  if (is.factor(u.vec)) {
    u_int <- as.integer(u.vec) - 1L
  } else {
    # assume it’s already 0/1
    u_int <- as.integer(u.vec)
  }
  
  ## Initialize a multicool object (avoiding factor):
  # 'multicool::initMC()' can accept an integer vector. 
  # Then 'nextPerm()' returns a permuted vector of the same integers.
  mc.object <- multicool::initMC(treatment.init)
  
  # Number of unique permutations:
  label_counts <- table(treatment.init)
  total_perms <- factorial(length(treatment.init)) / prod(factorial(label_counts))
  
  # We'll do partial sums for numerator/denominator
  numerator.summation <- 0
  denominator.summation <- 0
  
  # Prepare for iteration
  i <- 0L
  start.time <- Sys.time()
  
  while (i < total_perms) {
    # Generate next permutation of 'treatment.init'
    perm_treat <- nextPerm(mc.object)  
    # 'perm_treat' is an integer vector of length sum(treatment.margins),
    # containing values in {1,2} or {1,2,3}.
    
    ## 1) Compute the counts of how many "u=1" in each treatment category
    #    That is the same logic as "conf.table = table(perm_treat, u.vec)" then grabbing col=1.
    #    We do it manually:
    if (n_treat == 2L) {
      u1.counts <- c(0L, 0L)
      # for j in 1..N, if u_int[j]==1 => increment the relevant treatment
      # but we’ll do it more vectorized or with tapply for pure R, or a for-loop:
      for (idx in seq_along(perm_treat)) {
        if (u_int[idx] == 1L) {
          u1.counts[perm_treat[idx]] <- u1.counts[perm_treat[idx]] + 1L
        }
      }
    } else {
      u1.counts <- c(0L, 0L, 0L)
      for (idx in seq_along(perm_treat)) {
        if (u_int[idx] == 1L) {
          u1.counts[perm_treat[idx]] <- u1.counts[perm_treat[idx]] + 1L
        }
      }
    }
    
    # denominator.term = exp(gamma * sum(delta * u1.counts))
    # assume delta is numeric of length n_treat
    denominator.term <- exp(gamma * sum(delta * u1.counts))
    
    ## 2) Build the (treatment x outcome) table for 'transform.fun'
    # Instead of 'table(perm_treat, outcome)', do a small matrix fill:
    perm_table_mat <- matrix(0L, nrow = n_treat, ncol = n_outcome)
    for (idx in seq_along(perm_treat)) {
      # treat category is perm_treat[idx] (1-based),
      # outcome category is outcome[idx] (also 1-based).
      perm_table_mat[perm_treat[idx], outcome[idx]] <-
        perm_table_mat[perm_treat[idx], outcome[idx]] + 1L
    }
    
    # Evaluate the test statistic for this permutation:
    perm_stat <- transform.fun(perm_table_mat)
    
    # If perm_stat >= obs.stat => multiply denominator.term by 1, else 0
    numerator.term <- if (perm_stat >= obs.stat) denominator.term else 0
    
    numerator.summation   <- numerator.summation + numerator.term
    denominator.summation <- denominator.summation + denominator.term
    
    i <- i + 1L
  }
  
  end.time <- Sys.time()
  run.time <- difftime(end.time, start.time, units = "secs")
  
  list(
    permute.pvalue = numerator.summation / denominator.summation,
    run.time       = run.time,
    total_perms    = total_perms
  )
}



time_kernel_pvalue = function(transform.fun,obs.table,u_allocation,gamma, delta){
  start.time = Sys.time()
  obs.stat = transform.fun(obs.table)
  table.list = possible.table(threshold = obs.stat, table = obs.table,direction = "greater than", transform.fun = transform.fun)
  
  sen.result = exact.general.sen.IxJ(u_space = t(as.matrix(u_allocation)),table_space = table.list,gamma=gamma, delta = delta, obs.table = obs.table)
  
  end.time = Sys.time()
  run.time = difftime(end.time,start.time, units="secs")
  return(list(run.time = run.time,p.value = sen.result$max.prob))
  
}


exact_p_experiment = function(transform.fun,gamma,delta,obs.table,u.vec, run=10){
  outcome.margins = colSums(obs.table)
  treatment.margins = rowSums(obs.table)
  outcome = factor(c(rep(1,outcome.margins[1]),rep(2,outcome.margins[2]), rep(3,outcome.margins[3])),levels=1:3)
  conf.table = table(outcome,u.vec)
  u_allocation = rep(0, 3)
  if ("1" %in% colnames(conf.table)) {
    u_allocation[as.integer(rownames(conf.table))] = conf.table[,"1"]
  }
  permute.time.all.run = rep(NA,run)
  kernel.time.all.run = rep(NA,run)
  permute.p.all.run = rep(NA,run)
  kernel.p.all.run = rep(NA,run)
  for(i in 1:run){
    ## call two functions and return their p-values to verify it is the same except for numerical precision, and 
    ## return their times 
    permute.result = time_permute_pvalue(table = obs.table,treatment.margins = treatment.margins,outcome.margins = outcome.margins,u.vec=u.vec,transform.fun = transform.fun,gamma=gamma,delta = delta)
    
    kernel.result = time_kernel_pvalue(transform.fun = transform.fun, obs.table = obs.table, u_allocation = u_allocation, gamma=gamma, delta = delta)
    one.time.result = list(permute.time = permute.result$run.time, kernel.time = kernel.result$run.time,permute.p = permute.result$permute.pvalue, kernel.p = kernel.result$p.value, u_allocation = u_allocation)
    permute.time.all.run[i] = one.time.result$permute.time
    kernel.time.all.run[i] = one.time.result$kernel.time
    permute.p.all.run[i]   = one.time.result$permute.p
    kernel.p.all.run[i] = one.time.result$kernel.p
    
  }
  
  return(list(permute.time.all.run, kernel.time.all.run, 
              mean.permute.time = mean(permute.time.all.run), sd.permute.time = sd(permute.time.all.run),
              mean.kernel.time = mean(kernel.time.all.run), sd.kernel.time = sd(kernel.time.all.run),
              permute.p.all.run, kernel.p.all.run, u_allocation))
  
  
}


##############################################################################################################
#
#             Here conducts the experiment of comparing the speed finishing the exact computation
#             of p-value, we are going to try several tables, 3 (treatment) by 3 (outcome) 
#             tables with treatment margins (5,5,5), (6,6,6) with varying p-values, and 2(treatment)
#             by 3(outcome) tables. We illustrate two: (1) The exact computation based on 
#             kernel method is always faster than that with the permutation method (2) The time growth
#             to compute the exact p-value in terms of the treatment margins is way greater than the time
#             growth based on the kernel method.
#
##############################################################################################################

#### We recommend to clean your memory before the speed testing, using gc() before 
#### each batch of experiment



###### Experiment 1 - Based on a table with treatment margins (5,5,5), the first three rows  ################

#=======================
# First result
#=======================
sink(file = "output.txt", append = TRUE)

print("First result:")

gc()
alpha <- c(0,1,2)
beta <- c(0,1,2)
scores <- outer(alpha, beta)

transform.fun <- function(tb) {
  value <- sum(scores * tb)
  return(value)
}

gamma <- 1
delta <- c(0,1,1)
obs.table <- matrix(data = c(2,3,0,0,1,4,0,1,4), nrow = 3, byrow = TRUE)
u.vec <- factor(c(rep(0,12), rep(1,3)), levels = 0:1)

first.result <- exact_p_experiment(
  transform.fun = transform.fun, 
  gamma = gamma,
  delta = delta,
  obs.table = obs.table,
  u.vec = u.vec,
  run = 10
)

print(first.result)

sink()

#=======================
# Second result
#=======================
sink(file = "output.txt", append = TRUE)

print("Second result:")

u.vec <- factor(c(rep(0,8), rep(1,7)), levels = 0:1)

second.result <- exact_p_experiment(
  transform.fun = transform.fun,
  gamma = gamma,
  delta = delta,
  obs.table = obs.table,
  u.vec = u.vec,
  run = 10
)

print(second.result)

sink()

#=======================
# Third result
#=======================
sink(file = "output.txt", append = TRUE)

print("Third result:")

u.vec <- factor(c(rep(0,2), rep(1,13)), levels = 0:1)

third.result <- exact_p_experiment(
  transform.fun = transform.fun,
  gamma = gamma,
  delta = delta,
  obs.table = obs.table,
  u.vec = u.vec
)

print(third.result)

sink()




######  Experiment 2 - another table with treatment margins (5,5,5), different rejection region
#=======================
# Fourth result
#=======================
sink(file = "output.txt", append = TRUE)

print("Fourth result:")

gc()
alpha <- c(0,1,2)
beta <- c(0,1,2)
scores <- outer(alpha, beta)

transform.fun <- function(tb) {
  value <- sum(scores * tb)
  return(value)
}

gamma <- 1
delta <- c(0,1,1)
obs.table <- matrix(data = c(2,2,1,1,2,2,1,2,2), nrow = 3, byrow = TRUE)

u.vec <- factor(c(rep(0,13), rep(1,2)), levels = 0:1)

fourth.result <- exact_p_experiment(
  transform.fun = transform.fun,
  gamma = gamma,
  delta = delta,
  obs.table = obs.table,
  u.vec = u.vec
)

print(fourth.result)

sink()

#=======================
# Fifth result
#=======================
sink(file = "output.txt", append = TRUE)

print("Fifth result:")

u.vec <- factor(c(rep(0,7), rep(1,8)), levels = 0:1)

fifth.result <- exact_p_experiment(
  transform.fun = transform.fun,
  gamma = gamma,
  delta = delta,
  obs.table = obs.table,
  u.vec = u.vec
)

print(fifth.result)

sink()

#=======================
# Sixth result
#=======================
sink(file = "output.txt", append = TRUE)

print("Sixth result:")

u.vec <- factor(c(rep(0,4), rep(1,11)), levels = 0:1)

sixth.result <- exact_p_experiment(
  transform.fun = transform.fun,
  gamma = gamma,
  delta = delta,
  obs.table = obs.table,
  u.vec = u.vec
)

print(sixth.result)

sink()



########### Experiment 3 - A table with treatment margins (6,6,6) ##########################################

#=======================
# Seventh result
#=======================
sink(file = "output.txt", append = TRUE)

print("Seventh result:")

gc()
alpha <- c(0,1,2)
beta <- c(0,1,2)
scores <- outer(alpha, beta)

transform.fun <- function(tb) {
  value <- sum(scores * tb)
  return(value)
}

gamma <- 1
delta <- c(0,1,1)
obs.table <- matrix(data = c(3,2,1,0,2,4,0,1,5), nrow = 3, byrow = TRUE)

u.vec <- factor(c(rep(0,13), rep(1,5)), levels = 0:1)

seventh.result <- exact_p_experiment(
  transform.fun = transform.fun,
  gamma = gamma,
  delta = delta,
  obs.table = obs.table,
  u.vec = u.vec
)

print(seventh.result)

sink()

#=======================
# Eighth result
#=======================
sink(file = "output.txt", append = TRUE)

print("Eighth result:")

u.vec <- factor(c(rep(0,8), rep(1,10)), levels = 0:1)

eighth.result <- exact_p_experiment(
  transform.fun = transform.fun,
  gamma = gamma,
  delta = delta,
  obs.table = obs.table,
  u.vec = u.vec
)

print(eighth.result)

sink()

#=======================
# Ninth result
#=======================
sink(file = "output.txt", append = TRUE)

print("Ninth result:")

u.vec <- factor(c(rep(0,4), rep(1,14)), levels = 0:1)

ninth.result <- exact_p_experiment(
  transform.fun = transform.fun,
  gamma = gamma,
  delta = delta,
  obs.table = obs.table,
  u.vec = u.vec
)

print(ninth.result)

sink()





########### Experiment 4 - Another table with treatment margins (6,6,6) #######################################

#=======================
# Tenth result
#=======================
sink(file = "output.txt", append = TRUE)

print("Tenth result:")

gc()
alpha <- c(0,1,2)
beta <- c(0,1,2)
scores <- outer(alpha, beta)

transform.fun <- function(tb) {
  value <- sum(scores * tb)
  return(value)
}

gamma <- 1
delta <- c(0,1,1)
obs.table <- matrix(data = c(3,3,0,1,2,3,2,3,1), nrow = 3, byrow = TRUE)

u.vec <- factor(c(rep(0,14), rep(1,4)), levels = 0:1)

tenth.result <- exact_p_experiment(
  transform.fun = transform.fun,
  gamma = gamma,
  delta = delta,
  obs.table = obs.table,
  u.vec = u.vec
)

print(tenth.result)

sink()

#=======================
# Eleventh result
#=======================
sink(file = "output.txt", append = TRUE)

print("Eleventh result:")

u.vec <- factor(c(rep(0,11), rep(1,7)), levels = 0:1)

eleventh.result <- exact_p_experiment(
  transform.fun = transform.fun,
  gamma = gamma,
  delta = delta,
  obs.table = obs.table,
  u.vec = u.vec
)

print(eleventh.result)

sink()

#=======================
# Twelfth result
#=======================
sink(file = "output.txt", append = TRUE)

print("Twelfth result:")

u.vec <- factor(c(rep(0,6), rep(1,12)), levels = 0:1)

twelfth.result <- exact_p_experiment(
  transform.fun = transform.fun,
  gamma = gamma,
  delta = delta,
  obs.table = obs.table,
  u.vec = u.vec
)

print(twelfth.result)

sink()





####################### Experiment 5 - A 2 by 3 table ################################################## 

#=======================
# Thirteenth result
#=======================
sink(file = "output.txt", append = TRUE)

print("Thirteenth result:")

gc()
alpha <- c(0,1)
beta <- c(0,1,2)
scores <- outer(alpha, beta)

transform.fun <- function(tb) {
  value <- sum(scores * tb)
  return(value)
}

gamma <- 1
delta <- c(0,1)
obs.table <- matrix(data = c(4,6,0,1,3,6), nrow = 2, byrow = TRUE)

u.vec <- factor(c(rep(0,1), rep(1,19)), levels = 0:1)

thirteenth.result <- exact_p_experiment(
  transform.fun = transform.fun,
  gamma = gamma,
  delta = delta,
  obs.table = obs.table,
  u.vec = u.vec
)

print(thirteenth.result)

sink()

#=======================
# Fourteenth result
#=======================
sink(file = "output.txt", append = TRUE)

print("Fourteenth result:")

u.vec <- factor(c(rep(0,13), rep(1,7)), levels = 0:1)

fourteenth.result <- exact_p_experiment(
  transform.fun = transform.fun,
  gamma = gamma,
  delta = delta,
  obs.table = obs.table,
  u.vec = u.vec
)

print(fourteenth.result)

sink()





################ Experiment 6 - Another 2 by 3 table ########################################################

#=======================
# Fifteenth result
#=======================
sink(file = "output.txt", append = TRUE)

print("Fifteenth result:")

gc()
alpha <- c(0,1)
beta <- c(0,1,2)
scores <- outer(alpha, beta)

transform.fun <- function(tb) {
  value <- sum(scores * tb)
  return(value)
}

gamma <- 1
delta <- c(0,1)
obs.table <- matrix(data = c(2,4,4,2,2,6), nrow = 2, byrow = TRUE)

u.vec <- factor(c(rep(0,1), rep(1,19)), levels = 0:1)

fifteenth.result <- exact_p_experiment(
  transform.fun = transform.fun,
  gamma = gamma,
  delta = delta,
  obs.table = obs.table,
  u.vec = u.vec
)

print(fifteenth.result)

sink()

#=======================
# Sixteenth result
#=======================
sink(file = "output.txt", append = TRUE)

print("Sixteenth result:")

u.vec <- factor(c(rep(0,13), rep(1,7)), levels = 0:1)

sixteenth.result <- exact_p_experiment(
  transform.fun = transform.fun,
  gamma = gamma,
  delta = delta,
  obs.table = obs.table,
  u.vec = u.vec
)

print(sixteenth.result)

sink()



#######################################################################################################################################################
#
#           Sampling Based Method Simulation
#           The second part of this file discusses the sampling-based method approach
#           We develop two table-based sampling method, one uses the normalizing constant, called SIS
#           The other one does not use the known normalizing constant, so we call it snSIS (self-normalizing)
#           The final one directly sample from the permutation of treatment assignment space
#           The SIS_Algorithm is from Sampling for Conditional Inference on Contingency Tables,  in equation (6) 
#           draws contingency tables from an approximately hypergeometric distribution, in that each cell is drawn 
#           with a probability called Good's approximation: 
#           Note that the proposal distribution is not the same as the distribution of tables under the sensitivity analysis. 
#           However, this algorithm provides the proposal distribution (i.e, the probability of a table being drawn from this algorithm), 
#           and this distribution is similar to that under sensitivity analysis. Therefore, since we have an expression for each table, 
#           we could adapt their algorithm and device an importance sampling.
#
#####################################################################################################################################################



################## Experiment Functions for sampling-based method comparison, a fully-developed codes will be in the library not here ##############




######################################################################################################################################################
######################################################################################################################################################
importance.sampling.SIS = function(mc.iteration = 5000, table,row="treatment",u_allocation,gamma, delta, transform.fun){
  
  ## record time
  start.time = Sys.time()
  ## first obtain the obs.stat as the critical value 
  
  obs.stat = transform.fun(table)
  
  
  
  # Validate 'table' is a matrix or table object
  if (!is.matrix(table) && !inherits(table, "table")) {
    stop("'table' must be either a matrix or a table object.")
  }
  
  # If 'table' is a matrix, ensure it's numeric and two-dimensional
  if (is.matrix(table)) {
    if (!is.numeric(table)) {
      stop("If 'table' is a matrix, it must be numeric.")
    }
    if (length(dim(table)) != 2) {
      stop("'table' must be a two-dimensional matrix.")
    }
  }
  
  # If 'table' is a table object, ensure it's two-dimensional
  if (inherits(table, "table")) {
    if (length(dim(table)) != 2) {
      stop("'table' must be a two-dimensional table object.")
    }
  }
  
  # Validate 'row' argument
  if (!is.character(row) || length(row) != 1) {
    stop("'row' must be a single string, either 'outcome' or 'treatment'.")
  }
  
  if (!(row %in% c("outcome", "treatment"))) {
    stop("'row' must be either 'outcome' or 'treatment'.")
  }
  
  
  ## first detect the table size and the margins, rotate the table if the row is "outcome" because our paper
  ## set the rows as treatment as default
  N = sum(table)
  
  ## rotate table if needed
  
  if(row=="outcome"){
    table = t(table)
  }
  
  ## here detects the dimensions, I is the number of treatment levels, j is the number of outcome levels
  I = dim(table)[1]
  
  J = dim(table)[2]
  
  treatment_margins = rowSums(table)
  outcome_margins   = colSums(table)
  
  

  
  if((I==3)&(J==3)){
    n1 = treatment_margins[1]
    n2 = treatment_margins[2]
    n3 = treatment_margins[3]
    Us = sum(u_allocation)
    m1 = outcome_margins[1]
    m2 = outcome_margins[2]
    m3 = outcome_margins[3]
    u1 = u_allocation[1]
    u2 = u_allocation[2]
    u3 = u_allocation[3]
    
    
    
    ## a cumulant for \sum I(T >=c)*v(t)/h(t)/normalizing constant
    all.weight.summation = 0
    
    ## the normalizing constant 
    normalizing = denominator_three_treatment(n1=n1,n2=n2,n3=n3,N=N,gamma_delta=gamma*delta,Us=Us,shared_divisor = .Machine$double.xmax^0.9)
    
    ## here records the probability progression 
    mc.probability = c()
    
    mc = 1
    ## here starts the monte carlo estimation 
    while(mc < mc.iteration){
      ## sample a table from the sis algorithm 
      mc.table.result = sample_sis_one_table(row_s =treatment_margins,col_s = outcome_margins,seed=mc)
      inverse.h = mc.table.result$weight
      mc.table = mc.table.result$table
      if(transform.fun(mc.table)<obs.stat){
        ## do not update, namely, add zero
      }else{
        target.v = d_numerator_three_by_three(n11=mc.table[1,1],n12=mc.table[1,2],n21=mc.table[2,1],n22=mc.table[2,2],m1=m1,m2=m2,m3=m3,n1=n1,n2=n2,n3=n3,N=N,gamma_delta=gamma*delta,u1=u1,u2=u2,u3=u3,shared_divisor = .Machine$double.xmax^0.9)
        weight =inverse.h*target.v
        all.weight.summation = all.weight.summation+weight/normalizing
      }
      current.probability = all.weight.summation/mc
      mc.probability[mc] = current.probability
      
      ## update mc
      mc = mc + 1 
      
      
    }
    end.time = Sys.time()
    run.time <- difftime(end.time, start.time, units = "secs")
    return(list(probability.progression = mc.probability,final.prob.estimate = all.weight.summation/mc, run.time = run.time))
  }
  
}

#################################################################################
#################################################################################

### Self-normalizing importance sampling 

importance.sampling.snis = function(table, row = "treatment", mc.iteration = 5000,
                                    u_allocation, gamma, delta, transform.fun) {
  start.time = Sys.time()
  obs.stat = transform.fun(table)
  
  if (row == "outcome") {
    table = t(table)
  }
  
  N = sum(table)
  I = dim(table)[1]
  J = dim(table)[2]
  
  treatment_margins = rowSums(table)
  outcome_margins   = colSums(table)
  
  if ((I == 3) & (J == 3)) {
    n1 = treatment_margins[1]
    n2 = treatment_margins[2]
    n3 = treatment_margins[3]
    m1 = outcome_margins[1]
    m2 = outcome_margins[2]
    m3 = outcome_margins[3]
    u1 = u_allocation[1]
    u2 = u_allocation[2]
    u3 = u_allocation[3]
    
    numerator.summation = 0
    denominator.summation = 0
    mc.probability = numeric()
    
    mc = 1
    
    while (mc < mc.iteration) {
      mc.table.result = sample_sis_one_table(row_s = treatment_margins,
                                             col_s = outcome_margins, seed = mc)
      inverse.h = mc.table.result$weight
      mc.table = mc.table.result$table
      
      target.v = d_numerator_three_by_three(n11 = mc.table[1, 1], n12 = mc.table[1, 2],
                                            n21 = mc.table[2, 1], n22 = mc.table[2, 2],
                                            m1 = m1, m2 = m2, m3 = m3,
                                            n1 = n1, n2 = n2, n3 = n3,
                                            N = N, gamma_delta = gamma*delta,
                                            u1 = u1, u2 = u2, u3 = u3,
                                            shared_divisor = .Machine$double.xmax^0.9)
      
      weight = inverse.h * target.v
      
      indicator = ifelse(transform.fun(mc.table) >= obs.stat, 1, 0)
      
      numerator.summation = numerator.summation + indicator * weight
      denominator.summation = denominator.summation + weight
      
      current.probability = numerator.summation / denominator.summation
      mc.probability[mc] = current.probability
      
      mc = mc + 1
      
    }
    
    end.time = Sys.time()
    run.time <- difftime(end.time, start.time, units = "secs")
    return(list(
      probability.progression = mc.probability,
      final.prob.estimate = numerator.summation / denominator.summation,
      run.time = run.time,
      total.iteration = mc - 1
    ))
  }
}








################################################################################################################################################
################################################################################################################################################




permutation.sampling = function(table, treatment.margins, outcome.margins, u.vec, transform.fun, gamma, delta,mc.iteration=5000) {
  ## compute the observed test statistic
  obs.stat = transform.fun(table)
  
  ## define the multiset of treatment labels
  treatment.init = c(rep(1, treatment.margins[1]),
                     rep(2, treatment.margins[2]),
                     rep(3, treatment.margins[3]))
  n = length(treatment.init)
  
  ## define the fixed outcome vector
  outcome = factor(c(rep(1, outcome.margins[1]),
                     rep(2, outcome.margins[2]),
                     rep(3, outcome.margins[3])), levels = 1:3)
  
  u.vec <- factor(u.vec, levels = 0:1)
  
  numerator.summation = 0
  denominator.summation = 0
  mc.probability = c()
  
  mc = 1
  
  ## start the clock
  start.time = Sys.time()
  
  while (mc < mc.iteration) {
    ## sample a random permutation of treatment labels
    sample.treatment = factor(sample(treatment.init, size = n, replace = FALSE), levels = 1:3)
    
    ## compute weighted U=1 count per treatment group
    conf.table = table(sample.treatment, u.vec)
    u1.counts = rep(0, 3)
    if ("1" %in% colnames(conf.table)) {
      u1.counts[as.integer(rownames(conf.table))] = conf.table[,"1"]
    }
    
    denominator.term = exp(gamma * sum(delta * u1.counts))
    
    ## compute test statistic for this permutation
    perm.table = table(sample.treatment, outcome)
    numerator.term = ifelse(transform.fun(perm.table) >= obs.stat, 1, 0) * denominator.term
    
    numerator.summation = numerator.summation + numerator.term 
    denominator.summation = denominator.summation + denominator.term 
    
    
    current.probability = numerator.summation / denominator.summation
    mc.probability[mc] = current.probability
    
    ## update mc
    mc = mc + 1 
    
    
    
  }
  
  end.time = Sys.time()
  run.time <- difftime(end.time, start.time, units = "secs")
  
  return(list(
    probability.progression = mc.probability,
    final.prob.estimate = numerator.summation / denominator.summation,
    run.time = run.time
  ))
}


################################################################################################################################################
#################################################################################################################################################



## An experiment runner 
run_experiment = function(data.table, gamma,delta,u_allocation, u.vec, transform.fun,mc.iteration=10000){
  outcome.margins = colSums(data.table)
  treatment.margins = rowSums(data.table)
  
  ## total enumeration 
  start.time = Sys.time()
  table.list = possible.table(threshold = transform.fun(data.table),table=data.table, direction="greater than", transform.fun = transform.fun)
  
  exact.p = exact.general.sen.IxJ(u_space = t(as.matrix(u_allocation)), obs.table=data.table,table_space = table.list, gamma=gamma, delta = delta)$max.prob
  
  end.time = Sys.time()
  total.enumeration.time = difftime(end.time,start.time,units = "secs")
  
  #### Table-based sampling, algorithm 1 
  
  table.sample.SIS.result = importance.sampling.SIS(table=data.table,u_allocation=u_allocation,gamma=gamma, delta = delta,transform.fun = transform.fun,mc.iteration = mc.iteration)
  table.sample.SIS.time = table.sample.SIS.result$run.time
  table.sample.SIS.p=table.sample.SIS.result$final.prob.estimate
  
  
  table.sample.snSIS.result = importance.sampling.snis(table=data.table,u_allocation=u_allocation,gamma=gamma, delta = delta,transform.fun = transform.fun,mc.iteration = mc.iteration)
  table.sample.snSIS.time = table.sample.snSIS.result$run.time
  table.sample.snSIS.p=table.sample.snSIS.result$final.prob.estimate
  
  
  ### Permutation-treatment-based sampling 
  permutation.result = permutation.sampling(table=data.table,treatment.margins = treatment.margins,outcome.margins = outcome.margins,u.vec = u.vec,gamma=gamma,delta=delta,transform.fun = transform.fun, mc.iteration = mc.iteration)
  permutation.sample.time = permutation.result$run.time 
  permutation.p = permutation.result$final.prob.estimate
  
  return(list(data.table=data.table,total.enumeration.time=total.enumeration.time,
              table.sample.SIS.time = table.sample.SIS.time,
              table.sample.snSIS.time=table.sample.snSIS.time,
              permutation.sample.time = permutation.sample.time,
              exact.p = exact.p,
              table.sample.SIS.p=table.sample.SIS.p,
              table.sample.snSIS.p=table.sample.snSIS.p,
              permutation.p = permutation.p,
              table.SIS.p.progression = table.sample.SIS.result$probability.progression,
              table.snSIS.p.progression=table.sample.snSIS.result$probability.progression,
              permute.p.progression = permutation.result$probability.progression))
  
}


library(ggplot2)

plot_convergence <- function(result, skip.initial = 0, ylb = 0, yub = 1) {
  # Extract data
  table.SIS.p <- result$table.SIS.p.progression
  table.snSIS.p <- result$table.snSIS.p.progression
  permute.p <- result$permute.p.progression
  exact.p <- result$exact.p
  
  table.SIS.time <- result$table.sample.SIS.time
  table.snSIS.time <- result$table.sample.snSIS.time
  perm.time <- result$permutation.sample.time
  
  table.SIS.estimate <- result$table.sample.SIS.p
  table.snSIS.estimate <- result$table.sample.snSIS.p
  perm.estimate <- result$permutation.p
  
  # Build data frames
  df.SIS <- data.frame(Step = seq_along(table.SIS.p), Estimate = table.SIS.p, Method = "SIS")
  df.snSIS <- data.frame(Step = seq_along(table.snSIS.p), Estimate = table.snSIS.p, Method = "snSIS")
  df.permute <- data.frame(Step = seq_along(permute.p), Estimate = permute.p, Method = "PermTreat")
  
  # Combine and filter
  df <- rbind(df.SIS, df.snSIS, df.permute)
  df <- df[df$Step > skip.initial, ]
  df$Method <- factor(df$Method, levels = c("SIS", "snSIS", "PermTreat"))
  
  
  #subtitle.text <- paste0(
  # "Exact p = ", signif(exact.p, 4), "\n",
  #  "SIS: ", signif(table.SIS.estimate, 4), " (", round(as.numeric(table.SIS.time), 2), "s)\n",
  #"snSIS: ", signif(table.snSIS.estimate, 4), " (", round(as.numeric(table.snSIS.time), 2), "s)\n",
  #  "Perm: ", signif(perm.estimate, 4), " (", round(as.numeric(perm.time), 2), "s)"
  # )
  
  # Create plot
  p <- ggplot(df, aes(x = Step, y = Estimate, color = Method)) +
    geom_line(size = 0.7) +
    geom_hline(yintercept = exact.p, linetype = "dashed", color = "black") +
    labs(
      y = "Estimated p-value",
      x = "Monte Carlo Iteration"
    ) +
    ylim(ylb, yub) +   # <-- plus sign added here
    scale_color_manual(
      values = c("SIS" = "red", "snSIS" = "blue", "PermTreat" = "darkgreen"),
      labels = c(expression(hat(alpha)[SIS]), expression(hat(alpha)[snSIS]), expression(hat(alpha)[PermTreat]))
    ) +
    theme_minimal() +
    theme(
      legend.title = element_blank(),
      legend.position = "bottom",
      legend.spacing.x = unit(1, "cm"),  # Increase the value for more space
      plot.subtitle = element_text(size = 15, face = "italic", family = "Times New Roman"),
      axis.title = element_text(size = 15, family = "Times New Roman"),
      axis.text = element_text(size = 15, family = "Times New Roman"),
      legend.text = element_text(size = 15, family = "Times New Roman") # Legend text size
    )
  
  return(p)
}

#### The naming of the experiment is Experiment - Table - Corner #########################
################## Experiment 1-1 ##########################################################################

gc()
set.seed(100)

data.table = matrix(data=c(12,6,6,18,12,6,5,6,15),ncol=3)

data.table


delta = c(0,1,1)
gamma = 1
u_allocation = c(0,10,20)
outcome.margins = colSums(data.table)

u.vec = c(rep(0,outcome.margins[1]-u_allocation[1]),rep(1,u_allocation[1]),
          rep(0,outcome.margins[2]-u_allocation[2]),rep(1,u_allocation[2]),
          rep(0,outcome.margins[3]-u_allocation[3]),rep(1,u_allocation[3]))
alpha = c(0,1,2.5)
beta = c(0,1,2)
score = outer(alpha,beta)
transform.fun = function(tb){
  return(sum(score*tb))
}


Experiment.1.1.result = run_experiment(data.table=data.table,gamma=gamma,
                                       delta = delta,u_allocation = u_allocation,u.vec=u.vec,transform.fun = transform.fun)

p1.1= plot_convergence(Experiment.1.1.result, ylb = 0, yub=0.10)

p1.1

ggsave("plot_convergence_experiment_1_1.png", plot = p1.1, width = 8, height = 6, dpi = 1200)


experiment.1.1.result.no.progression <- Experiment.1.1.result[!grepl("progression", names(Experiment.1.1.result))]

experiment.1.1.result.no.progression






############################### Experiment 1 - 2 ############################### 


gc()
set.seed(100)

data.table = matrix(data=c(12,6,6,18,12,6,5,6,15),ncol=3)

data.table


delta = c(0,1,1)
gamma = 1
u_allocation = c(0,36,26)
outcome.margins = colSums(data.table)

u.vec = c(rep(0,outcome.margins[1]-u_allocation[1]),rep(1,u_allocation[1]),
          rep(0,outcome.margins[2]-u_allocation[2]),rep(1,u_allocation[2]),
          rep(0,outcome.margins[3]-u_allocation[3]),rep(1,u_allocation[3]))
alpha = c(0,1,2.5)
beta = c(0,1,2)
score = outer(alpha,beta)

Experiment.1.2.result = run_experiment(data.table=data.table,gamma=gamma,delta = delta,u_allocation = u_allocation,u.vec=u.vec,transform.fun = transform.fun)

p1.2= plot_convergence(Experiment.1.2.result, ylb = 0, yub=0.20)

p1.2

ggsave("plot_convergence_experiment_1_2.png", plot = p1.2, width = 8, height = 6, dpi = 1200)


experiment.1.2.result.no.progression <- Experiment.1.2.result[!grepl("progression", names(Experiment.1.2.result))]

experiment.1.2.result.no.progression



############################## Experiment 2 - 1 #################################

gc()
set.seed(100)

data.table = matrix(data=c(12,18,17,3,12,25,0,3,4),ncol=3)

data.table


delta = c(0,1,1)
gamma =0.5
u_allocation = c(0,30,7)
outcome.margins = colSums(data.table)

u.vec = c(rep(0,outcome.margins[1]-u_allocation[1]),rep(1,u_allocation[1]),
          rep(0,outcome.margins[2]-u_allocation[2]),rep(1,u_allocation[2]),
          rep(0,outcome.margins[3]-u_allocation[3]),rep(1,u_allocation[3]))
alpha = c(0,1,2.5)
beta = c(0,1,2)
score = outer(alpha,beta)
transform.fun = function(tb){
  return(sum(score*tb))
}




Experiment.2.1.result = run_experiment(data.table=data.table,gamma=gamma,delta = delta,u_allocation = u_allocation,u.vec=u.vec,transform.fun = transform.fun)



p2.1= plot_convergence(Experiment.2.1.result)

p2.1

# ggsave("plot_convergence_experiment_2_1.png", plot = p2.1, width = 8, height = 6, dpi = 300)


experiment.2.1.result.no.progression <- Experiment.2.1.result[!grepl("progression", names(Experiment.2.1.result))]

experiment.2.1.result.no.progression



######################### Experiment 2 - 2 ########################################


gc()
set.seed(100)

data.table = matrix(data=c(12,18,17,3,12,25,0,3,4),ncol=3)

data.table


delta = c(0,1,1)
gamma =0.5
u_allocation = c(0,40,7)
outcome.margins = colSums(data.table)

u.vec = c(rep(0,outcome.margins[1]-u_allocation[1]),rep(1,u_allocation[1]),
          rep(0,outcome.margins[2]-u_allocation[2]),rep(1,u_allocation[2]),
          rep(0,outcome.margins[3]-u_allocation[3]),rep(1,u_allocation[3]))
alpha = c(0,1,2.5)
beta = c(0,1,2)
score = outer(alpha,beta)
transform.fun = function(tb){
  return(sum(score*tb))
}




Experiment.2.2.result = run_experiment(data.table=data.table,gamma=gamma,delta = delta,u_allocation = u_allocation,u.vec=u.vec,transform.fun = transform.fun)



p2.2= plot_convergence(Experiment.2.2.result, ylb=0, yub=0.10)

p2.2

ggsave("plot_convergence_experiment_2_2.png", plot = p2.2, width = 8, height = 6, dpi = 1200)


experiment.2.2.result.no.progression <- Experiment.2.2.result[!grepl("progression", names(Experiment.2.2.result))]

experiment.2.2.result.no.progression




################## Experiment 2 - 3 with a larger gamma #####################################

gc()
set.seed(100)

data.table = matrix(data=c(12,18,17,3,12,25,0,3,4),ncol=3)

data.table


delta = c(0,1,1)
gamma =1
u_allocation = c(0,30,5)
outcome.margins = colSums(data.table)

u.vec = c(rep(0,outcome.margins[1]-u_allocation[1]),rep(1,u_allocation[1]),
          rep(0,outcome.margins[2]-u_allocation[2]),rep(1,u_allocation[2]),
          rep(0,outcome.margins[3]-u_allocation[3]),rep(1,u_allocation[3]))
alpha = c(0,1,2.5)
beta = c(0,1,2)
score = outer(alpha,beta)
transform.fun = function(tb){
  return(sum(score*tb))
}


Experiment.2.3.result = run_experiment(data.table=data.table,gamma=gamma,delta = delta,u_allocation = u_allocation,u.vec=u.vec,transform.fun = transform.fun)



p2.3= plot_convergence(Experiment.2.3.result,ylb=0, yub=0.15)

p2.3

ggsave("plot_convergence_experiment_2_3.png", plot = p2.3, width = 8, height = 6, dpi = 1200)


experiment.2.3.result.no.progression <- Experiment.2.3.result[!grepl("progression", names(Experiment.2.3.result))]

experiment.2.3.result.no.progression




################## Experiment 2 - 4 with a larger gamma #####################################

gc()
set.seed(100)

data.table = matrix(data=c(12,18,17,3,12,25,0,3,4),ncol=3)

data.table


delta = c(0,1,1)
gamma =1
u_allocation = c(0,40,7)
outcome.margins = colSums(data.table)

u.vec = c(rep(0,outcome.margins[1]-u_allocation[1]),rep(1,u_allocation[1]),
          rep(0,outcome.margins[2]-u_allocation[2]),rep(1,u_allocation[2]),
          rep(0,outcome.margins[3]-u_allocation[3]),rep(1,u_allocation[3]))
alpha = c(0,1,2.5)
beta = c(0,1,2)
score = outer(alpha,beta)
transform.fun = function(tb){
  return(sum(score*tb))
}



Experiment.2.4.result = run_experiment(data.table=data.table,gamma=gamma,delta = delta,u_allocation = u_allocation,u.vec=u.vec,transform.fun = transform.fun)



p2.4= plot_convergence(Experiment.2.4.result,ylb=0,yub=0.25)

p2.4

ggsave("plot_convergence_experiment_2_4.png", plot = p2.4, width = 8, height = 6, dpi = 1200)


experiment.2.4.result.no.progression <- Experiment.2.4.result[!grepl("progression", names(Experiment.2.4.result))]

experiment.2.4.result.no.progression









##################### Experiment 3 - 1 #########################################

gc()
set.seed(100)

data.table = matrix(data=c(10,29,20,8,11,24,1,3,6),ncol=3)

data.table


delta = c(0,1,1)
gamma =0.5
u_allocation = c(0,20,10)
outcome.margins = colSums(data.table)

u.vec = c(rep(0,outcome.margins[1]-u_allocation[1]),rep(1,u_allocation[1]),
          rep(0,outcome.margins[2]-u_allocation[2]),rep(1,u_allocation[2]),
          rep(0,outcome.margins[3]-u_allocation[3]),rep(1,u_allocation[3]))
alpha = c(0,1,2.5)
beta = c(0,1,2)
score = outer(alpha,beta)
transform.fun = function(tb){
  return(sum(score*tb))
}



Experiment.3.1.result = run_experiment(data.table=data.table,gamma=gamma,delta = delta,u_allocation = u_allocation,u.vec=u.vec,transform.fun = transform.fun)

p3.1 = plot_convergence(result = Experiment.3.1.result,ylb=0, yub=0.25)

p3.1

ggsave("plot_convergence_experiment_3_1.png", plot = p3.1, width = 8, height = 6, dpi = 1200)


experiment.3.1.result.no.progression <- Experiment.3.1.result[!grepl("progression", names(Experiment.3.1.result))]

experiment.3.1.result.no.progression





########################### Experiment 3 - 2 ##################################

gc()
set.seed(100)

data.table = matrix(data=c(10,29,20,8,11,24,1,3,6),ncol=3)

data.table


delta = c(0,1,1)
gamma =0.5
u_allocation = c(0,30,10)
outcome.margins = colSums(data.table)

u.vec = c(rep(0,outcome.margins[1]-u_allocation[1]),rep(1,u_allocation[1]),
          rep(0,outcome.margins[2]-u_allocation[2]),rep(1,u_allocation[2]),
          rep(0,outcome.margins[3]-u_allocation[3]),rep(1,u_allocation[3]))
alpha = c(0,1,2.5)
beta = c(0,1,2)
score = outer(alpha,beta)
transform.fun = function(tb){
  return(sum(score*tb))
}



Experiment.3.2.result = run_experiment(data.table=data.table,gamma=gamma,delta = delta,u_allocation = u_allocation,u.vec=u.vec,transform.fun = transform.fun)

p3.2 = plot_convergence(result = Experiment.3.2.result,ylb = 0, yub = 0.25)

p3.2
ggsave("plot_convergence_experiment_3_2.png", plot = p3.2, width = 8, height = 6, dpi = 1200)


experiment.3.2.result.no.progression <- Experiment.3.2.result[!grepl("progression", names(Experiment.3.2.result))]

experiment.3.2.result.no.progression



##################### Experiment 3 - 3 #########################################

gc()
set.seed(100)

data.table = matrix(data=c(10,29,20,8,11,24,1,3,6),ncol=3)

data.table


delta = c(0,1,1)
gamma =1
u_allocation = c(0,20,10)
outcome.margins = colSums(data.table)

u.vec = c(rep(0,outcome.margins[1]-u_allocation[1]),rep(1,u_allocation[1]),
          rep(0,outcome.margins[2]-u_allocation[2]),rep(1,u_allocation[2]),
          rep(0,outcome.margins[3]-u_allocation[3]),rep(1,u_allocation[3]))
alpha = c(0,1,2.5)
beta = c(0,1,2)



Experiment.3.3.result = run_experiment(data.table=data.table,gamma=gamma,delta = delta,u_allocation = u_allocation,u.vec=u.vec,transform.fun = transform.fun)

p3.3 = plot_convergence(result = Experiment.3.3.result)

p3.3

ggsave("plot_convergence_experiment_3_3.png", plot = p3.3, width = 8, height = 6, dpi = 300)


experiment.3.3.result.no.progression <- Experiment.3.3.result[!grepl("progression", names(Experiment.3.3.result))]

experiment.3.3.result.no.progression













########################################################################################
#                                                                                      #
#       This file presents simulations related to the paper:                           #
#       "Exact Sensitivity Analysis for Observational Studies of Contingency Tables"   #
#       Lemma 1, and Appendix E                                                        #
#                                                                                      #
#       This file shows that the maximizer of an ordinal test might not                #
#       live at the corners of [0,1]^N given some other sensitivity model when         #
#       the outcome and the treatment are both non-binary. However, the maximizer      #
#       of a sign-score test (special case of ordinal test when the outcome is binary) #
#       the maximizer lives at the corners of [0,1]^N for either generic bias model or #                                                                       #
#       exposure-dose sensitivity model. See also Theorem 3 and its proof in Appendix  #
#       B.                                                                             #
#                                                                                      #
########################################################################################


###########################################################
## Setups and write a general objective function          #
###########################################################

# Install and load the arrangements package if not already:
# install.packages("arrangements")
library(arrangements)

permute = unique(permutations(c(1,1,2,2,2,3)))

## consider a fixed response vector with two responses, we know the solution is 
## always at the corner, we just use this to check that this algorithm runs fine 
R = c(1,1,2,2,2,2)

## here prints out all possible critical value 
min.critical = 0
max.critical = 0

total.treatment = dim(permute)[1]
for(trt in 1:total.treatment){
  candidate  = permute[trt,]%*%R
  if(candidate>max.critical){
    max.critical = candidate
  }
  if(candidate<=min.critical){
    min.critical = candidate
  }
  
}
##print out the critical values 
min.critical
max.critical


sig.fun = function(R,permute, u, gamma, critical.value){
  ## compute the denominator 
  denom = 0
  numer = 0
  total.treatment = dim(permute)[1]
  for(trt in 1:total.treatment){
    denom = denom + as.numeric(exp(gamma*(permute[trt,]%*%u)))
    
    ## see if the indicator function is one
    indicator = ifelse(permute[trt,]%*%R>=critical.value,1,0)
    numer = numer + indicator*as.numeric(exp(gamma*(permute[trt,]%*%u)))
    
  }
  return(numer/denom)
  
}

obj_fun <- function(u, R, permute, gamma, critical.value) {
  # We want to MAXIMIZE sig.fun, 
  # but optim() MINIMIZES by default.
  # So we return the negative value.
  val <- sig.fun(R, permute, u, gamma, critical.value)
  return(- val)
}


##################################################################################################
##################################################################################################
####################### Start running the experiment #############################################




##################################################################################
# Example 1: When the outcome is binary, and suppose an exposure-dose sensitivity, 
#                  the maximizer is at the corner. 
#
#
################################################################################## 



# Suppose we have an 6-element u to be optimized, each constrained to [0, 1]:
u_init <- rep(0.5, 6)  # initial guess inside [0,1]

opt_result <- optim(
  par     = u_init,
  fn      = obj_fun,             # Our objective function
  R       = R,                   # Passed to obj_fun
  permute = permute,             # Passed to obj_fun
  gamma   = 2,                  # Passed to obj_fun
  critical.value = 19,           # Passed to obj_fun
  method  = "L-BFGS-B",          # Required for bounding
  lower   = rep(0, 6),           # Lower bound for each element of u
  upper   = rep(1, 6),           # Upper bound for each element of u
  control = list(
    trace  = 1,                  # Print iteration details
    REPORT = 1                   # Print info every 1 iteration
  )
)
# Check results, and verify that the entries are either zero or one 

opt_result$par




###############################################################################
#
#          Example 2:When the outcome is binary, the maximizer is always at the
#          corner given another critical value.
#
#################################################################################

# Suppose we have an 8-element u to be optimized, each constrained to [0, 1]:
u_init <- rep(0.5, 6)  # initial guess inside [0,1]

opt_result <- optim(
  par     = u_init,
  fn      = obj_fun,             # Our objective function
  R       = R,                   # Passed to obj_fun
  permute = permute,             # Passed to obj_fun
  gamma   = 3,                  # Passed to obj_fun
  critical.value = 20,           # Passed to obj_fun
  method  = "L-BFGS-B",          # Required for bounding
  lower   = rep(0, 6),           # Lower bound for each element of u
  upper   = rep(1, 6),           # Upper bound for each element of u
  control = list(
    trace  = 1,                  # Print iteration details
    REPORT = 1                   # Print info every 1 iteration
  )
)
# Check results, and make sure that each entry is either zero or one 
opt_result$par


##############################################################################################
####### Here we showed some counter examples that the maximizer is no longer at the boundary.  
##############################################################################################
#####

###############################################################################
##                                                                            #
##                                                                            #
##           Example 3: Non-binary and integer outcomes                       #
##                      The maximizer will no longer be at the                #
##                      the corner of [0,1]^N                                 #
##                                                                            #
##                                                                            #
###############################################################################

permute = unique(permutations(c(1,1,2,2,2,3)))

## consider a fixed response vector with three responses 
R = c(1,2,2,3,3,4)

## here prints out all possible critical value 
min.critical = 0
max.critical = 0

total.treatment = dim(permute)[1]
for(trt in 1:total.treatment){
  candidate  = permute[trt,]%*%R
  if(candidate>max.critical){
    max.critical = candidate
  }
  if(candidate<=min.critical){
    min.critical = candidate
  }
  
}
##print out the critical values 
min.critical
max.critical


# Suppose we have an 6-element u to be optimized, each constrained to [0, 1]:
u_init <- rep(0.5, 6)  # initial guess inside [0,1]

opt_result <- optim(
  par     = u_init,
  fn      = obj_fun,             # Our objective function
  R       = R,                   # Passed to obj_fun
  permute = permute,             # Passed to obj_fun
  gamma   = 2,                  # Passed to obj_fun
  critical.value = 31,           # Passed to obj_fun
  method  = "L-BFGS-B",          # Required for bounding
  lower   = rep(0, 6),           # Lower bound for each element of u
  upper   = rep(1, 6),           # Upper bound for each element of u
  control = list(
    trace  = 1,                  # Print iteration details
    REPORT = 1                   # Print info every 1 iteration
  )
)
# Check results
opt_result$par




################################################################################
#
#         Example 4, change the critical value and gamma
#
###############################################################################

permute = unique(permutations(c(1,1,2,2,2,3)))

## consider a fixed response vector with three responses 
R = c(1,2,2,3,3,4)

# Suppose we have an 6-element u to be optimized, each constrained to [0, 1]:
u_init <- rep(0.5, 6)  # initial guess inside [0,1]

opt_result <- optim(
  par     = u_init,
  fn      = obj_fun,             # Our objective function
  R       = R,                   # Passed to obj_fun
  permute = permute,             # Passed to obj_fun
  gamma   = 2.5,                  # Passed to obj_fun
  critical.value = 30,           # Passed to obj_fun
  method  = "L-BFGS-B",          # Required for bounding
  lower   = rep(0, 6),           # Lower bound for each element of u
  upper   = rep(1, 6),           # Upper bound for each element of u
  control = list(
    trace  = 1,                  # Print iteration details
    REPORT = 1                   # Print info every 1 iteration
  )
)
# Check results
opt_result$par





###############################################################################
#
#
#                   Example 5: non-binary and non-integer outcomes
#                  This is the example we include in the Appendix E
#
###############################################################################


# Install and load the arrangements package if not already:
# install.packages("arrangements")
# library(arrangements)

## four different treatments 
permute = unique(permutations(c(1,1,2,2,3,4)))

## four outcome levels, and suppose we assign them some scores 
## consider a fixed response vector with three responses 
R = c(1.4,1.4,2.1,2.1,3.5,4.7)

## here prints out all possible critical value 
min.critical = 0
max.critical = 0

total.treatment = dim(permute)[1]
for(trt in 1:total.treatment){
  candidate  = permute[trt,]%*%R
  if(candidate>max.critical){
    max.critical = candidate
  }
  if(candidate<=min.critical){
    min.critical = candidate
  }
  
}
##print out the critical values 
min.critical
max.critical

# Suppose we have an 8-element u to be optimized, each constrained to [0, 1]:
u_init <- rep(0.5, 6)  # initial guess inside [0,1]

opt_result <- optim(
  par     = u_init,
  fn      = obj_fun,             # Our objective function
  R       = R,                   # Passed to obj_fun
  permute = permute,             # Passed to obj_fun
  gamma   = 2,                  # Passed to obj_fun
  critical.value = 40,           # Passed to obj_fun
  method  = "L-BFGS-B",          # Required for bounding
  lower   = rep(0, 6),           # Lower bound for each element of u
  upper   = rep(1, 6),           # Upper bound for each element of u
  control = list(
    trace  = 1,                  # Print iteration details
    REPORT = 1                   # Print info every 1 iteration
  )
)
# Check results
opt_result$par


## the following verifies if the last example's found maximizer really dominates all of the $64$ corner solutions.  

# Generate all binary vectors of length 6
binary_vectors <- expand.grid(rep(list(c(0, 1)), 6))

binary_vectors_matrix <- as.matrix(binary_vectors)

## and the following evaluates the probabilities for the 64 binary vectors 
for(index in 1:64){
  binary.vec = binary_vectors_matrix[index,]
  p.value = obj_fun(u=binary.vec,R=R,permute=permute, gamma = 2, critical.value = 40)
  print(-p.value)
  if(-p.value>=0.098622){
    cat("index", index)
  }
}



## here check the p.value under the found maximizer again 
non.corner.p.value = -(obj_fun(u=opt_result$par, R=R, permute=permute, gamma=2, critical.value = 40))
cat("this is a p value not happening at the corner,\\
    but generates a pvalue larger than all the corners,", non.corner.p.value)




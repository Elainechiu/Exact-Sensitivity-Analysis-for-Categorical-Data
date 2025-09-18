############################################################################################################
#                                                                                                          #
# Exact Sensitivity Analysis for Observational Studies of Contingency Tables                               #
#                                                                                                          #
# This R file contains the implementation of methods described in the paper:                               #
# "Exact Sensitivity Analysis for Observational Studies of I×J Contingency Tables"                         #
#                                                                                                          #
# The code provides three complementary approaches for conducting sensitivity                              #
# analysis in observational studies with multiple treatments and multiple outcomes:                        #
#                                                                                                          #
# 1. EXACT METHODS: Compute exact p-values under the generic bias model                                    #
#    - Handles tables up to 5×5 dimensions                                                                 #
#    - Provides exact inference for both general test statistics and ordinal tests.                        #
#.   - See Section 4.                                                                                      #
#                                                                                                          #
# 2. SAMPLING-BASED METHODS: Monte Carlo approximation for larger tables                                   #
#    - Uses Sequential Importance Sampling (SIS) when exact computation is infeasible                      #
#    - Provides accurate approximations with computational efficiency                                      #
#    - See Section 4.                                                                                      #
#                                                                                                          #
# 3. NORMAL APPROXIMATION: Asymptotic methods for large sample sizes                                       #
#    - Based on multivariate normal approximations                                                         #
#    - Computationally efficient for very large tables                                                     #
#    - See Section 5.                                                                                      #
#                                                                                                          #
# The methods implement the generic bias model which allows for flexible                                   #
# sensitivity analysis in observational studies with multiple treatment levels,                            #
# accounting for potential unmeasured confounding.                                                         #
#                                                                                                          #
# Dependencies:                                                                                            #
# - Rcpp (for C++ integration)                                                                             #
# - SIS_methods.cpp (for table sampling)                                                                   #
# - sensitivityIxJ.cpp (for core computations)                                                             #
#                                                                                                          #
#                                                                                                          #
############################################################################################################





require(Rcpp)


###############################################################################
#
#        Auxiliary Functions 
#
#
###############################################################################


#' Generate All Possible Contingency Tables Exceeding or Not Exceeding a Threshold
#'
#' Given an observed contingency table (with fixed margins) and a user-defined
#' function that computes a test statistic, this function enumerates all
#' valid contingency tables that meet a particular criterion (i.e., those
#' with `T(table) >= threshold` or `T(table) <= threshold`).
#'
#' @param threshold A numeric value representing the cutoff for the test statistic.
#' @param table A matrix or table object specifying the observed contingency table.
#'   The function will preserve the row/column sums (margins) of this table when
#'   generating possible tables.
#' @param direction A character string, either \code{"greater than"} or \code{"less than"},
#'   indicating whether to return tables whose test statistic is \code{>= threshold}
#'   or \code{<= threshold}. Defaults to \code{c("greater than", "less than")} but
#'   you should pass in one value explicitly.
#' @param transform.fun A user-defined function of the form \code{f(tbl)}, which takes
#'   a contingency table (as a matrix or table) and returns a numeric test statistic.
#'
#' @details
#' The function systematically iterates over all valid ways to fill a contingency
#' table of dimension \code{I x J} such that row and column sums match the original
#' \code{table}'s margins. Only those tables that satisfy \code{transform.fun(tbl) >= threshold}
#' (for \code{direction = "greater than"}) or \code{transform.fun(tbl) <= threshold}
#' (for \code{direction = "less than"}) are returned.
#'
#' Currently, this function supports certain table dimensions (\code{I} up to 5
#' when \code{J=2}, etc.). If the dimension is not supported, it returns
#' \code{NULL} with a warning.
#'
#' @return A list of all contingency tables (each table is returned as a
#' \code{table} object) that meet the specified \code{threshold} criterion.
#' If no tables match or the dimension is not supported, \code{NULL} is returned.
#'
#' @examples
#' # Suppose we have a 3x3 table:
#' obs_table <- matrix(c(5, 3, 2, 6, 11, 7, 3, 0,3), ncol = 3, byrow = TRUE)
#' obs_table
#'
#' # Define a simple test statistic function (e.g., sum of diagonal)
#' diag_sum <- function(tbl) sum(diag(tbl))
#'
#' # Find all 3x3 tables with the same margins where the diagonal sum >= 19
#' result <- possible.table(threshold = 19, table = obs_table,
#'                         direction = "greater than", transform.fun = diag_sum)
#' result[[1]]         # Inspect the first matching table
#'

#' @export
possible.table <- function(threshold, table, direction = c("greater than", "less than"), transform.fun) {
  # Create an empty list to store valid tables
  table.list <- list()
  
  # Detect the margins of the tables
  row.margins <- rowSums(table)
  col.margins <- colSums(table)
  I <- dim(table)[1]
  J <- dim(table)[2]
  N <- sum(table)
  
  # Handle 'greater than' direction
  if(direction == "greater than") {
    if(I == 2 & J == 2) {
      n1 <- row.margins[1]
      n2 <- row.margins[2]
      m1 <- col.margins[1]
      m2 <- col.margins[2]
      
      n11_lb <- max(0, m1 + n1 - N)
      n11_ub <- min(m1, n1)
      
      for(n11 in n11_lb:n11_ub) {
        sampled.table <- as.table(matrix(data = c(n11, n1 - n11, m1 - n11, N - m1 - n1 + n11),
                                         ncol = 2, byrow = TRUE))
        if(all(sampled.table >= 0)) {
          if(transform.fun(sampled.table) >= threshold) {
            table.list <- append(table.list, list(sampled.table))
          }
        }
      }
      return(table.list)
      
    } else if(I == 2 & J == 3) {
      n1 <- row.margins[1]
      n2 <- row.margins[2]
      m1 <- col.margins[1]
      m2 <- col.margins[2]
      
      n11_lb <- max(0, m1 + n1 - N)
      n11_ub <- min(m1, n1)
      n12_lb <- max(0, m2 + n1 - N)
      n12_ub <- min(m2, n1)
      
      for(n11 in n11_lb:n11_ub) {
        for(n12 in n12_lb:n12_ub) {
          counts <- c(n11, n12, n1 - n11 - n12,
                      m1 - n11, m2 - n12, N - n1 - m1 - m2 + n11 + n12)
          sampled.table <- as.table(matrix(data = counts, ncol = 3, byrow = TRUE))
          if(all(sampled.table >= 0)) {
            if(transform.fun(sampled.table) >= threshold) {
              table.list <- append(table.list, list(sampled.table))
            }
          }
        }
      }
      return(table.list)
      
    } else if(I == 2 & J == 4) {
      n1 <- row.margins[1]
      n2 <- row.margins[2]
      m1 <- col.margins[1]
      m2 <- col.margins[2]
      m3 <- col.margins[3]
      m4 <- col.margins[4]
      
      n11_lb <- max(0, m1 + n1 - N)
      n11_ub <- min(m1, n1)
      n12_lb <- max(0, m2 + n1 - N)
      n12_ub <- min(m2, n1)
      n13_lb <- max(0, m3 + n1 - N)
      n13_ub <- min(m3, n1)
      
      for(n11 in n11_lb:n11_ub) {
        for(n12 in n12_lb:n12_ub) {
          for(n13 in n13_lb:n13_ub) {
            counts <- c(n11, n12, n13, n1 - n11 - n12 - n13,
                        m1 - n11, m2 - n12, m3 - n13, N - m1 - m2 - m3 - n1 + n11 + n12 + n13)
            sampled.table <- as.table(matrix(data = counts, ncol = 4, byrow = TRUE))
            if(all(sampled.table >= 0)) {
              if(transform.fun(sampled.table) >= threshold) {
                table.list <- append(table.list, list(sampled.table))
              }
            }
          }
        }
      }
      return(table.list)
      
    } else if(I == 2 & J == 5) {
      n1 <- row.margins[1]
      n2 <- row.margins[2]
      m1 <- col.margins[1]
      m2 <- col.margins[2]
      m3 <- col.margins[3]
      m4 <- col.margins[4]
      m5 <- col.margins[5]
      
      n11_lb <- max(0, m1 + n1 - N)
      n11_ub <- min(m1, n1)
      n12_lb <- max(0, m2 + n1 - N)
      n12_ub <- min(m2, n1)
      n13_lb <- max(0, m3 + n1 - N)
      n13_ub <- min(m3, n1)
      n14_lb <- max(0, m4 + n1 - N)
      n14_ub <- min(m4, n1)
      
      for(n11 in n11_lb:n11_ub) {
        for(n12 in n12_lb:n12_ub) {
          for(n13 in n13_lb:n13_ub) {
            for(n14 in n14_lb:n14_ub) {
              counts <- c(n11, n12, n13, n14, n1 - n11 - n12 - n13 - n14,
                          m1 - n11, m2 - n12, m3 - n13, m4 - n14,
                          N - m1 - m2 - m3 - m4 - n1 + n11 + n12 + n13 + n14)
              sampled.table <- as.table(matrix(data = counts, ncol = 5, byrow = TRUE))
              if(all(sampled.table >= 0)) {
                if(transform.fun(sampled.table) >= threshold) {
                  table.list <- append(table.list, list(sampled.table))
                }
              }
            }
          }
        }
      }
      return(table.list)
      
    } else if(I == 3 & J == 2) {
      n1 <- row.margins[1]
      n2 <- row.margins[2]
      n3 <- row.margins[3]
      m1 <- col.margins[1]
      m2 <- col.margins[2]
      
      n11_lb <- max(0, m1 + n1 - N)
      n11_ub <- min(m1, n1)
      n21_lb <- max(0, m1 + n2 - N)
      n21_ub <- min(m1, n2)
      
      for(n11 in n11_lb:n11_ub) {
        for(n21 in n21_lb:n21_ub) {
          counts <- c(n11, n1 - n11, n21, n2 - n21, m1 - n11 - n21, N - n1 - n2 + n21 - m1 + n11)
          sampled.table <- as.table(matrix(data = counts, ncol = 2, byrow = TRUE))
          if(all(sampled.table >= 0)) {
            if(transform.fun(sampled.table) >= threshold) {
              table.list <- append(table.list, list(sampled.table))
            }
          }
        }
      }
      return(table.list)
      
    } else if(I == 3 & J == 3) {
      n1 <- row.margins[1]
      n2 <- row.margins[2]
      n3 <- row.margins[3]
      m1 <- col.margins[1]
      m2 <- col.margins[2]
      m3 <- col.margins[3]
      
      n11_lb <- max(0, m1 + n1 - N)
      n11_ub <- min(m1, n1)
      n12_lb <- max(0, m2 + n1 - N)
      n12_ub <- min(m2, n1)
      n21_lb <- max(0, m1 + n2 - N)
      n21_ub <- min(m1, n2)
      n22_lb <- max(0, m2 + n2 - N)
      n22_ub <- min(m2, n2)
      
      for(n11 in n11_lb:n11_ub) {
        for(n12 in n12_lb:n12_ub) {
          for(n21 in n21_lb:n21_ub) {
            for(n22 in n22_lb:n22_ub) {
              counts <- c(n11, n12, n1 - n11 - n12, n21, n22, n2 - n21 - n22,
                          m1 - n11 - n21, m2 - n12 - n22, N - n1 - n2 - m1 + n11 + n21 - m2 + n12 + n22)
              sampled.table <- as.table(matrix(data = counts, ncol = 3, byrow = TRUE))
              if(all(sampled.table >= 0)) {
                if(transform.fun(sampled.table) >= threshold) {
                  table.list <- append(table.list, list(sampled.table))
                }
              }
            }
          }
        }
      }
      return(table.list)
      
    } else if(I == 3 & J == 4) {
      n1 <- row.margins[1]
      n2 <- row.margins[2]
      n3 <- row.margins[3]
      m1 <- col.margins[1]
      m2 <- col.margins[2]
      m3 <- col.margins[3]
      m4 <- col.margins[4]
      
      n11_lb <- max(0, m1 + n1 - N)
      n11_ub <- min(m1, n1)
      n12_lb <- max(0, m2 + n1 - N)
      n12_ub <- min(m2, n1)
      n13_lb <- max(0, m3 + n1 - N)
      n13_ub <- min(m3, n1)
      n21_lb <- max(0, m1 + n2 - N)
      n21_ub <- min(m1, n2)
      n22_lb <- max(0, m2 + n2 - N)
      n22_ub <- min(m2, n2)
      n23_lb <- max(0, m3 + n2 - N)
      n23_ub <- min(m3, n2)
      
      for(n11 in n11_lb:n11_ub) {
        for(n12 in n12_lb:n12_ub) {
          for(n13 in n13_lb:n13_ub) {
            for(n21 in n21_lb:n21_ub) {
              for(n22 in n22_lb:n22_ub) {
                for(n23 in n23_lb:n23_ub) {
                  counts <- c(n11, n12, n13, n1 - n11 - n12 - n13,
                              n21, n22, n23, n2 - n21 - n22 - n23,
                              m1 - n11 - n21, m2 - n12 - n22, m3 - n13 - n23, m4-(n1-n11-n12-n13)-(n2-n21-n22-n23))
                  sampled.table <- as.table(matrix(data = counts, ncol = 4, byrow = TRUE))
                  if(all(sampled.table >= 0)) {
                    if(transform.fun(sampled.table) >= threshold) {
                      table.list <- append(table.list, list(sampled.table))
                    }
                  }
                }
              }
            }
          }
        }
      }
      return(table.list)
      
    } else if(I == 3 & J == 5) {
      n1 <- row.margins[1]
      n2 <- row.margins[2]
      n3 <- row.margins[3]
      m1 <- col.margins[1]
      m2 <- col.margins[2]
      m3 <- col.margins[3]
      m4 <- col.margins[4]
      m5 <- col.margins[5]
      
      n11_lb <- max(0, m1 + n1 - N)
      n11_ub <- min(m1, n1)
      n12_lb <- max(0, m2 + n1 - N)
      n12_ub <- min(m2, n1)
      n13_lb <- max(0, m3 + n1 - N)
      n13_ub <- min(m3, n1)
      n14_lb <- max(0, m4 + n1 - N)
      n14_ub <- min(m4, n1)
      n21_lb <- max(0, m1 + n2 - N)
      n21_ub <- min(m1, n2)
      n22_lb <- max(0, m2 + n2 - N)
      n22_ub <- min(m2, n2)
      n23_lb <- max(0, m3 + n2 - N)
      n23_ub <- min(m3, n2)
      n24_lb <- max(0, m4 + n2 - N)
      n24_ub <- min(m4, n2)
      
      for(n11 in n11_lb:n11_ub) {
        for(n12 in n12_lb:n12_ub) {
          for(n13 in n13_lb:n13_ub) {
            for(n14 in n14_lb:n14_ub) {
              for(n21 in n21_lb:n21_ub) {
                for(n22 in n22_lb:n22_ub) {
                  for(n23 in n23_lb:n23_ub) {
                    for(n24 in n24_lb:n24_ub) {
                      counts <- c(n11, n12, n13, n14, n1 - n11 - n12 - n13 - n14,
                                  n21, n22, n23, n24, n2 - n21 - n22 - n23 - n24,
                                  m1 - n11 - n21, m2 - n12 - n22, m3 - n13 - n23,
                                  m4 - n14 - n24, m5 - n1 - n2 + n11 + n12 + n13 + n14 + n21 + n22 + n23 + n24)
                      sampled.table <- as.table(matrix(data = counts, ncol = 5, byrow = TRUE))
                      if(all(sampled.table >= 0)) {
                        if(transform.fun(sampled.table) >= threshold) {
                          table.list <- append(table.list, list(sampled.table))
                        }
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
      return(table.list)
      
    } else if(I == 4 & J == 2) {
      n1 <- row.margins[1]
      n2 <- row.margins[2]
      n3 <- row.margins[3]
      n4 <- row.margins[4]
      m1 <- col.margins[1]
      m2 <- col.margins[2]
      
      n11_lb <- max(0, m1 + n1 - N)
      n11_ub <- min(m1, n1)
      n21_lb <- max(0, m1 + n2 - N)
      n21_ub <- min(m1, n2)
      n31_lb <- max(0, m1 + n3 - N)
      n31_ub <- min(m1, n3)
      
      for(n11 in n11_lb:n11_ub) {
        for(n21 in n21_lb:n21_ub) {
          for(n31 in n31_lb:n31_ub) {
            counts <- c(n11, n1 - n11, n21, n2 - n21, n31, n3 - n31,
                        m1 - n11 - n21 - n31, n4 - m1 + n11 + n21 + n31)
            sampled.table <- as.table(matrix(data = counts, ncol = 2, byrow = TRUE))
            if(all(sampled.table >= 0)) {
              if(transform.fun(sampled.table) >= threshold) {
                table.list <- append(table.list, list(sampled.table))
              }
            }
          }
        }
      }
      return(table.list)
      
    } else if(I == 4 & J == 3) {
      n1 <- row.margins[1]
      n2 <- row.margins[2]
      n3 <- row.margins[3]
      n4 <- row.margins[4]
      m1 <- col.margins[1]
      m2 <- col.margins[2]
      m3 <- col.margins[3]
      
      n11_lb <- max(0, m1 + n1 - N)
      n11_ub <- min(m1, n1)
      n21_lb <- max(0, m1 + n2 - N)
      n21_ub <- min(m1, n2)
      n31_lb <- max(0, m1 + n3 - N)
      n31_ub <- min(m1, n3)
      n12_lb <- max(0, m2 + n1 - N)
      n12_ub <- min(m2, n1)
      n22_lb <- max(0, m2 + n2 - N)
      n22_ub <- min(m2, n2)
      n32_lb <- max(0, m2 + n3 - N)
      n32_ub <- min(m2, n3)
      
      for(n11 in n11_lb:n11_ub) {
        for(n21 in n21_lb:n21_ub) {
          for(n31 in n31_lb:n31_ub) {
            for(n12 in n12_lb:n12_ub) {
              for(n22 in n22_lb:n22_ub) {
                for(n32 in n32_lb:n32_ub) {
                  counts <- c(n11, n12, n1 - n11 - n12,
                              n21, n22, n2 - n21 - n22,
                              n31, n32, n3 - n31 - n32,
                              m1 - n11 - n21 - n31, m2 - n12 - n22 - n32,
                              n4 - m1 - m2 + n11 + n21 + n31 + n12 + n22 + n32)
                  sampled.table <- as.table(matrix(data = counts, ncol = 3, byrow = TRUE))
                  if(all(sampled.table >= 0)) {
                    if(transform.fun(sampled.table) >= threshold) {
                      table.list <- append(table.list, list(sampled.table))
                    }
                  }
                }
              }
            }
          }
        }
      }
      return(table.list)
      
    } else if(I == 5 & J == 2) {
      n1 <- row.margins[1]
      n2 <- row.margins[2]
      n3 <- row.margins[3]
      n4 <- row.margins[4]
      n5 <- row.margins[5]
      m1 <- col.margins[1]
      m2 <- col.margins[2]
      
      n11_lb <- max(0, m1 + n1 - N)
      n11_ub <- min(m1, n1)
      n21_lb <- max(0, m1 + n2 - N)
      n21_ub <- min(m1, n2)
      n31_lb <- max(0, m1 + n3 - N)
      n31_ub <- min(m1, n3)
      n41_lb <- max(0, m1 + n4 - N)
      n41_ub <- min(m1, n4)
      
      for(n11 in n11_lb:n11_ub) {
        for(n21 in n21_lb:n21_ub) {
          for(n31 in n31_lb:n31_ub) {
            for(n41 in n41_lb:n41_ub) {
              counts <- c(n11, n1 - n11, n21, n2 - n21,
                          n31, n3 - n31, n41, n4 - n41,
                          m1 - n11 - n21 - n31 - n41, n5 - m1 + n11 + n21 + n31 + n41)
              sampled.table <- as.table(matrix(data = counts, ncol = 2, byrow = TRUE))
              if(all(sampled.table >= 0)) {
                if(transform.fun(sampled.table) >= threshold) {
                  table.list <- append(table.list, list(sampled.table))
                }
              }
            }
          }
        }
      }
      return(table.list)
      
    } else {
      warning("Your table dimension is not supported in 'greater than'. Do you want to take a transpose?\n")
      return(NULL)
    }
  }
  
  # Handle 'less than' direction
  if(direction == "less than") {
    if(I == 2 & J == 2) {
      n1 <- row.margins[1]
      n2 <- row.margins[2]
      m1 <- col.margins[1]
      m2 <- col.margins[2]
      
      n11_lb <- max(0, m1 + n1 - N)
      n11_ub <- min(m1, n1)
      
      for(n11 in n11_lb:n11_ub) {
        sampled.table <- as.table(matrix(data = c(n11, n1 - n11, m1 - n11, N - m1 - n1 + n11),
                                         ncol = 2, byrow = TRUE))
        if(all(sampled.table >= 0)) {
          if(transform.fun(sampled.table) <= threshold) {
            table.list <- append(table.list, list(sampled.table))
          }
        }
      }
      return(table.list)
      
    } else if(I == 2 & J == 3) {
      n1 <- row.margins[1]
      n2 <- row.margins[2]
      m1 <- col.margins[1]
      m2 <- col.margins[2]
      
      n11_lb <- max(0, m1 + n1 - N)
      n11_ub <- min(m1, n1)
      n12_lb <- max(0, m2 + n1 - N)
      n12_ub <- min(m2, n1)
      
      for(n11 in n11_lb:n11_ub) {
        for(n12 in n12_lb:n12_ub) {
          counts <- c(n11, n12, n1 - n11 - n12,
                      m1 - n11, m2 - n12, N - n1 - m1 - m2 + n11 + n12)
          sampled.table <- as.table(matrix(data = counts, ncol = 3, byrow = TRUE))
          if(all(sampled.table >= 0)) {
            if(transform.fun(sampled.table) <= threshold) {
              table.list <- append(table.list, list(sampled.table))
            }
          }
        }
      }
      return(table.list)
      
    } else if(I == 2 & J == 4) {
      n1 <- row.margins[1]
      n2 <- row.margins[2]
      m1 <- col.margins[1]
      m2 <- col.margins[2]
      m3 <- col.margins[3]
      m4 <- col.margins[4]
      
      n11_lb <- max(0, m1 + n1 - N)
      n11_ub <- min(m1, n1)
      n12_lb <- max(0, m2 + n1 - N)
      n12_ub <- min(m2, n1)
      n13_lb <- max(0, m3 + n1 - N)
      n13_ub <- min(m3, n1)
      
      for(n11 in n11_lb:n11_ub) {
        for(n12 in n12_lb:n12_ub) {
          for(n13 in n13_lb:n13_ub) {
            counts <- c(n11, n12, n13, n1 - n11 - n12 - n13,
                        m1 - n11, m2 - n12, m3 - n13, N - m1 - m2 - m3 - n1 + n11 + n12 + n13)
            sampled.table <- as.table(matrix(data = counts, ncol = 4, byrow = TRUE))
            if(all(sampled.table >= 0)) {
              if(transform.fun(sampled.table) <= threshold) {
                table.list <- append(table.list, list(sampled.table))
              }
            }
          }
        }
      }
      return(table.list)
      
    } else if(I == 2 & J == 5) {
      n1 <- row.margins[1]
      n2 <- row.margins[2]
      m1 <- col.margins[1]
      m2 <- col.margins[2]
      m3 <- col.margins[3]
      m4 <- col.margins[4]
      m5 <- col.margins[5]
      
      n11_lb <- max(0, m1 + n1 - N)
      n11_ub <- min(m1, n1)
      n12_lb <- max(0, m2 + n1 - N)
      n12_ub <- min(m2, n1)
      n13_lb <- max(0, m3 + n1 - N)
      n13_ub <- min(m3, n1)
      n14_lb <- max(0, m4 + n1 - N)
      n14_ub <- min(m4, n1)
      
      for(n11 in n11_lb:n11_ub) {
        for(n12 in n12_lb:n12_ub) {
          for(n13 in n13_lb:n13_ub) {
            for(n14 in n14_lb:n14_ub) {
              counts <- c(n11, n12, n13, n14, n1 - n11 - n12 - n13 - n14,
                          m1 - n11, m2 - n12, m3 - n13, m4 - n14,
                          N - m1 - m2 - m3 - m4 - n1 + n11 + n12 + n13 + n14)
              sampled.table <- as.table(matrix(data = counts, ncol = 5, byrow = TRUE))
              if(all(sampled.table >= 0)) {
                if(transform.fun(sampled.table) <= threshold) {
                  table.list <- append(table.list, list(sampled.table))
                }
              }
            }
          }
        }
      }
      return(table.list)
      
    } else if(I == 3 & J == 2) {
      n1 <- row.margins[1]
      n2 <- row.margins[2]
      n3 <- row.margins[3]
      m1 <- col.margins[1]
      m2 <- col.margins[2]
      
      n11_lb <- max(0, m1 + n1 - N)
      n11_ub <- min(m1, n1)
      n21_lb <- max(0, m1 + n2 - N)
      n21_ub <- min(m1, n2)
      
      for(n11 in n11_lb:n11_ub) {
        for(n21 in n21_lb:n21_ub) {
          counts <- c(n11, n1 - n11, n21, n2 - n21, m1 - n11 - n21, N - n1 - n2 + n21 - m1 + n11)
          sampled.table <- as.table(matrix(data = counts, ncol = 2, byrow = TRUE))
          if(all(sampled.table >= 0)) {
            if(transform.fun(sampled.table) <= threshold) {
              table.list <- append(table.list, list(sampled.table))
            }
          }
        }
      }
      return(table.list)
      
    } else if(I == 3 & J == 3) {
      n1 <- row.margins[1]
      n2 <- row.margins[2]
      n3 <- row.margins[3]
      m1 <- col.margins[1]
      m2 <- col.margins[2]
      m3 <- col.margins[3]
      
      n11_lb <- max(0, m1 + n1 - N)
      n11_ub <- min(m1, n1)
      n12_lb <- max(0, m2 + n1 - N)
      n12_ub <- min(m2, n1)
      n21_lb <- max(0, m1 + n2 - N)
      n21_ub <- min(m1, n2)
      n22_lb <- max(0, m2 + n2 - N)
      n22_ub <- min(m2, n2)
      
      for(n11 in n11_lb:n11_ub) {
        for(n12 in n12_lb:n12_ub) {
          for(n21 in n21_lb:n21_ub) {
            for(n22 in n22_lb:n22_ub) {
              counts <- c(n11, n12, n1 - n11 - n12,
                          n21, n22, n2 - n21 - n22,
                          m1 - n11 - n21, m2 - n12 - n22, N - n1 - n2 - m1 + n11 + n21 - m2 + n12 + n22)
              sampled.table <- as.table(matrix(data = counts, ncol = 3, byrow = TRUE))
              if(all(sampled.table >= 0)) {
                if(transform.fun(sampled.table) <= threshold) {
                  table.list <- append(table.list, list(sampled.table))
                }
              }
            }
          }
        }
      }
      return(table.list)
      
    } else if(I == 3 & J == 4) {
      n1 <- row.margins[1]
      n2 <- row.margins[2]
      n3 <- row.margins[3]
      m1 <- col.margins[1]
      m2 <- col.margins[2]
      m3 <- col.margins[3]
      m4 <- col.margins[4]
      
      n11_lb <- max(0, m1 + n1 - N)
      n11_ub <- min(m1, n1)
      n12_lb <- max(0, m2 + n1 - N)
      n12_ub <- min(m2, n1)
      n13_lb <- max(0, m3 + n1 - N)
      n13_ub <- min(m3, n1)
      n21_lb <- max(0, m1 + n2 - N)
      n21_ub <- min(m1, n2)
      n22_lb <- max(0, m2 + n2 - N)
      n22_ub <- min(m2, n2)
      n23_lb <- max(0, m3 + n2 - N)
      n23_ub <- min(m3, n2)
      
      for(n11 in n11_lb:n11_ub) {
        for(n12 in n12_lb:n12_ub) {
          for(n13 in n13_lb:n13_ub) {
            for(n21 in n21_lb:n21_ub) {
              for(n22 in n22_lb:n22_ub) {
                for(n23 in n23_lb:n23_ub) {
                  counts <- c(n11, n12, n13, n1 - n11 - n12 - n13,
                              n21, n22, n23, n2 - n21 - n22 - n23,
                              m1 - n11 - n21, m2 - n12 - n22, m3 - n13 - n23,
                              m4 - n1 + n11 + n12+n13-n2+n21+n22+n23)
                  sampled.table <- as.table(matrix(data = counts, ncol = 4, byrow = TRUE))
                  if(all(sampled.table >= 0)) {
                    if(transform.fun(sampled.table) <= threshold) {
                      table.list <- append(table.list, list(sampled.table))
                    }
                  }
                }
              }
            }
          }
        }
      }
      return(table.list)
      
    } else if(I == 3 & J == 5) {
      n1 <- row.margins[1]
      n2 <- row.margins[2]
      n3 <- row.margins[3]
      m1 <- col.margins[1]
      m2 <- col.margins[2]
      m3 <- col.margins[3]
      m4 <- col.margins[4]
      m5 <- col.margins[5]
      
      n11_lb <- max(0, m1 + n1 - N)
      n11_ub <- min(m1, n1)
      n12_lb <- max(0, m2 + n1 - N)
      n12_ub <- min(m2, n1)
      n13_lb <- max(0, m3 + n1 - N)
      n13_ub <- min(m3, n1)
      n14_lb <- max(0, m4 + n1 - N)
      n14_ub <- min(m4, n1)
      n21_lb <- max(0, m1 + n2 - N)
      n21_ub <- min(m1, n2)
      n22_lb <- max(0, m2 + n2 - N)
      n22_ub <- min(m2, n2)
      n23_lb <- max(0, m3 + n2 - N)
      n23_ub <- min(m3, n2)
      n24_lb <- max(0, m4 + n2 - N)
      n24_ub <- min(m4, n2)
      
      for(n11 in n11_lb:n11_ub) {
        for(n12 in n12_lb:n12_ub) {
          for(n13 in n13_lb:n13_ub) {
            for(n14 in n14_lb:n14_ub) {
              for(n21 in n21_lb:n21_ub) {
                for(n22 in n22_lb:n22_ub) {
                  for(n23 in n23_lb:n23_ub) {
                    for(n24 in n24_lb:n24_ub) {
                      counts <- c(n11, n12, n13, n14, n1 - n11 - n12 - n13 - n14,
                                  n21, n22, n23, n24, n2 - n21 - n22 - n23 - n24,
                                  m1 - n11 - n21, m2 - n12 - n22, m3 - n13 - n23,
                                  m4 - n14 - n24, m5 - n1 - n2 + n11 + n12 + n13 + n14 + n21 + n22 + n23 + n24)
                      sampled.table <- as.table(matrix(data = counts, ncol = 5, byrow = TRUE))
                      if(all(sampled.table >= 0)) {
                        if(transform.fun(sampled.table) <= threshold) {
                          table.list <- append(table.list, list(sampled.table))
                        }
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
      return(table.list)
      
    } else if(I == 4 & J == 2) {
      n1 <- row.margins[1]
      n2 <- row.margins[2]
      n3 <- row.margins[3]
      n4 <- row.margins[4]
      m1 <- col.margins[1]
      m2 <- col.margins[2]
      
      n11_lb <- max(0, m1 + n1 - N)
      n11_ub <- min(m1, n1)
      n21_lb <- max(0, m1 + n2 - N)
      n21_ub <- min(m1, n2)
      n31_lb <- max(0, m1 + n3 - N)
      n31_ub <- min(m1, n3)
      
      for(n11 in n11_lb:n11_ub) {
        for(n21 in n21_lb:n21_ub) {
          for(n31 in n31_lb:n31_ub) {
            counts <- c(n11, n1 - n11, n21, n2 - n21, n31, n3 - n31,
                        m1 - n11 - n21 - n31, n4 - m1 + n11 + n21 + n31)
            sampled.table <- as.table(matrix(data = counts, ncol = 2, byrow = TRUE))
            if(all(sampled.table >= 0)) {
              if(transform.fun(sampled.table) <= threshold) {
                table.list <- append(table.list, list(sampled.table))
              }
            }
          }
        }
      }
      return(table.list)
      
    } else if(I == 4 & J == 3) {
      n1 <- row.margins[1]
      n2 <- row.margins[2]
      n3 <- row.margins[3]
      n4 <- row.margins[4]
      m1 <- col.margins[1]
      m2 <- col.margins[2]
      m3 <- col.margins[3]
      
      n11_lb <- max(0, m1 + n1 - N)
      n11_ub <- min(m1, n1)
      n21_lb <- max(0, m1 + n2 - N)
      n21_ub <- min(m1, n2)
      n31_lb <- max(0, m1 + n3 - N)
      n31_ub <- min(m1, n3)
      n12_lb <- max(0, m2 + n1 - N)
      n12_ub <- min(m2, n1)
      n22_lb <- max(0, m2 + n2 - N)
      n22_ub <- min(m2, n2)
      n32_lb <- max(0, m2 + n3 - N)
      n32_ub <- min(m2, n3)
      
      for(n11 in n11_lb:n11_ub) {
        for(n21 in n21_lb:n21_ub) {
          for(n31 in n31_lb:n31_ub) {
            for(n12 in n12_lb:n12_ub) {
              for(n22 in n22_lb:n22_ub) {
                for(n32 in n32_lb:n32_ub) {
                  counts <- c(n11, n12, n1 - n11 - n12,
                              n21, n22, n2 - n21 - n22,
                              n31, n32, n3 - n31 - n32,
                              m1 - n11 - n21 - n31, m2 - n12 - n22 - n32,
                              n4 - m1 - m2 + n11 + n21 + n31 + n12 + n22 + n32)
                  sampled.table <- as.table(matrix(data = counts, ncol = 3, byrow = TRUE))
                  if(all(sampled.table >= 0)) {
                    if(transform.fun(sampled.table) <= threshold) {
                      table.list <- append(table.list, list(sampled.table))
                    }
                  }
                }
              }
            }
          }
        }
      }
      return(table.list)
      
    } else {
      warning("Your table dimension is not supported in 'less than'. Do you want to take a transpose?\n")
      return(NULL)
    }
  }
  
  # If direction is neither 'greater than' nor 'less than'
  warning("Invalid direction. Use 'greater than' or 'less than'.")
  return(NULL)
}




#-------------------------------------------------------------------------------------------------------------



is_binary_vector <- function(x) {
  all(x %in% c(0, 1))
}


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
# set, making them suitable for smaller table dimensions where the             
# combinatorial space is manageable.                                           
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






#' Compute the exact Probability of a Single Table for the Generic Bias Model
#'
#' This function computes the probability of a single contingency table under the generic bias model.
#' It supports tables with varying dimensions, provided they meet specific constraints.

generic.I.by.J.sensitivity.point.probability = function(table,row="treatment",u_allocation,gamma, delta,shared_divisor=1000000){
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




  # Validate 'gamma_delta'
  if (!is_scalar_numeric(gamma)) {
    stop("'gamma must be a scalar.")
  }
  
  if(!is_binary_vector(delta)){
    stop("delta must be a binary vector")
  }

  if (length(delta) != I) {
    stop(paste0("'gamma_delta' must have the same length as the number of treatments (", I, ")."))
  }


  # Validate 'u_allocation'
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
    stop("Each element of 'u_allocation' must be less than or equal to the corresponding outcome margin.")
  }

  gamma_delta = gamma*delta

  # if I==2 and J ==2, then use the specialized two by two table computation
  if((I==2) &(J==2)){
    n1 = treatment_margins[1]
    n2 = treatment_margins[2]
    Us = sum(u_allocation)
    denominator = denominator_two_treatment(n1=n1,n2=n2,gamma_delta=gamma_delta,N=N,Us=Us,shared_divisor=shared_divisor)
    n11 = table[1,1]
    m1 = outcome_margins[1]
    m2 = outcome_margins[2]
    u1 = u_allocation[1]
    u2 = u_allocation[2]
    numerator = d_numerator_two_by_two(n11=n11,m1=m1,m2=m2,n1=n1,n2=n2,N=N,gamma_delta=gamma_delta,u1=u1,u2=u2,shared_divisor=shared_divisor)
    return(numerator/denominator)
  }else if((I==3)&(J==2)){
    ## this is the case of three treatment and two outcomes
    n1 = treatment_margins[1]
    n2 = treatment_margins[2]
    n3 = treatment_margins[3]
    Us = sum(u_allocation)
    denominator = denominator_three_treatment(n1=n1,n2=n2,n3=n3,gamma_delta=gamma_delta,N=N,Us=Us,shared_divisor=shared_divisor)
    n11 = table[1,1]
    n21 = table[2,1]
    m1 = outcome_margins[1]
    m2 = outcome_margins[2]
    u1 = u_allocation[1]
    u2 = u_allocation[2]
    numerator = d_numerator_three_by_two(n11=n11,n21=n21,m1=m1,m2=m2,n1=n1,n2=n2,n3=n3,N=N,gamma_delta=gamma_delta,u1=u1,u2=u2,shared_divisor = shared_divisor)
    return(numerator/denominator)
  }else if((I==4)&(J==2)){
    ## this is the case of four treatment and two outcomes
    n1 = treatment_margins[1]
    n2 = treatment_margins[2]
    n3 = treatment_margins[3]
    n4 = treatment_margins[4]
    Us = sum(u_allocation)
    denominator = denominator_four_treatment(n1=n1,n2=n2,n3=n3,n4=n4,gamma_delta=gamma_delta,N=N,Us=Us,shared_divisor=shared_divisor)
    n11 = table[1,1]
    n21 = table[2,1]
    n31 = table[3,1]
    m1 = outcome_margins[1]
    m2 = outcome_margins[2]
    u1 = u_allocation[1]
    u2 = u_allocation[2]
    numerator = d_numerator_four_by_two(n11=n11,n21=n21,n31=n31,m1=m1,m2=m2,n1=n1,n2=n2,n3=n3,n4=n4,N=N,gamma_delta=gamma_delta,u1=u1,u2=u2,shared_divisor = shared_divisor)
    return(numerator/denominator)

  }else if((I==5) &(J==2)){
    ## five treatment two outcome
    n1 = treatment_margins[1]
    n2 = treatment_margins[2]
    n3 = treatment_margins[3]
    n4 = treatment_margins[4]
    n5 = treatment_margins[5]
    Us = sum(u_allocation)
    denominator = denominator_five_treatment(n1=n1,n2=n2,n3=n3,n4=n4,n5=n5,gamma_delta=gamma_delta,N=N,Us=Us,shared_divisor=shared_divisor)
    n11 = table[1,1]
    n21 = table[2,1]
    n31 = table[3,1]
    n41 = table[4,1]
    m1 = outcome_margins[1]
    m2 = outcome_margins[2]
    u1 = u_allocation[1]
    u2 = u_allocation[2]
    numerator = d_numerator_five_by_two(n11=n11,n21=n21,n31=n31,n41=n41,m1=m1,m2=m2,n1=n1,n2=n2,n3=n3,n4=n4,n5=n5,N=N,gamma_delta=gamma_delta,u1=u1,u2=u2,shared_divisor = shared_divisor)
    return(numerator/denominator)

  }else if((I==2)&(J==3)){
    ## the case of two treatment and three outcome
    n1 = treatment_margins[1]
    n2 = treatment_margins[2]
    Us = sum(u_allocation)
    denominator = denominator_two_treatment(n1=n1,n2=n2,gamma_delta=gamma_delta,N=N,Us=Us,shared_divisor=shared_divisor)
    n11 = table[1,1]
    n12 = table[1,2]
    m1 = outcome_margins[1]
    m2 = outcome_margins[2]
    m3 = outcome_margins[3]
    u1 = u_allocation[1]
    u2 = u_allocation[2]
    u3 = u_allocation[3]
    numerator = d_numerator_two_by_three(n11=n11,n12=n12,m1=m1,m2=m2,m3=m3,n1=n1,n2=n2,N=N,gamma_delta=gamma_delta,u1=u1,u2=u2,u3=u3,shared_divisor = shared_divisor)
    return(numerator/denominator)


  }else if((I==3)&(J==3)){
    n1 = treatment_margins[1]
    n2 = treatment_margins[2]
    n3 = treatment_margins[3]
    Us = sum(u_allocation)
    denominator = denominator_three_treatment(n1=n1,n2=n2,n3=n3,gamma_delta=gamma_delta,N=N,Us=Us,shared_divisor=shared_divisor)
    n11 = table[1,1]
    n12 = table[1,2]
    n21 = table[2,1]
    n22 = table[2,2]
    m1 = outcome_margins[1]
    m2 = outcome_margins[2]
    m3 = outcome_margins[3]
    u1 = u_allocation[1]
    u2 = u_allocation[2]
    u3 = u_allocation[3]
    numerator = d_numerator_three_by_three(n11=n11,n12=n12,n21=n21,n22=n22,m1=m1,m2=m2,m3=m3,n1=n1,n2=n2,n3=n3,N=N,gamma_delta=gamma_delta,u1=u1,u2=u2,u3=u3,shared_divisor = shared_divisor)
    return(numerator/denominator)


  }else if((I==4)&(J==3)){
    ## four treatment and three outcome
    n1 = treatment_margins[1]
    n2 = treatment_margins[2]
    n3 = treatment_margins[3]
    n4 = treatment_margins[4]
    Us = sum(u_allocation)
    denominator = denominator_four_treatment(n1=n1,n2=n2,n3=n3,n4=n4,gamma_delta=gamma_delta,N=N,Us=Us,shared_divisor=shared_divisor)
    n11 = table[1,1]
    n12 = table[1,2]
    n21 = table[2,1]
    n22 = table[2,2]
    n31 = table[3,1]
    n32 = table[3,2]
    m1 = outcome_margins[1]
    m2 = outcome_margins[2]
    m3 = outcome_margins[3]
    u1 = u_allocation[1]
    u2 = u_allocation[2]
    u3 = u_allocation[3]
    numerator = d_numerator_four_by_three(n11=n11,n12=n12,n21=n21,n22=n22,n31=n31,n32=n32,m1=m1,m2=m2,m3=m3,n1=n1,n2=n2,n3=n3,n4=n4,N=N,gamma_delta=gamma_delta,u1=u1,u2=u2,u3=u3,shared_divisor = shared_divisor)
    return(numerator/denominator)



  }else if((I==2)&(J==4)){
    ## two treatment and four outcomes
    n1 = treatment_margins[1]
    n2 = treatment_margins[2]
    Us = sum(u_allocation)
    denominator = denominator_two_treatment(n1=n1,n2=n2,gamma_delta=gamma_delta,N=N,Us=Us,shared_divisor=shared_divisor)
    n11 = table[1,1]
    n12 = table[1,2]
    n13 = table[1,3]
    m1 = outcome_margins[1]
    m2 = outcome_margins[2]
    m3 = outcome_margins[3]
    m4 = outcome_margins[4]
    u1 = u_allocation[1]
    u2 = u_allocation[2]
    u3 = u_allocation[3]
    u4 = u_allocation[4]
    numerator = d_numerator_two_by_four(n11=n11,n12=n12,n13=n13,m1=m1,m2=m2,m3=m3,m4=m4,n1=n1,n2=n2,N=N,gamma_delta=gamma_delta,u1=u1,u2=u2,u3=u3,u4=u4,shared_divisor = shared_divisor)
    return(numerator/denominator)

  }else if((I==2)&(J==5)){
    ## two treatment and five outcomes
    n1 = treatment_margins[1]
    n2 = treatment_margins[2]
    Us = sum(u_allocation)
    denominator = denominator_two_treatment(n1=n1,n2=n2,gamma_delta=gamma_delta,N=N,Us=Us,shared_divisor=shared_divisor)
    n11 = table[1,1]
    n12 = table[1,2]
    n13 = table[1,3]
    n14 = table[1,4]
    m1 = outcome_margins[1]
    m2 = outcome_margins[2]
    m3 = outcome_margins[3]
    m4 = outcome_margins[4]
    m5 = outcome_margins[5]
    u1 = u_allocation[1]
    u2 = u_allocation[2]
    u3 = u_allocation[3]
    u4 = u_allocation[4]
    u5 = u_allocation[5]
    numerator = d_numerator_two_by_five(n11=n11,n12=n12,n13=n13,n14=n14, m1=m1,m2=m2,m3=m3,m4=m4,m5=m5,n1=n1,n2=n2,N=N,gamma_delta=gamma_delta,u1=u1,u2=u2,u3=u3,u4=u4,u5=u5,shared_divisor = shared_divisor)
    return(numerator/denominator)



  }else if((I==3)&(J==4)){
    ## three treatment and four outcomes
    n1 = treatment_margins[1]
    n2 = treatment_margins[2]
    n3 = treatment_margins[3]
    Us = sum(u_allocation)
    denominator = denominator_three_treatment(n1=n1,n2=n2,n3=n3,gamma_delta=gamma_delta,N=N,Us=Us,shared_divisor=shared_divisor)
    n11 = table[1,1]
    n12 = table[1,2]
    n13 = table[1,3]
    n21 = table[2,1]
    n22 = table[2,2]
    n23 = table[2,3]
    m1 = outcome_margins[1]
    m2 = outcome_margins[2]
    m3 = outcome_margins[3]
    m4 = outcome_margins[4]
    u1 = u_allocation[1]
    u2 = u_allocation[2]
    u3 = u_allocation[3]
    u4 = u_allocation[4]
    numerator = d_numerator_three_by_four(n11=n11,n12=n12,n13=n13,n21=n21,n22=n22,n23=n23,m1=m1,m2=m2,m3=m3,m4=m4,n1=n1,n2=n2,n3=n3,N=N,gamma_delta=gamma_delta,u1=u1,u2=u2,u3=u3,u4=u4,shared_divisor = shared_divisor)
    return(numerator/denominator)


  }else if((I==3)&(J==5)){
    n1 = treatment_margins[1]
    n2 = treatment_margins[2]
    n3 = treatment_margins[3]
    Us = sum(u_allocation)
    denominator = denominator_three_treatment(n1=n1,n2=n2,n3=n3,gamma_delta=gamma_delta,N=N,Us=Us,shared_divisor=shared_divisor)
    n11 = table[1,1]
    n12 = table[1,2]
    n13 = table[1,3]
    n14 = table[1,4]
    n21 = table[2,1]
    n22 = table[2,2]
    n23 = table[2,3]
    n24 = table[2,4]
    m1 = outcome_margins[1]
    m2 = outcome_margins[2]
    m3 = outcome_margins[3]
    m4 = outcome_margins[4]
    m5 = outcome_margins[5]
    u1 = u_allocation[1]
    u2 = u_allocation[2]
    u3 = u_allocation[3]
    u4 = u_allocation[4]
    u5 = u_allocation[5]
    numerator = d_numerator_three_by_five(n11=n11,n12=n12,n13=n13,n14=n14,n21=n21,n22=n22,n23=n23,n24=n24,m1=m1,m2=m2,m3=m3,m4=m4,m5=m5,n1=n1,n2=n2,n3=n3,N=N,gamma_delta=gamma_delta,u1=u1,u2=u2,u3=u3,u4=u4,u5=u5,shared_divisor = shared_divisor)
    return(numerator/denominator)

  }else{
    warning("sorry your table dimension not in the system!")
    return(NULL)
  }
}




#-------------------------------------------------------------------------------------------------------------------------------



#' Exact Sensitivity Analysis for General Test Statistics in I×J Tables
#'
#' This function computes exact p-values for sensitivity analysis in I×J contingency tables
#' under the generic bias model. It enumerates all tables in the reference set and
#' calculates the maximum p-value over the sensitivity parameter space (u allocations).
#' The function is designed for general permutation-invariant test statistics.
#'
#' @param u_space A matrix where each row represents a unique allocation of `u_i=1` across outcomes.
#'               If `NULL`, the function generates a default set of corner allocations based on 
#'               the number of outcomes. Defaults to `NULL`.
#' @param obs.table A matrix or table object representing the observed contingency table.
#' @param table_space A list of matrices or table objects representing the space of contingency 
#'                    tables to consider (typically all tables with test statistic >= observed).
#' @param gamma A scalar
#' @param delta A binary vector with no more than two unique values, corresponding to treatment levels.
#'              The length must match the number of treatments (rows of `obs.table` if `row = "treatment"`,
#'              or columns if `row = "outcome"`). Represents γ×δ in the generic bias model.
#' @param row A string indicating whether rows represent "outcome" or "treatment". 
#'            Must be either "outcome" or "treatment". Default is "treatment".
#' @param verbose A logical flag indicating whether to print progress messages showing the
#'               current maximizer and probability at each step. Default is `FALSE`.
#'
#' @return A list containing:
#' \describe{
#'   \item{rct.prob}{Probability under Randomized Controlled Trial (RCT) with `u_allocation` set to zero.}
#'   \item{max.prob}{Maximum probability found across all allocations in `u_space`.}
#'   \item{maximizer}{The `u_allocation` vector that yields `max.prob`.}
#'   \item{gamma}{Extracted γ value from the generic bias model.}
#'   \item{delta}{remind the practitioners of their delta}
#' }
#' 
#' @details
#' The function performs exact sensitivity analysis by:
#' 1. Enumerating all possible u allocations (or using provided u_space)
#' 2. Computing the p-value for each allocation by summing probabilities over table_space
#' 3. Finding the allocation that maximizes the p-value
#' 
#' For computational efficiency, the function only supports certain table dimensions.
#' If u_space is not provided, default corner allocations are generated for J <= 5.
#'
#' @examples
#' ## Example 1: 3×3 table with custom test statistic
#' # Create an observed table (example data)
#' obs.table <- matrix(c(5, 3, 2, 6, 11, 7, 3, 0, 3), ncol = 3, byrow = TRUE)
#' 
#' # Define a test statistic emphasizing certain cells
#' transform.fun <- function(tb){
#'   test.stat <- 4*tb[3,3] + 3*tb[2,3]
#'   return(test.stat)
#' }
#' obs.stat <- transform.fun(obs.table)
#' 
#' # Find the reference set (tables with test statistic >= observed)
#' table.set <- possible.table(threshold = obs.stat, table = obs.table,
#'                            direction = "greater than", transform.fun = transform.fun)
#'
#' # Perform sensitivity analysis
#' # This assumes longer hospital stays might indicate hidden bias (e.g., reluctance to leave)
#' sen.result <- exact.general.sen.IxJ(obs.table = obs.table,
#'                                    table_space = table.set, 
#'                                    gamma = 0.5
#'                                    delta = c(0, 1, 1))
#'
#' ## Example 2: Customizing corner allocations
#' # Create a 2×3 collapsed table
#' obs.table <- matrix(c(10, 15, 5, 8, 12, 10), ncol = 3, byrow = TRUE)
#' 
#' # Define test statistic with severity weights
#' transform.fun <- function(tb){
#'   test.stat <- 5*tb[1,1] + 2*tb[1,2]
#'   return(test.stat)
#' }
#' obs.stat <- transform.fun(obs.table)
#' 
#' # Get reference set
#' table.set <- possible.table(threshold = obs.stat, table = obs.table,
#'                            direction = "greater than", transform.fun = transform.fun)
#' 
#' # Create custom u-space
#' outcome.margins <- colSums(obs.table)
#' u.list <- possible.u.allocation(outcome_margins = outcome.margins, null.level = 3)
#' 
#' # Run sensitivity analysis with custom u-space
#' sen.result <- exact.general.sen.IxJ(obs.table = obs.table,
#'                                    table_space = table.set, 
#'                                    u_space = u.list,
#'                                    gamma = 0.3,
#'                                    delta = c(1, 0),
#'                                    verbose = TRUE)
#'
#' @seealso 
#' \code{\link{exact.score.sen.IxJ}} for score tests with monotone trends,
#' \code{\link{possible.table}} for generating reference sets,
#' \code{\link{possible.u.allocation}} for generating allocation spaces
#'
#' @export
exact.general.sen.IxJ = function(u_space = NULL, obs.table, table_space, gamma, delta, row = "treatment", verbose = FALSE) {
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
  obs.table = as.matrix(obs.table)
  
  # Validate 'row' argument
  if (!is.character(row) || length(row) != 1) {
    stop("'row' must be a single string, either 'outcome' or 'treatment'.")
  }
  
  if (!(row %in% c("outcome", "treatment"))) {
    stop("'row' must be either 'outcome' or 'treatment'.")
  }
  if(row == "treatment") {
    I = dim(obs.table)[1]
    J = dim(obs.table)[2]
    outcome_margins = colSums(obs.table)
    treatment_margins = rowSums(obs.table)
  }
  if(row == "outcome") {
    I = dim(obs.table)[2]
    J = dim(obs.table)[1]
    outcome_margins = rowSums(obs.table)
    treatment_margins = colSums(obs.table)
  }
  
  # Validate 'gamma'
  if (!is_scalar_numeric(gamma)) {
    stop("gamma must be a scalar.")
  }
  if(gamma < 0){
    stop("gamma needs to be nonnegative")
  }
  if(!is_binary_vector(delta)){
    stop("delta needs to be binary")
  }
  
  if (length(delta) != I) {
    stop(paste0("'delta' must have the same length as the number of treatments (", I, ")."))
  }
  

  
  # Validate 'u_allocation' if provided within 'u_space'
  if (!is.null(u_space)) {
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
    
    # 1) Construct the matrix that repeats 'outcome_margins' row-wise
    margin_mtx <- matrix(
      rep(outcome_margins, times = nrow(u_space)),
      nrow = nrow(u_space),
      byrow = TRUE
    )
    
    # 2) Compare element by element: (u_space > margin_mtx) gives a logical matrix
    check_mat <- (u_space > margin_mtx)
    
    # 3) If ANY element is TRUE, we have a violation.
    if (any(check_mat)) {
      # 'problem_idx' is a 2-column matrix of (row, col) positions where u_space > margin_mtx
      problem_idx <- which(check_mat, arr.ind = TRUE)
      
      # Extract just the row numbers with problems
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
  
  max.prob = 0
  maximizer = NULL
  
  
  # Create a helper function for creating safe sequences
  safe_seq <- function(n) if (n > 0) seq_len(n) else integer(0)
  
  # Construct u_space if not provided using the improved method
  if (is.null(u_space)) {
    u_space <- switch(as.character(J),
                      "2" = matrix(c(
                        0, 0,
                        outcome_margins[1], 0,
                        0, outcome_margins[2],
                        outcome_margins[1], outcome_margins[2]
                      ), ncol = 2, byrow = TRUE),
                      
                      "3" = as.matrix(rbind(
                        c(0, 0, 0),
                        expand.grid(u1 = 0, u2 = 0, u3 = safe_seq(outcome_margins[3])),
                        expand.grid(u1 = 0, u2 = safe_seq(outcome_margins[2]), u3 = 0),
                        expand.grid(u1 = safe_seq(outcome_margins[1]), u2 = 0, u3 = 0),
                        expand.grid(u1 = 0, u2 = safe_seq(outcome_margins[2]), u3 = outcome_margins[3]),
                        expand.grid(u1 = safe_seq(outcome_margins[1]), u2 = 0, u3 = outcome_margins[3]),
                        expand.grid(u1 = safe_seq(outcome_margins[1]), u2 = outcome_margins[2], u3 = 0),
                        expand.grid(u1 = safe_seq(outcome_margins[1]), u2 = outcome_margins[2], u3 = outcome_margins[3])
                      )),
                      
                      "4" = as.matrix(rbind(
                        c(0, 0, 0, 0),
                        expand.grid(u1 = 0, u2 = 0, u3 = 0, u4 = safe_seq(outcome_margins[4])),
                        expand.grid(u1 = 0, u2 = 0, u3 = safe_seq(outcome_margins[3]), u4 = 0),
                        expand.grid(u1 = 0, u2 = safe_seq(outcome_margins[2]), u3 = 0, u4 = 0),
                        expand.grid(u1 = safe_seq(outcome_margins[1]), u2 = 0, u3 = 0, u4 = 0),
                        # Single pairs of margins
                        expand.grid(u1 = 0, u2 = 0, u3 = safe_seq(outcome_margins[3]), u4 = outcome_margins[4]),
                        expand.grid(u1 = 0, u2 = safe_seq(outcome_margins[2]), u3 = 0, u4 = outcome_margins[4]),
                        expand.grid(u1 = safe_seq(outcome_margins[1]), u2 = 0, u3 = 0, u4 = outcome_margins[4]),
                        expand.grid(u1 = 0, u2 = safe_seq(outcome_margins[2]), u3 = outcome_margins[3], u4 = 0),
                        expand.grid(u1 = safe_seq(outcome_margins[1]), u2 = 0, u3 = outcome_margins[3], u4 = 0),
                        expand.grid(u1 = safe_seq(outcome_margins[1]), u2 = outcome_margins[2], u3 = 0, u4 = 0),
                        # Triplets of margins
                        expand.grid(u1 = 0, u2 = safe_seq(outcome_margins[2]), u3 = outcome_margins[3], u4 = outcome_margins[4]),
                        expand.grid(u1 = safe_seq(outcome_margins[1]), u2 = 0, u3 = outcome_margins[3], u4 = outcome_margins[4]),
                        expand.grid(u1 = safe_seq(outcome_margins[1]), u2 = outcome_margins[2], u3 = 0, u4 = outcome_margins[4]),
                        expand.grid(u1 = safe_seq(outcome_margins[1]), u2 = outcome_margins[2], u3 = outcome_margins[3], u4 = 0),
                        # All margins
                        expand.grid(u1 = safe_seq(outcome_margins[1]), u2 = outcome_margins[2], u3 = outcome_margins[3], u4 = outcome_margins[4])
                      )),
                      
                      "5" = as.matrix(rbind(
                        c(0, 0, 0, 0, 0),
                        # Generate single margin combinations
                        expand.grid(u1 = safe_seq(outcome_margins[1]), u2 = 0, u3 = 0, u4 = 0, u5 = 0),
                        expand.grid(u1 = 0, u2 = safe_seq(outcome_margins[2]), u3 = 0, u4 = 0, u5 = 0),
                        expand.grid(u1 = 0, u2 = 0, u3 = safe_seq(outcome_margins[3]), u4 = 0, u5 = 0),
                        expand.grid(u1 = 0, u2 = 0, u3 = 0, u4 = safe_seq(outcome_margins[4]), u5 = 0),
                        expand.grid(u1 = 0, u2 = 0, u3 = 0, u4 = 0, u5 = safe_seq(outcome_margins[5])),
                        # Generate key boundary points
                        expand.grid(u1 = 0, u2 = 0, u3 = 0, u4 = 0, u5 = outcome_margins[5]),
                        expand.grid(u1 = 0, u2 = 0, u3 = 0, u4 = outcome_margins[4], u5 = outcome_margins[5]),
                        expand.grid(u1 = 0, u2 = 0, u3 = outcome_margins[3], u4 = outcome_margins[4], u5 = outcome_margins[5]),
                        expand.grid(u1 = 0, u2 = outcome_margins[2], u3 = outcome_margins[3], u4 = outcome_margins[4], u5 = outcome_margins[5]),
                        expand.grid(u1 = outcome_margins[1], u2 = outcome_margins[2], u3 = outcome_margins[3], u4 = outcome_margins[4], u5 = outcome_margins[5])
                      )),
                      
                      # Default case: throw an error
                      stop("Tables with ", J, " outcome levels are not supported in the default u_space construction")
    )
    
    # Remove any duplicate rows that might have been created
    u_space <- unique(u_space)
  }
  
  ## search for all u in u_space
  u_length = nrow(u_space)
  for(u_index in 1:u_length){
    u_allocation = as.numeric(u_space[u_index,])
    
    ## compute the probability of summing all the tables considered
    table.sum.probability.given.u = 0
    table.length = length(table_space)
    for(table.index in 1:table.length){
      current_table <- table_space[[table.index]]
      
      # Define outcome and treatment margins for the current table
      if(row=="outcome"){
        current_outcome_margins <- rowSums(current_table)
        current_treatment_margins <- colSums(current_table)
      }
      if(row=="treatment"){
        current_outcome_margins <-colSums(current_table)
        current_treatment_margins<-rowSums(current_table)
      }
      # Verify that the current table's margins match the observed margins
      margins_match <- all(current_outcome_margins == outcome_margins) &&
        all(current_treatment_margins == treatment_margins)
      if(!margins_match){
        stop("some table in the table_space does not have the same margins as the observed table")
      }
      
      table.sum.probability.given.u = table.sum.probability.given.u + generic.I.by.J.sensitivity.point.probability(
        table=table_space[[table.index]],
        u_allocation=u_allocation,
        gamma = gamma,
        delta = delta,
        row=row
      )
    }
    
    if(verbose==TRUE){
      cat("current max prob", max.prob, "\n")
      cat("current maximizer", maximizer, "\n")
      cat("current u_allocation checked", u_allocation, "\n")
      cat("prob given current u_allocation",table.sum.probability.given.u,"\n")
    }
    if(table.sum.probability.given.u > max.prob){
      max.prob = table.sum.probability.given.u
      maximizer = u_allocation
    }
    
  }
  
  ## compute the probability in RCT to make sure the user really obtains the upper bound
  u_allocation = rep(0, length(outcome_margins))
  table.sum.probability.rct = 0
  table.length = length(table_space)
  for(table.index in 1:table.length){
    table.sum.probability.rct = table.sum.probability.rct + generic.I.by.J.sensitivity.point.probability(
      table=table_space[[table.index]],
      u_allocation=u_allocation,
      gamma = gamma,
      delta = delta,
      row=row
    )
  }
  
  return(list(
    rct.prob = table.sum.probability.rct, 
    max.prob = max.prob,
    maximizer = maximizer,
    gamma = gamma,
    delta = delta
  ))
}




#--------------------------------------------------------------------------------------------------------------------

#' Exact Sensitivity Analysis for Score Tests in I×J Tables
#'
#' This function computes exact p-values for score-based sensitivity analysis in I×J 
#' contingency tables under the generic bias model. It is specifically designed for 
#' ordinal data where both treatment levels and outcomes have a natural ordering,
#' and tests for trend using assigned scores.
#'
#' @param obs.table A matrix or table object representing the observed contingency table.
#' @param gamma A nonnegative scalar.
#' @param delta A binary vector with no more than two unique values, corresponding 
#'              to treatment levels. The length must match the number of treatments 
#'              (rows of `obs.table` if `row = "treatment"`, or columns if `row = "outcome"`). 
#'              Must have a monotone trend matching the treatment ordering.
#' @param row A string indicating whether rows represent "outcome" or "treatment". 
#'            Must be either "outcome" or "treatment". Default is "treatment".
#' @param treatment.scores A numeric vector specifying the scores for each treatment level.
#'                        Must have the same length as the number of treatments and 
#'                        exhibit a monotone trend.
#' @param outcome.scores A numeric vector specifying the scores for each outcome level.
#'                      Must have the same length as the number of outcomes and
#'                      exhibit a monotone trend.
#' @param verbose A logical flag indicating whether to print progress messages. 
#'               Default is `FALSE`.
#'
#' @return A list containing:
#' \describe{
#'   \item{rct.prob}{Probability under Randomized Controlled Trial (RCT) with `u_allocation` set to zero.}
#'   \item{max.prob}{Maximum probability found across all allocations in `u_space`.}
#'   \item{maximizer}{The `u_allocation` vector that yields `max.prob`.}
#'   \item{gamma}{Extracted γ value from the generic bias model.}
#'   \item{delta}{The delta vector}
#'   \item{obs.stat}{The observed test statistic based on the observed table and the assigned scores.}
#'   \item{obs.table}{The observed table.}
#' }
#' 
#' @details
#' The score test assumes both treatments and outcomes are ordinal with monotone trends.
#' The test statistic is computed as the sum of products of cell counts with their 
#' corresponding treatment and outcome scores. The function automatically generates an 
#' appropriate u-space for score tests, focusing on allocations that respect the 
#' ordinal nature of the data.
#'
#' @examples 
#' ## Example 1: Binary outcome table (2×2)
#' obs.table <- matrix(c(12, 18, 17, 3), ncol = 2, byrow = TRUE,
#'                    dimnames = list(treatment = c("control", "treated"),
#'                                   outcome = c("failure", "success")))
#' 
#' # Perform score-based sensitivity analysis
#' result_2x2 <- exact.score.sen.IxJ(obs.table = obs.table,
#'                                   gamma = 0.5,
#'                                  delta = c(0, 1),
#'                                  treatment.scores = c(0, 1),
#'                                  outcome.scores = c(0, 1))
#' result_2x2
#' 
#' ## Example 2: Three-level ordinal outcome (3×3)
#' obs.table <- matrix(c(12, 18, 17, 3, 12, 25, 0, 3, 4), 
#'                    ncol = 3, byrow = TRUE,
#'                    dimnames = list(treatment = c("low", "medium", "high"),
#'                                   outcome = c("poor", "fair", "good")))
#' 
#' # Test for trend with ordinal scores
#' result_3x3 <- exact.score.sen.IxJ(obs.table = obs.table,
#'                                   gamma = 0.5,
#'                                  delta = c(0, 1, 1),
#'                                  treatment.scores = c(0, 1, 2),
#'                                  outcome.scores = c(1, 2, 3))
#' result_3x3
#'
#' @seealso 
#' \code{\link{exact.general.sen.IxJ}} for general test statistics,
#' \code{\link{possible.table}} for generating reference sets
#'
#' @export
exact.score.sen.IxJ = function(obs.table, gamma, delta, row = "treatment", verbose = FALSE,
                               treatment.scores, outcome.scores) {
  # Helper: Safe sequence generator
  safe_seq <- function(n) {
    if (is.na(n) || !is.numeric(n) || n < 0) {
      return(integer(0))
    }
    if (n > 0) {
      return(seq_len(n))
    } else {
      return(integer(0))
    }
  }
  
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
  if(!gamma>=0){
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
  
  

  # Define test statistic function
  transform.fun <- if (row == "treatment") {
    function(tb) sum(tb * outer(treatment.scores, outcome.scores))
  } else {
    function(tb) sum(tb * outer(outcome.scores, treatment.scores))
  }
  
  # Calculate observed test statistic
  obs.stat <- transform.fun(obs.table)
  
  # Generate table space (all tables with test statistic >= observed)
  table_space <- possible.table(threshold = obs.stat, table = obs.table,
                                direction = "greater than", transform.fun = transform.fun)
  
  # Check if any tables were found
  if (length(table_space) == 0) {
    warning("No tables found in the reference set. Check your test statistic and table.")
    return(list(
      rct.prob = NA,
      max.prob = NA,
      maximizer = NA,
      gamma = gamma,
      delta = delta,
      obs.stat = obs.stat,
      obs.table = obs.table
    ))
  }
  
  # Construct u_space for score tests
  # For score tests with monotone trends, we focus on corner allocations
  # that respect the ordinal structure
  if (J == 2) {
    # For binary outcomes, check endpoints
    u_space <- matrix(c(0, 0,
                        0, outcome_margins[2]), 
                      ncol = 2, byrow = TRUE)
  } else if (J == 3) {
    # For 3 outcomes, build progressively from highest to lowest
    u_space_part_1 <- expand.grid(u1 = 0, u2 = 0, u3 = safe_seq(outcome_margins[3]))
    u_space_part_2 <- expand.grid(u1 = 0, u2 = safe_seq(outcome_margins[2]), u3 = outcome_margins[3])
    u_space_part_3 <- expand.grid(u1 = safe_seq(outcome_margins[1]), u2 = outcome_margins[2], u3 = outcome_margins[3])
    u_space <- as.matrix(rbind(u_space_part_1, u_space_part_2, u_space_part_3))
  } else if (J == 4) {
    u_space_part_1 <- expand.grid(u1 = 0, u2 = 0, u3 = 0, u4 = safe_seq(outcome_margins[4]))
    u_space_part_2 <- expand.grid(u1 = 0, u2 = 0, u3 = safe_seq(outcome_margins[3]), u4 = outcome_margins[4])
    u_space_part_3 <- expand.grid(u1 = 0, u2 = safe_seq(outcome_margins[2]), u3 = outcome_margins[3], u4 = outcome_margins[4])
    u_space_part_4 <- expand.grid(u1 = safe_seq(outcome_margins[1]), u2 = outcome_margins[2], u3 = outcome_margins[3], u4 = outcome_margins[4])
    u_space <- as.matrix(rbind(u_space_part_1, u_space_part_2, u_space_part_3, u_space_part_4))
  } else if (J == 5) {
    u_space_part_1 <- expand.grid(u1 = 0, u2 = 0, u3 = 0, u4 = 0, u5 = safe_seq(outcome_margins[5]))
    u_space_part_2 <- expand.grid(u1 = 0, u2 = 0, u3 = 0, u4 = safe_seq(outcome_margins[4]), u5 = outcome_margins[5])
    u_space_part_3 <- expand.grid(u1 = 0, u2 = 0, u3 = safe_seq(outcome_margins[3]), u4 = outcome_margins[4], u5 = outcome_margins[5])
    u_space_part_4 <- expand.grid(u1 = 0, u2 = safe_seq(outcome_margins[2]), u3 = outcome_margins[3], u4 = outcome_margins[4], u5 = outcome_margins[5])
    u_space_part_5 <- expand.grid(u1 = safe_seq(outcome_margins[1]), u2 = outcome_margins[2], u3 = outcome_margins[3], u4 = outcome_margins[4], u5 = outcome_margins[5])
    u_space <- as.matrix(rbind(u_space_part_1, u_space_part_2, u_space_part_3, u_space_part_4, u_space_part_5))
  } else {
    stop("Score test for tables with more than 5 outcomes is not implemented.")
  }
  
  # Search for maximizing u
  max.prob <- 0
  maximizer <- NULL
  
  # Also track probabilities for all u allocations (useful for debugging)
  all_probs <- numeric(nrow(u_space))
  
  for (u_index in seq_len(nrow(u_space))) {
    u_allocation <- as.numeric(u_space[u_index, ])
    
    # Compute probability sum over all tables in table_space
    table.sum.probability.given.u <- 0
    
    for (tb in table_space) {
      # Verify margins match (important for consistency)
      if (row == "outcome") {
        margins_match <- all(rowSums(tb) == outcome_margins) && 
          all(colSums(tb) == treatment_margins)
      } else {
        margins_match <- all(colSums(tb) == outcome_margins) && 
          all(rowSums(tb) == treatment_margins)
      }
      
      if (!margins_match) {
        stop("A table in 'table_space' has mismatched margins.")
      }
      
      # Compute probability for this table with current u_allocation
      prob <- generic.I.by.J.sensitivity.point.probability(
        table = tb, 
        u_allocation = u_allocation,
        gamma = gamma,
        delta = delta, 
        row = row
      )
      
      table.sum.probability.given.u <- table.sum.probability.given.u + prob
    }
    
    # Store probability
    all_probs[u_index] <- table.sum.probability.given.u
    
    # Update maximum if needed
    if (table.sum.probability.given.u > max.prob) {
      max.prob <- table.sum.probability.given.u
      maximizer <- u_allocation
    }
    
    if (verbose) {
      cat("Current u_allocation:", u_allocation, "\n")
      cat("Current probability:", table.sum.probability.given.u, "\n")
      cat("Current max probability:", max.prob, "\n")
      cat("Current maximizer:", maximizer, "\n\n")
    }
  }
  
  # Compute RCT (u = 0) probability for comparison
  u_allocation <- rep(0, J)
  table.sum.probability.rct <- 0
  
  for (tb in table_space) {
    prob <- generic.I.by.J.sensitivity.point.probability(
      table = tb, 
      u_allocation = u_allocation,
      gamma = gamma,
      delta = delta, 
      row = row
    )
    table.sum.probability.rct <- table.sum.probability.rct + prob
  }
  
  # Return results
  return(list(
    rct.prob = table.sum.probability.rct,
    max.prob = max.prob,
    maximizer = maximizer,
    gamma = gamma,
    delta = delta,
    obs.stat = obs.stat,
    obs.table = obs.table
  ))
}




################################################################################
#                                                                              
#                      SAMPLING-BASED COMPUTATION METHODS                      
#                                                                              
# This section implements sampling-based methods for sensitivity analysis in   
# I×J contingency tables under the generic bias model. These methods are       
# designed for larger tables where exact enumeration becomes computationally   
# prohibitive, providing Monte Carlo approximations to exact p-values.         
#                                                                              

# Main Functions:                                                              
# - compute_sampling_target_v()                                              
#    - Helper function to compute importance sampling weights                 
#    - Handles various table dimensions (I×J)                                  
#                                                                                
# -  single.u.sampling.general.IxJ()                                          
#    - Computes p-value estimate for a single u-allocation                    
#    - Uses importance sampling with SIS-generated tables                     
#                                                                                
# -  sampling.general.sen.IxJ()                                               
#    - Main interface for general test statistics                             
#    - Searches over u-space to find maximum p-value                          
#    - Returns Monte Carlo estimates with progression tracking                
#                                                                                
# -  sampling.score.sen.IxJ()                                                 
#    - Specialized for score tests with ordinal data                         
#    - Optimized u-space for monotone trends                                 
#    - Efficient for ordinal outcome scenarios                               
#                                                                                
#################################################################################


## A helper function to compute the required kernels 

compute_sampling_target_v <- function(I, J, mc.table, treatment_margins,
                                      outcome_margins, u_allocation,
                                      N, gamma_delta, shared_divisor) {
  if (I == 2 && J == 2) {
    return(d_numerator_two_by_two(
      n11 = mc.table[1, 1],
      m1 = outcome_margins[1], m2 = outcome_margins[2],
      n1 = treatment_margins[1], n2 = treatment_margins[2],
      N = N, gamma_delta = gamma_delta,
      u1 = u_allocation[1], u2 = u_allocation[2],
      shared_divisor = shared_divisor
    ))
  } else if (I == 2 && J == 3) {
    return(d_numerator_two_by_three(
      n11 = mc.table[1, 1], n12 = mc.table[1, 2],
      m1 = outcome_margins[1], m2 = outcome_margins[2], m3 = outcome_margins[3],
      n1 = treatment_margins[1], n2 = treatment_margins[2],
      N = N, gamma_delta = gamma_delta,
      u1 = u_allocation[1], u2 = u_allocation[2], u3 = u_allocation[3],
      shared_divisor = shared_divisor
    ))
  } else if (I == 2 && J == 4) {
    return(d_numerator_two_by_four(
      n11 = mc.table[1, 1], n12 = mc.table[1, 2], n13 = mc.table[1, 3],
      m1 = outcome_margins[1], m2 = outcome_margins[2],
      m3 = outcome_margins[3], m4 = outcome_margins[4],
      n1 = treatment_margins[1], n2 = treatment_margins[2],
      N = N, gamma_delta = gamma_delta,
      u1 = u_allocation[1], u2 = u_allocation[2],
      u3 = u_allocation[3], u4 = u_allocation[4],
      shared_divisor = shared_divisor
    ))
  } else if (I == 2 && J == 5) {
    return(d_numerator_two_by_five(
      n11 = mc.table[1, 1], n12 = mc.table[1, 2],
      n13 = mc.table[1, 3], n14 = mc.table[1, 4],
      m1 = outcome_margins[1], m2 = outcome_margins[2],
      m3 = outcome_margins[3], m4 = outcome_margins[4],
      m5 = outcome_margins[5],
      n1 = treatment_margins[1], n2 = treatment_margins[2],
      N = N, gamma_delta = gamma_delta,
      u1 = u_allocation[1], u2 = u_allocation[2],
      u3 = u_allocation[3], u4 = u_allocation[4], u5 = u_allocation[5],
      shared_divisor = shared_divisor
    ))
  } else if (I == 3 && J == 2) {
    return(d_numerator_three_by_two(
      n11 = mc.table[1, 1], n21 = mc.table[2, 1],
      m1 = outcome_margins[1], m2 = outcome_margins[2],
      n1 = treatment_margins[1], n2 = treatment_margins[2], n3 = treatment_margins[3],
      N = N, gamma_delta = gamma_delta,
      u1 = u_allocation[1], u2 = u_allocation[2],
      shared_divisor = shared_divisor
    ))
  } else if (I == 3 && J == 3) {
    return(d_numerator_three_by_three(
      n11 = mc.table[1, 1], n12 = mc.table[1, 2],
      n21 = mc.table[2, 1], n22 = mc.table[2, 2],
      m1 = outcome_margins[1], m2 = outcome_margins[2], m3 = outcome_margins[3],
      n1 = treatment_margins[1], n2 = treatment_margins[2], n3 = treatment_margins[3],
      N = N, gamma_delta = gamma_delta,
      u1 = u_allocation[1], u2 = u_allocation[2], u3 = u_allocation[3],
      shared_divisor = shared_divisor
    ))
  } else if (I == 3 && J == 4) {
    return(d_numerator_three_by_four(
      n11 = mc.table[1, 1], n12 = mc.table[1, 2], n13 = mc.table[1, 3],
      n21 = mc.table[2, 1], n22 = mc.table[2, 2], n23 = mc.table[2, 3],
      m1 = outcome_margins[1], m2 = outcome_margins[2],
      m3 = outcome_margins[3], m4 = outcome_margins[4],
      n1 = treatment_margins[1], n2 = treatment_margins[2], n3 = treatment_margins[3],
      N = N, gamma_delta = gamma_delta,
      u1 = u_allocation[1], u2 = u_allocation[2],
      u3 = u_allocation[3], u4 = u_allocation[4],
      shared_divisor = shared_divisor
    ))
  } else if (I == 3 && J == 5) {
    return(d_numerator_three_by_five(
      n11 = mc.table[1, 1], n12 = mc.table[1, 2],
      n13 = mc.table[1, 3], n14 = mc.table[1, 4],
      n21 = mc.table[2, 1], n22 = mc.table[2, 2],
      n23 = mc.table[2, 3], n24 = mc.table[2, 4],
      m1 = outcome_margins[1], m2 = outcome_margins[2],
      m3 = outcome_margins[3], m4 = outcome_margins[4], m5 = outcome_margins[5],
      n1 = treatment_margins[1], n2 = treatment_margins[2], n3 = treatment_margins[3],
      N = N, gamma_delta = gamma_delta,
      u1 = u_allocation[1], u2 = u_allocation[2],
      u3 = u_allocation[3], u4 = u_allocation[4], u5 = u_allocation[5],
      shared_divisor = shared_divisor
    ))
  } else if (I == 4 && J == 2) {
    return(d_numerator_four_by_two(
      n11 = mc.table[1, 1], n21 = mc.table[2, 1], n31 = mc.table[3, 1],
      m1 = outcome_margins[1], m2 = outcome_margins[2],
      n1 = treatment_margins[1], n2 = treatment_margins[2],
      n3 = treatment_margins[3], n4 = treatment_margins[4],
      N = N, gamma_delta = gamma_delta,
      u1 = u_allocation[1], u2 = u_allocation[2],
      shared_divisor = shared_divisor
    ))
  } else if (I == 4 && J == 3) {
    return(d_numerator_four_by_three(
      n11 = mc.table[1, 1], n12 = mc.table[1, 2],
      n21 = mc.table[2, 1], n22 = mc.table[2, 2],
      n31 = mc.table[3, 1], n32 = mc.table[3, 2],
      m1 = outcome_margins[1], m2 = outcome_margins[2], m3 = outcome_margins[3],
      n1 = treatment_margins[1], n2 = treatment_margins[2],
      n3 = treatment_margins[3], n4 = treatment_margins[4],
      N = N, gamma_delta = gamma_delta,
      u1 = u_allocation[1], u2 = u_allocation[2], u3 = u_allocation[3],
      shared_divisor = shared_divisor
    ))
  } else if (I == 5 && J == 2) {
    return(d_numerator_five_by_two(
      n11 = mc.table[1, 1], n21 = mc.table[2, 1], n31 = mc.table[3, 1], n41 = mc.table[4, 1],
      m1 = outcome_margins[1], m2 = outcome_margins[2],
      n1 = treatment_margins[1], n2 = treatment_margins[2],
      n3 = treatment_margins[3], n4 = treatment_margins[4], n5 = treatment_margins[5],
      N = N, gamma_delta = gamma_delta,
      u1 = u_allocation[1], u2 = u_allocation[2],
      shared_divisor = shared_divisor
    ))
  } else {
    stop("Sampling-based computation not implemented for this I × J.")
  }
}




#--------------------------------------------------------------------------------------



## A helper function to compute the upper bound on the p-value evaluated at a corner

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
                                           col_s = outcome_margins, seed = mc)
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


#-----------------------------------------------------------------------------------

#' Monte Carlo Score Test Sensitivity Analysis for I×J Tables
#'
#' This function implements a sampling-based approach to sensitivity analysis for 
#' score tests in I×J contingency tables under the generic bias model. It uses 
#' Sequential Importance Sampling (SIS) to approximate p-values when exact 
#' computation is infeasible due to large reference set
#'
#' @param obs.table A matrix or table object representing the observed contingency table.
#' @param gamma a nonnegative scalar
#' @param delta A binary vector with no more than two unique values, 
#'   corresponding to treatment levels. The length must match the number of 
#'   treatments. Must have a monotone trend for ordinal tests.
#' @param row A string indicating whether rows represent "outcome" or "treatment". 
#'   Must be either "outcome" or "treatment". Default is "treatment".
#' @param mc.iteration Integer specifying the number of Monte Carlo iterations 
#'   for each u-allocation. Higher values increase accuracy but require more 
#'   computation time. Default is 10000.
#' @param treatment.scores A numeric vector specifying scores for each treatment 
#'   level. Must have the same length as the number of treatments and exhibit 
#'   a monotone trend.
#' @param outcome.scores A numeric vector specifying scores for each outcome level.
#'   Must have the same length as the number of outcomes and exhibit a monotone trend.
#' @param verbose Logical flag indicating whether to print progress messages 
#'   during computation. Default is FALSE.
#'
#' @return A list containing:
#'   \item{rct.prob}{Estimated probability under RCT (all u-allocations zero)}
#'   \item{max.prob}{Maximum estimated probability across all u-allocations}
#'   \item{maximizer}{The u-allocation vector yielding max.prob}
#'   \item{obs.stat}{Observed test statistic value}
#'   \item{obs.table}{The input observed table}
#'   \item{gamma}{Extracted gamma value from the generic bias model}
#'   \item{delta}{delta vector}
#'
#' @details
#' The function uses importance sampling to estimate p-values for score tests
#' when both treatments and outcomes are ordinal. The score test statistic is
#' computed as the sum of products of cell counts with their corresponding
#' treatment and outcome scores.
#' 
#' The u-space is automatically constructed to respect the ordinal nature of
#' the data, focusing on allocations that assign bias to higher outcome levels
#' first (assuming an increasing trend in outcome scores).
#' 
#' Unlike exact methods, this function provides Monte Carlo estimates that
#' converge to true values as mc.iteration increases. Results may vary
#' slightly between runs unless the random seed is fixed.
#'
#' @examples
#' \dontrun{
#' # Binary outcome example
#' obs.table <- matrix(c(15, 10, 5, 20), ncol = 2, byrow = TRUE)
#' result <- sampling.score.sen.IxJ(obs.table = obs.table,
#'                                  gamma = 0.5,
#'                                  delta = c(0, 1),
#'                                  treatment.scores = c(0, 1),
#'                                  outcome.scores = c(0, 1),
#'                                  mc.iteration = 5000)
#' }
#'
#' @seealso 
#' \code{\link{exact.score.sen.IxJ}} for exact computation when feasible,
#' \code{\link{sampling.general.sen.IxJ}} for general test statistics
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
  
  
  # Define test statistic function
  transform.fun <- if (row == "treatment") {
    function(tb) sum(tb * outer(treatment.scores, outcome.scores))
  } else {
    function(tb) sum(tb * outer(outcome.scores, treatment.scores))
  }
  
  obs.stat <- transform.fun(obs.table)
  
  if (row == "outcome") {
    obs.table <- t(obs.table)
  }
  
  I <- nrow(obs.table)
  J <- ncol(obs.table)
  outcome_margins <- colSums(obs.table)
  treatment_margins <- rowSums(obs.table)
  
  # Construct u_space based on J
  safe_seq <- function(n) if (n > 0) seq_len(n) else integer(0)
  
  u_space <- switch(as.character(J),
                    "2" = matrix(data = c(0, as.numeric(outcome_margins[2])), ncol = 2),
                    "3" = as.matrix(rbind(
                      expand.grid(u1 = 0, u2 = 0, u3 = safe_seq(outcome_margins[3])),
                      expand.grid(u1 = 0, u2 = safe_seq(outcome_margins[2]), u3 = outcome_margins[3]),
                      expand.grid(u1 = safe_seq(outcome_margins[1]), u2 = outcome_margins[2], u3 = outcome_margins[3])
                    )),
                    "4" = as.matrix(rbind(
                      expand.grid(u1 = 0, u2 = 0, u3 = 0, u4 = safe_seq(outcome_margins[4])),
                      expand.grid(u1 = 0, u2 = 0, u3 = safe_seq(outcome_margins[3]), u4 = outcome_margins[4]),
                      expand.grid(u1 = 0, u2 = safe_seq(outcome_margins[2]), u3 = outcome_margins[3], u4 = outcome_margins[4]),
                      expand.grid(u1 = safe_seq(outcome_margins[1]), u2 = outcome_margins[2], u3 = outcome_margins[3], u4 = outcome_margins[4])
                    )),
                    "5" = as.matrix(rbind(
                      expand.grid(u1 = 0, u2 = 0, u3 = 0, u4 = 0, u5 = safe_seq(outcome_margins[5])),
                      expand.grid(u1 = 0, u2 = 0, u3 = 0, u4 = safe_seq(outcome_margins[4]), u5 = outcome_margins[5]),
                      expand.grid(u1 = 0, u2 = 0, u3 = safe_seq(outcome_margins[3]), u4 = outcome_margins[4], u5 = outcome_margins[5]),
                      expand.grid(u1 = 0, u2 = safe_seq(outcome_margins[2]), u3 = outcome_margins[3], u4 = outcome_margins[4], u5 = outcome_margins[5]),
                      expand.grid(u1 = safe_seq(outcome_margins[1]), u2 = outcome_margins[2], u3 = outcome_margins[3], u4 = outcome_margins[4], u5 = outcome_margins[5])
                    )),
                    stop("Table with this number of outcome levels not supported")
  )
  
  # Iterate through all candidate u vectors and evaluate via sampling
  max.prob <- 0
  maximizer <- NULL
  prob.vector <- numeric(nrow(u_space))
  
  for (idx in seq_len(nrow(u_space))) {
    u_allocation <- as.numeric(u_space[idx, ])
    result <- single.u.sampling.general.IxJ(
      table = obs.table,
      row = row,
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
    row = row,
    mc.iteration = mc.iteration,
    u_allocation = rep(0, J),
    gamma = gamma,
    delta = delta,
    transform.fun = transform.fun
  )
  
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



#' Monte Carlo Sensitivity Analysis for General Test Statistics in I×J Tables
#'
#' This function implements a sampling-based approach to sensitivity analysis for
#' general (permutation-invariant) test statistics in I×J contingency tables under
#' the generic bias model. It uses Sequential Importance Sampling (SIS) to
#' approximate p-values when exact enumeration is computationally infeasible.
#'
#' @param obs.table A matrix or table object representing the observed contingency table.
#' @param gamma a nonnegative scalar
#' @param delta A binary vector with no more than two unique values,
#'   corresponding to treatment levels. The length must match the number of
#'   treatments (rows if \code{row = "treatment"}, columns if \code{row = "outcome"}).
#' @param row A string indicating whether rows represent "outcome" or "treatment".
#'   Must be either "outcome" or "treatment". Default is "treatment".
#' @param mc.iteration Integer specifying the number of Monte Carlo iterations
#'   for each u-allocation. Higher values increase accuracy but require more
#'   computation time. Default is 5000.
#' @param transform.fun A user-defined function that takes a contingency table
#'   and returns a numeric test statistic. This function should be
#'   permutation-invariant within treatment groups.
#' @param verbose Logical flag indicating whether to print progress messages
#'   showing current u-allocation and estimated probabilities. Default is FALSE.
#' @param u_space Optional matrix where each row represents a candidate u-allocation.
#'   If NULL, a default set of corner allocations is generated. Each row should
#'   have length equal to the number of outcomes.
#'
#' @return A list containing:
#'   \item{rct.prob}{Estimated probability under RCT (all u-allocations zero)}
#'   \item{max.prob}{Maximum estimated probability across all u-allocations}
#'   \item{maximizer}{The u-allocation vector yielding max.prob}
#'   \item{obs.stat}{Observed test statistic value from transform.fun}
#'   \item{obs.table}{The input observed table}
#'   \item{gamma}{Extracted gamma value from the generic bias model}
#'   \item{delta}{The delta vector}
#'
#' @details
#' This function performs Monte Carlo sensitivity analysis for arbitrary test
#' statistics that are permutation-invariant within treatment groups. Unlike
#' the score test version, this function can handle any user-defined test
#' statistic through the \code{transform.fun} parameter.
#'
#' The algorithm:
#' 1. Generates tables using Sequential Importance Sampling (SIS)
#' 2. Evaluates the test statistic on each sampled table
#' 3. Estimates p-values using importance weights
#' 4. Searches over u-allocations to find the maximum p-value
#'
#' When \code{u_space} is not provided, the function generates default corner
#' allocations that explore extreme cases in the sensitivity parameter space.
#' This default is designed to be computationally efficient while capturing
#' key sensitivity scenarios.
#'
#' @examples
#' \dontrun{
#' # Example with custom test statistic emphasizing corner cells
#' obs.table <- matrix(c(10, 5, 8, 12), ncol = 2, byrow = TRUE)
#' 
#' # Define test statistic: sum of diagonal elements
#' diag_stat <- function(tb) {
#'   sum(diag(tb))
#' }
#' 
#' # Run sensitivity analysis
#' result <- sampling.general.sen.IxJ(
#'   obs.table = obs.table,
#'   gamma = 0.5,
#'   delta = c(0, 1),
#'   transform.fun = diag_stat,
#'   mc.iteration = 5000,
#'   verbose = TRUE
#' )
#' 
#' # Custom u-space example
#' custom_u <- matrix(c(0, 0,
#'                     5, 0,
#'                     0, 8,
#'                     5, 8), ncol = 2, byrow = TRUE)
#' 
#' result_custom <- sampling.general.sen.IxJ(
#'   obs.table = obs.table,
#'   gamma = 0.5,
#'   gamma_delta = c(0, 1),
#'   transform.fun = diag_stat,
#'   u_space = custom_u
#' )
#' }
#'
#' @seealso
#' \code{\link{exact.general.sen.IxJ}} for exact computation when feasible,
#' \code{\link{sampling.score.sen.IxJ}} for score tests with ordinal data,
#'
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
  
  # Apply the transformation function to get the observed statistic
  obs.stat <- transform.fun(obs.table)
  
  # Handle row orientation
  if (row == "outcome") {
    obs.table <- t(obs.table)
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
  
 
  
  # Helper function for safe sequences
  safe_seq <- function(n) if (n > 0) seq_len(n) else integer(0)
  
  # Construct u_space if not provided
  if (is.null(u_space)) {
    u_space <- switch(as.character(J),
                      "2" = matrix(data = c(0, 0, 0, outcome_margins[2]), ncol = 2, byrow = TRUE),
                      "3" = as.matrix(rbind(
                        expand.grid(u1 = 0, u2 = 0, u3 = safe_seq(outcome_margins[3])),
                        expand.grid(u1 = 0, u2 = safe_seq(outcome_margins[2]), u3 = outcome_margins[3]),
                        expand.grid(u1 = safe_seq(outcome_margins[1]), u2 = outcome_margins[2], u3 = outcome_margins[3])
                      )),
                      "4" = as.matrix(rbind(
                        expand.grid(u1 = 0, u2 = 0, u3 = 0, u4 = safe_seq(outcome_margins[4])),
                        expand.grid(u1 = 0, u2 = 0, u3 = safe_seq(outcome_margins[3]), u4 = outcome_margins[4]),
                        expand.grid(u1 = 0, u2 = safe_seq(outcome_margins[2]), u3 = outcome_margins[3], u4 = outcome_margins[4]),
                        expand.grid(u1 = safe_seq(outcome_margins[1]), u2 = outcome_margins[2], u3 = outcome_margins[3], u4 = outcome_margins[4])
                      )),
                      "5" = as.matrix(rbind(
                        expand.grid(u1 = 0, u2 = 0, u3 = 0, u4 = 0, u5 = safe_seq(outcome_margins[5])),
                        expand.grid(u1 = 0, u2 = 0, u3 = 0, u4 = safe_seq(outcome_margins[4]), u5 = outcome_margins[5]),
                        expand.grid(u1 = 0, u2 = 0, u3 = safe_seq(outcome_margins[3]), u4 = outcome_margins[4], u5 = outcome_margins[5]),
                        expand.grid(u1 = 0, u2 = safe_seq(outcome_margins[2]), u3 = outcome_margins[3], u4 = outcome_margins[4], u5 = outcome_margins[5]),
                        expand.grid(u1 = safe_seq(outcome_margins[1]), u2 = outcome_margins[2], u3 = outcome_margins[3], u4 = outcome_margins[4], u5 = outcome_margins[5])
                      )),
                      stop("Table with this number of outcome levels not supported")
    )
    
    # Remove any duplicate rows that might have been created
    u_space <- unique(u_space)
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
  max.prob <- 0
  maximizer <- NULL
  prob.vector <- numeric(nrow(u_space))
  
  for (idx in seq_len(nrow(u_space))) {
    u_allocation <- as.numeric(u_space[idx, ])
    
    # Call the sampling function for this u_allocation
    result <- single.u.sampling.general.IxJ(
      table = obs.table,
      row = row,
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
    row = row,
    mc.iteration = mc.iteration,
    u_allocation = rep(0, J),
    gamma = gamma,
    delta = delta,
    transform.fun = transform.fun
  )
  
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
# This section implements normal approximation methods for sensitivity         
# analysis in I×J contingency tables under the generic bias model. These      
# methods provide rapid, asymptotically valid approximations to exact         
# p-values, making them particularly useful for large tables where both       
# exact and sampling methods may be computationally intensive.                
#                                                                              

# Main Functions:                                                              
# - norm_single_u_allocation_p_value()                                       
#    - Computes p-value for a single u-allocation                            
#    - Calculates mean and variance under generic bias                       
#    - Returns z-score and associated p-value                                
#                                                                              
# - norm.sen.IxJ()                                                           
#    - Main interface for normal approximation sensitivity analysis           
#    - Searches over u-space to find maximum p-value                         
#    - Returns both RCT and maximum bias scenarios                           
#                                                                              
# Supporting Functions (C++ Backend):                                          
# - zr_moments_two_by_two()    : Moments for 2×2 tables                      
# - zr_moments_two_by_three()  : Moments for 2×3 tables                      
# - zr_moments_two_by_four()   : Moments for 2×4 tables                      
# - zr_moments_two_by_five()   : Moments for 2×5 tables                      
# - zr_moments_three_by_two()  : Moments for 3×2 tables                      
# - zr_moments_three_by_three(): Moments for 3×3 tables                      
# - zr_moments_four_by_two()   : Moments for 4×2 tables                      
#                                                                              

################################################################################




# Helper function to construct the A matrix based on the marginal constraints
# This function helps to recover the full covariance matrix by taking in
# the covariance matrix with only the top left free cells of a contingency table
# Some examples for the engineers. The free cells in a 3 by 3 tables are z1r1
# z1r2, z2r1, z2r2, and given those, the var_mat recovers the covariance matrix
# of nine cells. As in the main text, the operation is easier to apply by flattening
# the cell counts. In this particular example, the index mapping is:

# Flattened index <-> (i,j)
#  1                  (1,1)
#  2                  (1,2) 
#  3                  (1,3)
#  4                  (2,1)
#  5                  (2,2)
#  6                  (2,3)
#  7                  (3,1)
#  8                  (3,2)
#  9                  (3,3)
# therefore, the [9,3] entry of the following var_mat
#             gives the covariance between the cell counts [3,3] and [1,3]
#             in a 3 by 3 table, assuming that the c++ function fed the covariance 
#              matrix of the free cells already 

#### try for the 3 by 3 table 
### Sigma_gree means that for the free cells, we have 
### var(z1r1) = 1, var(z1r2) = 2, var(z2r1) = 0.5, var(z2r2) = 0.4
### cov(z1r1, z1r2) = -0.4, cov(z1r1, z2r1) = -0.3, cov(z1r1, z2r2) = 0.7
### cov(z1r2, z2r1) = 0.7, cov(z1r2, z2r2) = -0.5
### cov(z2r1, z2r2) = -0.15
#A = construct_A_matrix(I=3,J=3)
#Sigma_free = matrix(data=c(1,   -0.4, -0.3,  0.7,
#                           -0.4, 2,    0.7, -0.5,
#                           -0.3, 0.7,  0.5, -0.15,
#                           0.7, -0.5,  -0.15, 0.4 ),nrow=4)
#var_mat <- A %*% Sigma_free %*% t(A)

#var_mat





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


# Create full A vector taking the argument of the treatment.scores, and outcome.scores 
# Example for the engineers
# create_A_vector(I=2,J=2, treatment.scores = c(1,2), outcome.scores = c(1,3))
# create_A_vector(I=3,J=2, treatment.score = c(1,2,4), outcome.scores=c(1,4))
# Similar to the construct_A_matrix function, this function flattens the index
#   
# in the 2 by 2 case
# Flattened index <-> (i,j)
# 1                   (1,1)
# 2                   (1,2)
# 3                   (2,1)
# 4                   (2,2)
# therefore, the fourth index of the 2 by 2 score will be treatment.scores[2]*outcome.scores[2]

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


#' Compute the normal-approximation-based z-score and p-value for
#' a given 2x2, 2x3, 2x4, 2x5, 3x2, 4x2, or 3x3 contingency table.
#'
#' This function performs a score test for ordinal association in contingency tables
#' with fixed margins under a sensitivity analysis framework. It uses the test statistic
#' T = sum(A_{ij} * N_{ij}) where A_{ij} = alpha_i * beta_j (treatment score * outcome score)
#' and computes the mean and variance under the null hypothesis using the Poisson-binomial
#' approximation for the multivariate Fisher's noncentral hypergeometric distribution.
#'
#' @param table A matrix or table object representing the contingency table.
#'        Rows represent treatments, columns represent outcomes.
#' @param gamma A nonnegative scalar
#' @param delta A binary vector.
#'        Must have length equal to the number of treatments. Can contain at most two unique values.
#' @param u_allocation A numeric vector of unmeasured confounder allocations for each outcome category.
#'        Must have length equal to the number of outcomes. Each entry must be non-negative
#'        and not exceed the corresponding outcome margin.
#' @param row Character indicating whether rows represent "treatment" or "outcome".
#'        If "outcome", the table will be transposed. Default is "treatment".
#' @param treatment.scores A numeric vector of scores for treatments. Must be monotone (either
#'        increasing or decreasing). Higher scores typically indicate more intense treatments.
#' @param outcome.scores A numeric vector of scores for outcomes. Must be monotone (either
#'        increasing or decreasing). Higher scores typically indicate better outcomes.
#' @param shared_divisor A numeric value used for numerical stability in the Rcpp computations.
#'        Default is 1e6.
#' @param alternative Character specifying the alternative hypothesis.
#'        \itemize{
#'          \item \code{"two.sided"}: uses p-value = 2 * (1 - pnorm(abs(z)))
#'          \item \code{"greater than"}: uses p-value = 1 - pnorm(z)  
#'          \item \code{"less than"}: uses p-value = pnorm(z)
#'        }
#'
#' @return A list containing:
#'   \item{T_obs}{The observed test statistic}
#'   \item{mu_T}{The expected value of the test statistic under the null hypothesis}
#'   \item{var_T}{The variance of the test statistic under the null hypothesis}
#'   \item{z_score}{The standardized test statistic (T_obs - mu_T) / sqrt(var_T)}
#'   \item{p_value}{The p-value based on the normal approximation}
#'   \item{alternative}{The alternative hypothesis used}
#'   \item{treatment.scores}{The treatment scores used}
#'   \item{outcome.scores}{The outcome scores used}
#'
#' @details
#' The function implements a score test for ordinal association where the test statistic
#' is a weighted sum of the cell counts, with weights given by the product of treatment
#' and outcome scores. 
#' The function uses specialized Rcpp functions to compute the mean and variance-covariance
#' structure of the free cells (those with degrees of freedom), then recovers the full
#' mean vector and covariance matrix using the marginal constraints. The test statistic
#' and its moments are then computed using matrix operations.
#'
#' @examples
#' # 2x3 contingency table
#' table <- matrix(c(10, 20, 30, 15, 25, 10), nrow = 2, byrow = TRUE)
#' 
#' # Sensitivity parameters and u_allocation
#' gamma = 0.5
#' delta <- c(0, 1)
#' u_allocation <- c(5, 10, 8)
#' 
#' # Ordinal scores
#' treatment.scores <- c(0, 1)
#' outcome.scores <- c(0, 1, 2)
#' 
#' # Perform test
#' result <- norm_single_u_allocation_p_value(
#'   table = table,
#'   gamma = gamma,
#'   delta = delta,
#'   u_allocation = u_allocation,
#'   treatment.scores = treatment.scores,
#'   outcome.scores = outcome.scores,
#'   alternative = "greater than"
#' )
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
  
  # If gamma is 0, set u_allocation to all zeros of same length
  if(gamma == 0){
    u_allocation <- rep(0, length(u_allocation))
  }
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

#' Normal Approximation Sensitivity Analysis for I×J Tables
#'
#' This function implements normal approximation methods for sensitivity analysis
#' in I×J contingency tables under the generic bias model. It computes asymptotically
#' valid p-values for score test statistics based on the product of treatment and
#' outcome scores, providing rapid analysis for large tables.
#'
#' @param obs.table A matrix or table object representing the observed contingency table.
#' @param gamma a nonnegative scalar.
#' @param delta a binary vector
#'        to treatment levels. Its length must match the number of treatments
#'        (rows of \code{obs.table} if \code{row = "treatment"}, or columns if
#'        \code{row = "outcome"}).
#' @param row A string indicating whether rows represent "outcome" or "treatment". Must be
#'        either "outcome" or "treatment". Default is "treatment".
#' @param treatment.scores A numeric vector of scores for treatments. Must be monotone
#'        (either increasing or decreasing). Higher scores typically indicate more
#'        intense treatments. Length must equal the number of treatments.
#' @param outcome.scores A numeric vector of scores for outcomes. Must be monotone
#'        (either increasing or decreasing). Higher scores typically indicate better
#'        outcomes. Length must equal the number of outcomes.
#' @param shared_divisor Numeric value used for numerical stability in calculations.
#'        Default is 1e6.
#' @param alternative Character: "two.sided", "greater than", or "less than",
#'        specifying the alternative hypothesis. Default is "two.sided".
#' @param u_space A numeric matrix where each row is a candidate \code{u_allocation}.
#'        If \code{NULL} (default), corner allocations are generated automatically
#'        for tables with J ≤ 5 outcomes.
#' @param verbose Logical; if \code{TRUE}, prints progress messages including
#'        the current u-allocation and p-value at each step. Default is FALSE.
#'
#' @return A list containing:
#' \describe{
#'   \item{T_obs}{The observed test statistic value.}
#'   \item{RCT.mean}{Mean of the test statistic under RCT (γ = 1).}
#'   \item{max.mean}{Mean of the test statistic under the sensitivity model at maximizer.}
#'   \item{RCT.var}{Variance of the test statistic under RCT.}
#'   \item{max.var}{Variance of the test statistic under the sensitivity model at maximizer.}
#'   \item{RCT.prob}{P-value under RCT (no unmeasured confounding).}
#'   \item{max.prob}{Maximum p-value across all u-allocations (sensitivity bound).}
#'   \item{maximizer}{The u-allocation vector that yields max.prob.}
#'   \item{treatment.scores}{The treatment scores used in the analysis.}
#'   \item{outcome.scores}{The outcome scores used in the analysis.}
#' }
#'
#' @details
#' The function uses asymptotic normality of linear rank statistics to compute
#' p-values efficiently. For an I×J table, the test statistic is a weighted
#' sum of cell counts where weights are products of treatment and outcome scores:
#' T = sum(alpha_i * beta_j * Z_{ij}) across all cells.
#' 
#' The method computes:
#' \itemize{
#'   \item Mean and variance of the test statistic under the generic bias model
#'   \item Standardized z-scores assuming asymptotic normality
#'   \item P-values based on the specified alternative hypothesis
#' }
#' 
#' When \code{u_space} is not provided, the function automatically generates
#' corner allocations (all combinations of 0 and maximum values) which often
#' contain the worst-case scenarios for sensitivity analysis.
#'
#' @examples
#' # 2×3 table example with ordinal scores
#' obs.table <- matrix(c(10, 20, 30, 15, 25, 10), nrow = 2, byrow = TRUE)
#' treatment.scores <- c(0, 1)    # Control vs Treatment
#' outcome.scores <- c(0, 1, 2)   # Ordinal outcomes
#' 
#' result <- norm.score.sen.IxJ(obs.table = obs.table, 
#'                       gamma = 0.5,
#'                       delta = c(0, 1),
#'                       treatment.scores = treatment.scores,
#'                       outcome.scores = outcome.scores
#'                       )
#' 
#' # 3×4 table with custom scores
#' obs.table <- matrix(rnorm(12, 20, 5), nrow = 3)
#' treatment.scores <- c(0, 0.5, 1)  # Three treatment levels
#' outcome.scores <- c(0, 1, 2, 3)   # Four outcome levels
#' 
#' result <- norm.score.sen.IxJ(obs.table = obs.table,
#'                       gamma = 0.5,
#'                       gamma_delta = c(0, 0, 1),
#'                       treatment.scores = treatment.scores,
#'                       outcome.scores = outcome.scores
#'                       )
#'
#' @seealso
#' \code{\link{exact.general.sen.IxJ}} for exact methods,
#' \code{\link{sampling.general.sen.IxJ}} for Monte Carlo methods
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
  ## 5) If u_space is NULL, generate boundary allocations for J <= 5
  ## ----------------------------------------------------------------
  if (is.null(u_space)) {
    # Helper function for safe sequences
    safe_seq <- function(n) if (n > 0) seq_len(n) else integer(0)
    
    # Generate u_space based on number of outcomes
    u_space <- switch(as.character(J),
                      # Binary outcomes: corners only
                      "2" = matrix(c(
                        0, 0,
                        outcome_margins[1], 0,
                        0, outcome_margins[2],
                        outcome_margins[1], outcome_margins[2]
                      ), ncol = 2, byrow = TRUE),
                      
                      # Three outcomes: all corner combinations
                      "3" = {
                        u_grid <- expand.grid(
                          u1 = c(0, outcome_margins[1]),
                          u2 = c(0, outcome_margins[2]),
                          u3 = c(0, outcome_margins[3])
                        )
                        as.matrix(u_grid)
                      },
                      
                      # Four outcomes: all corner combinations  
                      "4" = {
                        u_grid <- expand.grid(
                          u1 = c(0, outcome_margins[1]),
                          u2 = c(0, outcome_margins[2]),
                          u3 = c(0, outcome_margins[3]),
                          u4 = c(0, outcome_margins[4])
                        )
                        as.matrix(u_grid)
                      },
                      
                      # Five outcomes: all corner combinations
                      "5" = {
                        u_grid <- expand.grid(
                          u1 = c(0, outcome_margins[1]),
                          u2 = c(0, outcome_margins[2]),
                          u3 = c(0, outcome_margins[3]),
                          u4 = c(0, outcome_margins[4]),
                          u5 = c(0, outcome_margins[5])
                        )
                        as.matrix(u_grid)
                      },
                      
                      # Error for unsupported dimensions
                      stop("Automatic boundary generation for J > 5 is not implemented.")
    )
    
    if (verbose) {
      cat("Generated default boundary allocations for J =", J, "\n")
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
  max.mean <- 0
  max.var <- 0
  max.prob <- 0
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
      treatment.scores    = treatment.scores,  # Now passing scores
      outcome.scores      = outcome.scores,    # Now passing scores
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
    treatment.scores    = treatment.scores,  # Now passing scores
    outcome.scores      = outcome.scores,    # Now passing scores
    shared_divisor      = shared_divisor
  )
  
  # Extract RCT results
  RCT.mean <- RCT.out$mu_T
  RCT.var  <- RCT.out$var_T
  RCT.prob <- RCT.out$p_value
  T_obs    <- RCT.out$T_obs
  
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






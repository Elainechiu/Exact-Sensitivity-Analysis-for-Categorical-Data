/* The helper sampling function for the paper 
 * "Exact Sensitivity Analysis for Observational Studies of Contingency Tables"
 * sampling-based kernel method discussed in Section 4. 
 * 
 * 
 * This files adapts the sequential Monte Carlo  
 * sampling method and its accompanying C codes described in the paper 
 * "Sampling for Conditional Inference on Contingency Tables, Eisinger and Chen, 2017"
 * and this paper provides accompanying codes on the website named "GoodCoding.c".
 * 
 * 
 * The following codes adapted from the original paper but not identical.
 * 
 * Two main functions are here:
 * 
 * estimate_tables_good_cpp takes the treatment and outcome margins and provides 
 * an estimate of the total number of tables satisfying the margin constraints.
 * 
 * 
 * sample_sis_good_tables samples a contingency table satisfying treatment and outcome 
 * margins constraints, and then outputs the contingency table and its the weight/probability
 * of being sampled from the sequential monte carlo algorithm. We will use this weight
 * information to develop the sampling-based kernel method, see the documentations and 
 * functions "sampling.general.sen.IxJ" and "sampling.score.sen.IxJ" in the file: sensitivityIxJ.R
 * for the computation of the upper bound on the p-value of general permutation-invariant tests 
 * and ordinal, sign-score tests under the generic bias sensitivity model of this approach.
 * 
 * 
*/


#include <Rcpp.h>
using namespace Rcpp;

#include <Rmath.h>  // Needed for R::runif

#include <vector>
#include <cmath>
#include <numeric>
#include <cstdlib>
#include <algorithm>  // For std::max, etc.

// [[Rcpp::plugins(cpp11)]]

// ----------------- Utility Functions -----------------

int Min2(int x1, int x2) {
  return (x1 <= x2) ? x1 : x2;
}
int Max2(int x1, int x2) {
  return (x1 >= x2) ? x1 : x2;
}

// Sum of vector from i+1 to m
int Sum1(const std::vector<int>& row, int m, int i) {
  int S = 0;
  for (int a = i + 1; a < m; a++) {
    S += row[a];
  }
  return S;
}

// Sum of int vector
int SumVec(const std::vector<int>& v) {
  int total = 0;
  for (auto x : v) total += x;
  return total;
}

// Mean of long double vector
long double MeanDouble(const std::vector<long double>& vec) {
  long double sum = 0.0L;
  for (auto x : vec) sum += x;
  return sum / vec.size();
}

// SD (population version)
long double SDDouble(const std::vector<long double>& vec) {
  long double mean = MeanDouble(vec);
  long double ss = 0.0L;
  for (auto x : vec) {
    long double diff = x - mean;
    ss += diff * diff;
  }
  return std::sqrt(ss / vec.size());
}

// Another SD function for sample-based approach
long double SDDouble2(const std::vector<long double>& vec) {
  long double mean = MeanDouble(vec);
  long double ss = 0.0L;
  for (auto x : vec) {
    long double diff = x - mean;
    ss += diff * diff;
  }
  return std::sqrt(ss / (vec.size() - 1));
}

// N Choose M
long double NChooseM(int n, int m) {
  if (m < 0 || n < 0 || n < m) return 0.0L;
  if (m == 0) return 1.0L;
  long double val = 1.0L;
  for (int i = 1; i <= m; i++) {
    val = val * (n - i + 1) / i;
  }
  return val;
}

// GoodCell
long double GoodCell(const std::vector<int>& row, int m,
                     const std::vector<int>& col, int n,
                     int i, int j) {
  int M = SumVec(row);
  long double S1 = NChooseM((n - j) + row[i - 1] - 1, row[i - 1]);
  long double S2 = NChooseM((m - i) + col[j - 1] - 1, col[j - 1]);
  long double S3 = NChooseM(M + m * (n - j) + (m - i) - 1, M);
  
  if (S3 == 0.0L) {
    return 0.0L;
  }
  return (S1 * S2) / S3;
}

// Weighted (multinomial) sampling
int WeightedSample(const std::vector<long double>& vec) {
  long double sum = 0.0L;
  for (auto val : vec) sum += val;
  if (sum <= 0.0L) {
    return (int)(R::runif(0.0, vec.size()));  // fallback uniform
  }
  std::vector<long double> cdf(vec.size());
  cdf[0] = vec[0] / sum;
  for (size_t i = 1; i < vec.size(); i++) {
    cdf[i] = cdf[i - 1] + (vec[i] / sum);
  }
  double u = R::runif(0.0, 1.0);
  for (size_t i = 0; i < cdf.size(); i++) {
    if (u <= cdf[i]) return (int)i;
  }
  return (int)cdf.size() - 1;
}

// ----------------- 1) Estimate # of tables -----------------

// [[Rcpp::export]]
List estimate_tables_good_cpp(std::vector<int> row_s,
                              std::vector<int> col_s,
                              int B = 10000,
                              int seed = 123) {
  srand(seed);
  
  int m = row_s.size();
  int n = col_s.size();
  
  std::vector<long double> w(B);
  
  // main SIS-G loop
  for (int b = 0; b < B; b++) {
    std::vector<int> row = row_s;
    std::vector<int> col = col_s;
    long double wv = 1.0L;
    
    for (int j = 0; j < n; j++) {
      for (int i = 0; i < m; i++) {
        int L = Max2(0, col[j] - Sum1(row, m, i));
        int U = Min2(col[j], row[i]);
        
        int range = U - L + 1;
        std::vector<long double> p(range);
        std::vector<int> val(range);
        
        // fill p
        for (int x = L; x <= U; x++) {
          std::vector<int> rownew = row;
          std::vector<int> colnew = col;
          rownew[i] -= x;
          colnew[j] -= x;
          long double tmp = GoodCell(rownew, m, colnew, n, i + 1, j + 1);
          tmp = std::max(tmp, 1e-300L); // avoid zero
          p[x - L] = tmp;
          val[x - L] = x;
        }
        // Weighted sample
        int idx = WeightedSample(p);
        int chosen = val[idx];
        // update row/col
        row[i] -= chosen;
        col[j] -= chosen;
        
        // replicate WeightedSample's normalization logic for p_chosen
        long double sum_p = 0.0L;
        for (auto &xx : p) sum_p += xx;
        long double p_chosen = p[idx] / sum_p;
        wv *= (1.0L / p_chosen);
      }
    }
    w[b] = wv;
  }
  
  // compute estimate, se, cv2
  long double mean_w = 0.0L;
  for (auto ww : w) mean_w += ww;
  mean_w /= B;
  
  long double sum_sq = 0.0L;
  for (auto ww : w) {
    long double diff = ww - mean_w;
    sum_sq += diff * diff;
  }
  long double sd_w = std::sqrt(sum_sq / (B - 1));
  long double se = sd_w / std::sqrt((long double)B);
  return List::create(
    _["estimate"]       = (long double)mean_w,
    _["standard_error"] = (long double)se,
    _["weights"]        = w
  );
}

// ----------------- 2) Sample B tables using SIS-G  -----------------

// [[Rcpp::export]]
List sample_sis_good_tables(std::vector<int> row_s,
                            std::vector<int> col_s,
                            int B = 1000,
                            int seed = 123) {
  srand(seed);
  
  int m = row_s.size();
  int n = col_s.size();
  
  // We'll store:
  // 1) A list of B integer matrices (the sampled tables)
  // 2) A numeric vector "weights" of length B
  List tableList(B);
  std::vector<long double> weights(B, 0.0L);
  
  for (int b = 0; b < B; b++) {
    // copy margins
    std::vector<int> row = row_s;
    std::vector<int> col = col_s;
    
    // create an R integer matrix for the final table
    IntegerMatrix currentTable(m, n);
    
    // SIS-G weight
    long double wv = 1.0L;
    
    // sample column by column
    for (int j = 0; j < n; j++) {
      for (int i = 0; i < m; i++) {
        int L = Max2(0, col[j] - Sum1(row, m, i));
        int U = Min2(col[j], row[i]);
        
        int range = U - L + 1;
        std::vector<long double> p(range);
        std::vector<int> val(range);
        
        for (int x = L; x <= U; x++) {
          std::vector<int> rowTest = row;
          std::vector<int> colTest = col;
          rowTest[i] -= x;
          colTest[j] -= x;
          long double tmp = GoodCell(rowTest, m, colTest, n, i + 1, j + 1);
          tmp = std::max(tmp, 1e-300L); // floor to avoid zero
          p[x - L] = tmp;
          val[x - L] = x;
        }
        
        // Weighted sample
        int idx = WeightedSample(p);
        int chosen = val[idx];
        
        // Fill the cell in currentTable
        currentTable(i, j) = chosen;
        
        // update row & col
        row[i] -= chosen;
        col[j] -= chosen;
        
        // replicate WeightedSample's normalization logic
        long double sumP = 0.0L;
        for (auto &aa : p) sumP += aa;
        long double p_chosen = p[idx] / sumP;
        
        // accumulate importance weight
        wv *= (1.0L / p_chosen);
      }
    }
    
    // store final table & weight
    tableList[b] = currentTable;
    weights[b] = wv;
  }
  
  return List::create(
    _["tables"]  = tableList,
    _["weights"] = weights
  );
}



// ----------------------------------------------------------------------------
//  sample_sis_one_table()
//    - Draw *one* table by SIS-G for the given row/col margins
//    - Return a List with exactly two elements:
//         $table  = the sampled 2D matrix (IntegerMatrix)
//         $weight = the importance weight (numeric scalar)
// ----------------------------------------------------------------------------
// [[Rcpp::export]]
List sample_sis_one_table(std::vector<int> row_s,
                          std::vector<int> col_s,
                          int seed = 123) {

  int m = row_s.size();
  int n = col_s.size();
  
  // create an R integer matrix for the final table
  IntegerMatrix currentTable(m, n);
  
  // SIS-G weight accumulator
  long double wv = 1.0L;
  
  // We'll sample column by column
  for (int j = 0; j < n; j++) {
    for (int i = 0; i < m; i++) {
      // figure out min (L) and max (U) possible fill for cell (i,j)
      int L = Max2(0, col_s[j] - Sum1(row_s, m, i));
      int U = Min2(col_s[j], row_s[i]);
      
      // Build the probability distribution "p" over possible values x = L..U
      int range = U - L + 1;
      std::vector<long double> p(range);
      std::vector<int> val(range);
      
      for (int x = L; x <= U; x++) {
        // try x in cell (i,j)
        // subtract from row_s[i], col_s[j]
        std::vector<int> rowTest = row_s;
        std::vector<int> colTest = col_s;
        rowTest[i] -= x;
        colTest[j] -= x;
        
        // Evaluate GoodCell(...) for the partial fill
        long double tmp = GoodCell(rowTest, m, colTest, n, i+1, j+1);
        // Avoid zero or negative
        tmp = std::max(tmp, 1e-300L);
        
        p[x - L]   = tmp;
        val[x - L] = x;
      }
      
      // WeightedSample(...) picks an index from p with probability ~ p[i]
      int idx = WeightedSample(p);
      int chosen = val[idx];
      
      // Fill the chosen count into the table
      currentTable(i,j) = chosen;
      
      // Update row & col margins
      row_s[i] -= chosen;
      col_s[j] -= chosen;
      
      // Recompute sum of p for normalization
      long double sumP = 0.0L;
      for (auto &pp : p) {
        sumP += pp;
      }
      long double p_chosen = p[idx] / sumP;
      
      // Multiply SIS weight
      wv *= (1.0L / p_chosen);
    }
  }
  
  // Return exactly one table + one weight
  return List::create(
    _["table"]  = currentTable,
    _["weight"] = (long double)wv
  );
}



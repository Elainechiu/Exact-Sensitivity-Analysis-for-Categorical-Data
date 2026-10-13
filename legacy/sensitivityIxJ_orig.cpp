/*
 * The cpp helper functions for the paper "Exact Sensitivity Analysis for 
 * Observational Studies of Contingency Tables", Section 4, exact method
 * and Section 5, normal approximation method. 
 * 
 * The helper functions include:
 * 
 * ** denominator_two_treatment, denominator_three_treatment, denominator_four_treatment,
 * ** denominator_five_treatment: these are the functions related to the exact p value 
 * ** denominator part, which is called kernel(q) in Appendix C, Lemma 5.
 * 
 * 
 * 
 * ** d_numerator_two_by_three,....d_numerator_three_by_two
 * ** Several functions are named as d_numerator_I_by_J computes the numerator part 
 * ** of the probability, which is called kernel(t,q) in Appendix C, Lemma 6. 
 * 
 * 
 * 
 * ** zu_moments_two_treatment_log, zu_moments_three_treatment_log,.... in general, 
 * ** zu_moments_I_treatment_log computes the first two moments of the vector (\sumI(Z=1,u=1),...\sumI(Z=I,u=1))
 * ** which forms the basis of efficient computation of the moments of cell counts of contingency
 * ** tables \sumI(Z=i,r=j) given the generic bias model. 
 * 
 * 
 * ** zr_moments_two_by_two, zr_moments_two_by_three...in general. zr_moments_I_by_three returns the moments 
 * ** of cell counts from an I by J contingency table, without directly using the probability expression in
 * ** Section 4. 
 * 
 * ** See the functions "exact.general.sen.IxJ", "exact.score.sen.IxJ" for the exact p value computation
 * ** and the function "norm.score.sen.IxJ" for the normally-approximated p value computation 
 * ** for the sensitivity analysis in the file "sensitivityIxJ.R"
 * 
 */








#include <Rcpp.h>
using namespace Rcpp;
#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <stdio.h>
#include <limits.h>
#include <float.h>




double Factorial(double nValue)
{
  double result = nValue;
  double result_next;
  double pc = nValue;
  do
  {
    result_next = result*(pc-1);
    result = result_next;
    pc--;
  }while(pc>2);
  nValue = result;
  return nValue;
}

// [[Rcpp::export]]
double choose_C(double n, double k)
{

  if(n<-0.5){
    return(0.0);
  }

  if(k<-0.5){
    return(0.0);
  }

  if (fabs(n - k) < 1e-7 || fabs(k)  < 1e-7) {
    return 1.0;
  }


  if(k>n){
    return(0.0);
  }



  if( fabs(k-1.0) < 1e-7 || fabs(k - (n-1)) < 1e-7)
  {
    return n;
  }


  return Factorial(n) /(Factorial(k)*Factorial((n - k)));
}



/*
 * Numerically stable version
 */
// Helper function to compute log combinations
double log_choose(int n, int k) {
  if (k < 0 || k > n) return -INFINITY; // log(0) = -Inf
  return lgamma(n + 1.0) - lgamma(k + 1.0) - lgamma(n - k + 1.0);
}


// Helper function to compute log-sum-exp
double compute_log_sum_exp(const std::vector<double>& log_terms) {
  if (log_terms.empty()) return -INFINITY; // log(0) = -Inf
  double max_log = -INFINITY;
  for (const auto& lt : log_terms) {
    if (lt > max_log) max_log = lt;
  }
  if (max_log == -INFINITY) return -INFINITY; // All terms are -Inf
  double sum = 0.0;
  for (const auto& lt : log_terms) {
    sum += std::exp(lt - max_log);
  }
  return max_log + std::log(sum);
}






/*
 * This function calculates the denominator of any table with two treatments, thus, will
 * be the denominator of two by two table, generic bias
 * the shared_divisor command is meant for scaling the denominator and the numerator at the same time
 * We always call the lower bound of something with the notation lb_XX, and upper bound ub_XX
 * We follow the convention in Rosenbaum's that the treatment margins is n_i, outcome margins are m_j
 */


// [[Rcpp::export]]
double denominator_two_treatment(
    int n1,
    int n2,
    int N,
    NumericVector gamma_delta,
    int Us,
    double shared_divisor)
{
  double gamma_val = gamma_delta[0] - gamma_delta[1];

  if (std::fabs(gamma_val) < 1e-14) {
    double val = choose_C(N, n1);
    return val / shared_divisor;
  }



  int lb_q1 = std::max(0,    n1 + Us - N);
  int ub_q1 = std::min(Us,   n1);

  double log_sum = -std::numeric_limits<double>::infinity(); // log(0)

  for (int q1 = lb_q1; q1 <= ub_q1; ++q1) {

    double lc1 = log_choose(Us, q1);
    double lc2 = log_choose(N - Us, n1 - q1);
    double log_term = lc1 + lc2 + gamma_val * q1;

    double max_val = std::max(log_sum, log_term);
    if (std::isinf(max_val)) {
      log_sum = (log_sum > log_term) ? log_sum : log_term;
    } else {
      log_sum = max_val + std::log(std::exp(log_sum - max_val)
                                     + std::exp(log_term - max_val));
    }
  }

  double convolution_sum = std::exp(log_sum) / shared_divisor;
  return convolution_sum;
}

/*
 * This function calculates the denominator of three treatments,generic bias
 */




// [[Rcpp::export]]
double denominator_three_treatment(
    int n1,
    int n2,
    int n3,
    int N,
    NumericVector gamma_delta,
    int Us,
    double shared_divisor)
{
  
  double gamma1 = gamma_delta[0] - gamma_delta[2];
  double gamma2 = gamma_delta[1] - gamma_delta[2];

  // 2) If both gamma1, gamma2 are ~0 => special case => direct binomial product
  if (std::fabs(gamma1) < 1e-14 && std::fabs(gamma2) < 1e-14) {
    double val = choose_C(N, n1) * choose_C(N - n1, n2);
    return val / shared_divisor;
  }

  double gamma_val = (std::fabs(gamma1) > 1e-14) ? gamma1
  : ((std::fabs(gamma2) > 1e-14) ? gamma2 : 0.0);

  // Define bounds for q, q1, q2
  int lb_q  = std::max(0, (n1 + n2 + Us) - N);
  int ub_q  = std::min(n1 + n2, Us);

  int lb_q1 = std::max(0, n1 + Us - N);
  int ub_q1 = std::min(Us, n1);

  int lb_q2 = std::max(0, n2 + Us - N);
  int ub_q2 = std::min(n2, Us);

  double log_sum = -std::numeric_limits<double>::infinity(); // log(0)

  for (int q = lb_q; q <= ub_q; ++q) {
    for (int q1 = lb_q1; q1 <= ub_q1; ++q1) {
      // q2 = q - q1
      int q2 = q - q1;
      if (q2 < lb_q2 || q2 > ub_q2) {
        continue; // invalid => skip
      }


      double log_exp_term = 0.0;
      if (std::fabs(gamma1) < 1e-14) {
        log_exp_term = gamma_val * (q - q1);
      } else if (std::fabs(gamma2) < 1e-14) {
        log_exp_term = gamma_val * (q - q2);
      } else {

        log_exp_term = gamma_val * q;
      }

      // 6) Calculate log of binomial coefficients
      double lc1 = log_choose(Us,            q1);
      double lc2 = log_choose(Us - q1,       q2);
      double lc3 = log_choose(N - Us,        (n1 - q1));
      double lc4 = log_choose((N - Us) - (n1 - q1), (n2 - q2));

      double log_term = lc1 + lc2 + lc3 + lc4 + log_exp_term;


      double max_val = std::max(log_sum, log_term);
      if (std::isinf(max_val)) {
        log_sum = (log_sum > log_term) ? log_sum : log_term;
      } else {
        log_sum = max_val + std::log(std::exp(log_sum - max_val)
                                       + std::exp(log_term - max_val));
      }
    }
  }

  double convolution_sum = std::exp(log_sum) / shared_divisor;
  return convolution_sum;
}



/*
 * This function calculates the denominator of four treatments, generic bias
 */



// [[Rcpp::export]]
double denominator_four_treatment(int n1,
                                                     int n2,
                                                     int n3,
                                                     int n4,
                                                     int N,
                                                     NumericVector gamma_delta,
                                                     int Us,
                                                     double shared_divisor)
{
  // Calculate gamma values
  double gamma1 = gamma_delta[0] - gamma_delta[3];
  double gamma2 = gamma_delta[1] - gamma_delta[3];
  double gamma3 = gamma_delta[2] - gamma_delta[3];

  // Check if all three are zero => special case => direct approach
  bool all_zero = (std::fabs(gamma1) < 1e-14
                     && std::fabs(gamma2) < 1e-14
                     && std::fabs(gamma3) < 1e-14);

                     if (all_zero) {
                       double val = choose_C(N, n1)
                       * choose_C(N - n1, n2)
                       * choose_C(N - n1 - n2, n3);
                       return val / shared_divisor;
                     }

                     double gamma_val = (std::fabs(gamma1) > 1e-14) ? gamma1
                     : ((std::fabs(gamma2) > 1e-14) ? gamma2 : gamma3);

                     // Bounds for q = q1 + q2 + q3
                     int lb_q  = std::max(0, (n1 + n2 + n3) + Us - N);
                     int ub_q  = std::min(n1 + n2 + n3, Us);

                     // Bounds for q1
                     int lb_q1 = std::max(0, n1 + Us - N);
                     int ub_q1 = std::min(Us, n1);

                     double log_sum = -std::numeric_limits<double>::infinity(); // log(0)

                     for (int q = lb_q; q <= ub_q; ++q) {
                       for (int q1 = lb_q1; q1 <= ub_q1; ++q1) {

                         // Bounds for q2
                         int lb_q2 = std::max(0, (n1 + n2 - q1) + Us - N);
                         int ub_q2 = std::min(n2, Us - q1);

                         for (int q2 = lb_q2; q2 <= ub_q2; ++q2) {
                           // q3 = q - q1 - q2
                           int lb_q3 = std::max(0, (n1 + n2 + n3 - q1 - q2) + Us - N);
                           int ub_q3 = std::min(n3, Us - q1 - q2);
                           int q3 = q - q1 - q2;

                           // Check if q3 is valid
                           if (q3 < lb_q3 || q3 > ub_q3) {
                             continue; // skip invalid
                           }

                           double log_exp_term = 0.0;

                           if (std::fabs(gamma1) < 1e-14 && std::fabs(gamma2) < 1e-14) {
                             log_exp_term = gamma_val * (q - q1 - q2);
                           } else if (std::fabs(gamma1) < 1e-14 && std::fabs(gamma3) < 1e-14) {
                             log_exp_term = gamma_val * (q - q1 - q3);
                           } else if (std::fabs(gamma2) < 1e-14 && std::fabs(gamma3) < 1e-14) {
                             log_exp_term = gamma_val * (q - q2 - q3);
                           } else if (std::fabs(gamma1) < 1e-14) {
                             log_exp_term = gamma_val * (q - q1);
                           } else if (std::fabs(gamma2) < 1e-14) {
                             log_exp_term = gamma_val * (q - q2);
                           } else if (std::fabs(gamma3) < 1e-14) {
                             log_exp_term = gamma_val * (q - q3);
                           } else {
                             log_exp_term = gamma_val * q;
                           }

                           double lc1 = log_choose(Us,              q1);
                           double lc2 = log_choose(Us - q1,         q2);
                           double lc3 = log_choose(Us - q1 - q2,    q3);
                           double lc4 = log_choose(N - Us,          (n1 - q1));
                           double lc5 = log_choose(N - Us - (n1 - q1), (n2 - q2));
                           double lc6 = log_choose(N - Us - (n1 - q1) - (n2 - q2), (n3 - q3));

                           double log_term = lc1 + lc2 + lc3 + lc4 + lc5 + lc6 + log_exp_term;

                           double max_val = std::max(log_sum, log_term);
                           if (std::isinf(max_val)) {
                             log_sum = (log_sum > log_term) ? log_sum : log_term;
                           } else {
                             log_sum = max_val + std::log(std::exp(log_sum - max_val)
                                                            + std::exp(log_term - max_val));
                           }
                         }
                       }
                     }

                     double convolution_sum = std::exp(log_sum) / shared_divisor;
                     return convolution_sum;
}



/*
 * the calculation of denominator five treatment, generic bias
 */


// [[Rcpp::export]]
double denominator_five_treatment(
    int n1,
    int n2,
    int n3,
    int n4,
    int n5,
    int N,
    NumericVector gamma_delta,
    int Us,
    double shared_divisor)
{
  // Calculate gamma values: gamma_i = gamma_delta[i] - gamma_delta[4]
  double gamma1 = gamma_delta[0] - gamma_delta[4];
  double gamma2 = gamma_delta[1] - gamma_delta[4];
  double gamma3 = gamma_delta[2] - gamma_delta[4];
  double gamma4 = gamma_delta[3] - gamma_delta[4];

  // Check if all are effectively zero
  bool all_zero = (std::fabs(gamma1) < 1e-14 &&
                   std::fabs(gamma2) < 1e-14 &&
                   std::fabs(gamma3) < 1e-14 &&
                   std::fabs(gamma4) < 1e-14);

  // If all zero, use direct binomial product
  if (all_zero) {
    double val = choose_C(N,    n1)
    * choose_C(N-n1, n2)
    * choose_C(N-n1-n2, n3)
    * choose_C(N-n1-n2-n3, n4);
    return val / shared_divisor;
  }

  double gamma_val = (std::fabs(gamma1) > 1e-14) ? gamma1
  : (std::fabs(gamma2) > 1e-14) ? gamma2
  : (std::fabs(gamma3) > 1e-14) ? gamma3
  : gamma4;

  double log_sum = -std::numeric_limits<double>::infinity();  // log(0)

  int lb_q = std::max(0, (n1 + n2 + n3 + n4) + Us - N);
  int ub_q = std::min(n1 + n2 + n3 + n4, Us);

  for (int q = lb_q; q <= ub_q; ++q) {

    int lb_q1 = std::max(0, n1 + Us - N);
    int ub_q1 = std::min(Us, n1);

    for (int q1 = lb_q1; q1 <= ub_q1; ++q1) {

      int lb_q2 = std::max(0, Us - n3 - n4 - n5 - q1);
      int ub_q2 = std::min(n2, Us - q1);

      for (int q2 = lb_q2; q2 <= ub_q2; ++q2) {

        int lb_q3 = std::max(0, Us - q1 - q2 - n4 - n5);
        int ub_q3 = std::min(n3, Us - q1 - q2);

        for (int q3 = lb_q3; q3 <= ub_q3; ++q3) {

          int lb_q4 = std::max(0, Us - q1 - q2 - q3 - n5);
          int ub_q4 = std::min(n4, Us - q1 - q2 - q3);
          int q4 = q - q1 - q2 - q3;

          if (q4 < lb_q4 || q4 > ub_q4) {
            continue;
          }

          double log_exp_term = 0.0;


          if (std::fabs(gamma1) < 1e-14 && std::fabs(gamma2) < 1e-14
                && std::fabs(gamma3) < 1e-14 && std::fabs(gamma4) < 1e-14) {
            log_exp_term = gamma_val * (q - q1 - q2 - q3 - q4);
          }
          else if (std::fabs(gamma1) < 1e-14 && std::fabs(gamma2) < 1e-14
                     && std::fabs(gamma3) < 1e-14) {
            log_exp_term = gamma_val * (q - q1 - q2 - q3);
          }
          else if (std::fabs(gamma1) < 1e-14 && std::fabs(gamma2) < 1e-14
                     && std::fabs(gamma4) < 1e-14) {
            log_exp_term = gamma_val * (q - q1 - q2 - q4);
          }
          else if (std::fabs(gamma1) < 1e-14 && std::fabs(gamma3) < 1e-14
                     && std::fabs(gamma4) < 1e-14) {
            log_exp_term = gamma_val * (q - q1 - q3 - q4);
          }
          else if (std::fabs(gamma2) < 1e-14 && std::fabs(gamma3) < 1e-14
                     && std::fabs(gamma4) < 1e-14) {
            log_exp_term = gamma_val * (q - q2 - q3 - q4);
          }
          else if (std::fabs(gamma1) < 1e-14 && std::fabs(gamma2) < 1e-14) {
            log_exp_term = gamma_val * (q - q1 - q2);
          }
          else if (std::fabs(gamma1) < 1e-14 && std::fabs(gamma3) < 1e-14) {
            log_exp_term = gamma_val * (q - q1 - q3);
          }
          else if (std::fabs(gamma1) < 1e-14 && std::fabs(gamma4) < 1e-14) {
            log_exp_term = gamma_val * (q - q1 - q4);
          }
          else if (std::fabs(gamma2) < 1e-14 && std::fabs(gamma3) < 1e-14) {
            log_exp_term = gamma_val * (q - q2 - q3);
          }
          else if (std::fabs(gamma2) < 1e-14 && std::fabs(gamma4) < 1e-14) {
            log_exp_term = gamma_val * (q - q2 - q4);
          }
          else if (std::fabs(gamma3) < 1e-14 && std::fabs(gamma4) < 1e-14) {
            log_exp_term = gamma_val * (q - q3 - q4);
          }
          else if (std::fabs(gamma1) < 1e-14) {
            log_exp_term = gamma_val * (q - q1);
          }
          else if (std::fabs(gamma2) < 1e-14) {
            log_exp_term = gamma_val * (q - q2);
          }
          else if (std::fabs(gamma3) < 1e-14) {
            log_exp_term = gamma_val * (q - q3);
          }
          else if (std::fabs(gamma4) < 1e-14) {
            log_exp_term = gamma_val * (q - q4);
          } else {
            log_exp_term = gamma_val * q;
          }



          double lc1 = log_choose(Us, q1);
          double lc2 = log_choose(Us - q1, q2);
          double lc3 = log_choose(Us - q1 - q2, q3);
          double lc4 = log_choose(Us - q1 - q2 - q3, q4);

          double lc5 = log_choose(N - Us,          n1 - q1);
          double lc6 = log_choose(N - Us - (n1 - q1), n2 - q2);
          double lc7 = log_choose(N - Us - (n1 - q1) - (n2 - q2), n3 - q3);
          double lc8 = log_choose(N - Us - (n1 - q1) - (n2 - q2) - (n3 - q3),
                                  n4 - q4);

          double log_term = lc1 + lc2 + lc3 + lc4
          + lc5 + lc6 + lc7 + lc8
          + log_exp_term;

          double max_val = std::max(log_sum, log_term);
          if (std::isinf(max_val)) {
            log_sum = (log_sum > log_term) ? log_sum : log_term;
          } else {
            log_sum = max_val + std::log(
              std::exp(log_sum - max_val) + std::exp(log_term - max_val)
            );
          }
        }
      }
    }
  }

  // 4) Final result
  double convolution_sum = std::exp(log_sum) / shared_divisor;
  return convolution_sum;
}

/*
 * here starts the writing of numerator, generic biase
 */


/*
 * This function computes the numerator of the probability for the joint distribution of the
 * cells from a two treatment two outcome table, the function only takes the support value n11
 * argument, name it d_numerator_two_by_two to follow the convention that d stands for density
 * or probability mass function, thus this function only concerns the probability mass of a point
 * and also given margins m1, m2, n1,n2 with constraints n1+ n2 = m1+m2 = N, under the specification
 * of an allocation u1,u2, again, to later on convert this to probability, make sure the margins
 * match also u1+u2 = Us by definition, and the shared_divisor need to be the same. Generic bias
 */



// [[Rcpp::export]]
double d_numerator_two_by_two(
    int n11,
    int n1,
    int n2,
    int m1,
    int m2,
    int N,
    NumericVector gamma_delta,
    int u1,
    int u2,
    double shared_divisor)
{
  if (n11 > std::min(n1, m1) || (n11 < std::max(0, n1 + m1 - N))) {
    return 0.0;
  }

  double gamma_val = gamma_delta[0] - gamma_delta[1];

  if (std::fabs(gamma_val) < 1e-14) {
    double val = choose_C(m1, n11) * choose_C(m2, n1 - n11);
    return val / shared_divisor;
  }


  int Us = u1 + u2;

  // Bounds for q1
  int lb_q1 = std::max(0,   n1 + Us - N);
  int ub_q1 = std::min(Us,  n1);

  double log_sum = -std::numeric_limits<double>::infinity(); // log(0)

  for (int q1 = lb_q1; q1 <= ub_q1; ++q1) {

    int lb_s11 = std::max({0, n11 + q1 - n1, n11 + u1 - m1, q1 - u2});
    int ub_s11 = std::min({n11, q1, u1, (m2 - u2) + n11 + q1 - n1});

    for (int s11 = lb_s11; s11 <= ub_s11; ++s11) {

      double lc1 = log_choose(u1,        s11);
      double lc2 = log_choose(m1 - u1,   n11 - s11);
      double lc3 = log_choose(u2,        q1 - s11);
      double lc4 = log_choose(m2 - u2,   (n1 - n11) - (q1 - s11));

      double lc_exp = gamma_val * q1;

      double log_term = lc1 + lc2 + lc3 + lc4 + lc_exp;

      double max_val = std::max(log_sum, log_term);
      if (std::isinf(max_val)) {
        log_sum = (log_sum > log_term) ? log_sum : log_term;
      } else {
        log_sum = max_val + std::log(std::exp(log_sum - max_val)
                                       + std::exp(log_term - max_val));
      }
    }
  }

  double convolution_sum = std::exp(log_sum) / shared_divisor;

  return convolution_sum;
}


/*
 * This function computes the numerator of two treatment by three outcome table,
 * To fully characterize a two by three table, we need to have (2-1)*(3-1) = 2 corners, which are n11 and n12
 * We have three allocations, u1, u2, and u3 for each Us loop
 * generic bias
 */


// [[Rcpp::export]]
double d_numerator_two_by_three(
    int n11,
    int n12,
    int n1,
    int n2,
    int m1,
    int m2,
    int m3,
    int N,
    NumericVector gamma_delta,
    int u1,
    int u2,
    int u3,
    double shared_divisor)
{
  // 1) Check if n11, n12 are in valid range...
  if (n11 > std::min(n1, m1) || (n11 < std::max(0, n1 + m1 - N))) {
    return 0.0;
  }
  if (n12 > std::min(n1, m2) || (n12 < std::max(0, n1 + m2 - N))) {
    return 0.0;
  }

  // 2) gamma
  double gamma_val = gamma_delta[0] - gamma_delta[1];

  // 3) If gamma == 0 => direct product
  if (std::fabs(gamma_val) < 1e-14) {
    double val = choose_C(m1, n11)
    * choose_C(m2, n12)
    * choose_C(m3, n1 - n11 - n12);
    return val / shared_divisor;
  }

  // total treated count
  int Us = u1 + u2 + u3;

  // Bounds
  int lb_q1 = std::max(0,   n1 + Us - N);
  int ub_q1 = std::min(Us,  n1);

  // Instead of a vector, maintain a running log-sum
  double log_sum = -std::numeric_limits<double>::infinity();  // log(0)

  for (int q1 = lb_q1; q1 <= ub_q1; ++q1) {
    int lb_s11 = std::max({0, n11 + u1 - m1});
    int ub_s11 = std::min({n11, u1});

    for (int s11 = lb_s11; s11 <= ub_s11; ++s11) {
      int lb_s12 = std::max({
        0,
        n12 + u2 - m2,
        q1 - s11 - u3,
        (q1 - n1 + n11 + n12 - s11)
      });
      int ub_s12 = std::min({
        n12,
        u2,
        q1 - s11,
        (m3 - u3) - (n1 - n11 - n12) + (q1 - s11)
      });

      for (int s12 = lb_s12; s12 <= ub_s12; ++s12) {
        // Log of the term
        double lc1 = log_choose(u1, s11);
        double lc2 = log_choose(m1 - u1, n11 - s11);
        double lc3 = log_choose(u2, s12);
        double lc4 = log_choose(m2 - u2, n12 - s12);
        double lc5 = log_choose(u3, (q1 - s11 - s12));
        double lc6 = log_choose((m3 - u3),
                                (n1 - n11 - n12) - (q1 - s11 - s12));

        double lc_exp = gamma_val * q1;

        double log_term = lc1 + lc2 + lc3 + lc4 + lc5 + lc6 + lc_exp;

        double max_val = std::max(log_sum, log_term);
        if (std::isinf(max_val)) {
          log_sum = (log_sum > log_term) ? log_sum : log_term;
        } else {
          log_sum = max_val + std::log(
            std::exp(log_sum - max_val) + std::exp(log_term - max_val)
          );
        }
      }
    }
  }

  // At the end, log_sum = log of the total.
  double convolution_sum = std::exp(log_sum) / shared_divisor;
  return convolution_sum;
}



/*
 * This function computes the numerator of two by four table, which is, four potential outcomes and two treatments
 * To fully characterize a four by two table, we need to have (4-1)*(2-1) = 3 corners, which are n11, n12, and n13
 * We have four allocations, u1, u2, and u3, and u4 for each Us loop
 * generic bias
 */


// [[Rcpp::export]]
double d_numerator_two_by_four(
    int n11, int n12, int n13,
    int n1,  int n2,
    int m1,  int m2,  int m3,  int m4,
    int N,
    NumericVector gamma_delta,
    int u1,  int u2,  int u3,  int u4,
    double shared_divisor)
{
  // 1) Validate n11, n12, n13 against ranges.
  if (n11 > std::min(n1, m1) || (n11 < std::max(0, n1 + m1 - N))) {
    return 0.0;
  }
  if (n12 > std::min(n1, m2) || (n12 < std::max(0, n1 + m2 - N))) {
    return 0.0;
  }
  if (n13 > std::min(n1, m3) || (n13 < std::max(0, n1 + m3 - N))) {
    return 0.0;
  }

  // 2) gamma = gamma_delta[0] - gamma_delta[1]
  double gamma_val = gamma_delta[0] - gamma_delta[1];

  // 3) If gamma is near zero => direct approach
  if (std::fabs(gamma_val) < 1e-14) {
    double val = choose_C(m1, n11)
    * choose_C(m2, n12)
    * choose_C(m3, n13)
    * choose_C(m4, n1 - n11 - n12 - n13);
    return val / shared_divisor;
  }

  int Us = u1 + u2 + u3 + u4;

  int lb_q1 = std::max(0, n1 + Us - N);
  int ub_q1 = std::min(Us, n1);

  double log_sum = -std::numeric_limits<double>::infinity(); // log(0)

  for (int q1 = lb_q1; q1 <= ub_q1; ++q1) {

    int lb_s11 = std::max({0, n11 + u1 - m1});
    int ub_s11 = std::min({n11, u1});

    for (int s11 = lb_s11; s11 <= ub_s11; ++s11) {

      int lb_s12 = std::max({0, n12 + u2 - m2});
      int ub_s12 = std::min({n12, u2});

      for (int s12 = lb_s12; s12 <= ub_s12; ++s12) {
        int lb_s13 = std::max({
          0,
          n13 + u3 - m3,
          q1 - s11 - s12 - u4,
          (q1 - n1 + n11 + n12 + n13 - s11 - s12)
        });
        int ub_s13 = std::min({
          n13,
          u3,
          (q1 - s11 - s12),
          (m4 - u4) - (n1 - (n11 + n12 + n13)) + (q1 - s11 - s12)
        });

        for (int s13 = lb_s13; s13 <= ub_s13; ++s13) {
          double lc1 = log_choose(u1,  s11);
          double lc2 = log_choose(m1 - u1,  n11 - s11);
          double lc3 = log_choose(u2,  s12);
          double lc4 = log_choose(m2 - u2,  n12 - s12);
          double lc5 = log_choose(u3,  s13);
          double lc6 = log_choose(m3 - u3,  n13 - s13);
          double lc7 = log_choose(u4, (q1 - s11 - s12 - s13));
          double lc8 = log_choose((m4 - u4),
                                  (n1 - n11 - n12 - n13)
                                    - (q1 - s11 - s12 - s13));

          double lc_exp = gamma_val * q1;

          double log_term = lc1 + lc2 + lc3 + lc4 + lc5 + lc6 + lc7 + lc8 + lc_exp;

          double max_val = std::max(log_sum, log_term);
          if (std::isinf(max_val)) {
            log_sum = (log_sum > log_term) ? log_sum : log_term;
          } else {
            log_sum = max_val + std::log(
              std::exp(log_sum - max_val) + std::exp(log_term - max_val)
            );
          }
        }
      }
    }
  }

  double convolution_sum = std::exp(log_sum) / shared_divisor;
  return convolution_sum;
}










/*
 * This function computes the numerator of two by five table, which is, five potential outcomes and two treatments
 * To fully characterize a two by five table, we need to have (5-1)*(2-1) = 4 corners, which are n11, n12, and n13,n14
 * We have four allocations, u1, u2, and u3, and u4, and u5 for each Us loop
 * generic bias
 */


// [[Rcpp::export]]
double d_numerator_two_by_five(
    int n11, int n12, int n13, int n14,
    int n1,  int n2,
    int m1,  int m2,  int m3,  int m4,  int m5,
    int N,
    NumericVector gamma_delta,
    int u1, int u2, int u3, int u4, int u5,
    double shared_divisor)
{
  // 1) Check bounds for n11, n12, n13, n14
  if (n11 > std::min(n1, m1) || (n11 < std::max(0, n1 + m1 - N))) {
    return 0.0;
  }
  if (n12 > std::min(n1, m2) || (n12 < std::max(0, n1 + m2 - N))) {
    return 0.0;
  }
  if (n13 > std::min(n1, m3) || (n13 < std::max(0, n1 + m3 - N))) {
    return 0.0;
  }
  if (n14 > std::min(n1, m4) || (n14 < std::max(0, n1 + m4 - N))) {
    return 0.0;
  }

  double gamma_val = gamma_delta[0] - gamma_delta[1];

  if (std::fabs(gamma_val) < 1e-14) {
    double val = choose_C(m1, n11)
    * choose_C(m2, n12)
    * choose_C(m3, n13)
    * choose_C(m4, n14)
    * choose_C(m5, n1 - n11 - n12 - n13 - n14);
    return val / shared_divisor;
  }

  double log_sum = -std::numeric_limits<double>::infinity(); // log(0)

  int Us = u1 + u2 + u3 + u4 + u5;

  int lb_q1 = std::max(0, n1 + Us - N);
  int ub_q1 = std::min(Us, n1);

  for (int q1 = lb_q1; q1 <= ub_q1; ++q1) {

    int lb_s11 = std::max({0, n11 + u1 - m1});
    int ub_s11 = std::min({n11, u1});

    for (int s11 = lb_s11; s11 <= ub_s11; ++s11) {

      int lb_s12 = std::max({0, n12 + u2 - m2});
      int ub_s12 = std::min({n12, u2});

      for (int s12 = lb_s12; s12 <= ub_s12; ++s12) {

        int lb_s13 = std::max({0, n13 + u3 - m3});
        int ub_s13 = std::min({n13, u3});

        for (int s13 = lb_s13; s13 <= ub_s13; ++s13) {


          int lb_s14 = std::max({
            0,
            n14 + u4 - m4,
            q1 - u5 - s11 - s12 - s13,
            (q1 - n1 + n11 + n12 + n13 + n14 - s11 - s12 - s13)
          });
          int ub_s14 = std::min({
            n14,
            u4,
            (q1 - s11 - s12 - s13),
            (m5 - u5) - (n1 - (n11 + n12 + n13 + n14)) + (q1 - s11 - s12 - s13)
          });

          for (int s14 = lb_s14; s14 <= ub_s14; ++s14) {

            double lc1 = log_choose(u1,  s11);
            double lc2 = log_choose(m1 - u1,  n11 - s11);
            double lc3 = log_choose(u2,  s12);
            double lc4 = log_choose(m2 - u2,  n12 - s12);
            double lc5 = log_choose(u3,  s13);
            double lc6 = log_choose(m3 - u3,  n13 - s13);
            double lc7 = log_choose(u4,  s14);
            double lc8 = log_choose(m4 - u4,  n14 - s14);
            double lc9 = log_choose(u5, (q1 - s11 - s12 - s13 - s14));
            double lc10 = log_choose(
              (m5 - u5),
              (n1 - n11 - n12 - n13 - n14) - (q1 - s11 - s12 - s13 - s14)
            );

            double lc_exp = gamma_val * q1;

            double log_term = (lc1 + lc2 + lc3 + lc4 + lc5
                                 + lc6 + lc7 + lc8 + lc9 + lc10
                                 + lc_exp);

                                 double max_val = std::max(log_sum, log_term);
                                 if (std::isinf(max_val)) {
                                   log_sum = (log_sum > log_term) ? log_sum : log_term;
                                 } else {
                                   log_sum = max_val + std::log(
                                     std::exp(log_sum - max_val) + std::exp(log_term - max_val)
                                   );
                                 }

          }
        }
      }
    }
  }

  double convolution_sum = std::exp(log_sum) / shared_divisor;
  return convolution_sum;
}




/*
 * this is a function computing the numerator of the three treatment two outcome table
 * natural constraint include, m1 + m2 + m3 = N = n1+n2
 * generic bias
 */


// [[Rcpp::export]]
double d_numerator_three_by_two(
    int n11, int n21,
    int n1,  int n2,
    int n3,
    int m1,  int m2,
    int N,
    NumericVector gamma_delta,
    int u1,  int u2,
    double shared_divisor)
{
  if (n11 > std::min(n1, m1) || n11 < std::max(0, n1 + m1 - N)) return 0.0;
  if (n21 > std::min(n2, m1) || n21 < std::max(0, n2 + m1 - N)) return 0.0;

  double gamma1 = gamma_delta[0] - gamma_delta[2];
  double gamma2 = gamma_delta[1] - gamma_delta[2];
  double gamma_val = (std::fabs(gamma1) > 1e-14)
    ? gamma1
  : ((std::fabs(gamma2) > 1e-14) ? gamma2 : 0.0);

  if (std::fabs(gamma1) < 1e-14 && std::fabs(gamma2) < 1e-14) {
    double val = choose_C(m1, n11)
    * choose_C(m1 - n11, n21)
    * choose_C(m2, n1 - n11)
    * choose_C(m2 - (n1 - n11), (n2 - n21));
    return val / shared_divisor;
  }

  int Us = u1 + u2;
  int lb_q  = std::max(0, (n1 + n2 + Us) - N);
  int ub_q  = std::min(n1 + n2, Us);
  int lb_q1 = std::max(0, n1 + Us - N);
  int ub_q1 = std::min(Us, n1);
  int lb_q2 = std::max(0, n2 + Us - N);
  int ub_q2 = std::min(n2, Us);

  double log_sum = -std::numeric_limits<double>::infinity(); // log(0)

  for (int q = lb_q; q <= ub_q; ++q) {
    for (int q1 = lb_q1; q1 <= ub_q1; ++q1) {
      int q2 = q - q1;
      if (q2 < lb_q2 || q2 > ub_q2) {
        continue;
      }

      int lb_s11 = std::max({0, n11 + q1 - n1, n11 + u1 - m1, q1 - u2});
      int ub_s11 = std::min({n11, q1, u1, (m2 - u2) + n11 + q1 - n1});

      for (int s11 = lb_s11; s11 <= ub_s11; ++s11) {

        int lb_s21 = std::max({
          0,
          n21 + q2 - n2,
          (n11 - s11) + n21 + u1 - m1,
          (q1 - s11) + (q2) - u2
        });
        int ub_s21 = std::min({
          n21,
          q2,
          (u1 - s11),
          (m2 + n11 + n21 - n1 - n2 - u2 + q1 + q2 - s11)
        });

        for (int s21 = lb_s21; s21 <= ub_s21; ++s21) {

          double log_exp_term = 0.0;
          if (std::fabs(gamma1) < 1e-14) {
            log_exp_term = gamma_val * (q - q1);
          } else if (std::fabs(gamma2) < 1e-14) {
            log_exp_term = gamma_val * (q - q2);
          } else {
            log_exp_term = gamma_val * q;
          }


          double lc1 = log_choose(u1, s11);
          double lc2 = log_choose(m1 - u1, (n11 - s11));
          double lc3 = log_choose(u2, (q1 - s11));
          double lc4 = log_choose((u1 - s11), s21);
          double lc5 = log_choose((m1 - u1) - (n11 - s11), (n21 - s21));
          double lc6 = log_choose((u2 - (q1 - s11)), (q2 - s21));
          double lc7 = log_choose(m2 - u2, (n1 - n11 - q1 + s11));
          double lc8 = log_choose(
            (m2 - u2) - (n1 - n11 - q1 + s11),
            (n2 - n21 - q2 + s21)
          );

          double log_term = lc1 + lc2 + lc3 + lc4 + lc5 + lc6 + lc7 + lc8 + log_exp_term;

          double max_val = std::max(log_sum, log_term);
          if (std::isinf(max_val)) {
            log_sum = (log_sum > log_term) ? log_sum : log_term;
          } else {
            log_sum = max_val + std::log(
              std::exp(log_sum - max_val) + std::exp(log_term - max_val)
            );
          }
        }
      }
    }
  }

  double convolution_sum = std::exp(log_sum) / shared_divisor;
  return convolution_sum;
}



/*
 * The calculation of the numerator of three by three table, generic bias
 */

// [[Rcpp::export]]
double d_numerator_three_by_three(
    int n11, int n21, int n12, int n22,
    int n1,  int n2,  int n3,
    int m1,  int m2,  int m3,
    int N,
    NumericVector gamma_delta,
    int u1, int u2, int u3,
    double shared_divisor)
{
  if (n11 > std::min(n1, m1) || n11 < std::max(0, n1 + m1 - N)) return 0.0;
  if (n21 > std::min(n2, m1) || n21 < std::max(0, n2 + m1 - N)) return 0.0;
  if (n12 > std::min(n1, m2) || n12 < std::max(0, n1 + m2 - N)) return 0.0;
  if (n22 > std::min(n2, m2) || n22 < std::max(0, n2 + m2 - N)) return 0.0;

  double gamma1 = gamma_delta[0] - gamma_delta[2];
  double gamma2 = gamma_delta[1] - gamma_delta[2];
  double gamma_val = (std::fabs(gamma1) > 1e-14)
    ? gamma1
  : ((std::fabs(gamma2) > 1e-14) ? gamma2 : 0.0);

  if (std::fabs(gamma1) < 1e-14 && std::fabs(gamma2) < 1e-14) {
    double val = choose_C(m1, n11)
    * choose_C(m1 - n11, n21)
    * choose_C(m2, n12)
    * choose_C(m2 - n12, n22)
    * choose_C(m3,            (n1 - n11 - n12))
    * choose_C(m3 - (n1 - n11 - n12), (n2 - n21 - n22));
    return val / shared_divisor;
  }

  int Us = u1 + u2 + u3;

  int lb_q  = std::max(0, (n1 + n2 + Us) - N);
  int ub_q  = std::min(n1 + n2, Us);
  int lb_q1 = std::max(0, n1 + Us - N);
  int ub_q1 = std::min(Us, n1);
  int lb_q2 = std::max(0, n2 + Us - N);
  int ub_q2 = std::min(n2, Us);

  double log_sum = -std::numeric_limits<double>::infinity(); // log(0)

  for (int q = lb_q; q <= ub_q; ++q) {
    for (int q1 = lb_q1; q1 <= ub_q1; ++q1) {
      int q2 = q - q1;
      if (q2 < lb_q2 || q2 > ub_q2) {
        continue;
      }
      int lb_s11 = std::max(0, u1 - m1 + n11);
      int ub_s11 = std::min(n11, u1);
      for (int s11 = lb_s11; s11 <= ub_s11; ++s11) {


        int lb_s21 = std::max(0, u1 - m1 + n11 - s11 + n21);
        int ub_s21 = std::min(n21, u1 - s11);
        for (int s21 = lb_s21; s21 <= ub_s21; ++s21) {

          int lb_s12 = std::max({0, u2 - m2 + n12, q1 - s11 - u3, q1 - s11 - n1 + n11 + n12});
          int ub_s12 = std::min({n12, u2, q1 - s11, m3 - u3 + q1 - n1 + n11 + n12 - s11});
          for (int s12 = lb_s12; s12 <= ub_s12; ++s12) {

            int lb_s22 = std::max({0, u2 - m2 + n12 + n22 - s12, q1 - s11 - u3 + q2 - s21 - s12, q2 - s21 - n2 + n21 + n22});
            int ub_s22 = std::min({n22, u2 - s12, q2 - s21, m3 - u3 - (n1 - n11 - n12 - q1 + s11 + s12) - (n2 - n21 - n22 - q2) - s21});
            for (int s22 = lb_s22; s22 <= ub_s22; ++s22) {


              double log_exp_term = 0.0;
              if (std::fabs(gamma1) < 1e-14) {
                log_exp_term = gamma_val * (q - q1);
              } else if (std::fabs(gamma2) < 1e-14) {
                log_exp_term = gamma_val * (q - q2);
              } else {
                log_exp_term = gamma_val * q;
              }

              double lc1  = log_choose(u1, s11);
              double lc2  = log_choose((m1 - u1), (n11 - s11));
              double lc3  = log_choose((u1 - s11), s21);
              double lc4  = log_choose((m1 - u1) - (n11 - s11), (n21 - s21));

              double lc5  = log_choose(u2, s12);
              double lc6  = log_choose((m2 - u2), (n12 - s12));
              double lc7  = log_choose((u2 - s12), s22);
              double lc8  = log_choose((m2 - u2) - (n12 - s12), (n22 - s22));

              double lc9  = log_choose(u3, (q1 - s11 - s12));
              double lc10 = log_choose(
                (u3 - (q1 - s11 - s12)),
                (q2 - s21 - s22)
              );

              double lc11 = log_choose(
                (m3 - u3),
                (n1 - n11 - n12) - (q1 - s11 - s12)
              );
              double lc12 = log_choose(
                (m3 - u3) - ((n1 - n11 - n12) - (q1 - s11 - s12)),
                ((n2 - n21 - n22) - (q2 - s21 - s22))
              );

              double log_term = (lc1 + lc2 + lc3 + lc4
                                   + lc5 + lc6 + lc7 + lc8
                                   + lc9 + lc10 + lc11 + lc12
                                   + log_exp_term);


                                   double max_val = std::max(log_sum, log_term);
                                   if (std::isinf(max_val)) {
                                     log_sum = (log_sum > log_term) ? log_sum : log_term;
                                   } else {
                                     log_sum = max_val + std::log(
                                       std::exp(log_sum - max_val) + std::exp(log_term - max_val)
                                     );
                                   }
            }
          }
        }
      }
    }
  }

  double convolution_sum = std::exp(log_sum) / shared_divisor;
  return convolution_sum;
}





/*
 * This function calculates the numerator of four potential outcomes and three treatments
 * generic bias
 */


// [[Rcpp::export]]
double d_numerator_three_by_four(int n11, int n21, int n12, int n22, int n13, int n23, int n1, int n2, int n3, int m1, int m2, int m3, int m4, int N, NumericVector gamma_delta, int u1, int u2, int u3, int u4, double shared_divisor) {
  if (n11 > std::min(n1, m1) || n11 < std::max(0, n1 + m1 - N)) return 0.0;
  if (n21 > std::min(n2, m1) || n21 < std::max(0, n2 + m1 - N)) return 0.0;
  if (n12 > std::min(n1, m2) || n12 < std::max(0, n1 + m2 - N)) return 0.0;
  if (n22 > std::min(n2, m2) || n22 < std::max(0, n2 + m2 - N)) return 0.0;
  if (n13 > std::min(n1, m3) || n13 < std::max(0, n1 + m3 - N)) return 0.0;
  if (n23 > std::min(n2, m3) || n23 < std::max(0, n2 + m3 - N)) return 0.0;

  int Us = u1 + u2 + u3 + u4;
  double gamma1 = gamma_delta[0] - gamma_delta[2];
  double gamma2 = gamma_delta[1] - gamma_delta[2];
  double gamma_val = (std::fabs(gamma1) >1e-14) ? gamma1 : (std::fabs(gamma2) > 1e-14 ? gamma2 : 0.0);
  if(std::fabs(gamma1)<1e-14 && std::fabs(gamma2)<1e-14){
    double val = choose_C(m1, n11)*choose_C(m2, n12)*choose_C(m3,n13)* choose_C(m1-n11, n21) *
      choose_C(m2- n12, n22)* choose_C(m3-n13, n23)*choose_C(m4, n1 - n11 - n12 - n13)*
      choose_C(m4- n1 + n11 + n12 + n13, n2 - n21 - n22 - n23);
    return(val/shared_divisor);
  }else{



    int lb_q = std::max(0, n1 + n2 + Us - N);
    int ub_q = std::min(n1 + n2, Us);
    int lb_q1 = std::max(0, n1 + Us - N);
    int ub_q1 = std::min(Us, n1);
    int lb_q2 = std::max(0, n2 + Us - N);
    int ub_q2 = std::min(n2, Us);

    double log_sum = -std::numeric_limits<double>::infinity(); // log(0)

    for (int q = lb_q; q <= ub_q; ++q) {
      for (int q1 = lb_q1; q1 <= ub_q1; ++q1) {
        int q2 = q - q1;
        if (q2 > ub_q2 || q2 < lb_q2) continue;

        int lb_s11 = std::max(0, u1 - m1 + n11);
        int ub_s11 = std::min(n11, u1);
        for (int s11 = lb_s11; s11 <= ub_s11; ++s11) {

          int lb_s12 = std::max(0, u2 - m2 + n12);
          int ub_s12 = std::min(n12, u2);
          for (int s12 = lb_s12; s12 <= ub_s12; ++s12) {

            int lb_s13 = std::max({0, u3 - m3 + n13, q1 - s11 - s12 - u4, q1 - s11 - s12 - n1 + n11 + n12 + n13});
            int ub_s13 = std::min({n13, u3, q1 - s11 - s12, m4 - u4 + q1 - n1 + n11 + n12 + n13 - s11 - s12});
            for (int s13 = lb_s13; s13 <= ub_s13; ++s13) {

              int lb_s21 = std::max({0, u1 - m1 + n21 + n11 - s11});
              int ub_s21 = std::min({n21, u1 - s11});
              for (int s21 = lb_s21; s21 <= ub_s21; ++s21) {

                int lb_s22 = std::max({0, u2 - m2 + n22 + n12 - s12});
                int ub_s22 = std::min({n22, u2 - s12});
                for (int s22 = lb_s22; s22 <= ub_s22; ++s22) {

                  int lb_s23 = std::max({0, u3 - m3 + n23 + n13 - s13, q2 - s21 - s22 - u4 + q1 - s11 - s12 - s13, q2 - n2 + n21 + n22 + n23 - s21 - s22});
                  int ub_s23 = std::min({n23, u3 - s13, q2 - s21 - s22, m4 - u4 - n1 + n11 + n12 + n13 + q1 - s11 - s12 - s13 - n2 + n21 + n22 + n23 + q2 - s21 - s22});
                  for (int s23 = lb_s23; s23 <= ub_s23; ++s23) {

                    double log_exp_term = 0.0;
                    if (std::fabs(gamma1) < 1e-14) {
                      log_exp_term = gamma_val * (q - q1);
                    } else if (std::fabs(gamma2) < 1e-14) {
                      log_exp_term = gamma_val * (q - q2);
                    } else {
                      log_exp_term = gamma_val*q;
                    }

                      double lc1 = log_choose(u1, s11);
                      double lc2 = log_choose(m1 - u1, n11 - s11);
                      double lc3 = log_choose(u2, s12);
                      double lc4 = log_choose(m2 - u2, n12 - s12);
                      double lc5 = log_choose(u3, s13);
                      double lc6 = log_choose(m3 - u3, n13 - s13);
                      double lc7 = log_choose(u1 - s11, s21);
                      double lc8 = log_choose(m1 - u1 - n11 + s11, n21 - s21);
                      double lc9 = log_choose(u2 - s12, s22);
                      double lc10 = log_choose(m2 - u2 - n12 + s12, n22 - s22);
                      double lc11 = log_choose(u3 - s13, s23);
                      double lc12 = log_choose(m3 - u3 - n13 + s13, n23 - s23);
                      double lc13 = log_choose(u4, q1 - s11 - s12 - s13);
                      double lc14 = log_choose(m4 - u4, n1 - n11 - n12 - n13 - q1 + s11 + s12 + s13);
                      double lc15 = log_choose(u4 - (q1 - s11 - s12 - s13), q2 - s21 - s22 - s23);
                      double lc16 = log_choose(m4 - u4 - n1 + n11 + n12 + n13 + q1 - s11 - s12 - s13, n2 - n21 - n22 - n23 - q2 + s21 + s22 + s23);

                      double log_term = (lc1+lc2+lc3+lc4+lc5+lc6+lc7+lc8+lc9+lc10+lc11+lc12+lc13+lc14+lc15+lc16+log_exp_term);


                      double max_val = std::max(log_sum, log_term);
                      if(std::isinf(max_val)){
                        log_sum = (log_sum>log_term)? log_sum : log_term;

                      }else{
                        log_sum = max_val + std::log(std::exp(log_sum-max_val)+std::exp(log_term-max_val));
                      }

                  }
                }
              }
            }
          }
        }
      }
    }

    double convolution_sum = std::exp(log_sum)/shared_divisor;
    return convolution_sum;
  }
}







/*
 * This function calculates the probability of three treatment five outcome table,
 * assuming treatment in rows and outcome in column
 * generic bias
 */

// [[Rcpp::export]]
double d_numerator_three_by_five(int n11, int n21, int n12, int n22, int n13, int n23, int n14, int n24, int n1, int n2, int n3, int m1, int m2, int m3, int m4, int m5, int N, NumericVector gamma_delta, int u1, int u2, int u3, int u4, int u5, double shared_divisor) {
  if (n11 > std::min(n1, m1) || n11 < std::max(0, n1 + m1 - N)) return 0.0;
  if (n21 > std::min(n2, m1) || n21 < std::max(0, n2 + m1 - N)) return 0.0;
  if (n12 > std::min(n1, m2) || n12 < std::max(0, n1 + m2 - N)) return 0.0;
  if (n22 > std::min(n2, m2) || n22 < std::max(0, n2 + m2 - N)) return 0.0;
  if (n13 > std::min(n1, m3) || n13 < std::max(0, n1 + m3 - N)) return 0.0;
  if (n23 > std::min(n2, m3) || n23 < std::max(0, n2 + m3 - N)) return 0.0;
  if (n14 > std::min(n1, m4) || n14 < std::max(0, n1 + m4 - N)) return 0.0;
  if (n24 > std::min(n2, m4) || n24 < std::max(0, n2 + m4 - N)) return 0.0;

  int Us = u1 + u2 + u3 + u4 + u5;
  double gamma1 = gamma_delta[0] - gamma_delta[2];
  double gamma2 = gamma_delta[1] - gamma_delta[2];
  double gamma_val = (std::fabs(gamma1) > 1e-14) ? gamma1 : (std::fabs(gamma2) > 1e-14 ? gamma2 : 0.0);
  if(std::fabs(gamma1)<1e-14 && std::fabs(gamma2)<1e-14){
    double val = choose_C(m1, n11)*choose_C(m2, n12)* choose_C(m3, n13)*
      choose_C(m4, n14)* choose_C(m1- n11, n21)*choose_C(m2- n12, n22)*
      choose_C(m3- n13, n23)*choose_C(m4- n14, n24)*
      choose_C(m5, n1 - n11 - n12 - n13 - n14) *
      choose_C(m5- n1 + n11 + n12 + n13 + n14, n2 - n21 - n22 - n23 - n24);
    return(val/shared_divisor);

  }else{
    int lb_q = std::max(0, n1 + n2 + Us - N);
    int ub_q = std::min(n1 + n2, Us);
    int lb_q1 = std::max(0, n1 + Us - N);
    int ub_q1 = std::min(Us, n1);
    int lb_q2 = std::max(0, n2 + Us - N);
    int ub_q2 = std::min(n2, Us);

    double log_sum = -std::numeric_limits<double>::infinity(); // log(0)

    for (int q = lb_q; q <= ub_q; ++q) {
      for (int q1 = lb_q1; q1 <= ub_q1; ++q1) {
        int q2 = q - q1;
        if (q2 > ub_q2 || q2 < lb_q2) continue;

        int lb_s11 = std::max(0, u1 - m1 + n11);
        int ub_s11 = std::min(n11, u1);
        for (int s11 = lb_s11; s11 <= ub_s11; ++s11) {

          int lb_s12 = std::max(0, u2 - m2 + n12);
          int ub_s12 = std::min(n12, u2);
          for (int s12 = lb_s12; s12 <= ub_s12; ++s12) {

            int lb_s13 = std::max({0, u3 - m3 + n13});
            int ub_s13 = std::min({n13, u3});
            for (int s13 = lb_s13; s13 <= ub_s13; ++s13) {

              int lb_s14 = std::max({0, n14 + u4 - m4, q1 - u5 - s11 - s12 - s13, q1 - n1 + n11 + n12 + n13 + n14 - s11 - s12 - s13});
              int ub_s14 = std::min({n14, u4, q1 - s11 - s12 - s13, m5 - u5 - n1 + n11 + n12 + n13 + n14 + q1 - s11 - s12 - s13});
              for (int s14 = lb_s14; s14 <= ub_s14; ++s14) {

                int lb_s21 = std::max({0, u1 - m1 + n21 + n11 - s11});
                int ub_s21 = std::min({n21, u1 - s11});
                for (int s21 = lb_s21; s21 <= ub_s21; ++s21) {

                  int lb_s22 = std::max({0, u2 - m2 + n22 + n12 - s12});
                  int ub_s22 = std::min({n22, u2 - s12});
                  for (int s22 = lb_s22; s22 <= ub_s22; ++s22) {

                    int lb_s23 = std::max({0, u3 - m3 + n23 + n13 - s13});
                    int ub_s23 = std::min({n23, u3 - s13});
                    for (int s23 = lb_s23; s23 <= ub_s23; ++s23) {

                      int lb_s24 = std::max({0, u4 - m4 + n24 + n14 - s14, q2 - s21 - s22 - s23 - u5 + q1 - s11 - s12 - s13 - s14, q2 - n2 + n21 + n22 + n23 + n24 - s21 - s22 - s23});
                      int ub_s24 = std::min({n24, u4 - s14, q2 - s21 - s22 - s23, m5 - u5 - n1 + n11 + n12 + n13 + n14 + q1 - s11 - s12 - s13 - s14 - n2 + n21 + n22 + n23 + n24 + q2 - s21 - s22 - s23});
                      for (int s24 = lb_s24; s24 <= ub_s24; ++s24) {

                        double log_exp_term;
                        if (std::fabs(gamma1) < 1e-14) {
                          log_exp_term = gamma_val * (q - q1);
                        } else if (std::fabs(gamma2) < 1e-14) {
                          log_exp_term = gamma_val * (q - q2);
                        } else {
                          log_exp_term = gamma_val * q;
                        }
                        double lc1 = log_choose(u1, s11);
                        double lc2 = log_choose(m1 - u1, n11 - s11);
                        double lc3 = log_choose(u2, s12);
                        double lc4 = log_choose(m2 - u2, n12 - s12);
                        double lc5 = log_choose(u3, s13);
                        double lc6 = log_choose(m3 - u3, n13 - s13);
                        double lc7 = log_choose(u4, s14);
                        double lc8 = log_choose(m4 - u4, n14 - s14);
                        double lc9 = log_choose(u1 - s11, s21);
                        double lc10 = log_choose(m1 - u1 - n11 + s11, n21 - s21);
                        double lc11 = log_choose(u2 - s12, s22);
                        double lc12 = log_choose(m2 - u2 - n12 + s12, n22 - s22);
                        double lc13 = log_choose(u3 - s13, s23);
                        double lc14 = log_choose(m3 - u3 - n13 + s13, n23 - s23);
                        double lc15 = log_choose(u4 - s14, s24);
                        double lc16 = log_choose(m4 - u4 - n14 + s14, n24 - s24);
                        double lc17 = log_choose(u5, q1 - s11 - s12 - s13 - s14);
                        double lc18 = log_choose(m5 - u5, n1 - n11 - n12 - n13 - n14 - q1 + s11 + s12 + s13 + s14);
                        double lc19 = log_choose(u5 - (q1 - s11 - s12 - s13 - s14), q2 - s21 - s22 - s23 - s24);
                        double lc20 = log_choose(m5 - u5 - n1 + n11 + n12 + n13 + n14 + q1 - s11 - s12 - s13 - s14, n2 - n21 - n22 - n23 - n24 - q2 + s21 + s22 + s23 + s24);

                        double log_term = log_term= (lc1+lc2+lc3+lc4+lc5+lc6+lc7+lc8+lc9+lc10+lc11+lc12+lc13+lc14+lc15+lc16+lc17 + lc18+lc19+lc20+log_exp_term);


                        double max_val = std::max(log_sum, log_term);
                        if(std::isinf(max_val)){
                          log_sum = (log_sum>log_term)? log_sum : log_term;

                        }else{
                          log_sum = max_val + std::log(std::exp(log_sum-max_val)+std::exp(log_term-max_val));
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
    }
    double convolution_sum = std::exp(log_sum)/shared_divisor;
    return convolution_sum;
  }
}



/*
 * This function assumes the treatment in row and the outcome in column and solves the problem of
 * four treatment and two outcome numerator
 * generic bias
 */


// [[Rcpp::export]]
double d_numerator_four_by_two(int n11, int n21, int n31, int n1, int n2, int n3, int n4, int m1, int m2, int N, NumericVector gamma_delta, int u1, int u2, double shared_divisor) {
  // Verify if n11, n12, and n13 are in the support; if not, return zero
  if (n11 > std::min(n1, m1) || n11 < std::max(0, n1 + m1 - N)) return 0.0;
  if (n21 > std::min(n2, m1) || n21 < std::max(0, n2 + m1 - N)) return 0.0;
  if (n31 > std::min(n3, m1) || n31 < std::max(0, n3 + m1 - N)) return 0.0;

  int Us = u1 + u2;

  double gamma1 = gamma_delta[0] - gamma_delta[3];
  double gamma2 = gamma_delta[1] - gamma_delta[3];
  double gamma3 = gamma_delta[2] - gamma_delta[3];
  double gamma_val =
    (std::fabs(gamma1) > 1e-14) ? gamma1
  : ((std::fabs(gamma2) > 1e-14) ? gamma2
       : ((std::fabs(gamma3) > 1e-14) ? gamma3 : 0.0));

  if(std::fabs(gamma1)<1e-14 && std::fabs(gamma2)<1e-14 && std::fabs(gamma3)<1e-14){
    double val = choose_C(m1, n11)*choose_C(m1-n11, n21)*
      choose_C(m1-n11-n21, n31)*choose_C(m2, n1 - n11)*choose_C(m2- (n1 - n11), n2-n21)*
      choose_C(m2 - (n1 - n11) - (n2 - n21), n3 - n31);
    return(val/shared_divisor);

  }
  int lb_q = std::max(0, n1 + n2 + n3 + Us - N);
  int ub_q = std::min(n1 + n2 + n3, Us);
  int lb_q1 = std::max(0, n1 + Us - N);
  int ub_q1 = std::min(Us, n1);

  double log_sum = -std::numeric_limits<double>::infinity(); // log(0)


  for (int q = lb_q; q <= ub_q; ++q) {
    for (int q1 = lb_q1; q1 <= ub_q1; ++q1) {
      int lb_q2 = std::max(0, n1 + n2 - q1 + Us - N);
      int ub_q2 = std::min(n2, Us - q1);

      for (int q2 = lb_q2; q2 <= ub_q2; ++q2) {
        int lb_q3 = std::max(0, n1 + n2 + n3 - q1 - q2 + Us - N);
        int ub_q3 = std::min(Us - q1 - q2, n3);
        int q3 = q - q1 - q2;

        if (q3 > ub_q3 || q3 < lb_q3) {
          continue;
        } else {
          double log_exp_term;
          if (std::fabs(gamma1) < 1e-14 && std::fabs(gamma2) < 1e-14 && std::fabs(gamma3) < 1e-14) {
            log_exp_term = gamma_val * (q - q1 - q2 - q3);
          } else if (std::fabs(gamma1) < 1e-14 && std::fabs(gamma2) < 1e-14) {
            log_exp_term = gamma_val * (q - q1 - q2);
          } else if (std::fabs(gamma1) < 1e-14 && std::fabs(gamma3) < 1e-14) {
            log_exp_term = gamma_val * (q - q1 - q3);
          } else if (std::fabs(gamma2) < 1e-14 && std::fabs(gamma3) < 1e-14) {
            log_exp_term = gamma_val * (q - q2 - q3);
          } else if (std::fabs(gamma1) <  1e-14) {
            log_exp_term = gamma_val * (q - q1);
          } else if (std::fabs(gamma2) < 1e-14) {
            log_exp_term = gamma_val * (q - q2);
          } else if (std::fabs(gamma3) < 1e-14) {
            log_exp_term = gamma_val * (q - q3);
          } else {
            log_exp_term = gamma_val * q;
          }

          int lb_s11 = std::max({0, n11 + q1 - n1, n11 + u1 - m1, q1 - u2});
          int ub_s11 = std::min({n11, q1, u1, m2 - u2 + n11 + q1 - n1});
          for (int s11 = lb_s11; s11 <= ub_s11; ++s11) {

            int lb_s21 = std::max({0, n21 + q2 - n2, n11 - s11 + n21 + u1 - m1, q1 - s11 + q2 - u2});
            int ub_s21 = std::min({n21, q2, u1 - s11, m2 + n11 + n21 - n1 - n2 - u2 + q1 + q2 - s11});

            for (int s21 = lb_s21; s21 <= ub_s21; ++s21) {

              int lb_s31 = std::max({0, n31 + q3 - n3, n11 - s11 + n21 - s21 + n31 + u1 - m1, q3 + q1 - s11 + q2 - s21 - u2});
              int ub_s31 = std::min({n31, q3, u1 - s11 - s21, m2 + n11 + n21 + n31 - n1 - n2 - n3 - u2 + q1 + q2 + q3 - s11 - s21});

              for (int s31 = lb_s31; s31 <= ub_s31; ++s31) {
                double lc1 = log_choose(u1, s11);
                double lc2 = log_choose(m1 - u1, n11 - s11);
                double lc3 = log_choose(u2, q1 - s11);
                double lc4 = log_choose(u1 - s11, s21);
                double lc5 = log_choose(u1 - s11 - s21, s31);
                double lc6 = log_choose((m1 - u1) - (n11 - s11), n21 - s21);
                double lc7 = log_choose((m1 - u1) - (n11 - s11) - (n21 - s21), n31 - s31);
                double lc8 = log_choose(u2 - (q1 - s11), q2 - s21);
                double lc9 = log_choose(u2 - (q1 - s11) - (q2 - s21), q3 - s31);
                double lc10 =log_choose(m2 - u2, n1 - n11 - q1 + s11);
                double lc11 =log_choose(m2 - u2 - (n1 - n11 - q1 + s11), n2 - n21 - q2 + s21);
                double lc12 = log_choose(m2 - u2 - (n1 - n11 - q1 + s11) - (n2 - n21 - q2 + s21), n3 - n31 - q3 + s31);

                double log_term = log_term= (lc1+lc2+lc3+lc4+lc5+lc6+lc7+lc8+lc9+lc10+lc11+lc12+log_exp_term);


                double max_val = std::max(log_sum, log_term);
                if(std::isinf(max_val)){
                  log_sum = (log_sum>log_term)? log_sum : log_term;

                }else{
                  log_sum = max_val + std::log(std::exp(log_sum-max_val)+std::exp(log_term-max_val));
                }

              }
            }
          }
        }
      }
    }
  }

  double convolution_sum = std::exp(log_sum)/shared_divisor;
  return convolution_sum;
}








/*
 * This function assumes the treatment in rows and outcome in columns and solve
 * the problem of four treatment and three outcome table probability numerator
 * generic bias
 */


// [[Rcpp::export]]
double d_numerator_four_by_three(int n11, int n21, int n31, int n12, int n22, int n32, int n1, int n2, int n3, int n4, int m1, int m2, int m3, int N, NumericVector gamma_delta, int u1, int u2, int u3, double shared_divisor) {
  // Verify if six corners are in the support; if not, return zero
  if (n11 > std::min(n1, m1) || n11 < std::max(0, n1 + m1 - N)) return 0.0;
  if (n21 > std::min(n2, m1) || n21 < std::max(0, n2 + m1 - N)) return 0.0;
  if (n31 > std::min(n3, m1) || n31 < std::max(0, n3 + m1 - N)) return 0.0;
  if (n12 > std::min(n1, m2) || n12 < std::max(0, n1 + m2 - N)) return 0.0;
  if (n22 > std::min(n2, m2) || n22 < std::max(0, n2 + m2 - N)) return 0.0;
  if (n32 > std::min(n3, m2) || n32 < std::max(0, n3 + m2 - N)) return 0.0;

  int Us = u1 + u2 + u3;

  // Calculate gamma values
  double gamma1 = gamma_delta[0] - gamma_delta[3];
  double gamma2 = gamma_delta[1] - gamma_delta[3];
  double gamma3 = gamma_delta[2] - gamma_delta[3];
  double gamma_val = (std::fabs(gamma1) > 1e-14) ? gamma1
  : ((std::fabs(gamma2) > 1e-14) ? gamma2
       : ((std::fabs(gamma3) > 1e-14) ? gamma3 : 0.0));


  if(std::fabs(gamma1) <1e-14 && std::fabs(gamma2)<1e-14 && std::fabs(gamma3)<1e-14){
      double val = choose_C(m1, n11)*choose_C(m1- n11, n21)*choose_C(m1- n11 - n21, n31)*
      choose_C(m2, n12)* choose_C(m2- n12, n22)* choose_C(m2- n12- n22, n32)*
      choose_C(m3, n1 - n11 - n12)*choose_C(m3 - n1 + n11 + n12, n2 - n21 - n22)*
      choose_C(m3-u3-n1+n11+n12-n2+n21+n22,n3-n31-n32);
    return(val/shared_divisor);
  }
  int lb_q = std::max(0, n1 + n2 + n3 + Us - N);
  int ub_q = std::min(n1 + n2 + n3, Us);
  int lb_q1 = std::max(0, n1 + Us - N);
  int ub_q1 = std::min(Us, n1);

  double log_sum = -std::numeric_limits<double>::infinity(); // log(0)



  for (int q = lb_q; q <= ub_q; ++q) {
    for (int q1 = lb_q1; q1 <= ub_q1; ++q1) {
      int lb_q2 = std::max(0, n1 + n2 - q1 + Us - N);
      int ub_q2 = std::min(n2, Us - q1);

      for (int q2 = lb_q2; q2 <= ub_q2; ++q2) {
        int lb_q3 = std::max(0, n1 + n2 + n3 - q1 - q2 + Us - N);
        int ub_q3 = std::min(Us - q1 - q2, n3);
        int q3 = q - q1 - q2;

        if (q3 > ub_q3 || q3 < lb_q3) {
          continue;
        } else {
          double log_exp_term;
          if (std::fabs(gamma1) < 1e-14 && std::fabs(gamma2) < 1e-14 && std::fabs(gamma3) < 1e-14) {
            log_exp_term = gamma_val * (q - q1 - q2 - q3);
          } else if (std::fabs(gamma1) < 1e-14 && std::fabs(gamma2) < 1e-14) {
            log_exp_term = gamma_val * (q - q1 - q2);
          } else if (std::fabs(gamma1) < 1e-14 && std::fabs(gamma3) < 1e-14) {
            log_exp_term = gamma_val * (q - q1 - q3);
          } else if (std::fabs(gamma2) < 1e-14 && std::fabs(gamma3) < 1e-14) {
            log_exp_term = gamma_val * (q - q2 - q3);
          } else if (std::fabs(gamma1) < 1e-14) {
            log_exp_term = gamma_val * (q - q1);
          } else if (std::fabs(gamma2) < 1e-14) {
            log_exp_term = gamma_val * (q - q2);
          } else if (std::fabs(gamma3) < 1e-14) {
            log_exp_term = gamma_val * (q - q3);
          } else {
            log_exp_term = gamma_val * q;
          }
          int lb_s11 = std::max({0, n11 + u1 - m1});
          int ub_s11 = std::min({n11, u1});
          for (int s11 = lb_s11; s11 <= ub_s11; ++s11) {
            int lb_s21 = std::max({0, n11 - s11 + n21 + u1 - m1});
            int ub_s21 = std::min({n21, u1 - s11});

            for (int s21 = lb_s21; s21 <= ub_s21; ++s21) {
              int lb_s31 = std::max({0, n11 - s11 + n21 - s21 + n31 + u1 - m1});
              int ub_s31 = std::min({n31, u1 - s11 - s21});

              for (int s31 = lb_s31; s31 <= ub_s31; ++s31) {
                int lb_s12 = std::max({0, n12 + u2 - m2, q1 - u3 - s11, q1 - n1 + n11 + n12 - s11});
                int ub_s12 = std::min({n12, u2, q1 - s11, m3 - n1 + n11 + n12 - u3 + q1 - s11});

                for (int s12 = lb_s12; s12 <= ub_s12; ++s12) {
                  int lb_s22 = std::max({0, n22 + n12 - s12 + u2 - m2, q2 - s21 + q1 - s11 - s12 - u3, q2 - n2 + n21 + n22 - s21});
                  int ub_s22 = std::min({n22, u2 - s12, q2 - s21, m3 - u3 - n1 + n11 + n12 + q1 - s11 - s12 + q2 - n2 + n21 + n22 - s21});

                  for (int s22 = lb_s22; s22 <= ub_s22; ++s22) {
                    int lb_s32 = std::max({0, n12 - s12 + n22 - s22 + n32 + u2 - m2, q1 - s11 - s12 + q2 - s21 - s22 + q3 - s31 - u3, q3 - n3 + n31 + n32 - s31});
                    int ub_s32 = std::min({n32, u2 - s12 - s22, q3 - s31, m3 - u3 - n1 + n11 + n12 + q1 - s11 - s12 - n2 + n21 + n22 + q2 - s21 - s22 - n3 + n31 + n32 + q3 - s31});

                    for (int s32 = lb_s32; s32 <= ub_s32; ++s32) {
                      double lc1 = log_choose(u1,s11);
                      double lc2 = log_choose(m1-u1,n11-s11);
                      double lc3 = log_choose(u1-s11,s21);
                      double lc4 = log_choose(m1-u1-n11+s11,n21-s21);
                      double lc5 = log_choose(u1-s11-s21,s31);
                      double lc6 = log_choose(m1-u1-n11+s11-n21+s21,n31-s31);
                      double lc7 = log_choose(u2,s12);
                      double lc8 = log_choose(m2-u2,n12-s12);
                      double lc9 = log_choose(u2-s12,s22);
                      double lc10 = log_choose(m2-u2-n12+s12,n22-s22);
                      double lc11 = log_choose(u2-s12-s22,s32);
                      double lc12 = log_choose(m2-u2-n12+s12-n22+s22,n32-s32);
                      double lc13 = log_choose(u3,q1-s11-s12);
                      double lc14 = log_choose(m3-u3,n1-n11-n12-q1+s11+s12);
                      double lc15 = log_choose(u3-q1+s11+s12,q2-s21-s22);
                      double lc16 = log_choose(m3-u3-n1+n11+n12+q1-s11-s12,n2-n21-n22-q2+s21+s22);
                      double lc17 = log_choose(u3-q1+s11+s12-q2+s21+s22,q3-s31-s32);
                      double lc18 = log_choose(m3-u3-n1+n11+n12+q1-s11-s12-n2+n21+n22+q2-s21-s22,n3-n31-n32-q3+s31+s32);

                      double log_term = log_term= (lc1+lc2+lc3+lc4+lc5+lc6+lc7+lc8+lc9+lc10+lc11+lc12+lc13+lc14+lc15+lc16+lc17+lc18+log_exp_term);


                      double max_val = std::max(log_sum, log_term);
                      if(std::isinf(max_val)){
                        log_sum = (log_sum>log_term)? log_sum : log_term;

                      }else{
                        log_sum = max_val + std::log(std::exp(log_sum-max_val)+std::exp(log_term-max_val));
                      }                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  double convolution_sum = std::exp(log_sum)/shared_divisor;
  return convolution_sum;
}




/*
 * This function computes the numerator of five treatment and two outcome table
 * assuming that treatment in rows and outcome in column
 * generic bias
 */

// [[Rcpp::export]]
double d_numerator_five_by_two(int n1, int n2, int n3, int n4, int n5, int n11, int n21, int n31, int n41, int m1, int m2, int N, NumericVector gamma_delta, int u1, int u2, double shared_divisor) {
  // Verify if n11, n12, n13, and n14 are in the support; if not, return zero
  if (n11 > std::min(n1, m1) || n11 < std::max(0, n1 + m1 - N)) return 0.0;
  if (n21 > std::min(n2, m1) || n21 < std::max(0, n2 + m1 - N)) return 0.0;
  if (n31 > std::min(n3, m1) || n31 < std::max(0, n3 + m1 - N)) return 0.0;
  if (n41 > std::min(n4, m1) || n41 < std::max(0, n4 + m1 - N)) return 0.0;

  int Us = u1 + u2;

  // Calculate gamma values
  double gamma1 = gamma_delta[0] - gamma_delta[4];
  double gamma2 = gamma_delta[1] - gamma_delta[4];
  double gamma3 = gamma_delta[2] - gamma_delta[4];
  double gamma4 = gamma_delta[3] - gamma_delta[4];
  double gamma_val= 0.0;
  if (std::fabs(gamma1) > 1e-14) {
    gamma_val = gamma1;
  } else if (std::fabs(gamma2) > 1e-14) {
    gamma_val = gamma2;
  } else if (std::fabs(gamma3) > 1e-14) {
    gamma_val = gamma3;
  } else if (std::fabs(gamma4) > 1e-14) {
    gamma_val = gamma4;
  } else {
    gamma_val = 0.0;
  }
  if(std::fabs(gamma1)<1e-14 && std::fabs(gamma2)<1e-14 && std::fabs(gamma3)<1e-14 && std::fabs(gamma4)<1e-14){
    double val = choose_C(m1, n11)*choose_C(m1-n11, n21) *choose_C(m1-n11-n21,n31)*
      choose_C(m1-n11- n21-n31, n41)*choose_C(m2, n1 - n11) *
      choose_C(m2- (n1 - n11), n2 - n21) *
      choose_C(m2- (n1 - n11) - (n2 - n21), n3 - n31) *
      choose_C(m2 - (n1 - n11) - (n2 - n21) - (n3 - n31), n4 - n41);
    return(val/shared_divisor);

  }
  int lb_q = std::max(0, n1 + n2 + n3 + n4 + Us - N);
  int ub_q = std::min(n1 + n2 + n3 + n4, Us);
  int lb_q1 = std::max(0, n1 + Us - N);
  int ub_q1 = std::min(Us, n1);

  double log_sum = -std::numeric_limits<double>::infinity(); // log(0)

  for (int q = lb_q; q <= ub_q; ++q) {
    for (int q1 = lb_q1; q1 <= ub_q1; ++q1) {
      int lb_q2 = std::max(0, Us - n3 - n4 - n5 - q1);
      int ub_q2 = std::min(n2, Us - q1);

      for (int q2 = lb_q2; q2 <= ub_q2; ++q2) {
        int lb_q3 = std::max(0, Us - q1 - q2 - n4 - n5);
        int ub_q3 = std::min(n3, Us - q1 - q2);

        for (int q3 = lb_q3; q3 <= ub_q3; ++q3) {
          int lb_q4 = std::max(0, Us - q1 - q2 - q3 - n5);
          int ub_q4 = std::min(n4, Us - q1 - q2 - q3);
          int q4 = q - q1 - q2 - q3;

          if (q4 < lb_q4 || q4 > ub_q4) {
            continue;
          } else {
            double log_exp_term;
            if (std::fabs(gamma1) < 1e-14 && std::fabs(gamma2) < 1e-14 && std::fabs(gamma3) < 1e-14 && std::fabs(gamma4) < 1e-14) {
              log_exp_term = gamma_val * (q - q1 - q2 - q3 - q4);
            } else if (std::fabs(gamma1) < 1e-14 && std::fabs(gamma2) < 1e-14 && std::fabs(gamma3) < 1e-14) {
              log_exp_term = gamma_val * (q - q1 - q2 - q3);
            } else if (std::fabs(gamma1) < 1e-14 && std::fabs(gamma2) < 1e-14 && std::fabs(gamma4) < 1e-14) {
              log_exp_term = gamma_val * (q - q1 - q2 - q4);
            } else if (std::fabs(gamma1) < 1e-14 && std::fabs(gamma3) < 1e-14 && std::fabs(gamma4) < 1e-14) {
              log_exp_term = gamma_val * (q - q1 - q3 - q4);
            } else if (std::fabs(gamma2) < 1e-14 && std::fabs(gamma3) < 1e-14 && std::fabs(gamma4) < 1e-14) {
              log_exp_term = gamma_val * (q - q2 - q3 - q4);
            } else if (std::fabs(gamma1) < 1e-14 && std::fabs(gamma2) < 1e-14) {
              log_exp_term = gamma_val * (q - q1 - q2);
            } else if (std::fabs(gamma1) < 1e-14 && std::fabs(gamma3) < 1e-14) {
              log_exp_term = gamma_val* (q - q1 - q3);
            } else if (std::fabs(gamma1) < 1e-14 && std::fabs(gamma4) < 1e-14) {
              log_exp_term = gamma_val * (q - q1 - q4);
            } else if (std::fabs(gamma2) < 1e-14 && std::fabs(gamma3) < 1e-14) {
              log_exp_term = gamma_val * (q - q2 - q3);
            } else if (std::fabs(gamma2) < 1e-14 && std::fabs(gamma4) < 1e-14) {
              log_exp_term = gamma_val * (q - q2 - q4);
            } else if (std::fabs(gamma3) < 1e-14 && std::fabs(gamma4) < 1e-14) {
              log_exp_term = gamma_val * (q - q3 - q4);
            } else if (std::fabs(gamma1) < 1e-14) {
              log_exp_term = gamma_val * (q - q1);
            } else if (std::fabs(gamma2) < 1e-14) {
              log_exp_term = gamma_val * (q - q2);
            } else if (std::fabs(gamma3) < 1e-14) {
              log_exp_term = gamma_val * (q - q3);
            } else if (std::fabs(gamma4) < 1e-14) {
              log_exp_term = gamma_val * (q - q4);
            } else {
              log_exp_term = gamma_val * q;
            }

            int lb_s11 = std::max({0, n11 + q1 - n1, n11 + u1 - m1, q1 - u2});
            int ub_s11 = std::min({n11, q1, u1, m2 - u2 + n11 + q1 - n1});
            for (int s11 = lb_s11; s11 <= ub_s11; ++s11) {

              int lb_s21 = std::max({0, n21 + q2 - n2, n11 - s11 + n21 + u1 - m1, q1 - u2 - s11 + q2});
              int ub_s21 = std::min({n21, q2, u1 - s11, m2 + n11 + n21 - n1 - n2 - u2 + q1 + q2 - s11});

              for (int s21 = lb_s21; s21 <= ub_s21; ++s21) {

                int lb_s31 = std::max({0, n31 + q3 - n3, n11 - s11 + n21 - s21 + n31 + u1 - m1, q1 - u2 - s11 + q2 - s21 + q3});
                int ub_s31 = std::min({n31, q3, u1 - s11 - s21, m2 + n11 + n21 + n31 - n1 - n2 - n3 - u2 + q1 + q2 + q3 - s11 - s21});

                for (int s31 = lb_s31; s31 <= ub_s31; ++s31) {

                  int lb_s41 = std::max({0, n41 + q4 - n4, n11 - s11 + n21 - s21 + n31 - s31 + n41 + u1 - m1, q1 - u2 - s11 + q2 - s21 + q3 - s31 + q4});
                  int ub_s41 = std::min({n41, q4, u1 - s11 - s21 - s31, m2 + n11 + n21 + n31 + n41 - n1 - n2 - n3 - n4 - u2 + q1 + q2 + q3 + q4 - s11 - s21 - s31});

                  for (int s41 = lb_s41; s41 <= ub_s41; ++s41) {


                      double lc1 = log_choose(u1, s11);
                      double lc2 = log_choose(m1 - u1, n11 - s11);
                      double lc3 = log_choose(u2, q1 - s11);
                      double lc4 = log_choose(u1 - s11, s21);
                      double lc5 = log_choose(u1 - s11 - s21, s31);
                      double lc6 = log_choose(u1 - s11 - s21 - s31, s41);
                      double lc7 = log_choose((m1 - u1) - (n11 - s11), n21 - s21);
                      double lc8 = log_choose((m1 - u1) - (n11 - s11) - (n21 - s21), n31 - s31);
                      double lc9 = log_choose((m1 - u1) - (n11 - s11) - (n21 - s21) - (n31 - s31), n41 - s41);
                      double lc10 =log_choose(u2 - (q1 - s11), q2 - s21);
                      double lc11 =log_choose(u2 - (q1 - s11) - (q2 - s21), q3 - s31);
                      double lc12 =log_choose(u2 - (q1 - s11) - (q2 - s21) - (q3 - s31), q4 - s41);
                      double lc13 =log_choose(m2 - u2, n1 - n11 - q1 + s11);
                      double lc14 =log_choose(m2 - u2 - (n1 - n11 - q1 + s11), n2 - n21 - q2 + s21);
                      double lc15 =log_choose(m2 - u2 - (n1 - n11 - q1 + s11) - (n2 - n21 - q2 + s21), n3 - n31 - q3 + s31);
                      double lc16 = log_choose(m2 - u2 - (n1 - n11 - q1 + s11) - (n2 - n21 - q2 + s21) - (n3 - n31 - q3 + s31), n4 - n41 - q4 + s41);

                      double log_term = log_term= (lc1+lc2+lc3+lc4+lc5+lc6+lc7+lc8+lc9+lc10+lc11+lc12+lc13+lc14+lc15+lc16+log_exp_term);


                      double max_val = std::max(log_sum, log_term);
                      if(std::isinf(max_val)){
                        log_sum = (log_sum>log_term)? log_sum : log_term;

                      }else{
                        log_sum = max_val + std::log(std::exp(log_sum-max_val)+std::exp(log_term-max_val));
                      }                  }
                }
              }
            }
          }
        }
      }
    }
  }
  double convolution_sum = std::exp(log_sum)/shared_divisor;
  return convolution_sum;
}





/*The following function assists the computation for normal approximation
 *
 */

/*
 * This function calculates the first two noncentral moment in a two treatment
 * setup, namely, E(Z_1^Tu) and E((Z_1^TU)^2)
 * The first returned value is the expectation, the second returned value is the
 * second non-central moment of Z1^Tu
 */





/*
 * This function calculates the first two moments of z1u
 * which will be used in covariance computation.
 * The returned vector is (E(z1u), E((z1u)^2))
 */
// [[Rcpp::export]]
NumericVector zu_moments_two_treatment_log(int n1, int n2, int N, NumericVector gamma_delta, int Us, double shared_divisor) {
  double gamma = gamma_delta[0] - gamma_delta[1];
  double z1u_first = 0.0;
  double z1u_second = 0.0;

  // Sanity check
  if (Us == 0 || gamma == 0.0) {
    z1u_first = 0.0;
    z1u_second = 0.0;
    return NumericVector::create(z1u_first, z1u_second);
  }

  if (Us == N) {
    z1u_first = static_cast<double>(n1);
    z1u_second = std::pow(static_cast<double>(n1), 2);
    return NumericVector::create(z1u_first, z1u_second);
  }

  // Define bounds for q1
  int lb_q1 = std::max(0, n1 + Us - N);
  int ub_q1 = std::min(Us, n1);

  // Vectors to store log-weighted terms for convolution and each moment
  std::vector<double> log_weights;
  std::vector<double> log_z1u_terms;
  std::vector<double> log_z1u_sq_terms;

  // Iterate through all valid q1 to compute log-weighted terms
  for (int q1 = lb_q1; q1 <= ub_q1; ++q1) {
    // Compute the exponent term
    double exponent = gamma * q1;

    // Compute log(kernel)
    double log_kernel = log_choose(Us, q1) + log_choose(N - Us, n1 - q1) - std::log(shared_divisor);

    // Compute log(weighted_kernel) = exponent + log_kernel
    double log_weighted = exponent + log_kernel;

    // Store the log-weighted term for the convolution sum
    log_weights.push_back(log_weighted);

    // Store log(q1 * w) = log(q1) + log_weighted, handle q1=0
    if (q1 > 0) {
      log_z1u_terms.push_back(std::log(static_cast<double>(q1)) + log_weighted);
      log_z1u_sq_terms.push_back(std::log(static_cast<double>(q1) * static_cast<double>(q1)) + log_weighted);
    }
  }

  // Compute log-sum-exp for convolution sum and each moment
  double log_convolution_sum = compute_log_sum_exp(log_weights);
  double log_z1u_sum = compute_log_sum_exp(log_z1u_terms);
  double log_z1u_sq_sum = compute_log_sum_exp(log_z1u_sq_terms);

  // Check if convolution sum is zero (log_convolution_sum == -Inf)
  if (log_convolution_sum == -INFINITY) {
    Rcpp::stop("Convolution sum is zero; division by zero detected.");
  }

  // Compute the moments by exponentiating the differences
  // Handle cases where log_z1u_sum or log_z1u_sq_sum might be -Inf
  z1u_first = (log_z1u_sum == -INFINITY) ? 0.0 : std::exp(log_z1u_sum - log_convolution_sum);
  z1u_second = (log_z1u_sq_sum == -INFINITY) ? 0.0 : std::exp(log_z1u_sq_sum - log_convolution_sum);

  return NumericVector::create(z1u_first, z1u_second);
}



/*
 * This function calculates the first two moments of z1u, z2u, and their product
 * which will be used in the covariance computation
 * the returned vector is (E(z1u), E(z2u), E((z1u)^2), E((z2u)^2), E(z1u*z2u))
 */
// [[Rcpp::export]]
NumericVector zu_moments_three_treatment_log(int n1, int n2, int n3, int N, NumericVector gamma_delta, int Us, double shared_divisor) {
  double gamma1 = gamma_delta[0] - gamma_delta[2];
  double gamma2 = gamma_delta[1] - gamma_delta[2];
  double gamma = (gamma1 != 0.0) ? gamma1 : (gamma2 != 0.0 ? gamma2 : 0.0);

  // Sanity check
  if (Us == 0 || gamma == 0.0) {
    return NumericVector::create(0.0, 0.0, 0.0, 0.0, 0.0);
  }

  if (Us == N) {
    return NumericVector::create(
      static_cast<double>(n1),
      static_cast<double>(n2),
      static_cast<double>(n1 * n1),
      static_cast<double>(n2 * n2),
      static_cast<double>(n1 * n2)
    );
  }

  int lb_q = std::max(0, n1 + n2 + Us - N);
  int ub_q = std::min(n1 + n2, Us);
  int lb_q1 = std::max(0, n1 + Us - N);
  int ub_q1 = std::min(Us, n1);
  int lb_q2 = std::max(0, n2 + Us - N);
  int ub_q2 = std::min(n2, Us);

  // Vectors to store log-weighted terms for convolution and each moment
  std::vector<double> log_weights;
  std::vector<double> log_z1u_terms;
  std::vector<double> log_z1u_sq_terms;
  std::vector<double> log_z2u_terms;
  std::vector<double> log_z2u_sq_terms;
  std::vector<double> log_z1u_z2u_terms;

  // Iterate through all valid q and q1 to compute log-weighted terms
  for (int q = lb_q; q <= ub_q; ++q) {
    for (int q1 = lb_q1; q1 <= ub_q1; ++q1) {
      int q2 = q - q1;
      if (q2 < lb_q2 || q2 > ub_q2) continue;

      // Compute the exponent term based on gamma1 and gamma2
      double exponent;
      if (gamma1 == 0.0) {
        exponent = gamma * (q - q1);
      } else if (gamma2 == 0.0) {
        exponent = gamma * (q - q2);
      } else {
        exponent = gamma * q;
      }

      // Compute log(exp_term * kernel) = exponent + log(kernel)
      double log_kernel = log_choose(Us, q1) + log_choose(Us - q1, q2) +
        log_choose(N - Us, n1 - q1) + log_choose((N - Us) - (n1 - q1), n2 - q2) -
        std::log(shared_divisor);

      double log_weighted = exponent + log_kernel;

      // Store the log-weighted term for the convolution sum
      log_weights.push_back(log_weighted);

      // Store log(q1 * w) = log(q1) + log_weighted, handle q1=0
      if (q1 > 0) {
        log_z1u_terms.push_back(std::log(static_cast<double>(q1)) + log_weighted);
        log_z1u_sq_terms.push_back(std::log(static_cast<double>(q1 * q1)) + log_weighted);
      }

      // Store log(q2 * w) = log(q2) + log_weighted, handle q2=0
      if (q2 > 0) {
        log_z2u_terms.push_back(std::log(static_cast<double>(q2)) + log_weighted);
        log_z2u_sq_terms.push_back(std::log(static_cast<double>(q2 * q2)) + log_weighted);
      }

      // Store log(q1 * q2 * w) = log(q1) + log(q2) + log_weighted, handle q1=0 or q2=0
      if (q1 > 0 && q2 > 0) {
        log_z1u_z2u_terms.push_back(std::log(static_cast<double>(q1)) + std::log(static_cast<double>(q2)) + log_weighted);
      }
    }
  }

  // Compute log-sum-exp for convolution sum and each moment
  double log_convolution_sum = compute_log_sum_exp(log_weights);
  double log_z1u_sum = compute_log_sum_exp(log_z1u_terms);
  double log_z1u_sq_sum = compute_log_sum_exp(log_z1u_sq_terms);
  double log_z2u_sum = compute_log_sum_exp(log_z2u_terms);
  double log_z2u_sq_sum = compute_log_sum_exp(log_z2u_sq_terms);
  double log_z1u_z2u_sum = compute_log_sum_exp(log_z1u_z2u_terms);

  // Check if convolution sum is zero (log_convolution_sum == -Inf)
  if (log_convolution_sum == -INFINITY) {
    Rcpp::stop("Convolution sum is zero; division by zero detected.");
  }

  // Compute the moments by exponentiating the differences
  double E_z1u = (log_z1u_sum == -INFINITY) ? 0.0 : std::exp(log_z1u_sum - log_convolution_sum);
  double E_z2u = (log_z2u_sum == -INFINITY) ? 0.0 : std::exp(log_z2u_sum - log_convolution_sum);
  double E_z1u_sq = (log_z1u_sq_sum == -INFINITY) ? 0.0 : std::exp(log_z1u_sq_sum - log_convolution_sum);
  double E_z2u_sq = (log_z2u_sq_sum == -INFINITY) ? 0.0 : std::exp(log_z2u_sq_sum - log_convolution_sum);
  double E_z1u_z2u = (log_z1u_z2u_sum == -INFINITY) ? 0.0 : std::exp(log_z1u_z2u_sum - log_convolution_sum);

  return NumericVector::create(E_z1u, E_z2u, E_z1u_sq, E_z2u_sq, E_z1u_z2u);
}





/*
 * This function obtains the moments for z1u, z2u, z3u, and their pairwise products.
 * The returned vector contains:
 * (E(z1u), E(z2u), E(z3u),
 *  E((z1u)^2), E((z2u)^2), E((z3u)^2),
 *  E(z1u*z2u), E(z1u*z3u), E(z2u*z3u))
 *  The output is the same as zu_moments_four_treatment
 *  But this version is more numerically stable
 */
// [[Rcpp::export]]
NumericVector zu_moments_four_treatment_log(int n1, int n2, int n3, int n4, int N, NumericVector gamma_delta, int Us, double shared_divisor) {

  // Calculate gamma values
  double gamma1 = gamma_delta[0] - gamma_delta[3];
  double gamma2 = gamma_delta[1] - gamma_delta[3];
  double gamma3 = gamma_delta[2] - gamma_delta[3];

  // Determine the primary gamma to use for the exponential term
  // Following the original logic: prioritize gamma1, then gamma2, then gamma3
  double gamma = (gamma1 != 0.0) ? gamma1 : ((gamma2 != 0.0) ? gamma2 : gamma3);

  // Initialize cumulants for first and second moments, and product moments
  // First moments
  std::vector<double> log_z1u_terms;
  std::vector<double> log_z2u_terms;
  std::vector<double> log_z3u_terms;

  // Second moments
  std::vector<double> log_z1u_sq_terms;
  std::vector<double> log_z2u_sq_terms;
  std::vector<double> log_z3u_sq_terms;

  // Product moments
  std::vector<double> log_z1u_z2u_terms;
  std::vector<double> log_z1u_z3u_terms;
  std::vector<double> log_z2u_z3u_terms;

  // Vector to store log-weighted terms for convolution sum
  std::vector<double> log_weights;

  // Sanity check
  if (Us == 0 || gamma == 0.0) {
    return NumericVector::create(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
  }

  if (Us == N) {
    // When Us == N, all units are treated; moments are deterministic
    double z1u_first = static_cast<double>(n1);
    double z2u_first = static_cast<double>(n2);
    double z3u_first = static_cast<double>(n3);

    double z1u_second = std::pow(static_cast<double>(n1), 2);
    double z2u_second = std::pow(static_cast<double>(n2), 2);
    double z3u_second = std::pow(static_cast<double>(n3), 2);

    double z1u_z2u = static_cast<double>(n1) * static_cast<double>(n2);
    double z1u_z3u = static_cast<double>(n1) * static_cast<double>(n3);
    double z2u_z3u = static_cast<double>(n2) * static_cast<double>(n3);

    return NumericVector::create(z1u_first, z2u_first, z3u_first,
                                 z1u_second, z2u_second, z3u_second,
                                 z1u_z2u, z1u_z3u, z2u_z3u);
  }

  // Define bounds for q1, q2, q3
  int lb_q1 = std::max(0, n1 + Us - N);
  int ub_q1 = std::min(Us, n1);

  for (int q1 = lb_q1; q1 <= ub_q1; ++q1) {
    // After selecting q1, define bounds for q2
    int lb_q2 = std::max(0, n2 + Us - N - q1);
    int ub_q2 = std::min(Us - q1, n2);

    for (int q2 = lb_q2; q2 <= ub_q2; ++q2) {
      // After selecting q1 and q2, define bounds for q3
      int lb_q3 = std::max(0, n3 + Us - N - q1 - q2);
      int ub_q3 = std::min(Us - q1 - q2, n3);

      for (int q3 = lb_q3; q3 <= ub_q3; ++q3) {
        // Compute q4 based on q1, q2, q3
        int q4 = Us - q1 - q2 - q3;

        // Ensure q4 is within valid bounds
        if (q4 < 0 || q4 > n4) continue;

        // Compute the exponent term: gamma1*q1 + gamma2*q2 + gamma3*q3
        double exponent = gamma1 * q1 + gamma2 * q2 + gamma3 * q3;

        // Compute log(kernel)
        double log_kernel = log_choose(Us, q1) +
          log_choose(Us - q1, q2) +
          log_choose(Us - q1 - q2, q3) +
          log_choose(N - Us, n1 - q1) +
          log_choose(N - Us - (n1 - q1), n2 - q2) +
          log_choose(N - Us - (n1 - q1) - (n2 - q2), n3 - q3) -
          std::log(shared_divisor);

        // Compute log(weighted_kernel) = exponent + log_kernel
        double log_weighted = exponent + log_kernel;

        // Store the log-weighted term for the convolution sum
        log_weights.push_back(log_weighted);

        // Store log terms for first moments
        if (q1 > 0) {
          log_z1u_terms.push_back(std::log(static_cast<double>(q1)) + log_weighted);
        }
        if (q2 > 0) {
          log_z2u_terms.push_back(std::log(static_cast<double>(q2)) + log_weighted);
        }
        if (q3 > 0) {
          log_z3u_terms.push_back(std::log(static_cast<double>(q3)) + log_weighted);
        }

        // Store log terms for second moments
        if (q1 > 0) {
          log_z1u_sq_terms.push_back(std::log(static_cast<double>(q1) * static_cast<double>(q1)) + log_weighted);
        }
        if (q2 > 0) {
          log_z2u_sq_terms.push_back(std::log(static_cast<double>(q2) * static_cast<double>(q2)) + log_weighted);
        }
        if (q3 > 0) {
          log_z3u_sq_terms.push_back(std::log(static_cast<double>(q3) * static_cast<double>(q3)) + log_weighted);
        }

        // Store log terms for product moments
        if (q1 > 0 && q2 > 0) {
          log_z1u_z2u_terms.push_back(std::log(static_cast<double>(q1)) + std::log(static_cast<double>(q2)) + log_weighted);
        }
        if (q1 > 0 && q3 > 0) {
          log_z1u_z3u_terms.push_back(std::log(static_cast<double>(q1)) + std::log(static_cast<double>(q3)) + log_weighted);
        }
        if (q2 > 0 && q3 > 0) {
          log_z2u_z3u_terms.push_back(std::log(static_cast<double>(q2)) + std::log(static_cast<double>(q3)) + log_weighted);
        }
      }
    }
  }

  // Compute log-sum-exp for convolution sum and each moment
  double log_convolution_sum = compute_log_sum_exp(log_weights);
  double log_z1u_sum = compute_log_sum_exp(log_z1u_terms);
  double log_z2u_sum = compute_log_sum_exp(log_z2u_terms);
  double log_z3u_sum = compute_log_sum_exp(log_z3u_terms);

  double log_z1u_sq_sum = compute_log_sum_exp(log_z1u_sq_terms);
  double log_z2u_sq_sum = compute_log_sum_exp(log_z2u_sq_terms);
  double log_z3u_sq_sum = compute_log_sum_exp(log_z3u_sq_terms);

  double log_z1u_z2u_sum = compute_log_sum_exp(log_z1u_z2u_terms);
  double log_z1u_z3u_sum = compute_log_sum_exp(log_z1u_z3u_terms);
  double log_z2u_z3u_sum = compute_log_sum_exp(log_z2u_z3u_terms);

  // Check if convolution sum is zero (log_convolution_sum == -Inf)
  if (log_convolution_sum == -INFINITY) {
    Rcpp::stop("Convolution sum is zero; division by zero detected.");
  }

  // Compute the moments by exponentiating the differences
  // Handle cases where log_sums might be -Inf
  double E_z1u = (log_z1u_sum == -INFINITY) ? 0.0 : std::exp(log_z1u_sum - log_convolution_sum);
  double E_z2u = (log_z2u_sum == -INFINITY) ? 0.0 : std::exp(log_z2u_sum - log_convolution_sum);
  double E_z3u = (log_z3u_sum == -INFINITY) ? 0.0 : std::exp(log_z3u_sum - log_convolution_sum);

  double E_z1u_sq = (log_z1u_sq_sum == -INFINITY) ? 0.0 : std::exp(log_z1u_sq_sum - log_convolution_sum);
  double E_z2u_sq = (log_z2u_sq_sum == -INFINITY) ? 0.0 : std::exp(log_z2u_sq_sum - log_convolution_sum);
  double E_z3u_sq = (log_z3u_sq_sum == -INFINITY) ? 0.0 : std::exp(log_z3u_sq_sum - log_convolution_sum);

  double E_z1u_z2u = (log_z1u_z2u_sum == -INFINITY) ? 0.0 : std::exp(log_z1u_z2u_sum - log_convolution_sum);
  double E_z1u_z3u = (log_z1u_z3u_sum == -INFINITY) ? 0.0 : std::exp(log_z1u_z3u_sum - log_convolution_sum);
  double E_z2u_z3u = (log_z2u_z3u_sum == -INFINITY) ? 0.0 : std::exp(log_z2u_z3u_sum - log_convolution_sum);

  return NumericVector::create(E_z1u, E_z2u, E_z3u,
                               E_z1u_sq, E_z2u_sq, E_z3u_sq,
                               E_z1u_z2u, E_z1u_z3u, E_z2u_z3u);
}




/*
 * This function computes first two central moments of z1^TR_1 for two by two table
 */


// [[Rcpp::export]]
NumericVector zr_moments_two_by_two(int n1, int n2, int m1, int m2, int N, NumericVector gamma_delta, int u1, int u2, double shared_divisor) {
  NumericVector z1u_moments = zu_moments_two_treatment_log(n1,n2,N,gamma_delta,u1+u2,shared_divisor);
  double z1u_first= z1u_moments[0];
  double z1u_second= z1u_moments[1];
  double Us = u1 + u2;
  //if Us = 0 or N, becomes the usual RCT case
  if(Us==0 || Us==N)
    return(NumericVector::create((double) m1*n1/N,(double) n1*n2*m1*m2/N/N/(N-1)));

  double q11 = u1/Us;
  double q10 = (m1-u1)/(N-Us);
  double z1r1_first = q11*z1u_first+q10*(n1-z1u_first);
  double w11, w10;
  if(Us==1)
    w11 = 0;
  else
    w11  = (double) u1/(Us-1)*(1-2*q11) + (double) Us/(Us-1)*q11*q11;

  if(Us==N-1)
    w10 = 0;
  else
    w10   = (double)(m1 - u1)/(N-Us-1)*(1-2*q10) + (double) (N-Us)/(N-Us-1)*pow(q10,2);
  double z1r1_var = (w11-w10)*z1u_first-(z1u_second)*(w11/Us+w10/(N-Us))
    +n1*(N-Us-n1+2*z1u_first)*w10/(N-Us)+(z1u_second-z1u_first*z1u_first)*(q11-q10)*(q11-q10);

  return(NumericVector::create(z1r1_first,z1r1_var));

}

//This function computes the moments of (E(z1r1),E(z1r2),var(z1r1),var(z1r2),cov(z1r1,z1r2))
// [[Rcpp::export]]
NumericVector zr_moments_two_by_three(int n1, int n2, int m1, int m2, int m3, int N, NumericVector gamma_delta, int u1, int u2, int u3, double shared_divisor) {
  NumericVector z1u_moments = zu_moments_two_treatment_log(n1,n2,N,gamma_delta,u1+u2+u3,shared_divisor);
  double z1u_first= z1u_moments[0];
  double z1u_second= z1u_moments[1];
  double Us = u1 + u2+u3;
  //if Us = 0 or N, becomes the usual RCT case
  if(Us==0 || Us==N)
    return(NumericVector::create((double) m1*n1/N,(double) m2*n1/N, (double) n1*(N-n1)*m1*(N-m1)/N/N/(N-1),(double) n1*(N-n1)*m2*(N-m2)/N/N/(N-1),(double) m1*m2*n1*(n1-N)/N/N/(N-1)));

  double q11 = u1/Us;
  double q10 = (m1-u1)/(N-Us);
  double q21 = u2/Us;
  double q20 = (m2-u2)/(N-Us);
  double z1r1_first = q11*z1u_first+q10*(n1-z1u_first);
  double z1r2_first = q21*z1u_first+q20*(n1-z1u_first);
  double w11, w10, w21, w20;
  if(Us==1)
  { w11 = 0;
    w21 = 0; }
  else
  {w11  = (double) u1/(Us-1)*(1-2*q11) + (double) Us/(Us-1)*q11*q11;
    w21  = (double) u2/(Us-1)*(1-2*q21) + (double) Us/(Us-1)*q21*q21;}
  if(Us==N-1)
  {w10 = 0;
    w20 = 0;
  }
  else
  {
    w10   = (double)(m1 - u1)/(N-Us-1)*(1-2*q10) + (double) (N-Us)/(N-Us-1)*pow(q10,2);
    w20   = (double)(m2 - u2)/(N-Us-1)*(1-2*q20) + (double) (N-Us)/(N-Us-1)*pow(q20,2);

  }
  double z1r1_var = (w11-w10)*z1u_first-(z1u_second)*(w11/Us+w10/(N-Us))
    +n1*(N-Us-n1+2*z1u_first)*w10/(N-Us)+(z1u_second-z1u_first*z1u_first)*(q11-q10)*(q11-q10);
  double z1r2_var = (w21-w20)*z1u_first-(z1u_second)*(w21/Us+w20/(N-Us))
    +n1*(N-Us-n1+2*z1u_first)*w20/(N-Us)+(z1u_second-z1u_first*z1u_first)*(q21-q20)*(q21-q20);
  double z1r1_z1r2_cov;
  if(Us>1 && Us <N-1){
    z1r1_z1r2_cov = q11*q21/(Us-1)*(z1u_second-Us*z1u_first)+q10*q20/(N-Us-1)*n1*(n1-N+Us)
    +q10*q20/(N-Us-1)*(z1u_second-(2*n1-N+Us)*z1u_first)
    +(q11*q21+q10*q20-q10*q21-q20*q11)*(z1u_second-z1u_first*z1u_first);
  }
  if(Us==1){
    z1r1_z1r2_cov = q10*q20/(N-Us-1)*n1*(n1-N+Us)
    +q10*q20/(N-Us-1)*(z1u_second-(2*n1-N+Us)*z1u_first)
    +(q11*q21+q10*q20-q10*q21-q20*q11)*(z1u_second-z1u_first*z1u_first);

  }
  if(Us==N-1){
    z1r1_z1r2_cov = q11*q21/(Us-1)*(z1u_second-Us*z1u_first)
    +(q11*q21+q10*q20-q10*q21-q20*q11)*(z1u_second-z1u_first*z1u_first);
  }


  return(NumericVector::create(z1r1_first,z1r2_first, z1r1_var, z1r2_var, z1r1_z1r2_cov));

}




//This function computes the moments of (E(z1r1),E(z1r2),E(z1r3), var(z1r1),var(z1r2), var(z1r3), cov(z1r1,z1r2), cov(z1r1, z1r3), cov(z1r2, z1r3))
// [[Rcpp::export]]
NumericVector zr_moments_two_by_four(int n1, int n2, int m1, int m2, int m3, int m4, int N, NumericVector gamma_delta, int u1, int u2, int u3, int u4, double shared_divisor) {
  NumericVector z1u_moments = zu_moments_two_treatment_log(n1,n2,N,gamma_delta,u1+u2+u3+u4,shared_divisor);
  double z1u_first= z1u_moments[0];
  double z1u_second= z1u_moments[1];
  double Us = u1 + u2+u3+u4;
  //if Us = 0 or N, becomes the usual RCT case
  if(Us==0 || Us==N)
    return(NumericVector::create((double) m1*n1/N,(double) m2*n1/N, (double) m3*n1/N, (double) n1*(N-n1)*m1*(N-m1)/N/N/(N-1), (double) n1*(N-n1)*m2*(N-m2)/N/N/(N-1),
                                 (double) n1*(N-n1)*m3*(N-m3)/N/N/(N-1), (double) m1*m2*n1*(n1-N)/N/N/(N-1),(double) m1*m3*n1*(n1-N)/N/N/(N-1),(double) m2*m3*n1*(n1-N)/N/N/(N-1)));

  double q11 = u1/Us;
  double q10 = (m1-u1)/(N-Us);
  double q21 = u2/Us;
  double q20 = (m2-u2)/(N-Us);
  double q31 = u3/Us;
  double q30 = (m3-u3)/(N-Us);
  double z1r1_first = q11*z1u_first+q10*(n1-z1u_first);
  double z1r2_first = q21*z1u_first+q20*(n1-z1u_first);
  double z1r3_first = q31*z1u_first+q30*(n1-z1u_first);
  double w11, w10, w21, w20, w31, w30;
  if(Us==1)
  { w11 = 0;
    w21 = 0;
    w31 = 0; }
  else
  {w11  = (double) u1/(Us-1)*(1-2*q11) + (double) Us/(Us-1)*q11*q11;
    w21  = (double) u2/(Us-1)*(1-2*q21) + (double) Us/(Us-1)*q21*q21;
    w31  = (double) u3/(Us-1)*(1-2*q31) + (double) Us/(Us-1)*q31*q31;}
  if(Us==N-1)
  {w10 = 0;
    w20 = 0;
    w30 = 0;
  }
  else
  {
    w10   = (double)(m1 - u1)/(N-Us-1)*(1-2*q10) + (double) (N-Us)/(N-Us-1)*pow(q10,2);
    w20   = (double)(m2 - u2)/(N-Us-1)*(1-2*q20) + (double) (N-Us)/(N-Us-1)*pow(q20,2);
    w30   = (double)(m3 - u3)/(N-Us-1)*(1-2*q30) + (double) (N-Us)/(N-Us-1)*pow(q30,2);


  }
  double z1r1_var = (w11-w10)*z1u_first-(z1u_second)*(w11/Us+w10/(N-Us))
    +n1*(N-Us-n1+2*z1u_first)*w10/(N-Us)+(z1u_second-z1u_first*z1u_first)*(q11-q10)*(q11-q10);
  double z1r2_var = (w21-w20)*z1u_first-(z1u_second)*(w21/Us+w20/(N-Us))
    +n1*(N-Us-n1+2*z1u_first)*w20/(N-Us)+(z1u_second-z1u_first*z1u_first)*(q21-q20)*(q21-q20);
  double z1r3_var = (w31-w30)*z1u_first-(z1u_second)*(w31/Us+w30/(N-Us))
    +n1*(N-Us-n1+2*z1u_first)*w30/(N-Us)+(z1u_second-z1u_first*z1u_first)*(q31-q30)*(q31-q30);
  double z1r1_z1r2_cov;
  double z1r1_z1r3_cov;
  double z1r2_z1r3_cov;
  if(Us>1 && Us <N-1){
    z1r1_z1r2_cov = q11*q21/(Us-1)*(z1u_second-Us*z1u_first)+q10*q20/(N-Us-1)*n1*(n1-N+Us)
    +q10*q20/(N-Us-1)*(z1u_second-(2*n1-N+Us)*z1u_first)
    +(q11*q21+q10*q20-q10*q21-q20*q11)*(z1u_second-z1u_first*z1u_first);

    z1r1_z1r3_cov = q11*q31/(Us-1)*(z1u_second-Us*z1u_first)+q10*q30/(N-Us-1)*n1*(n1-N+Us)
      +q10*q30/(N-Us-1)*(z1u_second-(2*n1-N+Us)*z1u_first)
      +(q11*q31+q10*q30-q10*q31-q30*q11)*(z1u_second-z1u_first*z1u_first);

      z1r2_z1r3_cov = q21*q31/(Us-1)*(z1u_second-Us*z1u_first)+q20*q30/(N-Us-1)*n1*(n1-N+Us)
        +q20*q30/(N-Us-1)*(z1u_second-(2*n1-N+Us)*z1u_first)
        +(q21*q31+q20*q30-q20*q31-q30*q21)*(z1u_second-z1u_first*z1u_first);
  }
  if(Us==1){
    z1r1_z1r2_cov = q10*q20/(N-Us-1)*n1*(n1-N+Us)
    +q10*q20/(N-Us-1)*(z1u_second-(2*n1-N+Us)*z1u_first)
    +(q11*q21+q10*q20-q10*q21-q20*q11)*(z1u_second-z1u_first*z1u_first);

    z1r1_z1r3_cov = q10*q30/(N-Us-1)*n1*(n1-N+Us)
      +q10*q30/(N-Us-1)*(z1u_second-(2*n1-N+Us)*z1u_first)
      +(q11*q31+q10*q30-q10*q31-q30*q11)*(z1u_second-z1u_first*z1u_first);

      z1r2_z1r3_cov = q20*q30/(N-Us-1)*n1*(n1-N+Us)
        +q20*q30/(N-Us-1)*(z1u_second-(2*n1-N+Us)*z1u_first)
        +(q21*q31+q20*q30-q20*q31-q30*q21)*(z1u_second-z1u_first*z1u_first);

  }
  if(Us==N-1){
    z1r1_z1r2_cov = q11*q21/(Us-1)*(z1u_second-Us*z1u_first)
    +(q11*q21+q10*q20-q10*q21-q20*q11)*(z1u_second-z1u_first*z1u_first);

    z1r1_z1r3_cov = q11*q31/(Us-1)*(z1u_second-Us*z1u_first)
      +(q11*q31+q10*q30-q10*q31-q30*q11)*(z1u_second-z1u_first*z1u_first);

    z1r2_z1r3_cov = q21*q31/(Us-1)*(z1u_second-Us*z1u_first)
      +(q21*q31+q20*q30-q20*q31-q30*q21)*(z1u_second-z1u_first*z1u_first);
  }


  return(NumericVector::create(z1r1_first,z1r2_first, z1r3_first, z1r1_var, z1r2_var, z1r3_var, z1r1_z1r2_cov, z1r1_z1r3_cov, z1r2_z1r3_cov));

}




//This function computes the moments of (E(z1r1),E(z1r2),E(z1r3), E(z1r4), var(z1r1),var(z1r2), var(z1r3),  var(z1r4), cov(z1r1,z1r2), cov(z1r1, z1r3), cov(z1r1, z1r4), cov(z1r2, z1r3), cov(z1r2,z1r4), cov(z1r3, z1r4))
// [[Rcpp::export]]
NumericVector zr_moments_two_by_five(int n1, int n2, int m1, int m2, int m3, int m4, int m5, int N, NumericVector gamma_delta, int u1, int u2, int u3, int u4,int u5, double shared_divisor) {
  NumericVector z1u_moments = zu_moments_two_treatment_log(n1,n2,N,gamma_delta,u1+u2+u3+u4+u5,shared_divisor);
  double z1u_first= z1u_moments[0];
  double z1u_second= z1u_moments[1];
  double Us = u1 + u2+u3+u4+u5;
  //if Us = 0 or N, becomes the usual RCT case
  if(Us==0 || Us==N)
    return(NumericVector::create((double) m1*n1/N,(double) m2*n1/N, (double) m3*n1/N,(double) m4*n1/N,(double) n1*(N-n1)*m1*(N-m1)/N/N/(N-1),(double) n1*(N-n1)*m2*(N-m2)/N/N/(N-1), (double) n1*(N-n1)*m3*(N-m3)/N/N/(N-1), (double) n1*(N-n1)*m4*(N-m4)/N/N/(N-1),
                                 (double) m1*m2*n1*(n1-N)/N/N/(N-1),(double) m1*m3*n1*(n1-N)/N/N/(N-1),(double) m1*m4*n1*(n1-N)/N/N/(N-1), (double) m2*m3*n1*(n1-N)/N/N/(N-1), (double) m2*m4*n1*(n1-N)/N/N/(N-1), (double) m3*m4*n1*(n1-N)/N/N/(N-1)));

  double q11 = u1/Us;
  double q10 = (m1-u1)/(N-Us);
  double q21 = u2/Us;
  double q20 = (m2-u2)/(N-Us);
  double q31 = u3/Us;
  double q30 = (m3-u3)/(N-Us);
  double q41 = u4/Us;
  double q40 = (m4-u4)/(N-Us);
  double z1r1_first = q11*z1u_first+q10*(n1-z1u_first);
  double z1r2_first = q21*z1u_first+q20*(n1-z1u_first);
  double z1r3_first = q31*z1u_first+q30*(n1-z1u_first);
  double z1r4_first = q41*z1u_first+q40*(n1-z1u_first);
  double w11, w10, w21, w20, w31, w30, w41, w40;
  if(Us==1)
  { w11 = 0;
    w21 = 0;
    w31 = 0;
    w41 = 0; }
  else
  {w11  = (double) u1/(Us-1)*(1-2*q11) + (double) Us/(Us-1)*q11*q11;
    w21  = (double) u2/(Us-1)*(1-2*q21) + (double) Us/(Us-1)*q21*q21;
    w31  = (double) u3/(Us-1)*(1-2*q31) + (double) Us/(Us-1)*q31*q31;
    w41  = (double) u4/(Us-1)*(1-2*q41) + (double) Us/(Us-1)*q41*q41;

  }
  if(Us==N-1)
  {w10 = 0;
    w20 = 0;
    w30 = 0;
    w40 = 0;
  }
  else
  {
    w10   = (double)(m1 - u1)/(N-Us-1)*(1-2*q10) + (double) (N-Us)/(N-Us-1)*pow(q10,2);
    w20   = (double)(m2 - u2)/(N-Us-1)*(1-2*q20) + (double) (N-Us)/(N-Us-1)*pow(q20,2);
    w30   = (double)(m3 - u3)/(N-Us-1)*(1-2*q30) + (double) (N-Us)/(N-Us-1)*pow(q30,2);
    w40   = (double)(m4 - u4)/(N-Us-1)*(1-2*q40) + (double) (N-Us)/(N-Us-1)*pow(q40,2);


  }
  double z1r1_var = (w11-w10)*z1u_first-(z1u_second)*(w11/Us+w10/(N-Us))
    +n1*(N-Us-n1+2*z1u_first)*w10/(N-Us)+(z1u_second-z1u_first*z1u_first)*(q11-q10)*(q11-q10);
  double z1r2_var = (w21-w20)*z1u_first-(z1u_second)*(w21/Us+w20/(N-Us))
    +n1*(N-Us-n1+2*z1u_first)*w20/(N-Us)+(z1u_second-z1u_first*z1u_first)*(q21-q20)*(q21-q20);
  double z1r3_var = (w31-w30)*z1u_first-(z1u_second)*(w31/Us+w30/(N-Us))
    +n1*(N-Us-n1+2*z1u_first)*w30/(N-Us)+(z1u_second-z1u_first*z1u_first)*(q31-q30)*(q31-q30);

  double z1r4_var = (w41-w40)*z1u_first-(z1u_second)*(w41/Us+w40/(N-Us))
    +n1*(N-Us-n1+2*z1u_first)*w40/(N-Us)+(z1u_second-z1u_first*z1u_first)*(q41-q40)*(q41-q40);

  double z1r1_z1r2_cov;
  double z1r1_z1r3_cov;
  double z1r1_z1r4_cov;
  double z1r2_z1r3_cov;
  double z1r2_z1r4_cov;
  double z1r3_z1r4_cov;

  if(Us>1 && Us <N-1){
    z1r1_z1r2_cov = q11*q21/(Us-1)*(z1u_second-Us*z1u_first)+q10*q20/(N-Us-1)*n1*(n1-N+Us)
    +q10*q20/(N-Us-1)*(z1u_second-(2*n1-N+Us)*z1u_first)
    +(q11*q21+q10*q20-q10*q21-q20*q11)*(z1u_second-z1u_first*z1u_first);

    z1r1_z1r3_cov = q11*q31/(Us-1)*(z1u_second-Us*z1u_first)+q10*q30/(N-Us-1)*n1*(n1-N+Us)
      +q10*q30/(N-Us-1)*(z1u_second-(2*n1-N+Us)*z1u_first)
      +(q11*q31+q10*q30-q10*q31-q30*q11)*(z1u_second-z1u_first*z1u_first);

      z1r1_z1r4_cov = q11*q41/(Us-1)*(z1u_second-Us*z1u_first)+q10*q40/(N-Us-1)*n1*(n1-N+Us)
        +q10*q40/(N-Us-1)*(z1u_second-(2*n1-N+Us)*z1u_first)
        +(q11*q41+q10*q40-q10*q41-q40*q11)*(z1u_second-z1u_first*z1u_first);

        z1r2_z1r3_cov = q21*q31/(Us-1)*(z1u_second-Us*z1u_first)+q20*q30/(N-Us-1)*n1*(n1-N+Us)
          +q20*q30/(N-Us-1)*(z1u_second-(2*n1-N+Us)*z1u_first)
          +(q21*q31+q20*q30-q20*q31-q30*q21)*(z1u_second-z1u_first*z1u_first);

          z1r2_z1r4_cov = q21*q41/(Us-1)*(z1u_second-Us*z1u_first)+q20*q40/(N-Us-1)*n1*(n1-N+Us)
            +q20*q40/(N-Us-1)*(z1u_second-(2*n1-N+Us)*z1u_first)
            +(q21*q41+q20*q40-q20*q41-q40*q21)*(z1u_second-z1u_first*z1u_first);

            z1r3_z1r4_cov = q31*q41/(Us-1)*(z1u_second-Us*z1u_first)+q30*q40/(N-Us-1)*n1*(n1-N+Us)
              +q30*q40/(N-Us-1)*(z1u_second-(2*n1-N+Us)*z1u_first)
              +(q31*q41+q30*q40-q30*q41-q40*q31)*(z1u_second-z1u_first*z1u_first);
  }
  if(Us==1){
    z1r1_z1r2_cov = q10*q20/(N-Us-1)*n1*(n1-N+Us)
    +q10*q20/(N-Us-1)*(z1u_second-(2*n1-N+Us)*z1u_first)
    +(q11*q21+q10*q20-q10*q21-q20*q11)*(z1u_second-z1u_first*z1u_first);

    z1r1_z1r3_cov = q10*q30/(N-Us-1)*n1*(n1-N+Us)
      +q10*q30/(N-Us-1)*(z1u_second-(2*n1-N+Us)*z1u_first)
      +(q11*q31+q10*q30-q10*q31-q30*q11)*(z1u_second-z1u_first*z1u_first);

      z1r1_z1r4_cov = q10*q40/(N-Us-1)*n1*(n1-N+Us)
        +q10*q40/(N-Us-1)*(z1u_second-(2*n1-N+Us)*z1u_first)
        +(q11*q41+q10*q40-q10*q41-q40*q11)*(z1u_second-z1u_first*z1u_first);

        z1r2_z1r3_cov = q20*q30/(N-Us-1)*n1*(n1-N+Us)
          +q20*q30/(N-Us-1)*(z1u_second-(2*n1-N+Us)*z1u_first)
          +(q21*q31+q20*q30-q20*q31-q30*q21)*(z1u_second-z1u_first*z1u_first);

          z1r2_z1r4_cov = q20*q40/(N-Us-1)*n1*(n1-N+Us)
            +q20*q40/(N-Us-1)*(z1u_second-(2*n1-N+Us)*z1u_first)
            +(q21*q41+q20*q40-q20*q41-q40*q21)*(z1u_second-z1u_first*z1u_first);

            z1r3_z1r4_cov = q30*q40/(N-Us-1)*n1*(n1-N+Us)
              +q30*q40/(N-Us-1)*(z1u_second-(2*n1-N+Us)*z1u_first)
              +(q31*q41+q30*q40-q30*q41-q40*q31)*(z1u_second-z1u_first*z1u_first);

  }
  if(Us==N-1){
    z1r1_z1r2_cov = q11*q21/(Us-1)*(z1u_second-Us*z1u_first)
    +(q11*q21+q10*q20-q10*q21-q20*q11)*(z1u_second-z1u_first*z1u_first);

    z1r1_z1r3_cov = q11*q31/(Us-1)*(z1u_second-Us*z1u_first)
      +(q11*q31+q10*q30-q10*q31-q30*q11)*(z1u_second-z1u_first*z1u_first);

    z1r1_z1r4_cov = q11*q41/(Us-1)*(z1u_second-Us*z1u_first)
      +(q11*q41+q10*q40-q10*q41-q40*q11)*(z1u_second-z1u_first*z1u_first);

    z1r2_z1r3_cov = q21*q31/(Us-1)*(z1u_second-Us*z1u_first)
      +(q21*q31+q20*q30-q20*q31-q30*q21)*(z1u_second-z1u_first*z1u_first);

    z1r2_z1r4_cov = q21*q41/(Us-1)*(z1u_second-Us*z1u_first)
      +(q21*q41+q20*q40-q20*q41-q40*q21)*(z1u_second-z1u_first*z1u_first);

    z1r3_z1r4_cov = q31*q41/(Us-1)*(z1u_second-Us*z1u_first)
      +(q31*q41+q30*q40-q30*q41-q40*q31)*(z1u_second-z1u_first*z1u_first);
  }


  return(NumericVector::create(z1r1_first,z1r2_first, z1r3_first, z1r4_first, z1r1_var, z1r2_var, z1r3_var, z1r4_var, z1r1_z1r2_cov, z1r1_z1r3_cov, z1r1_z1r4_cov, z1r2_z1r3_cov, z1r2_z1r4_cov, z1r3_z1r4_cov));

}



/*
 * This function computes the moments of z1^TR_1 and z2^Tr1 for three by two table
 */


// [[Rcpp::export]]
NumericVector zr_moments_three_by_two(int n1, int n2, int n3, int m1, int m2, int N, NumericVector gamma_delta, int u1, int u2, double shared_divisor) {
  NumericVector zu_moments = zu_moments_three_treatment_log(n1,n2,n3,N,gamma_delta,u1+u2,shared_divisor);
  double z1u_first= zu_moments[0];
  double z2u_first= zu_moments[1];
  double z1u_second=zu_moments[2];
  double z2u_second=zu_moments[3];
  double z1u_z2u_product = zu_moments[4];
  double Us = u1 + u2;
  //if Us = 0 or N, becomes the usual RCT case
  if(Us==0 || Us==N)
    return(NumericVector::create((double) m1*n1/N,(double) m1*n2/N,
                                 (double) n1*(N-n1)*m1*(N-m1)/N/N/(N-1),(double) n2*(N-n2)*m1*(N-m1)/N/N/(N-1),
                                 (double) m1*(m1-N)*n1*n2/N/N/(N-1)));

  double q11 = u1/Us;
  double q10 = (m1-u1)/(N-Us);
  double z1r1_first = q11*z1u_first+q10*(n1-z1u_first);
  double z2r1_first = q11*z2u_first+q10*(n2-z2u_first);
  double w11, w10;
  if(Us==1)
    w11 = 0;
  else
    w11  = (double) u1/(Us-1)*(1-2*q11) + (double) Us/(Us-1)*q11*q11;

  if(Us==N-1)
    w10 = 0;
  else
    w10   = (double)(m1 - u1)/(N-Us-1)*(1-2*q10) + (double) (N-Us)/(N-Us-1)*pow(q10,2);
  double z1r1_var = (w11-w10)*z1u_first-(z1u_second)*(w11/Us+w10/(N-Us))
    +n1*(N-Us-n1+2*z1u_first)*w10/(N-Us)+(z1u_second-z1u_first*z1u_first)*(q11-q10)*(q11-q10);
  double z2r1_var = (w11-w10)*z2u_first-(z2u_second)*(w11/Us+w10/(N-Us))
    +n2*(N-Us-n2+2*z2u_first)*w10/(N-Us)+(z2u_second-z2u_first*z2u_first)*(q11-q10)*(q11-q10);

  double z1r1_z2r1_cov;
  if(Us>1 && Us<N-1){
    z1r1_z2r1_cov = (q11*q11+q10*q10-2*q11*q10)*(z1u_z2u_product-z1u_first*z2u_first)
    +u1*(u1-Us)/Us/Us/(Us-1)*z1u_z2u_product
    +(m1-u1)*(m1-u1-N+Us)/(N-Us)/(N-Us)/(N-Us-1)*(z1u_z2u_product+n1*n2)
    -(m1-u1)*(m1-u1-N+Us)/(N-Us)/(N-Us)/(N-Us-1)*(n1*z2u_first+n2*z1u_first);

  }
  if(Us==1){
    z1r1_z2r1_cov = (q11*q11+q10*q10-2*q11*q10)*(z1u_z2u_product-z1u_first*z2u_first)
    +(m1-u1)*(m1-u1-N+Us)/(N-Us)/(N-Us)/(N-Us-1)*(z1u_z2u_product+n1*n2)
    -(m1-u1)*(m1-u1-N+Us)/(N-Us)/(N-Us)/(N-Us-1)*(n1*z2u_first+n2*z1u_first);


  }
  if(Us==N-1){
    z1r1_z2r1_cov = (q11*q11+q10*q10-2*q11*q10)*(z1u_z2u_product-z1u_first*z2u_first)
    +u1*(u1-Us)/Us/Us/(Us-1)*z1u_z2u_product;

  }


  return(NumericVector::create(z1r1_first,z2r1_first,z1r1_var,z2r1_var,z1r1_z2r1_cov));

}





//This function computes the moments of (E(z1r1),E(z1r2),E(z2r1),E(z2r2), var(z1r1),var(z1r2), var(z2r1), var(z2r2),
//cov(z1r1,z1r2), cov(z1r1,z2r1), cov(z1r1,z2r2),cov(z1r2,z2r1),cov(z1r2,z2r2),cov(z2r1,z2r2))
// [[Rcpp::export]]
NumericVector zr_moments_three_by_three(int n1, int n2, int n3, int m1, int m2, int m3, int N, NumericVector gamma_delta, int u1, int u2, int u3, double shared_divisor) {
  NumericVector zu_moments = zu_moments_three_treatment_log(n1,n2,n3, N,gamma_delta,u1+u2+u3,shared_divisor);
  double z1u_first= zu_moments[0];
  double z2u_first = zu_moments[1];
  double z1u_second = zu_moments[2];
  double z2u_second = zu_moments[3];
  double z1u_z2u_product = zu_moments[4];
  double Us = u1 + u2+u3;
  //if Us = 0 or N, becomes the usual RCT case, have not done
  if(Us==0 || Us==N)
    return(NumericVector::create((double) m1*n1/N,(double)m2*n1/N,(double) m1*n2/N, (double) m2*n2/N,
                                 (double) n1*(N-n1)*m1*(N-m1)/N/N/(N-1),(double) n1*(N-n1)*m2*(N-m2)/N/N/(N-1),(double) n2*(N-n2)*m1*(N-m1)/N/N/(N-1),(double) n2*(N-n2)*m2*(N-m2)/N/N/(N-1),
                                 (double) m1*m2*n1*(n1-N)/N/N/(N-1), (double) n1*n2*m1*(m1-N)/N/N/(N-1),(double) m1*m2*n1*n2/N/N/(N-1),(double) m1*m2*n1*n2/N/N/(N-1), (double) n1*n2*m2*(m2-N)/N/N/(N-1),(double) m1*m2*n2*(n2-N)/N/N/(N-1)));

  double q11 = u1/Us;
  double q10 = (m1-u1)/(N-Us);
  double q21 = u2/Us;
  double q20 = (m2-u2)/(N-Us);
  double z1r1_first = q11*z1u_first+q10*(n1-z1u_first);
  double z1r2_first = q21*z1u_first+q20*(n1-z1u_first);
  double z2r1_first = q11*z2u_first+q10*(n2-z2u_first);
  double z2r2_first = q21*z2u_first+q20*(n2-z2u_first);
  double w11, w10, w21, w20;
  if(Us==1)
  { w11 = 0;
    w21 = 0;
  }
  else
  {w11  = (double) u1/(Us-1)*(1-2*q11) + (double) Us/(Us-1)*q11*q11;
    w21  = (double) u2/(Us-1)*(1-2*q21) + (double) Us/(Us-1)*q21*q21;}
  if(Us==N-1)
  {w10 = 0;
    w20 = 0;
  }
  else
  {
    w10   = (double)(m1 - u1)/(N-Us-1)*(1-2*q10) + (double) (N-Us)/(N-Us-1)*pow(q10,2);
    w20   = (double)(m2 - u2)/(N-Us-1)*(1-2*q20) + (double) (N-Us)/(N-Us-1)*pow(q20,2);

  }
  double z1r1_var = (w11-w10)*z1u_first-(z1u_second)*(w11/Us+w10/(N-Us))
    +n1*(N-Us-n1+2*z1u_first)*w10/(N-Us)+(z1u_second-z1u_first*z1u_first)*(q11-q10)*(q11-q10);
  double z1r2_var = (w21-w20)*z1u_first-(z1u_second)*(w21/Us+w20/(N-Us))
    +n1*(N-Us-n1+2*z1u_first)*w20/(N-Us)+(z1u_second-z1u_first*z1u_first)*(q21-q20)*(q21-q20);
  double z2r1_var = (w11-w10)*z2u_first-(z2u_second)*(w11/Us+w10/(N-Us))
    +n2*(N-Us-n2+2*z2u_first)*w10/(N-Us)+(z2u_second-z2u_first*z2u_first)*(q11-q10)*(q11-q10);
  double z2r2_var = (w21-w20)*z2u_first-(z2u_second)*(w21/Us+w20/(N-Us))
    +n2*(N-Us-n2+2*z2u_first)*w20/(N-Us)+(z2u_second-z2u_first*z2u_first)*(q21-q20)*(q21-q20);
  double z1r1_z1r2_cov, z1r1_z2r1_cov, z1r1_z2r2_cov, z1r2_z2r1_cov, z1r2_z2r2_cov, z2r1_z2r2_cov;
  if(Us>1 && Us <N-1){
    z1r1_z1r2_cov = q11*q21/(Us-1)*(z1u_second-Us*z1u_first)+q10*q20/(N-Us-1)*n1*(n1-N+Us)
    +q10*q20/(N-Us-1)*(z1u_second-(2*n1-N+Us)*z1u_first)
    +(q11*q21+q10*q20-q10*q21-q20*q11)*(z1u_second-z1u_first*z1u_first);

    z2r1_z2r2_cov = q11*q21/(Us-1)*(z2u_second-Us*z2u_first)+q10*q20/(N-Us-1)*n2*(n2-N+Us)
      +q10*q20/(N-Us-1)*(z2u_second-(2*n2-N+Us)*z2u_first)
      +(q11*q21+q10*q20-q10*q21-q20*q11)*(z2u_second-z2u_first*z2u_first);


      z1r1_z2r1_cov = (q11*q11+q10*q10-2*q11*q10)*(z1u_z2u_product-z1u_first*z2u_first)
        +u1*(u1-Us)/Us/Us/(Us-1)*z1u_z2u_product
      +(m1-u1)*(m1-u1-N+Us)/(N-Us)/(N-Us)/(N-Us-1)*(z1u_z2u_product+n1*n2)
      -(m1-u1)*(m1-u1-N+Us)/(N-Us)/(N-Us)/(N-Us-1)*(n1*z2u_first+n2*z1u_first);

      z1r2_z2r2_cov = (q21*q21+q20*q20-2*q21*q20)*(z1u_z2u_product-z1u_first*z2u_first)
        +u2*(u2-Us)/Us/Us/(Us-1)*z1u_z2u_product
      +(m2-u2)*(m2-u2-N+Us)/(N-Us)/(N-Us)/(N-Us-1)*(z1u_z2u_product+n1*n2)
      -(m2-u2)*(m2-u2-N+Us)/(N-Us)/(N-Us)/(N-Us-1)*(n1*z2u_first+n2*z1u_first);

      z1r1_z2r2_cov = z1r2_z2r1_cov = (q11*q21+q10*q20-q21*q10-q11*q20)*(z1u_z2u_product-z1u_first*z2u_first)
        +q11*q21/(Us-1)*z1u_z2u_product
      +(q10*q20)/(N-Us-1)*(z1u_z2u_product+n1*n2)
      -q20*q10/(N-Us-1)*(n1*z2u_first+n2*z1u_first);



  }
  if(Us==1){
    z1r1_z1r2_cov = q10*q20/(N-Us-1)*n1*(n1-N+Us)
    +q10*q20/(N-Us-1)*(z1u_second-(2*n1-N+Us)*z1u_first)
    +(q11*q21+q10*q20-q10*q21-q20*q11)*(z1u_second-z1u_first*z1u_first);

    z2r1_z2r2_cov = q10*q20/(N-Us-1)*n2*(n2-N+Us)
      +q10*q20/(N-Us-1)*(z2u_second-(2*n2-N+Us)*z2u_first)
      +(q11*q21+q10*q20-q10*q21-q20*q11)*(z2u_second-z2u_first*z2u_first);


      z1r1_z2r1_cov = (q11*q11+q10*q10-2*q11*q10)*(z1u_z2u_product-z1u_first*z2u_first)
        +(m1-u1)*(m1-u1-N+Us)/(N-Us)/(N-Us)/(N-Us-1)*(z1u_z2u_product+n1*n2)
        -(m1-u1)*(m1-u1-N+Us)/(N-Us)/(N-Us)/(N-Us-1)*(n1*z2u_first+n2*z1u_first);

        z1r2_z2r2_cov = (q21*q21+q20*q20-2*q21*q20)*(z1u_z2u_product-z1u_first*z2u_first)
          +(m2-u2)*(m2-u2-N+Us)/(N-Us)/(N-Us)/(N-Us-1)*(z1u_z2u_product+n1*n2)
          -(m2-u2)*(m2-u2-N+Us)/(N-Us)/(N-Us)/(N-Us-1)*(n1*z2u_first+n2*z1u_first);

          z1r1_z2r2_cov = z1r2_z2r1_cov = (q11*q21+q10*q20-q21*q10-q11*q20)*(z1u_z2u_product-z1u_first*z2u_first)
            +(q10*q20)/(N-Us-1)*(z1u_z2u_product+n1*n2)
            -q20*q10/(N-Us-1)*(n1*z2u_first+n2*z1u_first);



  }
  if(Us==N-1){
    z1r1_z1r2_cov = q11*q21/(Us-1)*(z1u_second-Us*z1u_first)
    +(q11*q21+q10*q20-q10*q21-q20*q11)*(z1u_second-z1u_first*z1u_first);

    z2r1_z2r2_cov = q11*q21/(Us-1)*(z2u_second-Us*z2u_first)
      +(q11*q21+q10*q20-q10*q21-q20*q11)*(z2u_second-z2u_first*z2u_first);


    z1r1_z2r1_cov = (q11*q11+q10*q10-2*q11*q10)*(z1u_z2u_product-z1u_first*z2u_first)
      +u1*(u1-Us)/Us/Us/(Us-1)*z1u_z2u_product;

    z1r2_z2r2_cov = (q21*q21+q20*q20-2*q21*q20)*(z1u_z2u_product-z1u_first*z2u_first)
      +u2*(u2-Us)/Us/Us/(Us-1)*z1u_z2u_product;


    z1r1_z2r2_cov = z1r2_z2r1_cov = (q11*q21+q10*q20-q21*q10-q11*q20)*(z1u_z2u_product-z1u_first*z2u_first)
      +q11*q21/(Us-1)*z1u_z2u_product;


  }


  return(NumericVector::create(z1r1_first,z1r2_first, z2r1_first, z2r2_first,
                               z1r1_var, z1r2_var, z2r1_var,z2r2_var,
                               z1r1_z1r2_cov, z1r1_z2r1_cov,
                               z1r1_z2r2_cov, z1r2_z2r1_cov,
                               z1r2_z2r2_cov, z2r1_z2r2_cov));

}

/*
 * This function computes the mean and covariance matrix components for four (treatment) and two (outcome) table
 */
// [[Rcpp::export]]
NumericVector zr_moments_four_by_two(int n1, int n2, int n3, int n4, int m1, int m2, int N, NumericVector gamma_delta, int u1, int u2, double shared_divisor) {
  NumericVector zu_moments = zu_moments_four_treatment_log(n1,n2,n3,n4,N,gamma_delta,u1+u2,shared_divisor);
  double z1u_first= zu_moments[0];
  double z2u_first= zu_moments[1];
  double z3u_first= zu_moments[2];
  double z1u_second=zu_moments[3];
  double z2u_second=zu_moments[4];
  double z3u_second=zu_moments[5];
  double z1u_z2u_product=zu_moments[6];
  double z1u_z3u_product=zu_moments[7];
  double z2u_z3u_product=zu_moments[8];
  double Us = u1 + u2;
  //if Us = 0 or N, becomes the usual RCT case
  if(Us==0 || Us==N)
    return(NumericVector::create((double) m1*n1/N,(double) m1*n2/N,(double) m1*n3/N,
                                 (double) n1*(N-n1)*m1*(N-m1)/N/N/(N-1),(double) n2*(N-n2)*m1*(N-m1)/N/N/(N-1), (double) n3*(N-n3)*m1*(N-m1)/N/N/(N-1),
                                 (double) m1*(m1-N)*n1*n2/N/N/(N-1),(double) m1*(m1-N)*n1*n3/N/N/(N-1),(double) m1*(m1-N)*n2*n3/N/N/(N-1)));

  double q11 = u1/Us;
  double q10 = (m1-u1)/(N-Us);
  double z1r1_first = q11*z1u_first+q10*(n1-z1u_first);
  double z2r1_first = q11*z2u_first+q10*(n2-z2u_first);
  double z3r1_first = q11*z3u_first+q10*(n3-z3u_first);
  double w11, w10;
  if(Us==1)
    w11 = 0;
  else
    w11  = (double) u1/(Us-1)*(1-2*q11) + (double) Us/(Us-1)*q11*q11;

  if(Us==N-1)
    w10 = 0;
  else
    w10   = (double)(m1 - u1)/(N-Us-1)*(1-2*q10) + (double) (N-Us)/(N-Us-1)*pow(q10,2);
  double z1r1_var = (w11-w10)*z1u_first-(z1u_second)*(w11/Us+w10/(N-Us))
    +n1*(N-Us-n1+2*z1u_first)*w10/(N-Us)+(z1u_second-z1u_first*z1u_first)*(q11-q10)*(q11-q10);
  double z2r1_var = (w11-w10)*z2u_first-(z2u_second)*(w11/Us+w10/(N-Us))
    +n2*(N-Us-n2+2*z2u_first)*w10/(N-Us)+(z2u_second-z2u_first*z2u_first)*(q11-q10)*(q11-q10);
  double z3r1_var = (w11-w10)*z3u_first-(z3u_second)*(w11/Us+w10/(N-Us))
    +n3*(N-Us-n3+2*z3u_first)*w10/(N-Us)+(z3u_second-z3u_first*z3u_first)*(q11-q10)*(q11-q10);

  double z1r1_z2r1_cov, z1r1_z3r1_cov, z2r1_z3r1_cov;
  if(Us>1 && Us<N-1){
    z1r1_z2r1_cov = (q11*q11+q10*q10-2*q11*q10)*(z1u_z2u_product-z1u_first*z2u_first)
    +u1*(u1-Us)/Us/Us/(Us-1)*z1u_z2u_product
    +(m1-u1)*(m1-u1-N+Us)/(N-Us)/(N-Us)/(N-Us-1)*(z1u_z2u_product+n1*n2)
    -(m1-u1)*(m1-u1-N+Us)/(N-Us)/(N-Us)/(N-Us-1)*(n1*z2u_first+n2*z1u_first);

    z1r1_z3r1_cov = (q11*q11+q10*q10-2*q11*q10)*(z1u_z3u_product-z1u_first*z3u_first)
      +u1*(u1-Us)/Us/Us/(Us-1)*z1u_z3u_product
    +(m1-u1)*(m1-u1-N+Us)/(N-Us)/(N-Us)/(N-Us-1)*(z1u_z3u_product+n1*n3)
    -(m1-u1)*(m1-u1-N+Us)/(N-Us)/(N-Us)/(N-Us-1)*(n1*z3u_first+n3*z1u_first);

    z2r1_z3r1_cov = (q11*q11+q10*q10-2*q11*q10)*(z2u_z3u_product-z2u_first*z3u_first)
      +u1*(u1-Us)/Us/Us/(Us-1)*z2u_z3u_product
    +(m1-u1)*(m1-u1-N+Us)/(N-Us)/(N-Us)/(N-Us-1)*(z2u_z3u_product+n2*n3)
    -(m1-u1)*(m1-u1-N+Us)/(N-Us)/(N-Us)/(N-Us-1)*(n2*z3u_first+n3*z2u_first);

  }
  if(Us==1){
    z1r1_z2r1_cov = (q11*q11+q10*q10-2*q11*q10)*(z1u_z2u_product-z1u_first*z2u_first)
    +(m1-u1)*(m1-u1-N+Us)/(N-Us)/(N-Us)/(N-Us-1)*(z1u_z2u_product+n1*n2)
    -(m1-u1)*(m1-u1-N+Us)/(N-Us)/(N-Us)/(N-Us-1)*(n1*z2u_first+n2*z1u_first);

    z1r1_z3r1_cov = (q11*q11+q10*q10-2*q11*q10)*(z1u_z3u_product-z1u_first*z3u_first)
      +(m1-u1)*(m1-u1-N+Us)/(N-Us)/(N-Us)/(N-Us-1)*(z1u_z3u_product+n1*n3)
      -(m1-u1)*(m1-u1-N+Us)/(N-Us)/(N-Us)/(N-Us-1)*(n1*z3u_first+n3*z1u_first);

      z2r1_z3r1_cov = (q11*q11+q10*q10-2*q11*q10)*(z2u_z3u_product-z2u_first*z3u_first)
        +(m1-u1)*(m1-u1-N+Us)/(N-Us)/(N-Us)/(N-Us-1)*(z2u_z3u_product+n2*n3)
        -(m1-u1)*(m1-u1-N+Us)/(N-Us)/(N-Us)/(N-Us-1)*(n2*z3u_first+n3*z2u_first);




  }
  if(Us==N-1){
    z1r1_z2r1_cov = (q11*q11+q10*q10-2*q11*q10)*(z1u_z2u_product-z1u_first*z2u_first)
    +u1*(u1-Us)/Us/Us/(Us-1)*z1u_z2u_product;

    z1r1_z3r1_cov = (q11*q11+q10*q10-2*q11*q10)*(z1u_z3u_product-z1u_first*z3u_first)
      +u1*(u1-Us)/Us/Us/(Us-1)*z1u_z3u_product;

    z2r1_z3r1_cov = (q11*q11+q10*q10-2*q11*q10)*(z2u_z3u_product-z2u_first*z3u_first)
      +u1*(u1-Us)/Us/Us/(Us-1)*z2u_z3u_product;

  }


  return(NumericVector::create(z1r1_first,z2r1_first,z3r1_first,z1r1_var,z2r1_var,z3r1_var,z1r1_z2r1_cov,z1r1_z3r1_cov, z2r1_z3r1_cov));

}








/*
 * The following functions are kept for designer's numerically checking, will not be implemented
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 */





// [[Rcpp::export]]
NumericVector zu_moments_two_treatment(int n1, int n2, int N, NumericVector gamma_delta, int Us, double shared_divisor) {
  double gamma = gamma_delta[0] - gamma_delta[1];
  double z1u_first = 0.0;
  double z1u_second = 0.0;
  double convolution_sum = 0.0;

  // Sanity check
  if (Us == 0 || gamma == 0) {
    z1u_first = 0.0;
    z1u_second = 0.0;
    return NumericVector::create(z1u_first, z1u_second);
  }

  if (Us == N) {
    z1u_first = n1;
    z1u_second = std::pow(n1, 2);
    return NumericVector::create(z1u_first, z1u_second);
  }

  // General case
  int lb_q1 = std::max(0, n1 + Us - N);
  int ub_q1 = std::min(Us, n1);
  double z1u_first_cumulant = 0.0;
  double z1u_second_cumulant = 0.0;

  for (int q1 = lb_q1; q1 <= ub_q1; ++q1) {
    // Precompute the kernel
    double weighted_kernel = choose_C(Us, q1) * choose_C(N - Us, n1 - q1) / shared_divisor * exp(gamma * q1);

    convolution_sum += weighted_kernel;
    z1u_first_cumulant += q1 * weighted_kernel;
    z1u_second_cumulant += std::pow(q1, 2) * weighted_kernel;
  }

  // Check for division by zero
  if (convolution_sum == 0.0) {
    stop("Convolution sum is zero; division by zero detected.");
  }

  // Normalize cumulants by the convolution sum
  z1u_first = z1u_first_cumulant / convolution_sum;
  z1u_second = z1u_second_cumulant / convolution_sum;

  return NumericVector::create(z1u_first, z1u_second);
}


/*
 * This function calculates the first two moments of z1u, z2u, and their product
 * which will be used in the covariance computation
 * the returned vector is (E(z1u),E(z2u),E((z1u)^2),E((z2u)^2), E(z1u*z2u))
 * This one is in theory the same as zu_moments_three_treatment_log,
 * yet, the other one uses log sum trick to give numerical stability but takes
 * longer to compute
 */
// [[Rcpp::export]]
NumericVector zu_moments_three_treatment(int n1, int n2, int n3, int N, NumericVector gamma_delta, int Us, double shared_divisor) {
  double gamma1 = gamma_delta[0] - gamma_delta[2];
  double gamma2 = gamma_delta[1] - gamma_delta[2];
  double gamma = (gamma1 != 0) ? gamma1 : (gamma2 != 0 ? gamma2 : 0);

  double convolution_sum = 0.0;

  // Sanity check
  if (Us == 0 || gamma == 0) {
    return NumericVector::create(0.0, 0.0, 0.0, 0.0, 0.0);
  }

  if (Us == N) {
    return NumericVector::create(n1, n2, std::pow(n1, 2), std::pow(n2, 2), n1 * n2);
  }

  int lb_q = std::max(0, n1 + n2 + Us - N);
  int ub_q = std::min(n1 + n2, Us);
  int lb_q1 = std::max(0, n1 + Us - N);
  int ub_q1 = std::min(Us, n1);
  int lb_q2 = std::max(0, n2 + Us - N);
  int ub_q2 = std::min(n2, Us);

  double z1u_first_cumulant = 0.0, z1u_second_cumulant = 0.0, z2u_first_cumulant = 0.0;
  double z2u_second_cumulant = 0.0, z1u_z2u_product_cumulant = 0.0;

  for (int q = lb_q; q <= ub_q; ++q) {
    for (int q1 = lb_q1; q1 <= ub_q1; ++q1) {
      int q2 = q - q1;
      if (q2 > ub_q2 || q2 < lb_q2) {
        continue;
      }

      double exp_term = (gamma1 == 0) ? exp(gamma * (q - q1)) :
        (gamma2 == 0) ? exp(gamma * (q - q2)) : exp(gamma * q);

      double kernel = choose_C(Us, q1) * choose_C(Us - q1, q2) *
        choose_C(N - Us, n1 - q1) * choose_C((N - Us) - (n1 - q1), n2 - q2) / shared_divisor;

      double weighted_kernel = exp_term * kernel;

      convolution_sum += weighted_kernel;
      z1u_first_cumulant += q1 * weighted_kernel;
      z1u_second_cumulant += std::pow(q1, 2) * weighted_kernel;
      z2u_first_cumulant += q2 * weighted_kernel;
      z2u_second_cumulant += std::pow(q2, 2) * weighted_kernel;
      z1u_z2u_product_cumulant += q1 * q2 * weighted_kernel;
    }
  }

  if (convolution_sum == 0.0) stop("Convolution sum is zero; division by zero detected.");

  return NumericVector::create(
    z1u_first_cumulant / convolution_sum,
    z2u_first_cumulant / convolution_sum,
    z1u_second_cumulant / convolution_sum,
    z2u_second_cumulant / convolution_sum,
    z1u_z2u_product_cumulant / convolution_sum
  );
}



/*
 * This function obtains the moments for z1u z2u z3u, the returned vector
 * contains (E(z1u),E(z2u),E(z3u),E((z1u)^2),E((z2u)^2),E((z3u)^2),E(z1u*z2u),E(z1u*z3u),E(z2u*z3u))
 */

// [[Rcpp::export]]
NumericVector zu_moments_four_treatment(int n1, int n2, int n3, int n4, int N, NumericVector gamma_delta, int Us, double shared_divisor) {
  // Calculate gamma values
  double gamma1 = gamma_delta[0] - gamma_delta[3];
  double gamma2 = gamma_delta[1] - gamma_delta[3];
  double gamma3 = gamma_delta[2] - gamma_delta[3];
  double gamma = (gamma1 != 0) ? gamma1 : (gamma2 != 0 ? gamma2 : gamma3);

  double convolution_sum = 0.0;
  // Sanity check
  if (Us == 0 || gamma == 0) {
    return NumericVector::create(0.0, 0.0, 0.0, 0.0, 0.0,0.0,0.0,0.0,0.0);
  }

  if (Us == N) {
    return NumericVector::create(n1, n2,n3,std::pow(n1, 2), std::pow(n2, 2),std::pow(n3,2), n1 * n2, n1*n3, n2*n3);
  }

  //general case, define the cumulant

  double z1u_first_cumulant = 0.0, z1u_second_cumulant = 0.0, z2u_first_cumulant = 0.0, z2u_second_cumulant = 0.0,
    z3u_first_cumulant = 0.0, z3u_second_cumulant = 0.0;
  double z1u_z2u_product_cumulant = 0.0, z1u_z3u_product_cumulant = 0.0, z2u_z3u_product_cumulant = 0.0;

  // Define bounds for q, q1, q2, and q3
  int lb_q = std::max(0, n1 + n2 + n3 + Us - N);
  int ub_q = std::min(n1 + n2 + n3, Us);
  int lb_q1 = std::max(0, n1 + Us - N);
  int ub_q1 = std::min(Us, n1);

  for (int q = lb_q; q <= ub_q; ++q) {
    for (int q1 = lb_q1; q1 <= ub_q1; ++q1) {
      int lb_q2 = std::max(0, n1 + n2 - q1 + Us - N);
      int ub_q2 = std::min(n2, Us - q1);

      for (int q2 = lb_q2; q2 <= ub_q2; ++q2) {
        int lb_q3 = std::max(0, n1 + n2 + n3 - q1 - q2 + Us - N);
        int ub_q3 = std::min(Us - q1 - q2, n3);
        int q3 = q - q1 - q2;

        if (q3 > ub_q3 || q3 < lb_q3) {
          continue;
        } else {
          // Determine the exponential term based on which gamma values are zero
          double exp_term;
          if (gamma1 == 0 && gamma2 == 0 && gamma3 == 0) {
            exp_term = exp(gamma * (q - q1 - q2 - q3));
          } else if (gamma1 == 0 && gamma2 == 0) {
            exp_term = exp(gamma * (q - q1 - q2));
          } else if (gamma1 == 0 && gamma3 == 0) {
            exp_term = exp(gamma * (q - q1 - q3));
          } else if (gamma2 == 0 && gamma3 == 0) {
            exp_term = exp(gamma * (q - q2 - q3));
          } else if (gamma1 == 0) {
            exp_term = exp(gamma * (q - q1));
          } else if (gamma2 == 0) {
            exp_term = exp(gamma * (q - q2));
          } else if (gamma3 == 0) {
            exp_term = exp(gamma * (q - q3));
          } else {
            exp_term = exp(gamma * q);
          }
          double weighted_kernel = exp_term / shared_divisor *
            choose_C(Us, q1) * choose_C(Us - q1, q2) * choose_C(Us - q1 - q2, q3) *
            choose_C(N - Us, n1 - q1) * choose_C(N - Us - n1 + q1, n2 - q2) *
            choose_C(N - Us - n1 + q1 - n2 + q2, n3 - q3);
          // Accumulate the convolution sum
          convolution_sum += weighted_kernel;
          z1u_first_cumulant += q1 * weighted_kernel;
          z1u_second_cumulant += std::pow(q1, 2) * weighted_kernel;
          z2u_first_cumulant += q2 * weighted_kernel;
          z2u_second_cumulant += std::pow(q2, 2) * weighted_kernel;
          z3u_first_cumulant += q3 * weighted_kernel;
          z3u_second_cumulant += std::pow(q3, 2) * weighted_kernel;
          z1u_z2u_product_cumulant += q1 * q2 * weighted_kernel;
          z1u_z3u_product_cumulant += q1 * q3 * weighted_kernel;
          z2u_z3u_product_cumulant += q2 * q3 * weighted_kernel;



        }
      }
    }
  }

  if (convolution_sum == 0.0) stop("Convolution sum is zero; division by zero detected.");

  return NumericVector::create(
    z1u_first_cumulant / convolution_sum,
    z2u_first_cumulant / convolution_sum,
    z3u_first_cumulant / convolution_sum,
    z1u_second_cumulant / convolution_sum,
    z2u_second_cumulant / convolution_sum,
    z3u_second_cumulant / convolution_sum,
    z1u_z2u_product_cumulant / convolution_sum,
    z1u_z3u_product_cumulant / convolution_sum,
    z2u_z3u_product_cumulant / convolution_sum

  );
}
















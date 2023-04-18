// https://figshare.com/articles/C_code_needed_for_running_R_code/4836503
// Simple multinomial fitting.
#include <TMB.hpp>


template<class Type>
Type objective_function<Type>::operator() ()
{

  Type nll;           // Define variable that holds the return value (neg. log. lik)
  nll = Type(0.0);    // Assign value 1.2; a cast is needed.

  DATA_MATRIX(Y); // observations (binary matrix)
  int d = Y.cols(); // number of features
  int n = Y.rows(); // number of samples
  DATA_MATRIX(x); // matrix of covariates for fixed effects
  PARAMETER_MATRIX(beta); // coefficients for the fixed effects
  matrix<Type> theta_prime(n,d); // The logratios of each signature in a sample being zero

  theta_prime = x * beta;

  for(int l=0; l<n; l++){ // binomial draws, for each observation individually
    for(int j=0; j<d; j++){
      Type num_zeros_lth_row = Y(l,j); // number of zeros in observation
      Type p_l = (exp(theta_prime(l,j))) / (1 + exp(theta_prime(l,j))); // softmax transformation; pick first probability
      nll -= log(p_l)*(num_zeros_lth_row) + log(1-p_l)*(1-num_zeros_lth_row); // binomial log-likelihood likelihood
   }
  }
  return nll;

}

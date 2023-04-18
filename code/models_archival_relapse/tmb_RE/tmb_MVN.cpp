// https://figshare.com/articles/C_code_needed_for_running_R_code/4836503
// Simple multinomial fitting.
#include <TMB.hpp>


template<class Type>
Type objective_function<Type>::operator() ()
{

  Type nll;           // Define variable that holds the return value (neg. log. lik)
  nll = Type(0.0);    // Assign value 1.2; a cast is needed.

  DATA_MATRIX(Y); // MVN observations
  int d = Y.cols(); // number of features
  int n = Y.rows(); // 2*n samples
  PARAMETER_VECTOR(logs_sd_RE); // RE: variances (scaling)
  PARAMETER_VECTOR(cov_RE); // RE: covariances
  Type d_min1 = d-1;

  using namespace density;
    
  UNSTRUCTURED_CORR_t<Type> nll_mvn(cov_RE);
  for(int i=0;i<n;i++){ // likelihood for the random effects (multivariate normal)
    nll += VECSCALE_t(nll_mvn, exp(logs_sd_RE))(Y.row(i));
  }
      
  return nll;
    



}

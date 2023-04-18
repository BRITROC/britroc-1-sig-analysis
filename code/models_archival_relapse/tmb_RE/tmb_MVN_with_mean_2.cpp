// https://figshare.com/articles/C_code_needed_for_running_R_code/4836503
// Simple multinomial fitting.
#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() ()
{

  Type nll;           // Define variable that holds the return value (neg. log. lik)
  nll = Type(0.0);    // Assign value 1.2; a cast is needed.

  DATA_MATRIX(Y); // MVN observations
  DATA_MATRIX(x); // matrix of covariates for fixed effects
  int d = Y.cols(); // number of features
  int n = Y.rows(); // numnber of samples
  PARAMETER_VECTOR(logs_sd_RE); // RE: variances (scaling)
  PARAMETER_VECTOR(cov_RE); // RE: covariances
  PARAMETER_MATRIX(beta); // means

  using namespace density;

  matrix<Type> means(n,d);
  means = x * beta;
    
  UNSTRUCTURED_CORR_t<Type> nll_mvn(cov_RE);
  for(int i=0;i<n;i++){ // likelihood for the random effects (multivariate normal)
    vector<Type> subtract(d);
    subtract = vector<Type>(Y.row(i)) - vector<Type>(means.row(i));
    nll += VECSCALE(nll_mvn, exp(logs_sd_RE))(subtract);
  }
  return nll;

}

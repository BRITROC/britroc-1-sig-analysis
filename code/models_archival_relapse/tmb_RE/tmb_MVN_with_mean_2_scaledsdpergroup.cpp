// adding a dirichlet step for possible differential precision
// note we're only changing the variances, not the covariances
#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() ()
{

  Type nll;           // Define variable that holds the return value (neg. log. lik)
  nll = Type(0.0);    // Assign value 1.2; a cast is needed.

  DATA_MATRIX(Y); // MVN observations
  DATA_MATRIX(x); // matrix of covariates for fixed effects
  DATA_MATRIX(lambda_accessory_mat); // matrix to get a length 2*n vector lambda from a length 2 vector lambda
  int d = Y.cols(); // number of features
  int n = Y.rows(); // number of samples
  PARAMETER_VECTOR(logs_sd_RE); // RE: variances (scaling)
  PARAMETER_VECTOR(cov_RE); // RE: covariances
  PARAMETER_MATRIX(beta); // means
  PARAMETER(log_lambda); // scaling factor for differential precision

  using namespace density;

  matrix<Type> means(n,d);
  means = x * beta;

  matrix<Type> logs_sd_RE_i(2,d);
  logs_sd_RE_i.row(0) = logs_sd_RE;
  vector<Type> log_lambda_rep(d);
  log_lambda_rep.fill(log_lambda);
  logs_sd_RE_i.row(1) = logs_sd_RE + log_lambda_rep; // add some value to increase or decrease standard deviations of MVN

  matrix<Type> log_lambda_mat(n,d);
  log_lambda_mat = lambda_accessory_mat * logs_sd_RE_i; // get a vector of lambdas of length 2*n

  UNSTRUCTURED_CORR_t<Type> nll_mvn(cov_RE);
  for(int i=0;i<n;i++){ // likelihood for the random effects (multivariate normal)
    vector<Type> subtract(d);
    subtract = vector<Type>(Y.row(i)) - vector<Type>(means.row(i));
    // int group = x.row(i)(1);
    // nll += VECSCALE(nll_mvn, exp(logs_sd_RE_i.row(group)))(subtract); // get the appropriate standard deviations for the group
    nll += VECSCALE(nll_mvn, exp(vector<Type>(log_lambda_mat.row(i))))(subtract); 
  }
  return nll;

}

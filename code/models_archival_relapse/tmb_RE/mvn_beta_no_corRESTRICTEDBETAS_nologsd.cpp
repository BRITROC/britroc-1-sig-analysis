//
//  same as mvn_beta_no_cor.cpp, but with zero beta slope

#include <TMB.hpp>
template<class Type>
Type objective_function<Type>::operator() ()
{

  Type nll;           // Define variable that holds the return value (neg. log. lik)
  nll = Type(0.0);    // Assign value 0.0; a cast is needed.

  DATA_MATRIX(Y); // observations (multivariate gaussian)
  PARAMETER_VECTOR(beta); // coefficients for the intercept and change between conditions
  DATA_MATRIX(x); // matrix of covariates for fixed effects

  int n = Y.rows();
  int d = Y.cols();

  matrix<Type> beta_mat(2,d);
  beta_mat.fill(0); // slopes will be zero
  beta_mat.row(0) = beta; // we are estimating the intercepts

  matrix<Type> mu(n,d); // 
  mu = x * beta_mat;

  using namespace density;
  vector<Type> vec1(d);
  vec1.fill(1.0);
  // multivariate normal
  matrix<Type> cov(d,d);
  cov.fill(0);
  cov.diagonal() = vec1;
  MVNORM_t<Type> nll_mvn(cov); // no correlations
  for(int i=0;i<n;i++){
      nll += nll_mvn(Y.row(i)-mu.row(i));
  }

  return nll;

}


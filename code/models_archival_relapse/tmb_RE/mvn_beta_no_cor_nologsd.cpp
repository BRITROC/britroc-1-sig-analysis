//
//  based on mvn.cpp (and tmb_MVN.cpp, which are the same)
// here we hadd a mean value

#include <TMB.hpp>
template<class Type>
Type objective_function<Type>::operator() ()
{

  Type nll;           // Define variable that holds the return value (neg. log. lik)
  nll = Type(0.0);    // Assign value 0.0; a cast is needed.

  DATA_MATRIX(Y); // observations (multivariate gaussian)
  PARAMETER_MATRIX(beta); // coefficients for the intercept and change between conditions
  DATA_MATRIX(x); // matrix of covariates for fixed effects

  int n = Y.rows();
  int d = Y.cols();

  matrix<Type> mu(n,d); // 
  mu = x * beta;

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

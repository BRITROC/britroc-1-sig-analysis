//
//  mvn.cpp
//  
//
//  Created by Lena Morrill on 30/07/2021.
//
// Simple MVN fitting, to see if the covariance matrix that it gives is PSD.
// same as tmb_MVN.cpp
#include <TMB.hpp>
template<class Type>
Type objective_function<Type>::operator() ()
{

  Type nll;           // Define variable that holds the return value (neg. log. lik)
  nll = Type(0.0);    // Assign value 0.0; a cast is needed.

  DATA_MATRIX(Y); // observations (multivariate gaussian)
  PARAMETER_VECTOR(cov_par_RE); // RE. vector with non-redundant entries of the covariance matrix for the random effects variances and covariances
  PARAMETER_VECTOR(logs_sd_RE); // RE. vector with scaling factors for the matrix of covariances of RE
//  int d = Y.cols();
//  int d_min1 = d - 1;
  int n = Y.rows();

  // multivariate normal
  using namespace density;
  UNSTRUCTURED_CORR_t<Type> nll_mvn(cov_par_RE);
for(int i=0;i<n;i++){
      nll += VECSCALE_t(nll_mvn, exp(logs_sd_RE))(Y.row(i));
    }


    report(nll_mvn.cov());
  return nll;

}

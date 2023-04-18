// https://figshare.com/articles/C_code_needed_for_running_R_code/4836503
// Simple multinomial fitting.
#include <TMB.hpp>

// template<class Type>
// Type LN_lik( vector<Type> Y_row, vector<Type> theta_prime_l, int d)
// {
// 	using namespace density;
//     Type yald = log(Y_row[d-1]);
//     vector<Type> y_alr(d-1);
//     for(int i=0; i<(d-1); i++){
//       y_alr[i] = log(Y_row[i]) - yald;
//     }
//     matrix<Type> Sigma_diag(d-1,d-1);
//     vector<Type> ones(d-1);
//     ones.fill(1);
//     Sigma_diag.diagonal() = ones;  // Fill the diagonal
//     MVNORM_t<Type> nll_mvn_2(Sigma_diag);
//     return nll_mvn_2(y_alr-theta_prime_l);
// }


template<class Type>
Type objective_function<Type>::operator() ()
{

  Type nll;           // Define variable that holds the return value (neg. log. lik)
  nll = Type(0.0);    // Assign value 1.2; a cast is needed.

  DATA_MATRIX(Y); // observations (count matrix)
  int d = Y.cols(); // number of features
  int n = Y.rows(); // number of observations: all samples in all groups * number of categories
  DATA_INTEGER(num_individuals); // number of different indivdiduals (each with a random effect)
  DATA_INTEGER(num_sigs); // number of categories or signatures
  DATA_MATRIX(x); // matrix of covariates for fixed effects
  DATA_MATRIX(z); // matrix for random effects
  PARAMETER_MATRIX(u_large); // coefficients for the random effects. Even though it is defined as matrix (for TMB matrix multiplication), it is a vector
  DATA_MATRIX(V); // matrix of 0 and 1 for the signatures
  PARAMETER_MATRIX(shared_intercept); // shared intercept by all entries
  PARAMETER_VECTOR(logs_sd_RE); // RE
  PARAMETER_MATRIX(mu_beta); // scaling factor for the diagonal entries of the matrix for signatures
  PARAMETER_VECTOR(logs_sd_cov_beta); // non-diagonal and non-redundant entries of beta (indicating correlations between signatures)
  PARAMETER_VECTOR(logs_var_cov_beta); // non-diagonal and non-redundant entries of the covariance matrix of beta (indicating correlations between signatures)
  int d_min1 = d - 1;

  // for covariance matrix for random effects
  using namespace density;
  matrix<Type> theta_prime(n,d_min1); // The probabilities of each event (in ALR)
  matrix<Type> theta(n,d); // The probabilities of each event
  vector<Type> Q(n); // The probabilities of each event (marginal of exp of ALR)
  matrix<Type> Sigma(d_min1,d_min1);
  Sigma.fill(0);                   // Fill the whole matrix
  Sigma.diagonal() = exp(logs_sd_RE);  // Multiply diagonal by 10 to positive definite Sigma
  MVNORM_t<Type> nll_mvn(Sigma);


  for(int i=0;i<num_individuals;i++){ // likelihood for the random effects (multivariate normal)
    nll += VECSCALE_t(nll_mvn, exp(logs_sd_RE))(u_large.row(i));
  }

  for(int i=0;i<num_individuals;i++){ // likelihood for the random effects (multivariate normal)
    nll += nll_mvn(u_large.row(i));
  }

  // // third matrix: covariance matrix between signatures
  matrix<Type> realisation_beta(num_sigs, 1);
  UNSTRUCTURED_CORR_t<Type> mvn_beta_cov(logs_var_cov_beta);
  nll += VECSCALE_t(mvn_beta_cov, exp(logs_sd_cov_beta))(vector<Type>(realisation_beta-mu_beta));

  theta_prime = shared_intercept + x * realisation_beta + z * u_large;

  for(int l=0;l<n;l++){
    Q(l) = 0;
    for(int j=0;j<d_min1;j++){
      Q(l) += exp(theta_prime(l,j));
    }
  }

  for(int l=0;l<n;l++){
    for(int j=0;j<d_min1;j++){
      theta(l,j) = exp(theta_prime(l,j))/(Type(1.0)+Q(l));
    }
    theta(l,d_min1) = Type(1.0)/(Type(1.0)+Q(l));
  }

  for(int l=0; l<n; l++){ // Multinomial draws
    vector<Type> lth_row = Y.row(l);
    vector<Type> theta_l = theta.row(l);
    nll -= dmultinom(lth_row, theta_l, true);
  }


  return nll;

}

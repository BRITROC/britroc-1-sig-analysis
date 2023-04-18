// // https://figshare.com/articles/C_code_needed_for_running_R_code/4836503
// // Simple multinomial fitting.
// #include <TMB.hpp>

Incorrect, because the entries of cov_RE_part are a part of L, not of the covariance matrix itself

// template<class Type>
// Type objective_function<Type>::operator() ()
// {

//   Type nll;           // Define variable that holds the return value (neg. log. lik)
//   nll = Type(0.0);    // Assign value 1.2; a cast is needed.

//   DATA_MATRIX(Y); // observations (count matrix)
//   int d = Y.cols(); // number of features
//   int n = Y.rows(); // 2*n samples
//   DATA_INTEGER(num_individuals); // number of different indivdiduals (each with a random effect)
//   DATA_MATRIX(x); // matrix of covariates for fixed effects
//   DATA_MATRIX(z); // matrix for random effects
//   PARAMETER_MATRIX(beta); // coefficients for the fixed effects
//   PARAMETER_MATRIX(u_large); // coefficients for the random effects. Even though it is defined as matrix (for TMB matrix multiplication), it is a vector
//   PARAMETER_VECTOR(logs_sd_RE); // RE: variances (scaling)
//   PARAMETER(cov_RE_part); // RE: covariances, but only partial (see "If using a sparse covariance matrix..." section)
//   matrix<Type> theta_prime(n,d); // The logratios of each sample adn signature being zero

//   using namespace density;


//   // Using a sparse covariance matrix
//   int len_cov = (d*d-d)/2;
//   vector<Type> cov_RE(len_cov);
//   cov_RE.fill(0);
//   cov_RE[4] = cov_RE_part;

//   // [1] correlated random effects
//   // for covariance matrix for random effects
  
//   UNSTRUCTURED_CORR_t<Type> nll_mvn(cov_RE);
//   for(int i=0;i<num_individuals;i++){ // likelihood for the random effects (multivariate normal)
//     nll += VECSCALE_t(nll_mvn, exp(logs_sd_RE))(u_large.row(i));
//   }


//   theta_prime = x * beta + z * u_large;

//   for(int l=0; l<n; l++){ // binomial draws, for each observation individually
//     for(int j=0; j<d; j++){
//       Type num_zeros_lth_row = Y(l,j); // number of zeros in observation
//       Type p_l = (exp(theta_prime(l,j))) / (1 + exp(theta_prime(l,j))); // softmax transformation; pick first probability
//       nll -= log(p_l)*(num_zeros_lth_row) + log(1-p_l)*(1-num_zeros_lth_row); // binomial log-likelihood likelihood
//    }
//   }
//   return nll;

// }

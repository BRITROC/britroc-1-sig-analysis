I don't see how this can work
// //
// //  tmb_correlated_multinom_only_intercepts.cpp
// //  
// //
// //  Created by Lena Morrill on 03/12/2020.
// //

// #include <stdio.h>
// #include <TMB.hpp>

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
//   PARAMETER_VECTOR(logs_sd_RE); // RE: variances (scaling)
//   PARAMETER_VECTOR(cov_RE); // RE: covariances
//   matrix<Type> theta_prime(n,d); // The logratios of each sample adn signature being zero
//   matrix<Type> u_large(num_individuals,d); // The logratios of each sample adn signature being zero
//   Type d_min1 = d-1;


//   using namespace density;
    
//   // [1] correlated random effects
//   // for covariance matrix for random effects
//   UNSTRUCTURED_CORR_t<Type> nll_mvn(cov_RE);
// //  vector<Type> logs_sd_RE(d);
// //  logs_sd_RE.fill(0);
//   for(int i=0;i<num_individuals;i++){ // likelihood for the random effects (multivariate normal)
// //        nll += MVNORM(nll_mvn)(u_large.row(i));
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

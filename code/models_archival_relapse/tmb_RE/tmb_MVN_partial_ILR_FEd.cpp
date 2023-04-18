INCORRECT: CHECK tmb_MVN_partial_ILR_FEe INSTEAD

// // Partial MVM fit for partial ILR
// // we are estimating all u_large! this model cannot converge

// #include <TMB.hpp>


// template<class Type>
// Type objective_function<Type>::operator() ()
// {

//   Type nll;           // Define variable that holds the return value (neg. log. lik)
//   nll = Type(0.0);    // Assign value 1.2; a cast is needed.

//   DATA_MATRIX(Y); // observations (ILR)
//   int d = Y.cols(); // number of features - 1, i.e. the number of columns in Y
//   int n = Y.rows(); // the numbre of observations or samples, i.e. the number of rows in Y
//   DATA_MATRIX(x); // matrix of covariates for fixed effects
//   PARAMETER_MATRIX(beta); // coefficients for the fixed effects
//   PARAMETER_MATRIX(u_large); // coefficients for the random effects, it's a matrix
//   PARAMETER_VECTOR(logs_sd_RE); // Random effects: variances (scaling)
//   PARAMETER_VECTOR(cov_RE); // Random effects: covariances
//   matrix<Type> mu(n,d); // 

//   using namespace density;

//   // correlated random effects: full covariance matrix for random effects
//   mu = x * beta;

//   UNSTRUCTURED_CORR_t<Type> nll_mvn(cov_RE);
//   for(int i=0;i<n;i++){ // likelihood for the random effects (multivariate normal)
//     nll += VECSCALE_t(nll_mvn, exp(logs_sd_RE))(u_large.row(i)-mu.row(i));
//   }

//   for(int l=0; l<n; l++){

//     // compute how many zero exposures there are
//     int count_zeros = 0;
//     for(int loop_zeros=0; loop_zeros<d; loop_zeros++){
//       if(Y(l,loop_zeros) == 0){
//         count_zeros += 1;
//       }
//     }
//     // save in vector nonzero_idx_vector the indices of the non-zero exposures
//     int count_nonzeros = 0;
//     int count_nonzero = d-count_zeros;
//     vector <int> nonzero_idx_vector(count_nonzero);
//     for(int loop_zeros=0; loop_zeros<d; loop_zeros++){
//       if(Y(l,loop_zeros) != 0){
//         nonzero_idx_vector(count_nonzeros) = loop_zeros;
//         count_nonzeros += 1;
//       }
//     }


//     // for the non-zero exposures, compute the subcomposition for mu
//     vector<Type> mu_subcomposition_l(count_nonzero);
//     for(int loop_nonzeros=0; loop_nonzeros<count_nonzero; loop_nonzeros++){
//       mu_subcomposition_l(loop_nonzeros) = mu.row(l)(nonzero_idx_vector(loop_nonzeros))+u_large.row(l)(nonzero_idx_vector(loop_nonzeros));
//     }
//     // create the same for the observed data
//     vector<Type> ilr_subcomposition_l(count_nonzero);
//     for(int loop_nonzeros=0; loop_nonzeros<count_nonzero; loop_nonzeros++){
//       ilr_subcomposition_l(loop_nonzeros) = Y.row(l)(nonzero_idx_vector(loop_nonzeros));
//     }

//     // create a sigma of the appropriate size
//     matrix<Type> Sigma_ones(count_nonzero,count_nonzero);
//     Sigma_ones.fill(0);                   // Fill the whole matrix
//     vector<Type> vec1(count_nonzero);
//     vec1.fill(1.0);
//     Sigma_ones.diagonal() = vec1; // diagonal matrix
//     MVNORM_t<Type> N_basic(Sigma_ones);   // N_basic is now a Distribution  

//     nll += N_basic(ilr_subcomposition_l-mu_subcomposition_l);
  
//   }
//   return nll;

// }

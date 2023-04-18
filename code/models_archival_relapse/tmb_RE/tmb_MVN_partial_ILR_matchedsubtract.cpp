// see matchedsubtractB instead

// #include <TMB.hpp>


// template<class Type>
// Type objective_function<Type>::operator() ()
// {

//   Type nll;           // Define variable that holds the return value (neg. log. lik)
//   nll = Type(0.0);    // Assign value 1.2; a cast is needed.

//   DATA_MATRIX(Y); // observations (ILR)
//   int d = Y.cols(); // number of features - 1, i.e. the number of columns in Y
//   int n = Y.rows(); // the numbre of observations or samples, i.e. the number of rows in Y
//   // int num_individuals = Y.rows();
//   DATA_INTEGER(num_individuals); // number of different indivdiduals (each with a random effect)
//   DATA_IVECTOR(group_1_idx); 
//   DATA_IVECTOR(group_2_idx); 
//   DATA_MATRIX(x); // matrix of covariates for fixed effects
//   // DATA_MATRIX(z); // matrix for random effects
//   PARAMETER_MATRIX(beta); // coefficients for the fixed effects
//   PARAMETER_VECTOR(logs_sd_RE); // Random effects: variances (scaling)
//   // PARAMETER_VECTOR(cov_RE); // Random effects: covariances
//   matrix<Type> mu(n,d); // 

//   using namespace density;

//   mu = x * beta; // x is only 1s and beta only the second row

//   for(int i=0;i<num_individuals;i++){ 

//     // compute how many zero exposures there are
//     int count_zeros = 0;
//     for(int loop_zeros=0; loop_zeros<d; loop_zeros++){
//       if(Y(group_1_idx(i),loop_zeros) == 0){
//         count_zeros += 1;
//       }
//     }
//     // save in vector nonzero_idx_vector the indices of the non-zero exposures
//     int count_nonzeros = 0;
//     int count_nonzero = d-count_zeros;
//     vector <int> nonzero_idx_vector(count_nonzero);
//     for(int loop_zeros=0; loop_zeros<d; loop_zeros++){
//       if(Y(group_1_idx(i),loop_zeros) != 0){
//         nonzero_idx_vector(count_nonzeros) = loop_zeros;
//         count_nonzeros += 1;
//       }
//     }

//     // for the non-zero exposures, compute the subcomposition for mu
//     vector<Type> mu_subcomposition_l(count_nonzero);
//     for(int loop_nonzeros=0; loop_nonzeros<count_nonzero; loop_nonzeros++){
//       mu_subcomposition_l(loop_nonzeros) = mu.row(group_1_idx(i))(nonzero_idx_vector(loop_nonzeros));
//     }
//     // create the same for the observed data
//     vector<Type> ilr_subcomposition_l(count_nonzero);
//     for(int loop_nonzeros=0; loop_nonzeros<count_nonzero; loop_nonzeros++){
//       ilr_subcomposition_l(loop_nonzeros) = Y.row(group_1_idx(i))(nonzero_idx_vector(loop_nonzeros));
//     }

//     vector<Type> logsd_subcomposition_l(count_nonzero);
//     for(int loop_nonzeros=0; loop_nonzeros<count_nonzero; loop_nonzeros++){
//       logsd_subcomposition_l(loop_nonzeros) = logs_sd_RE(nonzero_idx_vector(loop_nonzeros));
//     }

//   /* ------------- */
//     count_zeros = 0;
//     for(int loop_zeros=0; loop_zeros<d; loop_zeros++){
//       if(Y(group_2_idx(i),loop_zeros) == 0){
//         count_zeros += 1;
//       }
//     }
//     // save in vector nonzero_idx_vector the indices of the non-zero exposures
//     count_nonzeros = 0;
//     count_nonzero = d-count_zeros;
//     for(int loop_zeros=0; loop_zeros<d; loop_zeros++){
//       if(Y(group_2_idx(i),loop_zeros) != 0){
//         nonzero_idx_vector(count_nonzeros) = loop_zeros;
//         count_nonzeros += 1;
//       }
//     }

//     vector<Type> ilr_subcomposition_l2(count_nonzero);
//     for(int loop_nonzeros=0; loop_nonzeros<count_nonzero; loop_nonzeros++){
//       ilr_subcomposition_l2(loop_nonzeros) = Y.row(group_2_idx(i))(nonzero_idx_vector(loop_nonzeros));
//     }

//   // /* ------------- */

//   matrix<Type> Sigma_ones(count_nonzero,count_nonzero);
//   Sigma_ones.fill(0);                   // Fill the whole matrix
//   vector<Type> vec1(count_nonzero);
//   Sigma_ones.diagonal() = exp(logsd_subcomposition_l); // diagonal matrix with values from the sd of the subcomposition
//   MVNORM_t<Type> nll_mvn(Sigma_ones);   // N_basic is now a Distribution  

//   nll += nll_mvn((ilr_subcomposition_l2- ilr_subcomposition_l- mu_subcomposition_l));
  
//   }
//   return nll;

// }

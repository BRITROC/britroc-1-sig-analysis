// // Partial MVM fit for partial ILR

// #include <TMB.hpp>


// template<class Type>
// Type objective_function<Type>::operator() ()
// {

//   Type nll;           // Define variable that holds the return value (neg. log. lik)
//   nll = Type(0.0);    // Assign value 1.2; a cast is needed.

//   DATA_MATRIX(Y); // observations (ILR)
//   int d = Y.cols(); // number of features - 1, i.e. the number of columns in Y
//   int n = Y.rows(); // the numbre of observations or samples, i.e. the number of rows in Y
//   DATA_INTEGER(num_individuals); // number of different indivdiduals (each with a random effect)
//   DATA_MATRIX(x); // matrix of covariates for fixed effects
//   DATA_MATRIX(z); // matrix for random effects
//   PARAMETER_MATRIX(beta); // coefficients for the fixed effects
//   PARAMETER_MATRIX(u_large); // coefficients for the random effects, it's a matrix
//   PARAMETER_VECTOR(logs_sd_RE); // Random effects: variances (scaling)
//   matrix<Type> mu(n,d); // 

//   using namespace density;

//   // correlated random effects: full covariance matrix for random effects
//   matrix<Type> Sigma_ones(d,d);
//   Sigma_ones.fill(0);                   // Fill the whole matrix
//   vector<Type> vec1(d);
//   vec1.fill(1.0);
//   Sigma_ones.diagonal() = vec1; // diagonal matrix
//   MVNORM_t<Type> N_basic(Sigma_ones);   // N_basic is now a Distribution  

//   UNSTRUCTURED_CORR_t<Type> nll_mvn(Sigma_ones);
//   for(int i=0;i<num_individuals;i++){ // likelihood for the random effects (multivariate normal)
//     nll += VECSCALE_t(nll_mvn, exp(logs_sd_RE))(u_large.row(i));
//   }

//   mu = x * beta + z * u_large;

//   for(int l=0; l<n; l++){
//     // ILR are multivariate normal draws
//     nll += N_basic(Y.row(l)-mu.row(l));
//     // nll += N_basic(ilr_subcomposition-mu_subcomposition);
//   }
//   return nll;

// }

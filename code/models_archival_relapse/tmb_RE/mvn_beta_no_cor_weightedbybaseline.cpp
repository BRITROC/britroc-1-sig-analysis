// //
// //  coped from mvn_beta_no_cor.cpp but weighting by the exposure of the baseline

total nonsense

// #include <TMB.hpp>
// template<class Type>
// Type objective_function<Type>::operator() ()
// {

//   Type nll;           // Define variable that holds the return value (neg. log. lik)
//   nll = Type(0.0);    // Assign value 0.0; a cast is needed.

//   DATA_MATRIX(Y); // observations (multivariate gaussian)
//   DATA_VECTOR(exposures_baseline); // observations (multivariate gaussian)
//   PARAMETER_VECTOR(logs_sd); // vector with scaling factors for the matrix of covariances (here full of zeros; uncorrelayed)
//   PARAMETER_MATRIX(beta); // coefficients for the intercept and change between conditions
//   DATA_MATRIX(x); // matrix of covariates for fixed effects

//   int n = Y.rows();
//   int d = Y.cols();

//   matrix<Type> mu(n,d); // 
//   mu = x * beta;

//   using namespace density;
//   vector<Type> vec1(d);
//   vec1.fill(1.0);
//   // multivariate normal
//   matrix<Type> cov(d,d);
//   cov.fill(0);
//   cov.diagonal() = vec1;
//   MVNORM_t<Type> nll_mvn(cov); // no correlations
//   for(int i=0;i<n;i++){
//       nll += VECSCALE_t(nll_mvn, exp(logs_sd))(Y.row(i)-mu.row(i))/exposures_baseline(i);
//   }

//   return nll;

// }

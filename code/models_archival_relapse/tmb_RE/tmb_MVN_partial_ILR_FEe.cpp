// Partial MVM fit for partial ILR, only for fixed effects
// correlated categories
// modified version of tmb_MVN_partial_ILR_Feb.cpp, and in substitution of tmb_MVN_partial_ILR_Fed.cpp, which wanted to estimate all u_large, which is not possible

#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() ()
{

  Type nll;           // Define variable that holds the return value (neg. log. lik)
  nll = Type(0.0);    // Assign value 0.0; a cast is needed.

  DATA_MATRIX(Y); // observations (ILR)
  int d = Y.cols(); // number of features - 1, i.e. the number of columns in Y
  int n = Y.rows(); // the numbre of observations or samples, i.e. the number of rows in Y
  DATA_MATRIX(x); // matrix of covariates for fixed effects
  PARAMETER_MATRIX(beta); // coefficients for the fixed effects
  PARAMETER_VECTOR(logsd); // log standard deviations for the log-ratios of categories
  // ------ cannot get estimates ------//
  PARAMETER_VECTOR(Lcov); // covariances
  // ------ cannot get estimates ------//
  matrix<Type> mu(n,d); // 

  using namespace density;

  mu = x * beta;

  // ------ cannot get estimates ------//
  // create a sigma of the whole set of signatures
  UNSTRUCTURED_CORR_t<Type> nll_full(Lcov);
  matrix<Type> cov_all(d,d);
  cov_all = nll_full.cov();
  // ------ cannot get estimates ------//

  for(int l=0; l<n; l++){

    // compute how many zero exposures there are
    int count_zeros = 0;
    for(int loop_zeros=0; loop_zeros<d; loop_zeros++){
      if(Y(l,loop_zeros) == 0){
        count_zeros += 1;
      }
    }
    // save in vector nonzero_idx_vector the indices of the non-zero exposures
    int count_nonzeros = 0;
    int count_nonzero = d-count_zeros;
    vector <int> nonzero_idx_vector(count_nonzero);
    for(int loop_zeros=0; loop_zeros<d; loop_zeros++){
      if(Y(l,loop_zeros) != 0){
        nonzero_idx_vector(count_nonzeros) = loop_zeros;
        count_nonzeros += 1;
      }
    }


    // for the non-zero exposures, compute the subcomposition for mu
    vector<Type> mu_subcomposition_l(count_nonzero);
    for(int loop_nonzeros=0; loop_nonzeros<count_nonzero; loop_nonzeros++){
      mu_subcomposition_l(loop_nonzeros) = mu.row(l)(nonzero_idx_vector(loop_nonzeros));
    }
    // create the same for the observed data
    vector<Type> ilr_subcomposition_l(count_nonzero);
    for(int loop_nonzeros=0; loop_nonzeros<count_nonzero; loop_nonzeros++){
      ilr_subcomposition_l(loop_nonzeros) = Y.row(l)(nonzero_idx_vector(loop_nonzeros));
    }
    // and the same for the standard deviations
    vector<Type> logsd_subcomposition_l(count_nonzero);
    for(int loop_nonzeros=0; loop_nonzeros<count_nonzero; loop_nonzeros++){
      logsd_subcomposition_l(loop_nonzeros) = logsd(nonzero_idx_vector(loop_nonzeros));
    }



    // ------ cannot get estimates ------//

    // // now get a subset of the covariance matrix, Sigma_subset
    matrix<Type> Sigma_subset(count_nonzero,count_nonzero);
    // // Eigen::Matrix<Type,Dynamic,Dynamic> Sigma_subset(count_nonzero,count_nonzero);
    for(int loop_nonzerosA=0; loop_nonzerosA<count_nonzero; loop_nonzerosA++){
      for(int loop_nonzerosB=0; loop_nonzerosB<count_nonzero; loop_nonzerosB++){
        Sigma_subset(loop_nonzerosA,loop_nonzerosB) = cov_all(nonzero_idx_vector(loop_nonzerosA),nonzero_idx_vector(loop_nonzerosB));
      }
    }

    // convert Sigma_subset to L_subset
    Eigen::LDLT<Eigen::Matrix<Type,Dynamic,Dynamic> > ldlt(Sigma_subset);
    matrix<Type> L_subset_mat = ldlt.matrixL();

    // // triangular matrix L to vector of theta entries in the lower diagonal, by row
    int length_count_nonzero_theta = (pow(count_nonzero,2) - count_nonzero)/2;
    vector<Type> L_subset(length_count_nonzero_theta);
    int count_nonzeros_it = 0;
    for(int loop_nonzerosrow=1; loop_nonzerosrow<count_nonzero; loop_nonzerosrow++){ // starting at the second column
      for(int loop_nonzeroscol=0; loop_nonzeroscol<loop_nonzerosrow; loop_nonzeroscol++){ // only lower triangle
        L_subset(count_nonzeros_it) = L_subset_mat(loop_nonzerosrow,loop_nonzeroscol);
        count_nonzeros_it += 1;
      }
    }

    // // use L_subset
    // // ILR are multivariate normal draws
    UNSTRUCTURED_CORR_t<Type> N_basic(L_subset); /// what it should be  
    nll += VECSCALE_t(N_basic, exp(logsd_subcomposition_l))(ilr_subcomposition_l-mu_subcomposition_l);
    // ------ cannot get estimates ------//


    // ------ for debugging only ------//
    // matrix<Type> Sigma_ones(count_nonzero,count_nonzero);     // create a sigma of the appropriate size
    // Sigma_ones.fill(0);                   // Fill the whole matrix
    // Sigma_ones.diagonal() = exp(logsd_subcomposition_l); // diagonal matrix with values from the sd of the subcomposition
    // MVNORM_t<Type> N_basic(Sigma_ones);   // N_basic is now a Distribution  
    // nll += N_basic(ilr_subcomposition_l-mu_subcomposition_l);
    // ------ for debugging only ------//
  
  }
  return nll;

}

#include <TMB.hpp>                                // Links in the TMB libraries
//#include <fenv.h>
#include <Eigen/Sparse>

template<class Type>
Type objective_function<Type>::operator() ()
{
  using namespace Eigen;
  using namespace tmbutils;
  using namespace density;
  
  // Bernoulli GLMM with random effects for smooth and iid terms
  // UPDATED to allow multiple smooth terms and random effects params
  // DATA
  DATA_VECTOR(y); // Response
  // eta = XF*betaF + XR*br + AU
  DATA_MATRIX(XF); // Unpenalized smooth term design matrix
  DATA_MATRIX(XR); // Penalized smooth term design matrix
  // PARAMETERS
  PARAMETER_VECTOR(betaF); // Unpenalized smooth term regression coefficients
  PARAMETER_VECTOR(bR); // Penalized smooth term regression coefficients
  PARAMETER_VECTOR(logsmoothing); // Vector of log smoothing parameters
  
  // DIMENSIONS
   // Input the dimensions as data, integer types
  DATA_INTEGER(p); // Number of smooth terms
  DATA_IVECTOR(r); // Vector of dimension of each smooth term. sum(r) = dim(bR), length(r) = p
  int bdim = r.sum();
  
  // TRANSFORMATIONS
  // Create the vector of smoothing params
  // vector<Type> reprec(logresd.size());
  vector<Type> smoothprec(logsmoothing.size());
  // for (int i=0;i<reprec.size();i++) reprec(i) = exp(-2.0*logresd(i));
  for (int i=0;i<smoothprec.size();i++) smoothprec(i) = exp(logsmoothing(i));

    
  // Linear predictor
  vector<Type> eta = XF*betaF + XR*bR;
  int n = y.size();
  
  // NEGATIVE log likelihood
  Type loglik = 0;
  for (int i=0;i<n;i++) loglik -= y(i)*eta(i) - log(1 + exp(eta(i)));
  
  // NEGATIVE log priors
  double pi = 3.141592653589793115998;
  // bR
  double Rdub = double(bdim);
  vector<double> rdub(p);
  for (int j=0;j<p;j++) rdub(j) = double(r(j));
  loglik -= -0.5*Rdub*log(2.0*pi);
  int tmpidx2 = 0;
  for (int j=0;j<p;j++) {
    loglik -=0.5*rdub(j)*logsmoothing(j);
    for (int l=0;l<r(j);l++) {
      loglik -= -0.5*smoothprec(j) * bR(tmpidx2)*bR(tmpidx2);
      tmpidx2++;
    }
  }
  
  return loglik; // Actually minus loglik
}

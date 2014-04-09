//// Implemented by Pu Wang http://cs.gmu.edu/~pwang7/gibbsBLasso.html
#include <Rmath.h>
#include <RcppArmadillo.h>
#include "utilarmadillo.hpp"
using namespace Rcpp ;
using namespace arma ;

RcppExport SEXP call_Bayes_Lasso( SEXP y_R, SEXP X_R, SEXP burn_R, SEXP thin_R, SEXP nmc_R );

void Bayes_Lasso( const vec y, const mat x, 
                  int burn, int thin, int nmc,
                  int initialize, 
                  vec betaini, double sigma2ini, vec invTau2ini, double lambdaini,
                  mat &betaSamples,    vec &sigma2Samples,  
                  mat &invTau2Samples, vec &lambdaSamples  );
#ifndef _BAYES_GLASSO_H
#define _BAYES_GLASSO_H

#include <RcppArmadillo.h>
using namespace Rcpp ;
using namespace arma ;
//RcppExport SEXP testBayesGLasso_Columnwise( SEXP S_R, SEXP Sig_R, SEXP C_R, SEXP n, SEXP burnin, SEXP nmc );
void BayesGLasso_Columnwise( mat S, size_t n,  mat Sig, mat C, size_t burnin, size_t nmc,
                            double a_lambda, double b_lambda,
                      cube& Sig_save, cube& C_save, vec& lambda_save  );

#endif

#ifndef _BAYES_AL_H
#define _BAYES_AL_H

#include <RcppArmadillo.h>

using namespace Rcpp ;
using namespace arma ; 

RcppExport SEXP call_Bayes_AL(SEXP y_R, SEXP X_R, SEXP b_R, SEXP delta_R, SEXP tau_R, 
                     SEXP burn_R, SEXP thin_R, SEXP nmc_R);

int Bayes_AL( vec y, mat X, int burn, int thin, int nmc,
               vec& b, double& delta, double& tau,
               mat& bSample, vec& deltaSam, vec& tauSam,
               bool update_b0, bool update_delta, bool update_tau
               ) ;
#endif

// # library(lars)
// # data(diabetes)
// # x = scale(diabetes$x)
// # y = scale(diabetes$y)

// # res=gibbsBLasso(x,y)


//  testGBLasso <- function(N=10000){
//  library(lars)
//  data(diabetes)
//  x = scale(diabetes$x)
//  y = scale(diabetes$y)
//  res =  .Call( "call_Bayes_AL", 
//          as.matrix(y),
//          as.matrix(x),
//          as.double(rnorm(dim(x)[2])),
//          as.double(1),
//          as.double(1),
//          as.integer(0),
//          as.integer(1),
//          as.integer(N),
//          PACKAGE = "testpackage" )
//  return(res)
//  }
//  #source("/Users/masanaoyajima/Documents/Rpackage/BayesLasso/gibbsBLasso.R")
//  #> res
//  # [1] -0.005121731 -0.147221767  0.321549292  0.199366391 -0.414064172
//  # [6]  0.237102499  0.025667757  0.095410286  0.437027590  0.041545572

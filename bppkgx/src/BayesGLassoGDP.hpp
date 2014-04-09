#ifndef _BAYES_GLASSO_GDP_H
#define _BAYES_GLASSO_GDP_H

#include <RcppArmadillo.h>
using namespace Rcpp ;
using namespace arma ;
//RcppExport SEXP testBayesGLassoGDP( SEXP S_R, SEXP ns_R, SEXP Sig_R, SEXP C_R, SEXP burnin, SEXP nmc );
void BayesGLassoGDP( mat S, size_t n,
                     mat Sig, mat C,  
                     size_t burnin, size_t nmc,
                     double s, double t,
                     cube& Sig_save, cube& C_save, mat& lambda_save  );
inline vec rand_ig( vec theta, vec chi );

#endif


// # R code to test the function

// bglasso <- function(S_R, Sig_R,   C_R,     n, burnin,  nmc, adaptive = FALSE ){
//     if( adaptive == FALSE ){
//         result = .Call( "testBayesGLasso_Columnwise", 
//             as.matrix(     S_R), 
//             as.matrix(   Sig_R), 
//             as.matrix(     C_R), 
//             as.integer(      n), 
//             as.integer( burnin), 
//             as.integer(    nmc),
//             PACKAGE = "test" )
//     }else{
//             result = .Call( "testBayesGLassoGDP", 
//             as.matrix(     S_R), 
//             as.matrix(   Sig_R), 
//             as.matrix(     C_R), 
//             as.integer(      n), 
//             as.integer( burnin), 
//             as.integer(    nmc),
//             PACKAGE = "test" )
//     }
//     return(result)
// }

// test<-function(type=6, adaptive = FALSE  ){
//     p = 30
//     n = 50

// if(type==1){
//     #### AR(1) case
//     # SigTrue = toeplitz(0.7.^[0:p-1]);
//     # CTrue = solve(SigTrue);
// } else if(type==2) {
//     #### AR(2) case 
//     # CTrue = toeplitz([1,0.5,0.25,zeros(1,p-3)]);
//     # SigTrue = solve(CTrue);
// } else if(type==3) {

//     #### Block case
//     SigTrue = diag(p);
//     SigTrue[1:(p/2),1:(p/2)] = 0.5*matrix(1,p/2,p/2)+(1-0.5)*diag(p/2);
//     SigTrue[(p/2+1):p,(p/2+1):p] = 0.5*matrix(1,p/2,p/2)+(1-0.5)*diag(p/2);

// } else if(type==4) {
//     ### Star case
//     CTrue = diag(p); 
//     CTrue[1,2:p] = 0.1; 
//     CTrue[2:p,1] = 0.1;
//     SigTrue = solve(CTrue);

// } else if(type==5) {
//     ### Circle case
//     # SigTrue = solve(toeplitz([2,1,zeros(1,p-3),0.9]));
//     # CTrue = solve(SigTrue); 
// } else if(type==6) {
//     ### Full case
//     CTrue = matrix(1,p,p)+diag(p);
//     SigTrue = solve(CTrue);
// }
// library(mvtnorm)
// Y      = rmvnorm( n, matrix( 0, p, 1 ), SigTrue );
// S      = t(Y)%*%(Y)
// Sig    = S/n
// C      = solve(Sig)
// burnin = 10000  
// nmc    = 10000
// res    = bglasso(S, Sig,   C,     n, burnin,  nmc, adaptive = adaptive )
// return(res)
// }


// # fitBGlasso=test(adaptive=TRUE)
// #    CTrue = matrix(1,p,p)+diag(p);
// #    SigTrue = solve(CTrue);
// # Sig   = fitBGlasso$Sig
// # dS    = dim(Sig)
// # #library(test); a=test()
// # pdf("SigmaAdaptive3.pdf",width=60,heigh=55)
// # density.plot.array(Sig,SigTrue)
// # #mtext( expression(paste("posterior mean of ", beta) ),
// # #      NORTH<-3, line=1, adj=0.5, cex=4 , outer=TRUE)
// # dev.off()


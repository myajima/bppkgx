#include "BayesGLassoGDP.hpp"
#include "BayesGLasso_Columnwise.hpp"
#include "utilarmadillo.hpp"
// SEXP testBayesGLasso_Columnwise( SEXP S_R, SEXP Sig_R, SEXP C_R, SEXP n, SEXP burnin, SEXP nmc )
// {
//   mat S    = Rcpp::as<mat>(   S_R );
//   mat Sig  = Rcpp::as<mat>( Sig_R );
//   mat C    = Rcpp::as<mat>(   C_R );
//   int NS   = INTEGER(      n )[0];  // No. of samples 
//   int BURN = INTEGER( burnin )[0];  // No. of samples 
//   int NMC  = INTEGER(    nmc )[0];  // No. of samples 
//   cube Sig_save;
//   cube C_save; 
//   vec lambda_save; 
//   BayesGLasso_Columnwise( S,  NS, Sig, C, BURN, NMC, Sig_save, C_save, lambda_save  );
//     return  Rcpp::List::create(
//                               Rcpp::Named(   "lambda" ) = wrap( lambda_save ),
//                               Rcpp::Named(      "Sig" ) = wrap(    Sig_save ),
//                               Rcpp::Named(        "C" ) = wrap(      C_save )               
//             ); 
// }

void BayesGLasso_Columnwise( mat S, size_t n, mat Sig, mat C,  size_t burnin, size_t nmc,
                              double a_lambda, double b_lambda,
                              cube& Sig_save, cube& C_save, vec& lambda_save  )
{
// Efficient Bayesian Graphical Lasso MCMC sampler using data-augmented
// block (column-wise) Gibbs sampler
//Input:
//     S = Y'*Y : sample covariance matrix * n
//     n: sample size
//     lambda:   |C|^(n/2) exp(-1/2 (SC) -  lambda/2 ||C||_1 );
//     Sig,C: initial Covariance and precision matrix C = inv(Sig);
//     burnin, nmc : number of MCMC burnins and saved samples
//     lambda ~ Ga(a_lambda,b_lambda)

//Output:
//     Sig_save: p by p by nmc matrices of saved posterior samples of covariance
//     C_save: p by p by nmc matrices of saved posterior samples of precision
//     lambda: 1 by nmc vector of saved posterior samples of lambda

//  Ref:  Wang 2012 Bayesain Graphical Lasso and Efficient Posterior
//  Computation

//  Written by Hao Wang & U of South Carolina
// double a_lambda = 1.0; 
// double b_lambda = 0.1;
size_t p   = S.n_rows; //size(S,1); 
uvec upperind = trimat_index( p, true );
uvec lowerind = trimat_index( p, true );


C_save     = zeros<cube>( p, p, nmc );// zeros(p,p,nmc); 
Sig_save   = C_save;
lambda_save = zeros<vec>( nmc ); //lambda_save = zeros(1,nmc);
mat tau = zeros<mat>( p, p );

umat ind_noi_all = zeros<umat>( p - 1 , p );
for( size_t j = 0; j < p; j++ ){
  int idx = 0;
  for( size_t i = 0; i < p; i++ ){ //1:p
    if( i != j ){
      ind_noi_all( idx++, j ) = i;
    }
  }
}//end

double apost = a_lambda + p * ( p + 1 )/2; 

size_t total_iter = burnin + nmc;
for( size_t iter=0; iter< total_iter; iter++  )  {
            
    //if( iter%1000 ==0 ){
        //fprintf('iter = %d \n',iter);
    //}//end
    // %%% Sample lambda 
    double bpost = b_lambda + accu(abs(C))/2;    
    double lambda = Rf_rgamma(apost,1/bpost);//gamrnd(apost,1/bpost,1);
    //%% sample tau off-diagonal        
    
    vec Cadjust = abs( C.elem( upperind ) ) ;
    uvec belowidx= find( Cadjust < 1.0e-6 );
    Cadjust.elem(belowidx) = ones<vec>(belowidx.n_elem)*1.0e-6;     
    vec lambda_prime = ones<mat>(Cadjust.n_elem,1) * pow( lambda, 2 );  
    //vec mu_prime(1); 
    vec mu_prime   = lambda / Cadjust;
    uvec aboveidx= find( mu_prime > 1.0e12 );
    mu_prime.elem( aboveidx ) = ones<vec>( aboveidx.n_elem ) * 1.0e12;
    vec tau_temp =  1.0/rand_ig( mu_prime, lambda_prime );
    tau.elem(upperind) = tau_temp;
    tau.elem(lowerind) = tau_temp;
    //%% sample Sig and C = inv(Sig)   
      uvec iidx(1);       
    for( size_t i=0; i< p; i++ ){
        iidx(0) = i;
        uvec ind_noi = ind_noi_all.col(i);//(:,i);
        vec tau_temp = tau.submat( ind_noi,iidx );

        mat Sig11 = Sig.submat( ind_noi, ind_noi ); 
        mat Sig12 = Sig.submat( ind_noi,    iidx );
        mat invC11 = Sig11 - Sig12 * Sig12.t() / Sig( i, i );
        mat Ci = ( S( i, i ) + lambda )* invC11 + diagmat( 1.0 / tau_temp );
        mat Ci_chol = chol( Ci );
        
        mat mu_i = solve(-Ci, S.submat( ind_noi, iidx ) );
        mat beta = mu_i + solve( Ci_chol, randn<mat>( p - 1, 1 ) );
        C.submat( ind_noi, iidx ) = beta;
        C.submat( iidx, ind_noi ) = beta.t();
        double gam = Rf_rgamma( n / 2 + 1,  2.0/(S(i,i) + lambda) );//gamrnd( n/2+1, (S(i,i)+lambda_ii)\2);
        C(i,i) = as_scalar(gam + beta.t() * invC11 * beta);
        mat invC11beta = invC11 * beta;
        Sig.submat( ind_noi, ind_noi ) = invC11 + invC11beta * invC11beta.t() / gam;
        Sig12 = -invC11beta / gam;
        Sig.submat(ind_noi, iidx ) = Sig12;
        Sig.submat(iidx, ind_noi ) = Sig12.t();
        Sig( i, i ) = 1.0 / gam;
        
    }//end
    if( iter >= burnin  )  {       
      Sig_save.slice( iter-burnin ) = Sig; 
      C_save.slice(   iter-burnin ) = C;
      lambda_save(    iter-burnin ) = lambda;
    }//end

  }//end

}
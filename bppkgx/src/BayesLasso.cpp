#include "BayesLasso.hpp"


RcppExport SEXP call_Bayes_Lasso( SEXP y_R, SEXP X_R, SEXP burn_R, SEXP thin_R, SEXP nmc_R )
{
   vec y    = Rcpp::as<vec>( y_R );
   mat X  = Rcpp::as<mat>(   X_R );
   int nmc  = INTEGER(    nmc_R )[0];
   int burn = INTEGER(   burn_R )[0];
   int thin = INTEGER(   thin_R )[0];
   mat betaSamples;  
   vec sigma2Samples;  
   mat invTau2Samples; 
   vec lambdaSamples;  
   GetRNGstate();
   vec betaini;
   double sigmaini;
   vec invTau2ini;
   double lambdaini;
   Bayes_Lasso( y, X, burn, thin, nmc, 1, betaini, sigmaini, invTau2ini, lambdaini,
                betaSamples, sigma2Samples, invTau2Samples, lambdaSamples  );
    // vec test1(1000000),test2(1000000);
    // vec a = ones<vec>(1000000);
    // vec b = 3.0*ones<vec>(1000000);
    // double start1 = std::clock();
    // std::cout << "test1 start " <<start1 << std::endl;
    // test1 = rinvgaussian( a, b );
    // std::cout << "test1 end " << std::clock()-start1 << std::endl;
    // //std::cout << "test1 diff " << std::clock() << std::endl;
    // double start2 = std::clock();
    // std::cout << "test2 start " << start2 << std::endl;
    // for (int i=0; i<1000000; i++ ) {
    //     test2(i) = inverseGaussian( a[i], b[i]);
    // }
    // std::cout << "test2 end " << std::clock()-start2 << std::endl;
   PutRNGstate();


   return  Rcpp::List::create(
                          //Rcpp::Named(        "test1" ) = wrap( test1 ),
                          //Rcpp::Named(     "test2" ) = wrap(  test2 ), 
                          Rcpp::Named(          "b" ) = wrap(    betaSamples ),
                          Rcpp::Named(      "sigma" ) = wrap(  sigma2Samples ),
                          Rcpp::Named(        "tau" ) = wrap( invTau2Samples ),
                          Rcpp::Named(     "lambda" ) = wrap(  lambdaSamples )             
        ); 
}

void Bayes_Lasso( const vec y, const mat x, 
                  int burn, int thin, int nmc,
                  int initialize, 
                  vec betaini, double sigma2ini, vec invTau2ini, double lambdaini,
                  mat &betaSamples,    vec &sigma2Samples,  
                  mat &invTau2Samples, vec &lambdaSamples  
                ) 
{

    int n     = x.n_rows;
    int m     = x.n_cols;
    mat XtX   = x.t() * x;   //Time saving
    vec xy    = x.t() * y;
    int r     = 1;
    double delta   = 0.1; // 1.78
    betaSamples    = zeros<mat>( nmc, m );
    sigma2Samples  = zeros<vec>( nmc );
    invTau2Samples = zeros<mat>( nmc, m );
    lambdaSamples  = zeros<vec>( nmc );
    vec beta;      
    double sigma2; 
    vec invTau2;   
    double lambda; 

    if( initialize == 1 ){
        beta    = solve( XtX + eye( m, m ), xy );
        vec residue = y - x * beta;
        sigma2  = dot( residue, residue) / n;
        invTau2 = 1 / (beta % beta);
        lambda  = m * sqrt(sigma2) / sum(abs(beta));
    } else {
        beta    = betaini;
        sigma2  = sigma2ini;
        invTau2 = invTau2ini;
        lambda  = lambdaini;
    }

    int idx = 0;
    int tot_iter = nmc * thin + burn;
    for( int k = 0; k < tot_iter; k++ ){
        // sample beta
        mat invD   = diagmat( invTau2 );
        mat invA   = inv( XtX + invD );
        mat bmean  = invA * xy;
        mat varcov = sigma2 * invA;
        beta = MVNORM( 0, bmean, varcov );

        // sample sigma2
        double shape   = 0.5 * ( n + m - 1.0 ); 
        vec residue    = ( y - x * beta );
        double scale   = 0.5 * as_scalar( ( dot( residue, residue ) + beta.t() * invD * beta ) );
        sigma2  = 1 / Rf_rgamma( shape, scale );

        // sample tau2
        vec muPrime        = sqrt( pow( lambda, 2 ) * sigma2 / pow( beta, 2 ) ); 
        double lambdaPrime = pow( lambda, 2 );
        for( int i = 0; i < m; i++ ) {
            invTau2( i ) = rinvgaussian( muPrime[i], lambdaPrime );
        }
        uvec belowidx= find( invTau2 < 1.0e-6 );
        invTau2.elem(belowidx) = ones<vec>(belowidx.n_elem)*1.0e-6;  
        // update lambda
        shape  = 0.5 * ( r + m );
        scale  = delta + 0.5 * sum( 1 / invTau2 );
        lambda = Rf_rgamma( shape, scale );

        // save
        if( ( k >= burn ) && ( k - burn ) % thin == 0 ){
            betaSamples.row(    idx ) = beta.t();
            sigma2Samples(      idx ) = sigma2;
            invTau2Samples.row( idx ) = invTau2.t();
            lambdaSamples(      idx ) = lambda;
            idx += 1;
        }
    }
}

// RcppExport SEXP testBLasso( SEXP y_R, SEXP X_R, SEXP burn_R, SEXP thin_R, SEXP nmc_R )
// {
//    vec y    = Rcpp::as<vec>( y_R );
//    mat X  = Rcpp::as<mat>(   X_R );
//    int nmc  = INTEGER(    nmc_R )[0];
//    int burn = INTEGER(   burn_R )[0];
//    int thin = INTEGER(   thin_R )[0];
//    mat betaSamples;  
//    vec sigma2Samples;  
//    mat invTau2Samples; 
//    vec lambdaSamples;  
//    GetRNGstate();
//    Bayes_Lasso( y, X, burn, thin, nmc, betaSamples, sigma2Samples, invTau2Samples, lambdaSamples  );
//     // vec test1(1000000),test2(1000000);
//     // vec a = ones<vec>(1000000);
//     // vec b = 3.0*ones<vec>(1000000);
//     // double start1 = std::clock();
//     // std::cout << "test1 start " <<start1 << std::endl;
//     // test1 = rinvgaussian( a, b );
//     // std::cout << "test1 end " << std::clock()-start1 << std::endl;
//     // //std::cout << "test1 diff " << std::clock() << std::endl;
//     // double start2 = std::clock();
//     // std::cout << "test2 start " << start2 << std::endl;
//     // for (int i=0; i<1000000; i++ ) {
//     //     test2(i) = inverseGaussian( a[i], b[i]);
//     // }
//     // std::cout << "test2 end " << std::clock()-start2 << std::endl;
//    PutRNGstate();


//    return  Rcpp::List::create(
//                           Rcpp::Named(        "test1" ) = wrap( test1 ),
//                           Rcpp::Named(     "test2" ) = wrap(  test2 ), 
//                           Rcpp::Named(       "beta" ) = wrap(    betaSamples ),
//                           Rcpp::Named(      "sigma" ) = wrap(  sigma2Samples ),
//                           Rcpp::Named(        "tau" ) = wrap( invTau2Samples ),
//                           Rcpp::Named(     "lambda" ) = wrap(  lambdaSamples )             
//         ); 
// }
//library(lars)
//data(diabetes)
//x = scale(diabetes$x)
//y = scale(diabetes$y)

//res=gibbsBLasso(x,y)


// testGBLasso <- function(){
// library(lars)
// data(diabetes)
// x = scale(diabetes$x)
// y = scale(diabetes$y)
// res =  .Call( "testBLasso", 
//         as.matrix(y),
//         as.matrix(x),
//         as.integer(10),
//         as.integer(2),
//         as.integer(10000),
//         PACKAGE = "test" )
// return(res)
// }
// #source("/Users/masanaoyajima/Documents/Rpackage/BayesLasso/gibbsBLasso.R")
// #> res
// # [1] -0.005121731 -0.147221767  0.321549292  0.199366391 -0.414064172
// # [6]  0.237102499  0.025667757  0.095410286  0.437027590  0.041545572

// slower version of inverse gaussian
// vec rinvgaussian( vec mu, vec lambda ) 
// {
//     //use.n = if ((length.n <- length(n)) > 1) 
//     //    length.n
//     //else if (!is.Numeric(n, integer.valued = TRUE, allowable.length = 1, 
//     //    positive = TRUE)) 
//     //    stop("bad input for argument 'n'")
//     //else n
//     //mu    = rep(mu, len = use.n)
//     //lambda = rep(lambda, len = use.n)
//     int n    = mu.n_elem;
//     vec u    = runif( n );
//     vec Z    = pow(rnorm(n),2);
//     vec phi  = lambda/mu;
//     vec y1   = 1 - 0.5 * (sqrt(pow( Z, 2 ) + 4 * phi % Z) - Z )/phi;
//     uvec idx = find((1 + y1) % u > 1 );
//     vec ans  = y1;
//     ans.elem(idx) = 1 / y1.elem( idx );
//     //ans = mu * ifelse((1 + y1) * u > 1, 1/y1, y1)
//     //ans[mu <= 0] = NaN
//     //ans[lambda <= 0] = NaN
//     return( ans );
// }
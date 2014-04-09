#include "BayesGLassoGDP.hpp"
#include "utilarmadillo.hpp"

void BayesGLassoGDP( mat S, size_t n, mat Sig, mat C, size_t burnin, size_t nmc, double s, double t,
                     cube& Sig_save, cube& C_save, mat& lambda_save  )
{
// Efficient Bayesian Adaptive Graphical Lasso (or generalized double Pareto)  MCMC sampler using data-augmented
// block(column) Gibbs sampler
//
//     y \sim N(0, C^{-1}), Y = (y_1,...,y_n), S = Y*Y';     p(Y) \propto |C|^(n/2) exp(-1/2 (SC) )
//     C_jj  \sim Exp(lamba_ii/2);                                p(C_jj) \prop \lambda_ii exp(- \lambda_ii/2 C_ii)
//     C_{ij} \sim DE(lambda_{ij} ),                      p(C_{ij} \propto exp(- lambda_{ij} |C_{ij}| )
//     lambda_{ij} \sim Ga(s,t).                   

//  This adaptive glasso prior is equivalent to generalized double Pareto: 

//     y \sim N(0, C^{-1}), Y = (y_1,...,y_n), S = Y*Y';  |C|^(n/2) exp(-1/2 (SC) )
//     C_ii   \sim  Exp(lamba_ii/2 );    
//     C_{ij} \sim  GDP( xi = t/s, alpha = s) where x \sim GDP(xi,alpha) has density
//     f(x) = 1/(2 xi)(1+ |x|/(alpha xi))^{-(1+alpha)}   


//Input:
//     S = Y'*Y : sample covariance matrix * n
//     n: sample size
//     Sig,C: initial Covariance and precision matrix C = inv(Sig);
//     burnin, nmc : number of MCMC burnins and saved samples

//Output:
//     Sig_save: p by p by nmc matrices of saved posterior samples of covariance
//     C_save: p by p by nmc matrices of saved posterior samples of precision
//     lambda: 1 by nmc vector of saved posterior samples of lambda

//  Ref:  Wang 2012 Bayesain Graphical Lasso and Efficient Posterior
//  Computation

//  Written by Hao Wang & U of South Carolina

size_t p = S.n_rows; //size(S,1); 
uvec upperind = trimat_index( p, true );
uvec lowerind = trimat_index( p, true );

C_save      = zeros<cube>( p, p, nmc );// zeros(p,p,nmc); 
Sig_save    = C_save;
lambda_save = zeros<mat>( nmc, p*( p - 1 )/ 2 );

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

//% Hyperparameter of \lambda_{ij} ~ Ga(s,t)
//double s = 1.0e-2;  
//double t = 1.0e-6;

double lambda_ii = 1.0;
vec lambda( p*( p - 1 )/ 2 );
size_t total_iter = burnin + nmc;
for( size_t iter = 0; iter< total_iter; iter++ ){

  vec Cadjust   = abs( C.elem( upperind ) ) ;
  uvec belowidx = find( Cadjust < 1.0e-12 );
  Cadjust.elem(belowidx) = ones<vec>( belowidx.n_elem )*1.0e-12;
  //%% sample lambda off-diagonal
  double s_post = 1.0 + s;  
  vec t_post = Cadjust + t;

  for( size_t k = 0; k < t_post.n_elem; k++ ){
    lambda(k) = Rf_rgamma( s_post, 1.0/t_post( k ) ); //lambda = gamrnd( s_post,1.0/t_post ); 
  }
                
  //%% sample tau off-diagonal        
  vec lambda_prime = pow( lambda, 2 );  
  vec mu_prime = lambda / Cadjust;
  uvec aboveidx= find( mu_prime > 1.0e12 );
  mu_prime.elem( aboveidx ) = ones<vec>( aboveidx.n_elem ) * 1.0e12;
  vec tau_temp =  1.0/rand_ig( mu_prime, lambda_prime );
  uvec taubelowidx= find( tau_temp < 1.0e-12 );
  tau_temp.elem(taubelowidx) = ones<vec>( taubelowidx.n_elem )*1.0e-12;
  //tau = offdiag * tau_temp;
  tau.elem(upperind) = tau_temp;
  tau.elem(lowerind) = tau_temp;
  //%% sample Sig and C = inv(Sig)    
  uvec iidx(1);    
  for( size_t i = 0; i < p; i++ ){ // 1:p
    iidx(0) = i;
    uvec ind_noi = ind_noi_all.col(i);//(:,i);
    vec tau_temp = tau.submat( ind_noi,iidx );
    mat Sig11 = Sig.submat( ind_noi, ind_noi ); 
    mat Sig12 = Sig.submat( ind_noi,    iidx );
    
    mat invC11 = Sig11 - Sig12 * Sig12.t() / Sig( i, i );
    mat Ci = ( S( i, i ) + lambda_ii )* invC11 + diagmat( 1.0 / tau_temp );
        Ci = ( Ci + Ci.t() ) / 2.0; 
    mat Ci_chol = chol( Ci );
      
    mat mu_i = solve( -Ci, S.submat( ind_noi, iidx ) );
    mat beta = mu_i + solve( Ci_chol, randn<mat>( p - 1, 1 ) );
    C.submat( ind_noi, iidx ) = beta;
    C.submat( iidx, ind_noi ) = beta.t();
    double gam = Rf_rgamma( n / 2 + 1,  2.0/( S( i, i ) + lambda_ii ) );//gamrnd( n/2+1, (S(i,i)+lambda_ii)\2);
    C( i, i ) = as_scalar( gam + beta.t() * invC11 * beta );
    
    //% Below updating Covariance matrix according to one-column change of precision matrix
    mat invC11beta = invC11 * beta;
    Sig.submat( ind_noi, ind_noi ) = invC11 + invC11beta * invC11beta.t() / gam;
    Sig12 = -invC11beta / gam;
    Sig.submat( ind_noi, iidx ) = Sig12;
    Sig.submat( iidx, ind_noi ) = Sig12.t();
    Sig( i, i ) = 1.0 / gam;
  }//end
      
  if( iter >= burnin ){           
    Sig_save.slice(  iter-burnin ) = Sig; 
    C_save.slice(    iter-burnin ) = C;
    lambda_save.row( iter-burnin ) = lambda.t();
  }//end
}
}


vec rand_ig( vec theta, vec chi ){

  //
  // START ig HELP START inversegauss HELP  START invgauss HELP
  // THE INVERSE GAUSSIAN DISTRIBUTION
  //
  // The Inverse Gaussian distribution is left skewed distribution whose
  // location is set by the mean with the profile determined by the
  // scale factor.  The random variable can take a value between zero and
  // infinity.  The skewness increases rapidly with decreasing values of
  // the scale parameter.
  //
  //
  // pdf(y) = sqrt(chi/(2*pi*y^3)) * exp(-chi./(2*y).*(y/theta-1).^2);
  // cdf(y) = normcdf(sqrt(chi./y).*(y/theta-1)) + ...
  //            exp(2*chi/theta)*normcdf(sqrt(chi./y).*(-y/theta-1));
  //
  //   where  normcdf(x) = 0.5*(1+erf(y/sqrt(2))); is the standard normal CDF
  //         
  // Mean     = theta;
  // Variance = theta^3/chi;
  // Skewness = sqrt(9*theta/chi);
  // Kurtosis = 15*mean/scale;
  // Mode = theta/(2*chi)*(sqrt(9*theta^2+4*chi^2)-3*theta);
  //
  // PARAMETERS:
  //  theta - location; (theta>0)
  //  chi - scale; (chi>0)
  //
  // SUPPORT:
  //  y,  y>0
  //
  // CLASS:
  //   Continuous skewed distribution
  //
  // NOTES:
  //   1. There are several alternate forms for the PDF, 
  //      some of which have more than two parameters
  //   2. The Inverse Gaussian distribution is often called the Inverse Normal
  //   3. Wald distribution is a special case of The Inverse Gaussian distribution
  //      where the mean is a constant with the value one.
  //   4. The Inverse Gaussian distribution is a special case of The Generalized
  //        Hyperbolic Distribution
  //
  // USAGE:
  //   randraw('ig', [theta, chi], sampleSize) - generate sampleSize number
  //         of variates from the Inverse Gaussian distribution with 
  //         parameters theta and chi;
  //   randraw('ig') - help for the Inverse Gaussian distribution;
  //
  // EXAMPLES:
  //  1.   y = randraw('ig', [0.1, 1], [1 1e5]);
  //  2.   y = randraw('ig', [3.2, 10], 1, 1e5);
  //  3.   y = randraw('ig', [100.2, 6], 1e5 );
  //  4.   y = randraw('ig', [10, 10.5], [1e5 1] );
  //  5.   randraw('ig');
  // 
  // SEE ALSO:
  //   WALD distribution
  // END ig HELP END inversegauss HELP  END invgauss HELP 
  
  // Method:
  //
  // There is an efficient procedure that utilizes a transformation
  // yielding two roots.
  // If Y is Inverse Gauss random variable, then following to [1]
  // we can write:
  // V = chi*(Y-theta)^2/(Y*theta^2) ~ Chi-Square(1),
  //
  // i.e. V is distributed as a chi-square random variable with
  // one degree of freedom.
  // So it can be simply generated by taking a square of a
  // standard normal random number.
  // Solving this equation for Y yields two roots:
  //
  // y1 = theta + 0.5*theta/chi * ( theta*V - sqrt(4*theta*chi*V + ...
  //      theta^2*V.^2) );
  // and
  // y2 = theta^2/y1;
  //
  // In [2] showed that  Y can be simulated by choosing y1 with probability
  // theta/(theta+y1) and y2 with probability 1-theta/(theta+y1)
  //h
  // References:
  // [1] Shuster, J. (1968). On the Inverse Gaussian Distribution Function,
  //         Journal of the American Statistical Association 63: 1514-1516.
  //
  // [2] Michael, J.R., Schucany, W.R. and Haas, R.W. (1976).
  //     Generating Random Variates Using Transformations with Multiple Roots,
  //     The American Statistician 30: 88-90.
  //
  //

  //assert(theta.n_elem ==  chi.n_elem);
  size_t sampleSize = theta.n_elem; //= max( length( theta ), length( chi ) );
  vec chisq1 = pow( randn<mat>( sampleSize, 1 ), 2 );//
  vec y = theta + 0.5 * ( theta /chi % (theta % chisq1 - sqrt( 4 * theta % chi % chisq1 + pow( theta, 2 ) % pow( chisq1, 2 ) ) ));
  uvec l = find( randu<mat>( sampleSize, 1 ) >= theta/( theta + y ) );
  y.elem( l ) = pow( theta.elem( l ), 2 ) / y.elem(l);
  return( y ); 
}



// SEXP testBayesGLassoGDP( SEXP S_R, SEXP ns_R, SEXP Sig_R, SEXP C_R, SEXP burnin, SEXP nmc )
// {
//   mat S    = Rcpp::as<mat>(   S_R );
//   mat Sig  = Rcpp::as<mat>( Sig_R );
//   mat C    = Rcpp::as<mat>(   C_R );
//   int NS   = INTEGER(   ns_R )[0];  // No. of samples 
//   int BURN = INTEGER( burnin )[0];  // No. of samples 
//   int NMC  = INTEGER(    nmc )[0];  // No. of samples 
//   cube Sig_save;
//   cube C_save; 
//   mat lambda_save; 
//   BayesGLassoGDP( S, Sig, C, NS, BURN, NMC, Sig_save, C_save, lambda_save  );
//     return  Rcpp::List::create(
//                               Rcpp::Named(   "lambda" ) = wrap( lambda_save ),
//                               Rcpp::Named(      "Sig" ) = wrap(    Sig_save ),
//                               Rcpp::Named(        "C" ) = wrap(      C_save )               
//             ); 
// }

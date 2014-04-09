#include "utilarmadillo.hpp"
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <algorithm>
//#include <stdio.h>
#include <stdarg.h>
//#include <omp.h>
//#include <assert.h>


// Inverse Gamma
double dinvgamma( double x, double shape, double scale, int give_log ) {
  if( shape <= 0 ) Rf_error( "shape must be > 0" );
  if( scale <= 0 ) Rf_error( "scale must be > 0" );
  double alpha = shape, beta = scale;
  double log_density = alpha * log( beta ) - Rf_lgammafn( alpha ) - ( alpha + 1 ) * log( x ) - ( beta/x );
  return( give_log ? log_density : pow( M_E, log_density ) );
}

double rinvgamma( double shape, double scale ) {
  if( shape <= 0 ) Rf_error( "shape must be > 0" );
  if( scale <= 0 ) Rf_error( "scale must be > 0" );
  //NumericVector samp = rgamma( 1, shape, scale );
  //return( samp(0)  );
  return( 1 / Rf_rgamma( shape, 1 / scale ) );
}
 
// Inverse Chi-square
double dinvchisq( double x, double df, double scale, int give_log )
{
  if( df <= 0 ) Rf_error( "df must be > 0" );
  double nu = df / 2;
  if( x <= 0 ){
    return( give_log ? R_NegInf : 0 );
  }
  double log_density = nu * log( nu ) - log( Rf_lgammafn( nu ) ) + nu * log( scale ) - ( nu + 1 ) * log( x ) - ( nu * scale / x );
  return( give_log ? log_density : pow( M_E, log_density ) );
}

double rinvchisq( double df, double scale )
{
  if( df <= 0 ) Rf_error( "df must be > 0" );
  if( scale <= 0 ) Rf_error( "scale must be > 0" );
  //NumericVector csq = rchisq( 1, df );
  //return( ( df * scale ) / csq(0)  ); 
  return( ( df * scale ) / Rf_rchisq( df )  ); 
  
}

double dbetabinom( double x, int n, double a, double b, int give_log ) {
  //assert( a > 0 && b > 0 );
  //assert( n > 0 && n > x && x >=0 );
  if( a <= 0 ) Rf_error( "a must be > 0" );
  if( b <= 0 ) Rf_error( "b must be > 0" );
  if( n <= 0 ) Rf_error( "n must be > 0" );
  if( x <  0 ) Rf_error( "x must be >= 0" );
  if( n <= x ) Rf_error( "n must be > x" );
  return( give_log?
  Rf_lchoose( n, x ) + Rf_lbeta( x + a, n - x + b ) - Rf_lbeta( a, b ) :
  Rf_choose( n, x ) * Rf_beta( x + a, n - x + b ) / Rf_beta( a, b ) );
}

double rbetabinom( double n, double alpha, double beta ) {
  //NumericVector prob = rbeta( 1, alpha, beta );
  //NumericVector samp = rbinom( 1, n, prob(0) );
  //return( samp(0) );
  double prob = Rf_rbeta( alpha, beta );
  return( Rf_rbinom( n, prob ) );
}


double rinvgaussian( double mu, double lambda ) {
    double v    = Rf_rnorm( 0, 1 );   // sample from a normal distribution with a mean of 0 and 1 standard deviation
    double y    = v * v;
    double mumu = mu * mu;
    double x    = mu + ( mumu * y )/( 2 * lambda ) - ( mu / ( 2 * lambda ) ) * sqrt( 4 * mu * lambda * y + mumu * y * y );
    return( Rf_runif( 0, 1 ) <= ( mu )/( mu + x ) ? x :  mumu / x );
}

double rand_truncated_normal_rej_below( const double mu, const double sig, const double below){
  //% Generate one sample from N(mu,sig^2)1_{x>below}
  //% Ref: 
  //%     Proposition 2.3  C. P. Robert 'Simulation of truncated normal variables'
  double mu_neg = (below - mu)/sig;
  double x;
  if( mu_neg < 0){
      x = Rf_rnorm(0,1);
      while( x < mu_neg){
        x = Rf_rnorm(0,1);
      }
  }else{
      double alpha = (mu_neg + sqrt(pow(mu_neg,2)+4))/2; //% Optimal alpha
      x = Rf_rexp(1/alpha) + mu_neg;
      while( log( Rf_runif( 0, 1 ) ) > ( - pow( ( x - alpha ), 2 ) / 2 ) ){
        x = Rf_rexp(1/alpha) + mu_neg;
      }
  }
  x = x * sig + mu;
  return(x);
}
double rand_truncated_normal_rej_above( const double mu, const double sig, const double above ){
//% Generate one sample from N(mu,sig^2)1_{x < above}
//% Ref: 
//%     Proposition 2.3  C. P. Robert 'Simulation of truncated normal variables'
  double below  = -above;
  double negmu  = -mu;
  double x = rand_truncated_normal_rej_below( negmu, sig, below);
  x = -x;
  return( x );
}


//# Sample from truncated normal N(mu,sig^2)1_{lower<r<upper}
double rand_tn(const double mu,const double sig, const double lower,const double upper){

    NumericVector x1 = NumericVector( 1, ( lower-mu ) / sig );
    NumericVector x2 = NumericVector( 1, ( upper-mu ) / sig );
    NumericVector p1 = pnorm( x1, 0.0, 1.0 );
    NumericVector p2 = pnorm( x2, 0.0, 1.0 );
    double x;
    if( p2(0) < 1e-30 ){
        x = upper;
    } else {
        double u = ( p2(0) - p1(0) ) * Rf_runif( 0, 1 ) + p1(0);
        u = Rf_fmin2( Rf_fmax2( 1e-30, u ), 1-1e-16 );
        NumericVector u2 = NumericVector( 1, u );
        NumericVector x4 = qnorm( u2, 0.0, 1.0 );
        x = x4(0) * sig + mu;
    }   
    return(x);
}
// Subroutine for rand_truncated_normal_rej_twosided
double oneside(const double mu, const double sig, const double below, const double above, const double mu_neg,const double mu_pos)
{
  double x;
  if( (mu_pos-mu_neg) > 2 * sqrt(exp(1))/(mu_neg+sqrt(pow(mu_neg,2)+4))*exp((pow(mu_neg,2)-mu_neg*sqrt(pow(mu_neg,2)+4))/4) ){
    x = rand_truncated_normal_rej_below( mu, sig, below );
    while( x > above ){
      x = rand_truncated_normal_rej_below( mu, sig, below );
    }
  }else{
    x = Rf_runif(0,1)*(mu_pos-mu_neg)+mu_neg;
    while( log(Rf_runif(0,1)) > (pow(mu_neg,2) - pow(x,2))/2 ){
      x = Rf_runif(0,1)*(mu_pos-mu_neg)+mu_neg;            
    }
    x = x*sig + mu;
  }
  return(x);
}

double rand_truncated_normal_rej_twosided(const double mu, const double sig, const double below, const double above){
  //% Generate one sample from N(mu,sig^2)1_{above > x > below}
  //% Ref: 
  //%     Proposition 2.3  C. P. Robert 'Simulation of truncated normal variables'
  double mu_neg = (below - mu)/sig;
  double mu_pos = (above - mu)/sig;
  double x;
  if( mu_neg>0 && mu_pos>0 ){
    x = oneside(mu, sig, below, above, mu_neg, mu_pos );
  } else if( mu_neg<0 && mu_pos<0  ) {
    double aa = mu_neg;
    mu_neg = -mu_pos;
    mu_pos = -aa;
    double negmu = -mu;
    double negabove = -below;
    double negbelow = -above;
    x = oneside(negmu, sig, negbelow, negabove, mu_neg, mu_pos );
    x = -x;
  } else {
    x = rand_tn( mu, sig, below, above ); 
  } 
  return(x);
}
vec rtnorm( const int n, const double mu, const double sigma, const double lower, const double upper )
{
    NumericVector lv = NumericVector( 1, lower );
    NumericVector uv = NumericVector( 1, upper );
    NumericVector lo = pnorm( lv, mu, sigma );
    NumericVector hi = pnorm( uv, mu, sigma );
    NumericVector vs = runif( n, lo(0), hi(0) );
    NumericVector re = qnorm( vs, mu, sigma );
    vec result(re.begin(), re.size(), false );
    return(result);
}
 

vec rtgamma( const int n, const double shape, const double scale, const double lower, const double upper )
{
    NumericVector lv = NumericVector( 1, lower );
    NumericVector uv = NumericVector( 1, upper );
    NumericVector lo = pgamma( lv, shape, scale );
    NumericVector hi = pgamma( uv, shape, scale );
    NumericVector vs = runif( n, lo(0), hi(0) );
    NumericVector re = qgamma( vs ,shape, scale );
    vec result(re.begin(), re.size(), false );
    return(result);
}

uvec trimat_index( int p, bool upper ){
    mat m0    = ones<mat>( p, p );
    m0.diag() = zeros<vec>( p );
    uvec idx;
    //lower triangular index
    if( upper == false ){
        idx = find( trimatl( m0 ) );
    } else{ // upper triangular index
        idx = find( trimatu( m0 ) );
    }
    return(idx);
}




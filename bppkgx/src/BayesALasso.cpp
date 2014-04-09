#include "BayesALasso.hpp"

extern"C"{ 
#include "arms.h"
}

SEXP call_Bayes_AL( SEXP y_R, SEXP X_R, SEXP b_R, SEXP delta_R, SEXP tau_R, 
                    SEXP burn_R, SEXP thin_R, SEXP nmc_R )
{

   vec y    = Rcpp::as<vec>( y_R );
   mat X    = Rcpp::as<mat>( X_R );
   vec bini = Rcpp::as<vec>( b_R );
   double deltaini = REAL( delta_R )[0];
   double tauini   = REAL( tau_R )[0];
   int THIN = INTEGER(   thin_R )[0];  // No. of samples 
   int BURN = INTEGER(   burn_R )[0];  // No. of samples 
   int NMC  = INTEGER(    nmc_R )[0];  // No. of samples 
   mat bSample;
   vec deltaSam, tauSam;

   Bayes_AL( y, X, BURN, THIN, NMC,
             bini, deltaini, tauini,
             bSample, deltaSam, tauSam,
             true, true, true ); 

   return  Rcpp::List::create(
                               Rcpp::Named(          "b" ) = wrap(     bSample ),
                               Rcpp::Named(      "delta" ) = wrap(    deltaSam ),
                               Rcpp::Named(        "tau" ) = wrap(      tauSam )               
             ); 
 }
 /*********************************************************************
  *
  * ld_bj
  *
  * log density function for b_j used in Bayesian Adaptive Lasso
  *
  *********************************************************************/
  struct bj_parm {
    double bj_bar;
    double sigma_j2;
    double kappa_j;
  };
 
  double ld_bj( double x, void *bj_parameters )
  {
    struct bj_parm *d;
    double y, tmp;

    //d   = bj_parm;
    d  = reinterpret_cast<bj_parm*>( bj_parameters );
    tmp = x - d->bj_bar;
    y   = -tmp * tmp / ( 2.0 * d->sigma_j2 ) - fabs( x ) / d->kappa_j;
    return y;
  }

 /********************************************************************
  *
  * ld_delta
  *
  * log density function for delta used in Bayesian AL or Bayesan t
  *
  ********************************************************************/
  struct delta_parm {
    double tau;
    double sum_log_v; // sum_{j=0}^{p-1} log(sigma_j^2)
    int p;
  };
 
  double ld_delta( double x, void *delta_parameters )
  {
    struct delta_parm *d;
    double y;
 
    d = reinterpret_cast<delta_parm*>( delta_parameters );
    y = x * ( d->p * log( d-> tau ) - d->sum_log_v ) - d->p * Rf_lgammafn( x );
    return y;
  }

 /*********************************************************************
  *
  * bayes_AL
  *
  * The Bayesian version of adaptive Lasso
  *
  *********************************************************************/
int Bayes_AL( vec y, mat X, int burn, int thin, int nmc,
               vec& b, double& delta, double& tau,
               mat& bSample, vec& deltaSam, vec& tauSam,
               bool update_b0, bool update_delta, bool update_tau
               ) 
{
  Rprintf("0\n");
  int n = y.n_elem; 
  int p = X.n_cols;
  struct delta_parm delta_data;
  delta_data.p = p;
  struct bj_parm bj_data;
  /* variables for Adpative Rejection Metroplois Sampling (arms) */
  int err, neval=0, ninit = 4, npoint = 50, nsamp=1, ncent = 4;
  double binit[10];
  double bl = -100.0, br = 100.0;
  double xinit[10]={ 0.01, 2.0, 10.0, 20.0 }, xl = 0.0, xr = 500.0;
  double xcent[10], qcent[10] = {5., 30., 70., 95.};
  double convex = 1.0;
  int dometrop  = 0;
  double xprev  = 0.01;
  double b0_bar, bj_bar, s0, sj,bj_t, sj_t, bj;
  /* allocate memory */
   vec v    =  vec( p ); //(double *)calloc(p, sizeof(double));
   vec res  =  vec( n ); //(double *)calloc(n, sizeof(double));
   vec Xj2  =  vec( p );  //(double *)calloc(p, sizeof(double));
  ivec Jset = ivec( p ); //(int *)callo
  bSample  = mat( nmc, p ); 
  deltaSam = vec( nmc ); 
  tauSam   = vec( nmc );
  double xij;
  for( int j = 0; j < p; j++ ){
    Xj2( j ) = 0.0;
    for( int i = 0; i < n; i++ ){
      xij = X( i, j );
      Xj2( j ) += xij * xij;
    }
  }
  /**
   * step 1. Initialization
   */
  double b0 = 0.0;
  double v0 = 0.01;

  v = ones<vec>(p)*0.01;

  int tot_iter = burn + thin * nmc;
  // Set random seed ---------------------------------------------------------------
  //GetRNGstate();  

  for( int j = 0; j < p; j++ ){
    Jset(j) = j;
  }
   int k = 0;
   for( int w = 0; w < tot_iter; w++ ){    
     /**
      * step 2. Update b0
      */
     res    = y - X*b;
     if( update_b0 ){
        b0_bar = mean( res );
        s0     = v0 / n;
        b0     = Rf_rnorm( b0_bar, sqrt( s0 ) );
        res    = res-b0;
     }
     /**
      * step 3.1 choose the order of updating b[j]
      */

     // if(*b_update_order > 1){
     //   get_b_order(*b_update_order, Jset, w, n, p, y, X, b);
     // }  


     /**
      * step 3.2 Update b[j]
      */
     for( int j = 0; j < p; j++ ){
       vec Xj = X.col( j );
       res = res + Xj * b( j );
       bj_bar = sum( Xj % res )/ Xj2( j );
       sj   = v0/Xj2( j );
       bj_data.bj_bar   = bj_bar;
       bj_data.sigma_j2 = sj;
       bj_data.kappa_j  = v( j );
       xprev            = b( j );
       if( bj_bar > sj/v( j ) ){
         bj_t = bj_bar - sj/v( j );
       }else if( bj_bar < -sj/v( j ) ){
         bj_t = bj_bar + sj/v( j );
       }else{
         bj_t = 0.0;
       }
       sj_t = sqrt( sj );
       binit[0] = bj_t - 2.0*sj_t;
       binit[1] = bj_t - 0.5*sj_t;
       binit[2] = bj_t + 0.5*sj_t;
       binit[3] = bj_t + 2.0*sj_t;
       err = arms( binit,ninit,&bl,&br,ld_bj,&bj_data,&convex,npoint,
                    dometrop,&xprev,&bj,nsamp,qcent,xcent,ncent,&neval);
       if( neval > 1000 ){
         Rf_warning( "In total %d evaluations is used in ARMS\n", neval );
       }
       if( err > 0 ){
         Rf_warning("error in arms, error code = %d, b[j] is not updated\n", err);
         bj = b(j);
       }else{
         b(j) = bj;
       }
       /* add the effect of Xj back into the residual */
       res = res - Xj*bj;
     }

     /**
      * step 4 updatedte v0, i.e. sigma_0^2
      */
     v0 = dot(res,res)/Rf_rchisq(n);

  
     /**
      * step 5 Update v[j], i.e. kappa_j
      */
     v = 2.0 * ( abs( b ) + tau ) / Rf_rchisq( 2 + 2 * delta );
     
     /**
     * step 6 Update tau, if needed
     */
     if( update_tau ){
       tau = Rf_rgamma( p * delta, 1.0/ sum( 1.0 / v ) );
     }

     /**
     * step 7 Update delta, if needed
     */ 
     if( update_delta ){
         delta_data.tau = tau;
         delta_data.sum_log_v = sum( log( v ) );
         xprev = delta;
         err = arms( xinit, ninit, &xl, &xr, ld_delta, &delta_data, &convex,
                     npoint, dometrop, &xprev, &delta, nsamp, qcent, xcent, ncent, &neval );
         if( err > 0 ){
           Rf_error("error in arms, error code = %d\n",err);
         }
         for( int j = 0; j < ninit; j++ ){
           xinit[j] = xcent[j];
         }
     }

    /*****************************************/
    /*  last step Record b when neccesary    */
    /*****************************************/
    if( ( w >= burn ) && ( w - burn ) % thin == 0 ){
         bSample.row( k ) = b.t();
         deltaSam( k ) = delta;
         tauSam( k ) = tau;
         k += 1;
     }
   }
   //PutRNGstate();
   return( err?0:err );
 }
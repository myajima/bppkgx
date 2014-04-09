#ifndef UTILR_HEADER_20120507
#define UTILR_HEADER_20120507

/*
 *  utilR.h
 *
 *  Created by Masanao Yajima on 2012/05/07.
 *  Copyright 2012 __MyCompanyName__. All rights reserved.
 *
 */
#include <RcppArmadillo.h>
#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <Rdefines.h>
//#include <RInside.h>                            // for the embedded R via RInside

using namespace Rcpp;
using namespace arma;

extern "C"{
  SEXP as_list( char**, int, ... );
}

// Evaluate R expressions --------------------------------------------------------------------
SEXP mkpointer(double *x, int n);
//SEXP matrixpointer(double **x, int nrow, int ncol);
SEXP mkint(int i);

//mat fev(SEXP expr, SEXP rho, int cK, int cT, rowvec x, int i);
/* Evaluate R Expression taking as arguments *par and i ------------------------- */

template< typename T, template <typename> class ARMA_VECTOR_TYPE >
mat fev( SEXP expr, SEXP rho, int clK, int clT, ARMA_VECTOR_TYPE<T> x, int i ) //double *x, int n, int i )
{
   SEXP     r_fit;
   PROTECT( r_fit = Rf_allocMatrix( REALSXP, clK, clT ) );   // allocate SEXP matrix
   Rf_defineVar( Rf_install( "param" ), mkpointer( x.memptr(), x.n_elem ), rho ); // assign param
   Rf_defineVar( Rf_install(     "i" ), mkint( i ),                        rho ); // assign i
   r_fit = Rf_eval( expr, rho );                             // Evaluates the R expression          
   //NumericMatrix res( r_fit );                               // Create RCpp matrix 
   //mat result( res.begin(), res.nrow(), res.ncol(), false ); // Create Armadillo matrix
   mat result  = Rcpp::as<mat>( r_fit );
   UNPROTECT( 1 );                                           // Unprotect r_fit
   return( result );
}

template< typename T, template <typename> class ARMA_MATRIX_TYPE >
cube fev_parallel( SEXP expr, SEXP rho, SEXP cluster, int clK, int clT, int clN, ARMA_MATRIX_TYPE<T> x, int i ) //double *x, int n, int i )
{
     Rprintf("1\n");
   SEXP     r_fit;
   PROTECT( r_fit = Rf_allocMatrix( REALSXP, clK, clT ) );   // allocate SEXP matrix
   SEXP ans;
      Rprintf("2\n");
   PROTECT( ans = Rf_allocVector( REALSXP, x.n_rows * x.n_cols ) );
      Rprintf("31\n");
   for( unsigned int i = 0; i < x.n_rows; i++ ) {
      for( unsigned int j = 0; j < x.n_cols; j++ ){ 
         REAL( ans )[i+ j*x.n_rows] = x(i,j); 
      }
   }
         Rprintf("8\n");
   Rf_defineVar( Rf_install( "param" ), ans, rho ); // assign param
     Rprintf("4\n");
   Rf_defineVar( Rf_install( "cl" ), cluster, rho ); // assign param
   //Rf_defineVar( Rf_install(     "i" ), mkint( i ),                        rho ); // assign i
        Rprintf("5\n");
   r_fit = Rf_eval( expr, rho );                             // Evaluates the R expression          
   //NumericMatrix res( r_fit );                               // Create RCpp matrix 
   //mat result( res.begin(), res.nrow(), res.ncol(), false ); // Create Armadillo matrix
   //cube result  = Rcpp::as<cube>( r_fit );
        Rprintf("6\n");
   NumericVector vecArray( r_fit );
   cube result( vecArray.begin(), clT, clK, clN, false );
        Rprintf("7\n");
   UNPROTECT( 2 );                                           // Unprotect r_fit
   return( result );
}

#endif // UTILR_HEADER_20120507






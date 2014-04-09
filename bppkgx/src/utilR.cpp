#include "utilR.hpp"

/***********************************************************************************/
/* Evaluating R Expressions                                                        */
/***********************************************************************************/
     
/* Return an R SEXP that points to x --------------------------------------------- */
SEXP mkpointer( double *x, int n )
{
  int  i;
  SEXP ans;
  PROTECT( ans = Rf_allocVector( REALSXP, n ) );
  for( i = 0; i < n; i++ ){ REAL( ans )[i] = x[i]; }
  UNPROTECT( 1 );
  return( ans );
}
// extern "C"{
// /* Return an R SEXP that points to x --------------------------------------------- */
// SEXP matrixpointer( double **x, int nrow, int ncol )
// {
//   int  i, j;
//   SEXP ans;
//   PROTECT( ans = Rf_allocMatrix( REALSXP, nrow, ncol ) );
//     for( i = 0; i < nrow; i++ )for( j = 0; j < ncol; j++ ){ REAL( ans )[i + j*ncol] = x[i][j]; }
//   UNPROTECT( 1 );
//   return( ans );
// }
// }


/* Return an R SEXP that points to i --------------------------------------------- */
SEXP mkint( int i ) 
{
  SEXP ans;
  PROTECT( ans = Rf_allocVector( INTSXP, 1 ) );
  INTEGER( ans )[0] = i + 1;
  UNPROTECT( 1 );
  return( ans );
}




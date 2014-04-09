#ifndef UTILARMADILLO_HEADER
#define UTILARMADILLO_HEADER
#include <RcppArmadillo.h>

//#include "armadillo.h"
using namespace Rcpp;
using namespace arma;

#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <R_ext/Lapack.h>
//#include <math.h>
#include <algorithm>


// distribution functions
// Inverse Gamma
double dinvgamma( double, double, double, int );
double rinvgamma( double, double );

// Inverse Chi-square
double dinvchisq( double, double, double, int );
double rinvchisq( double, double );

double rinvgaussian( double, double );

// beta binomial
double dbetabinom( double, int, double, double, int );
double rbetabinom( double, double, double );
double rand_truncated_normal_rej_below( const double mu, const double sig, const double below );
double rand_truncated_normal_rej_above( const double mu, const double sig, const double above );
double rand_tn( const double mu,const double sig, const double lower,const double upper );

uvec trimat_index( int p, bool upper );

//  Takes a matrix and returns its cholesky
//  Note: the matrix you pass this function is overwritten with the result.
template< class Matrix >
Matrix La_Chol( Matrix X )
{
  int i,j;
  char uplo = 'U';
  int info;
  Matrix cX(X);
  int p = X.n_rows;
  //----Use Lapack-----
  F77_NAME( dpotrf )( &uplo,&p,cX.memptr(), &p, &info);
  //-------------------

  //---- The lapack function above is messy, it doesn't zero out the lower diagonal -----
  for(i = 0; i < p; i++){
      for(j = 0; j < i; j++){
        cX(i,j)=0.0;
      }
  }
  //-------------------------------------------------------------------------------
  return(cX);
}

//uvec sequence( int start, int end );
//template< typename T, template <typename> class ARMA_VECTOR_TYPE >
template< class Vector >
Vector seq( int start, int end )
{
  size_t i;
  size_t len = abs( end - start ) + 1;
  int sng = ( end - start )/abs( end - start );
  Vector  result( len );
  for( i = 0; i < len; i++ ){
   result( i ) = start + i * sng;
  }
  return( result );
}
//imat gen_rand_graph( const int p, const double prob );
/*********************************************************************************/
/* Generate Random Graph                                                         */   
/*********************************************************************************/ 
//template< typename T, template <typename> class ARMA_MATRIX_TYPE >
template< class Matrix >
Matrix gen_rand_graph( const int p, const double prob )
{
  int i, j;
  Matrix A( p, p );
  A.eye();
  for( i = 0; i < ( p - 1 ); i++ ){
    for( j = ( i + 1 ); j < p; j++ ){
      int a = Rf_rbinom( 1, prob );
      A( i, j ) = a;
      A( j, i ) = a;
    }
  }
  return(A);
}  //gen.rand.A


 
/******************************************************************************/ 
/*         multinomial                                                        */
/******************************************************************************/ 
/* returns draw from multinomial with cell prob pr  
  ------------------------------------------   
  Value: 
    x:   vector of indices indicating draws 
         x in [0,..n_cell-1] 
  ------------------------------------------   
  Args: 
    n_draw: number of draws 
    n_cell: number of cells 
    pr:     n_cell vector of cell probs 
            (not necessarily standardized) 
    x:      n_draw vector of indices 
            (OUTPUT) 
 ********************************************/ 
template< typename T, template <typename> class ARMA_VECTOR_TYPE >
ivec multinomial( int n_draw, int n_cell, ARMA_VECTOR_TYPE<T> pr ) 
{ 
  double uj; 
  int i,j; 
  ivec x;
  ARMA_VECTOR_TYPE<T> cum_p = cumsum( pr ); 
  for( j = 0; j < n_draw; j++ ){ 
    uj = runif( 0, 1 ) * cum_p( n_cell - 1 ); 
    for( i = 0; ( (uj > cum_p( i ) ) & ( i < n_cell ) ); i++ ); 
    x( j ) = i; 
  } 
  return( x ); 
} 
 
template< typename T, template <typename> class ARMA_VECTOR_TYPE >
ivec multinomial_table(int n_draw, int n_cell, ARMA_VECTOR_TYPE<T> pr ) 
{ /* same as multinomial, but returns counts instead of indices */ 
  double uj; 
  int i,j; 
  ivec ct=zeros<ivec>(n_cell);
  ARMA_VECTOR_TYPE<T>  cum_p = cumsum( pr ); 
  for( j = 0; j < n_draw; j++ ){ 
    uj = runif( 0, 1 ) * cum_p( n_cell - 1 ); 
    for( i = 0; ( (uj > cum_p(i) ) & ( i < n_cell ) ); i++ ){ 
          ct(i) += 1; 
    }
  } 
  return(ct);
} 

double rand_truncated_normal_rej_twosided(const double mu, const double sig, const double below, const double above);

vec rtgamma( const int n, const double shape, const double scale, const double lower, const double upper );
vec rtnorm( const int n, const double mu, const double sigma, const double lower, const double upper );

template< typename T, template <typename> class ARMA_VECTOR_TYPE >
void vec_append( ARMA_VECTOR_TYPE<T> & target, const T element )
{
    target.resize( target.n_elem + 1 );
    target( target.n_elem - 1 ) = element;
}

// armadillo vector has just one template type parameter
template< typename T, template <typename> class ARMA_VECTOR_TYPE >
ARMA_VECTOR_TYPE<T> vec_intersection( const ARMA_VECTOR_TYPE<T> first, const ARMA_VECTOR_TYPE<T> second )
{
    std::vector<T> output ;
    std::set_intersection( first.begin(), first.end(), second.begin(), second.end(),
                           std::back_inserter(output) ) ;
    std::reverse( output.begin(), output.end() ) ;   
    ARMA_VECTOR_TYPE<T> result = conv_to< ARMA_VECTOR_TYPE<T> >::from(output);
    return result ;
}

// armadillo vector has just one template type parameter
template< typename T, template <typename> class ARMA_VECTOR_TYPE >
ARMA_VECTOR_TYPE<T> vec_union( const ARMA_VECTOR_TYPE<T> first, const ARMA_VECTOR_TYPE<T> second )
{
    std::vector<T> output ;
    std::set_union( first.begin(), first.end(), second.begin(), second.end(),
                           std::back_inserter(output) ) ;
    std::reverse( output.begin(), output.end() ) ;   
    ARMA_VECTOR_TYPE<T> result = conv_to< ARMA_VECTOR_TYPE<T> >::from(output);
    return result ;
}

// armadillo vector has just one template type parameter
// in first but not second
template< typename T, template <typename> class ARMA_VECTOR_TYPE >
ARMA_VECTOR_TYPE<T> vec_difference( const ARMA_VECTOR_TYPE<T> first, const ARMA_VECTOR_TYPE<T> second )
{
    std::vector<T> output ;
    std::set_difference( first.begin(), first.end(), second.begin(), second.end(),
                           std::back_inserter(output) ) ;
    std::reverse( output.begin(), output.end() ) ;   
    ARMA_VECTOR_TYPE<T> result = conv_to< ARMA_VECTOR_TYPE<T> >::from(output);
    return result ;
}

// armadillo vector has just one template type parameter
template< typename T, template <typename> class ARMA_VECTOR_TYPE >
ARMA_VECTOR_TYPE<T> vec_symdifference( const ARMA_VECTOR_TYPE<T> first, const ARMA_VECTOR_TYPE<T> second )
{
    std::vector<T> output ;
    std::set_symmetric_difference( first.begin(), first.end(), 
                                   second.begin(), second.end(),
                                   std::back_inserter(output) ) ;
    std::reverse( output.begin(), output.end() ) ;   
    ARMA_VECTOR_TYPE<T> result = conv_to< ARMA_VECTOR_TYPE<T> >::from(output);
    return result ;
}

template< typename T, template <typename> class ARMA_VECTOR_TYPE >
ARMA_VECTOR_TYPE<T>  vec_merge( const ARMA_VECTOR_TYPE<T> first, const ARMA_VECTOR_TYPE<T> second )
{
    ARMA_VECTOR_TYPE<T> result( first.n_elem + second.n_elem );
    std::copy( first.begin(),  first.end(),  result.begin() );
    std::copy( second.begin(), second.end(), result.begin() + first.n_elem );
    return( result );
}
template< typename T, template <typename> class ARMA_MATRIX_TYPE >
int  joint_diag( const ARMA_MATRIX_TYPE<T> A, 
                 const ARMA_MATRIX_TYPE<T> B,
                       ARMA_MATRIX_TYPE<T>& Da, 
                       ARMA_MATRIX_TYPE<T>& Db,
                       ARMA_MATRIX_TYPE<T>& V
                )
{
  mat Ra  = chol(A);
  mat iRa = inv( trimatu(Ra) );
  vec eigval;
  mat eigvec;
  mat iRaBiRa=iRa.t() * B * iRa;
  eig_sym( eigval, eigvec, iRaBiRa ); 
  V = iRa * eigvec;
  Db.diag = eigval;
  Da = V.t() * A * V;   
  return 0;  
}
template< typename T, template <typename> class ARMA_FIELD_TYPE, template <typename> class ARMA_MATRIX_TYPE >
int  approx_joint_diag( const ARMA_FIELD_TYPE< ARMA_MATRIX_TYPE<T> > A, 
                                      ARMA_FIELD_TYPE< ARMA_MATRIX_TYPE<T> >& D,
                                      ARMA_MATRIX_TYPE<T>& V )
{
  double eps = 1e-06;
  int maxiter = 100;
  
  int k = A.n_elems;   //k=kp/p
  ARMA_MATRIX_TYPE<T> X = A( 0 );
  for( int i = 1; i < k; i++ ){
      X = join_cols( X, A( i ) );   
  }
  int p = X.n_cols;
  int kp= X.n_rows;

  V = ARMA_MATRIX_TYPE<T> ( p, p );
  V.eye(); 
  uvec iIdx( 1 );
  uvec jIdx( 1 );
  uvec Ii(k);
  uvec Ij(k);
  int encore = 1;
  int iter = 0;
  while( encore == 1 ){
    iter = iter + 1;
    encore = 0;
    for( int i = 0; i < ( p - 1 ); i++ ){
      for( int j = ( i + 1 ); j < p; j++ ){
        //Ii<-seq(i,kp,p) //Ii<-seq(i,kp,p)
        //Ij<-seq(j,kp,p) //Ij<-seq(j,kp,p)
        for( int h = 0; h < k; h++ ){
            Ii(h) = i + h*p;
            Ij(h) = j + h*p;
        }
        iIdx(0)=i;
        jIdx(0)=j;
        ARMA_MATRIX_TYPE<T> g1 = X(Ii,iIdx)-X(Ij,jIdx); //g1<-X[Ii,i]-X[Ij,j]
        ARMA_MATRIX_TYPE<T> g2 = X(Ij,iIdx)+X(Ii,jIdx); //g2<-X[Ij,i]+X[Ii,j]

        ARMA_MATRIX_TYPE<T> g  = join_rows( g1, g2 );  //g<-cbind(g1,g2)
        ARMA_MATRIX_TYPE<T> gg = g.t()* g;   //gg<-t(g)%*%g
        T ton = gg(1,1)-gg(2,2);
        T toff= gg(1,2)+gg(2,1);
        T theta = 0.5 * atan2( toff, ton + sqrt( ton * ton + toff * toff ) );
        
        T cos_theta = cos( theta );
        T sin_theta = sin( theta );
        
        if( abs( sin_theta ) > eps ){
          encore = 1;
          
          ARMA_MATRIX_TYPE<T> Mi = X(Ii, span::all); // Mi <- X[Ii,]
          ARMA_MATRIX_TYPE<T> Mj = X(Ij, span::all); // Mj <- X[Ij,]
          
          X(Ii, span::all)=cos_theta * Mi + sin_theta * Mj; //X[Ii,]<- cos.theta * Mi + sin.theta * Mj
          X(Ij, span::all)=cos_theta * Mj - sin_theta * Mi; //X[Ij,]<- cos.theta * Mj - sin.theta * Mi
          
          colvec col_i = X.col(i); // col.i <- X[,i]
          colvec col_j = X.col(i); //X[,j]
          
          X.col(i) = cos_theta * col_i + sin_theta * col_j;
          X.col(j) = cos_theta * col_j - sin_theta * col_i;
          
          rowvec temp = V.row(i);
          V.row(i) = cos_theta * V.row(i) + sin_theta * V.row(j);
          V.row(j) = cos_theta * V.row(j) - sin_theta * temp;
        }
      }
    }
  if( iter >= maxiter ) Rf_error("maxiter reached without convergence");    
  }
  D = ARMA_FIELD_TYPE< ARMA_MATRIX_TYPE<T> >( k );
  for( int i = 0; i < k; i++ ){
      D(i) = V * A( i )* V.t(); 
  }
  V = V.t();
  return 0;
}


template< class Vector >
Vector SampleNoReplace( int k, int n )
{
  int i, j;
  Vector y( k ), x( n );
  for ( i = 0; i < n; i++ )
      x( i ) = i;
  for ( i = 0; i < k; i++ ) {
      j = conv_to< int >::from( (n * randu(1))); //n * randu(1);
      y( i ) = x( j ) + 1;
      x( j ) = x( --n );
  }
  return( y );
}

template< class Vector >
int ProbSample1( Vector& p )
{
  int accepted = 0, n = p.size();
  vec result(1);
  //std::cout << "p: " <<p<< std::endl;
  while( accepted == 0 ){
    result = SampleNoReplace<vec>( 1, n );
    //std::cout << "result: " << result<< std::endl;
    if( randu(1) < n * p( result( 0 ) - 1 ) ){ accepted = 1; }
  }
  return( result( 0 ) );
}

//------------------------------------------------------------------------------
// Check cycle in DAG
template < class Matrix >
int checkCycle( Matrix& A, int numEdge )
{

  Matrix Al(A);
  int l, cycle = 0; 
  int p = A.n_rows;
  int maxIter = Rf_imin2( p, numEdge );
  for( l = 0; l < ( maxIter - 1 ); l++ ){
    Al = Al * A;
    int traA = trace( Al );
    if( traA != 0 ) {
      cycle = 1; 
      break;
    }
  }
  return( cycle );
}


template <class Matrix>
void qr_solve( Matrix x, Matrix y, Matrix coef)
/* Translation of the R function qr.solve into pure C
   NB We have to transpose the matrices since the ordering of an array is different in Fortran
   NB2 We have to copy x to avoid it being overwritten.
*/
{
    int i, info = 0, rank, n, p;
    int n1 = x.n_rows;
    double tol = 1.0E-7;
    n = n1;
    p = n1;
    double pivot[n1];
    for( i = 0; i < n1; i++ ){
      pivot[i] = i+1;
    }
    double qraux[n1];
    double work[2*n1];
    F77_CALL(dqrdc2)( x.memptr(), &n, &n, &p, &tol, &rank, qraux, pivot, work );

    if (rank != p){
      //error("Singular matrix in qr_solve\n");
    }
    F77_CALL(dqrcf)(x.memptr(), &n, &rank, qraux, y.memptr(), &n, coef.memptr(), &info);

}


/**************************************************************************************/
/* Function: DMVNORM()                                                                */
/**************************************************************************************/
template <class Matrix, class Vector>
double DMVNORM( const Vector& x, const Vector& mu, const Matrix& Sigma, int prec, int give_log )
{
  arma_debug_check( x.n_elem != mu.n_elem, "x and mu must have the same size");
  arma_debug_check( (Sigma.is_square() == false), "sigma is not square" );
  arma_debug_check( Sigma.n_rows != mu.n_elem, "sigma and mu are not compatible");
  double log_density = 0;
  int d = x.n_elem;
  double Const= - M_LN_SQRT_2PI * d; //M_LN_SQRT_2PI
  Vector xc = x - mu;
  double val;
  double sign;
  if( prec ) { // if precision matrix 
    log_det( val, sign, Sigma );
    log_density = Const + 0.5 * ( val + log( sign ) ) - as_scalar( xc.t() * Sigma * xc ) / 2;
  } else {     // if covariance  matrix
    log_det(val, sign, Sigma);
    log_density = Const - 0.5 * ( val + log( sign ) ) - as_scalar( xc.t() * inv( Sigma ) * xc ) / 2;
  }
  return( give_log ? log_density : R_pow( M_E, log_density ) );
}
// template <class Matrix, class Vector>
// double DMVNORM( const Vector& x, const Vector& mu, const Matrix& sigma, int prec, int give_log )
// {
//   arma_debug_check( x.n_elem != mu.n_elem, "x and mu must have the same size");
//   arma_debug_check( (sigma.is_square() == false), "sigma is not square" );
//   arma_debug_check( sigma.n_rows != mu.n_elem, "sigma and mu are not compatible");
//   double log_density = 0;
//   int d = x.n_elem;
//   double Const= - M_LN_SQRT_2PI * d; //M_LN_SQRT_2PI
//   Vector xc = x - mu;
//   Matrix inv_sigma = inv( sigma );
//   double val;
//   double sign;
//   if( prec ) { // if precision matrix 
//     log_det(val, sign, inv_sigma);
//     log_density = Const - ( val + sign ) / 2 - ( xc * sigma * xc ) / 2;
//   } else {     // if covariance  matrix
//     log_det(val, sign, sigma);
//     log_density = Const - ( val + sign ) / 2 - ( xc * inv_sigma * xc ) / 2;
//   }
//   return( give_log ? log_density : R_pow( M_E, log_density ) );
// }

/**************************************************************************************/
/* Function: MVNORM()                                                                 */
/**************************************************************************************/
template <class Matrix, class Vector>
Vector MVNORM( int C, const Vector& Mvec, const Matrix& S )
{
  // Dimension check
  arma_debug_check( S.n_rows != Mvec.n_elem, "sigma and mu are not compatible");
  arma_debug_check( S.n_cols != Mvec.n_elem, "sigma and mu are not compatible");
  // Local Memory ---------------------------------------------------------------
  int D = Mvec.n_elem;
  // Generate random standard normals ------------------
  Vector z = randn(D);
  Vector Y;
  // If s is a concentration matrix ---------------------------------------------
  if( C == 1 ){
    // Cholesky UtU = S -> returns Upper triangula matrix U
    mat U = chol(S); 
    // Add variance then add mean-------------------------
    Y = solve( trimatu(U), z ) + Mvec;
  }
  // If s is a covariance matrix ------------------------------------------------   
  if( C == 0 ){
    // Cholesky decomp --------------------------------------------------------
    mat U = chol(S); 
    // Y = Lz + mean
    Y = trans(U) * z + Mvec;
  }
  return( Y );
} // END MVNORM -------------------------------------------------------------------

template <class Matrix, class Vector>
Vector MVNORM_SVD( int C, const Vector& Mvec, const Matrix& S )
{
  // Dimension check
  arma_debug_check( S.n_rows != Mvec.n_elem, "sigma and mu are not compatible");
  arma_debug_check( S.n_cols != Mvec.n_elem, "sigma and mu are not compatible");
  // Generate random standard normals ------------------
  Vector z = randn(Mvec.n_elem);
  Matrix U;
  Vector s;
  Matrix V;
  svd(U,s,V,S);
  Vector Y = z * U * diagmat( sqrt( s ) ) * trans( V ) + Mvec;
  return( Y );
} // END MVNORM_SVD -------------------------------------------------------------------


/*********************************************************************************/
/* Matrix Normal using Cholesky:                    M_norm()                     */
/*********************************************************************************/
template <class Matrix>
Matrix M_norm( Matrix mu, Matrix Vrow, Matrix Vcol )
{
   int nrow = mu.n_rows;
   int ncol = mu.n_cols;
   // Gaussian Simulations ---------------------------------------------------
  Matrix Z =randn(nrow, ncol);
   // Cholesky Decompositions ------------------------------------------------
  mat VR = chol(Vrow); 
  mat VC = chol(Vcol); 
  Matrix sample = trans(VR) * Z * trans(VC)+ mu;
  return( sample );
}// END M_nrom -----------------------------------------------------------------
/*********************************************************************************/
/* Matrix Normal without Cholesky:                  M_C_norm()                   */
/* CVrow is chol( Vrow ) and CVcol is chol( Vcol )                               */
/*********************************************************************************/
template <class Matrix>
Matrix M_C_norm( Matrix mu, Matrix CVrow, Matrix CVcol )
{
   int nrow = mu.n_rows;
   int ncol = mu.n_cols;
   // Gaussian Simulations ---------------------------------------------------
  Matrix Z =randn(nrow, ncol);
  Matrix sample = trans(CVrow) * Z * trans(CVcol)+ mu;
  return( sample );
}// END M_nrom -----------------------------------------------------------------


/*********************************************************************************/
/* Multivariate Banded Gaussian Process            r_gauss()                     */   
/*********************************************************************************/ 
template <class Matrix, class Vector>
Vector triang_system(Matrix L, Vector x, int n, int lower)
{
   int i,j;
   double sum;
   Vector solution;
   if(lower == 1){ // Lower (forward) ---------------------
    solution(0) = x(0)/L(0,0);
    for(i=1; i<n; i++){
      sum = 0.0;
      for(j=0; j<i; j++)  sum = sum + L(i,j)*solution(j);
      solution(i) = (x(i) - sum)/L(i,i);     
    }
   } 
   if(lower == 0){ // Upper (backward) --------------------
    solution((n-1)) = x((n-1))/L((n-1),(n-1));
    for(i=(n-2); i>=0; i--) {
      sum =0.0; 
      for(j=i+1;j<=n;j++) sum = sum + L(i,j)*solution(j);
      solution(i) = (x(i) - sum)/L(i,i);
    }
   }
   return(solution);
}

template <class Matrix, class Vector>
Vector r_gauss(Matrix P, Vector XY, int n){
   int i,j;
   // Cholesky Decomposition -------------------------------- 
   vec z = randn(n);// z  Norm(0,I) 
   //choldc(P, n, p);                           // P  Lower
   mat U = chol(P);                             // U  Upper
   //transpose_matrix(P, Pt, n, n);             // Pt Upper
   // Add Covariance Structure ------------------------------
   // Solve Pt B = z  (backward)
   //Vector sample = triang_system( U, z, n, 0);
   Vector sample = solve( trimatu(U), z );
   // Calculate Posterior Mean ------------------------------
   // Solve P v = XY  (forward)
   //Vector v = triang_system( P, XY, n, 1 );
   Vector v = solve( trimatl(trans(U)), XY );
   // Solve Pt m = v  (backward)
   //Vector m = triang_system( U, v, n, 0 );
   Vector m = solve( trimatu(U), v );
   // Add mean ----------------------------------------------
   sample = sample + m;
   return( sample );
}
/*********************************************************************************/
/* Wishart Distribution                              rwish()                     */   
/*********************************************************************************/ 
// template <class Matrix>
// Matrix rwish( int v, Matrix S) {
//   //if ( S.n_rows <> S.n_cols ) {
//   //  Rf_error(message="S not square in rwish().\n")
//   //}
//   // if ( v < S.n_rows ) {
//   //   Rf_error("v is less than the dimension of S in rwish().\n")
//   // }
//   int i, j;
//   int p     = S.n_rows;
//   //Matrix CC = chol( S );
//   Matrix CC = La_Chol(S); 
//   Matrix Z  = zeros<Matrix>( p, p );
//   for( i = 0; i < p; i++ ){
//     Z( i, i ) = sqrt( Rf_rchisq( v - i ) );
//   }
//   for( i = 0; i < p; i++ ){ 
//     for( j = i + 1; j < p; j++ )
//     //Z( i, span( i + 1, p - 1 ) ) = mat( rnorm( p - i - 1, 0, 1), 1, p - i );
//       Z( i, j ) = Rf_rnorm( 0, 1 ); //mat( rnorm( p - i - 1, 0, 1), 1, p - i );
//   } 
//   Matrix ZCC = Z * CC;
//   return( ZCC.t() * ZCC );
// }

//  Standard gaussian stuff.  Takes a Cov matrix Full and returns the conditional
//  Cov matrix for a subset clique_ID.
template< typename T, template <typename> class ARMA_MATRIX_TYPE >
ARMA_MATRIX_TYPE<T> get_cond_matrix( const int p_cl, 
                                      const uvec clique_ID, 
                                      const uvec V_ID, 
                                      const ARMA_MATRIX_TYPE<T> Full)
{
  //int i, j;
  //int p_V = Full.n_rows - p_cl;
  //ARMA_MATRIX_TYPE<T1> cV(p_cl, p_V);
 // ARMA_MATRIX_TYPE<T1> VV(p_V,p_V);

  //----  Start: Form Matrix Objects --------
  //-----  Make the two cross parts  ------
  //  printf("Here Full\n");
  //  Full->print();
  // for(i = 0; i < p_cl; i++){
  //     for(j = 0;j < p_V; j++){
  //       cV(i,j) = Full(clique_ID(i),V_ID(j) );
  //     }
  // }
  ARMA_MATRIX_TYPE<T> cV = Full( clique_ID, V_ID );
  //----------------------------------------------
  
  //---Now the inverse of the remaining elements----
  // for(i = 0; i < p_V; i++){
  //     for(j = 0; j < p_V; j++){
  //       VV(i,j)= Full( V_ID(i),V_ID(j) );
  //     }
  // }

  ARMA_MATRIX_TYPE<T> VV = Full( V_ID, V_ID );
  //-------------------------------------------------
  ARMA_MATRIX_TYPE<T> Tinv = cholbacksolve( VV );
  ARMA_MATRIX_TYPE<T> cV_Tinv = cV * Tinv;
  ARMA_MATRIX_TYPE<T> Result  = cV_Tinv * cV_Tinv.t();
  //---------------------------------------------
  return( Result );
}//get_cond_matrix

//Cholesky and backsolve
template< typename T, template <typename> class ARMA_MATRIX_TYPE >
ARMA_MATRIX_TYPE<T> cholbacksolve( ARMA_MATRIX_TYPE<T> M )
{
  ARMA_MATRIX_TYPE<T> Mc = chol(M);
  int p = Mc.n_rows;
  char side = 'L';
  char up = 'U';
  char transa = 'N';
  char diag = 'N';
  double alpha = 1.0;
  ARMA_MATRIX_TYPE<T> B( p, p );
  B.eye();
  F77_NAME( dtrsm )( &side,&up,&transa,&diag,&p,&p,&alpha,Mc.memptr(),&p,B.memptr(),&p );
  return( B );
}
// //-------- Samples a Wishart Variate (Implies Full Graph) ---------
// // Note that I use the Barlett decomposition to achieve this
template< typename T, template <typename> class ARMA_MATRIX_TYPE >
ARMA_MATRIX_TYPE<T> rwish( int delta, ARMA_MATRIX_TYPE<T> D)
{
  int i, j;
  int p = D.n_rows;

  ARMA_MATRIX_TYPE<T> D_inv = inv( D );
  ARMA_MATRIX_TYPE<T> Tmat = La_Chol( D_inv );
  ARMA_MATRIX_TYPE<T> psi( p, p );

  // ***  Sample Values in Psi   ******
  //------Diagaonal Elements----------
  for (i = 0; i < p; i++) 
    psi( i, i ) = sqrt( Rf_rchisq( p - i - 1 + delta ) );
  //----------------------------------
  
  //------Offdiagonal in G------------
  for( i = 0; i < p-1; i++ ){
    for( j = i + 1; j < p; j++ ){
      psi( i, j ) = Rf_rnorm( 0, 1 );
    }
  }
  //----------------------------------
  
  //******** End Sampling  *************

  //----  Now Complete things  ----------
  ARMA_MATRIX_TYPE<T> psiT = psi * Tmat;
  ARMA_MATRIX_TYPE<T> K    = psiT.t()*psiT;
  //-------------------------------------
  return(K);
}
/*********************************************************************************/
/* Moralize Ancestral Matrix:                    Moralize(...)                   */
/* Input : Parent list                                                           */
/* Output: Returns an adjacency matrix with entries:                             */
/*  -1,1,2 for negative/positive/unsigned dependence and                         */
/*   3 for unsigned moral edge                                                   */
/*********************************************************************************/
template <class Matrix, class Vector>
void moralize(Vector npg, Matrix Anc, Matrix Adj, int G) 
{
  int i, ii, jj, j, l;

  Matrix M = zeros<Matrix>(G,G);
  // Construct Ancestral Matrix -----------------------------------------------  
  for(i=0; i<G; i++){ 
    if(npg(i) > 0){
      for(j=0; j<npg(i); j++){
        jj = abs(Anc(i,j));
        if(Anc(i,j) < 0 ){
          M(i,jj) = -1;
        } else { 
          M(i,jj) = 1; 
        }
        Adj(i,jj) = Adj(jj,i) = M(i,jj);
      }//(j)
    }//(if npg[i] > 0)
    Adj(i,i) = M(i,i) = 0;                    // Make sure diag(M) is 0
  }//(i)

  // Add unsigned edges -----------------------------------------------
  for(i=0; i<(G-1); i++){
    for(j=(i+1); j<G; j++){
      if(Adj(i,j) < 2){                          // if not already recorded
        if( M(j,i)*M(i,j) < 0) Adj(j,i) = Adj(i,j) = 2;
      }
    } // j 
  } // i  
  // Marry the parents -----------------------------------------------------
  for(i=0; i<G; i++){
    if(npg(i) > 1){
      for(j=0; j<(npg(i)-1); j++){
        for(l=(j+1); l<npg(i); l++){
          ii = abs(Anc(i,j)); // Anc j 
          jj = abs(Anc(i,l)); // Anc j + 1
          if((abs(Adj(jj,ii)) < 1) & (ii != jj)) Adj(jj,ii) = Adj(ii,jj) = 3;
        } // l 
      } // j
    }// (npg[i] > 1)
    Adj(i,i) = 0;                          // Make sure diag(Adj) is 0
  } // i
}  

#endif // UTILARMADILLO_HEADER
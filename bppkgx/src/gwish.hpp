#ifndef _GWISH_HPP
#define _GWISH_HPP
//--------------------------------------------------------
// gwish.Rcpp:  This is a collection of functions that manipulate graph objects
//   according to the G-Wishart distribution.
//   Note: this depends on the newgraph.Rcpp library, and the matrix.Rcpp library
//
//-----------------------------------------------------------

//Alex Lenkoski  alex.lenkoski@uni-heidelberg.de
#include <RcppArmadillo.h>

//#include "armadillo.h"
using namespace Rcpp;
using namespace arma;

#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include "utilarmadillo.hpp"
#include "utilR.hpp"
#include "cliques.hpp"
#include "SymChol.hpp"

//  Returns the complement of a set clique_ID and stores it in V_ID
//template< class Vector >
//uvec get_complementary_set( int p, uvec clique_ID );

//const double PI = std::atan(1.0)*4;

//#include "cliques.Rcpp"

// void gwish_dmh_update_all( mat K, imat G, mat delta, int n, 
//                            mat D_prior, mat D_post,
//                            mat& K_result, imat& G_result );
template< class Matrix, class iMatrix >
void gwish_dmh_update_all( Matrix K, iMatrix G, int delta, int n, 
                           Matrix D_prior, Matrix D_post,
                           Matrix& K_result, iMatrix& G_result
                           )
{
  int p  = K.n_rows;// <- dim(K)[1]
  uvec e( 2 );
  uvec ind = seq<uvec>( 0, p-1 );
  // Rprintf("p=%d\n",p);
  // Rprintf("delta+n=%d\n",delta + n);
  // G.print("G=\n");
  // D_post.print("D_post=\n");
  // K.print("K=\n");

  Matrix K_prop = gwish_blgibbs_iterate_new_R( p, G, delta + n, D_post, K );
  // K_prop.print("K_prop=\n");
  //##  print("Updating Edges")
  for( int i =0; i < ( p - 1 ); i++ ){
      for( int j = ( i+ 1 ); j < p; j++ ){
          e( 0 )   = i;
          e( 1 )   = j;
          uvec less_e = sort( vec_difference( ind, e ) ); //ind[-e]
          iMatrix G_prop( G );
          G_prop( e(0), e(1) ) = 1 - G_prop( e(0), e(1) );
          G_prop( e(1), e(0) ) = 1 - G_prop( e(1), e(0) );
          
          //## 2. permute
          uvec tt( p );
          tt.subvec(  0, p-3  ) = less_e;
          tt.subvec( p-2, p-1 ) = e;
          //tt.print("tt=\n");
          G_prop     = G_prop( tt, tt );
          Matrix D_post_new = D_post( tt, tt );
          K_prop     = K_prop( tt, tt );
          Matrix phi;
          // bool state = chol( phi, K_prop );
          // if(state==false) {
          //   Rf_error("Choleski error");
          // }
          phi = La_Chol( K_prop );
          //Rprintf("state=%d\n",state);
          //phi.print("phi=\n");
          //## 3. accept/reject update
          //D_post_new.print("D_post_new=\n");
          double w = gwish_dmh_logH( D_post_new, phi );
            //phi.print("phi=\n");
            //D_post_new.print("D_post_new=\n");

          // Rprintf("w=%f \n",w);
          int include_edge;
          w = 1 / ( exp( w ) + 1 );
          // Rprintf(",w=%f ",w);
          if( Rf_runif( 0, 1 ) < w ){
            include_edge = 1;
          }else{
            include_edge = 0;
          }
          //Rprintf("include_edge=%d ",include_edge);
          //Rprintf("G_prop( p-2, p-1 )=%d ",G_prop( p-2, p-1 ));

          if( include_edge == G_prop( p-2, p-1 ) ){
            //Matrix phi_aux = gwish_dmh_R( p, G_prop, phi, delta ); 
            //for(i in 1:1) phi_aux = gwish_dmh_R( p, G_prop, phi_aux, delta );
            Matrix K_aux = gwish_blgibbs_iterate_new_R( p, G_prop, delta, D_prior, K_prop );
            Matrix phi_aux;
            // state = chol( phi_aux, K_aux );
            // if(state==false) {
            //   Rf_error("Choleski error");
            // }
            phi_aux = La_Chol(K_aux);
            double w2 =  ( 2 * G_prop( p-2 ,p-1 ) - 1 ) * gwish_dmh_logR( phi_aux );
            //phi_aux.print("phi_aux=\n");
            // Rprintf("w2=%f \n",w2);

            if( !( log( Rf_runif( 0, 1 ) ) < w2 ) ){
              //Rprintf(" flip1\n");
              G_prop( p-2, p-1 ) = 1 - G_prop( p-2, p-1 );
              G_prop( p-1, p-2 ) = 1 - G_prop( p-1, p-2 );
            }
          } else {
             //Rprintf(" flip2\n");
            G_prop( p-2, p-1 ) = include_edge;
            G_prop( p-1, p-2 ) = include_edge;
          }

          //##-------- Resample -----------------------------------------
          phi( p-1, p-1 ) = sqrt( Rf_rgamma( ( delta + n ) / 2,  2 / D_post_new( p-1, p-1 ) )   );
          if( G_prop( p - 2, p - 1 ) == 1 ){
              double mu1 = - phi( p-2, p-2 ) * D_post_new( p-2, p-1 ) / D_post_new( p-1, p-1 );
              phi( p - 2, p - 1 ) = Rf_rnorm( mu1, sqrt( 1/D_post_new( p-1, p-1 ) ) );
            }else{
              phi( p - 2, p - 1 ) = - 1/phi( p - 2, p - 2 ) * sum( phi(span( 0, p-3 ), p-2 ) % phi(span( 0, p-3 ), p-1 ) );
            }
          //##---------------------------------------------------------
          Matrix K_tilde  = phi.t() * phi;
          G(tt,tt)     = G_prop;
          K_prop(tt,tt)= K_tilde;
        }
    }
  G_result = G;
  //   Rprintf("p=%d\n",p);
  // Rprintf("delta+n=%d\n",delta + n);
  // G.print("G=\n");
  // D_post.print("D_post=\n");
  // K_prop.print("K_prop=\n");
  K_result = gwish_blgibbs_iterate_new_R( p, G, delta + n, D_post, K_prop );
  // K_result.print("K_result=\n");
  //return(state);
}
//mat gwish_dmh_R( int p, imat G, mat Phi, int delta);
template< class Matrix, class iMatrix >
Matrix gwish_dmh_R( int p, iMatrix G, Matrix Phi, int delta )
{
    SymChol S = new SymChol_CLASS( G );
    Matrix Psi( Phi );
    Matrix Phi_R( Psi );
    int i,j;
    vec nu( p );
    double temp1, temp2;
    int t_i, t_j;
    double alpha = 0;

    for( i = 0; i < ( p - 1 ); i++ ){
      nu( i ) = 0;
      for( j = ( i + 1 ); j < p; j++ ){
        nu( i ) += G( i, j );
      }
    }
    nu( p - 1 ) = 0;
    for( i = 0; i < p; i++ ){ Psi( i, i ) = sqrt( Rf_rchisq( delta + nu( i ) ) );}
    for( i = 0; i < ( p - 1 ); i++ ){
      for( j = ( i + 1 ); j < p; j++ ){
        if( G( i, j ) ){
          Psi( i, j ) = Rf_rnorm( 0, 1 );
        }
      }
    }
    S->CompletePhi(Psi);
    for( i = 0; i < S->n_diff; i++ ){
      t_i = S->diff[i]->d_i;
      t_j = S->diff[i]->d_j;
      temp1 = Psi( t_i, t_j );
      temp2 = Phi( t_i, t_j );
      alpha += -0.5 * ( temp1 * temp1 - temp2 * temp2 );
      if( alpha < -10 ) break;//trust me, it won't recover.
    }
    //    printf("alpha: %f\n",alpha);
    if( log( Rf_runif( 0, 1 ) ) < alpha ){
      Phi_R = Psi;
    }
    return( Phi_R );
}


template< class Matrix >
double gwish_dmh_logH ( Matrix D_post, Matrix phi ){
    int p         = phi.n_cols; //p      <- dim(phi)[2]
    double my_PI  =3.141592653589793238462;
    double phi_0  = -1/phi( p-2, p-2 ) * accu( phi( span( 0, p-3 ), p-2 ) % phi( span(0, p - 3 ), p-1 ) );
    double A      = phi(p-2,p-2) * D_post(p-2,p-1) / D_post(p-1,p-1);
    double result = -0.5 * ( D_post( p-1, p-1 ) * pow( phi_0 + A, 2 ));
    result = result - log( phi( p-2, p-2 ) );
    result = result + 0.5 * log( D_post( p-1, p-1 ) ) - 0.5 * log( 2 * my_PI );
    return(result);
}
template< class Matrix >
double gwish_dmh_logR( const Matrix phi )
{
    double result=0;
    double my_PI  =3.141592653589793238462;
    int p  = phi.n_cols;
    double phi_0  = -1/phi( p - 2, p - 2 ) * accu( phi( span(0, (p-3) ), p-2 ) % phi( span(0, (p-3) ), p-1 ));
    result = -0.5 * pow( phi_0, 2 );
    result = result - log( phi( p - 2, p - 2 ) );
    result = result - 0.5 * log( 2 * my_PI );
    return( result );
}


template< class Matrix >
Matrix gwish_blgibbs_iterate_new_R( int p, 
                                     imat A,
                                     int delta,
                                     Matrix D,
                                     Matrix K
                                    ) //(int *p_R, int* A_R, int* delta, double* D_R, double* K_R )
{
    A.diag() = zeros<ivec>( p );
    //int ee = p * (p - 1) / 2;
    imat cl( p * p, p );
    uvec cl_dims = zeros<uvec>( p * p );
    int n_cl = FindCliques( cl, cl_dims, A );
    Matrix K_update = gwish_blgibbs_iterate( cl, cl_dims, n_cl, delta, D, K );
    return K_update;
}

//  Returns the complement of a set clique_ID and stores it in V_ID
template< class Vector >
Vector get_complementary_set( size_t p, Vector clique_ID )
{
  size_t i,j,k;
  int temp_in_clique;

  k = 0;
  size_t p_clique =  clique_ID.n_elem;
  Vector V_ID( p - p_clique );
  for( i = 0; i < p; i++ ){
    temp_in_clique = 0;
    for( j = 0; j < p_clique; j++ ){
      if( i == clique_ID( j ) ){
        temp_in_clique = 1;
      }
    }
    if( temp_in_clique == 0 ) {
      V_ID( k ) = i;
      k++;
    }
  }
  return( V_ID );
}

//Performs one iteration of the Block Gibbs Sampler

template< typename Tmd, template <typename> class ARMA_MATRIX_TYPE,
          typename Tvi, template <typename> class ARMA_VECTOR_TYPE>
ARMA_MATRIX_TYPE<Tmd> gwish_blgibbs_iterate( imat cl,  
                                             ARMA_VECTOR_TYPE<Tvi> cl_dims, 
                                             int n_cl, 
                                             int delta, 
                                             ARMA_MATRIX_TYPE<Tmd> D, 
                                             ARMA_MATRIX_TYPE<Tmd> K )
{
  size_t p = D.n_rows;
  size_t i, p_clique;
  uvec clique_ID, V_ID;
  int clique_num;
  //uvec clique_ID;
  ARMA_MATRIX_TYPE<Tmd> submatrix_D;
  ARMA_MATRIX_TYPE<Tmd> submatrix_cc;
  ARMA_MATRIX_TYPE<Tmd> submatrix_cond;

  //------------ Loop through cliques and update accordingly -------
  //  K->print();
  //  util_print_vec_int(cl_dims, p * (p - 1)/2);
  for( clique_num = 0; clique_num < n_cl; clique_num++ ){

      //*******Just recording Indices***********
      p_clique = cl_dims( clique_num );
      //      printf("p: %d, pcl: %d, clnum: %d\n",p,p_clique,clique_num);
      //clique_ID = new int[p_clique];
      uvec clique_ID( p_clique );
      for( i = 0; i < p_clique; i++){ 
        clique_ID( i ) = cl( clique_num, i ); 
      }
      //      util_print_vec_int(clique_ID, p_clique);
      //      util_print_vec_int(V_ID, p - p_clique);
      //******************************************
      
      //******  Start: Making Submatrices  ************
      
      //---First make the part specific to the clique-----
      submatrix_D  = D( clique_ID, clique_ID ); //D->sub_mat(p_clique, clique_ID);
      submatrix_cc = rwish( delta, submatrix_D );
      //--------------------------------------------------
      
      //-----   Now form the part that's based ------
      //----- on the entries outside of the clique --
      if( p_clique < p & p_clique > 1 ){
        //V_ID = new int[p - p_clique];
        V_ID = get_complementary_set<uvec>( p, clique_ID );
        submatrix_cond = get_cond_matrix( p_clique, clique_ID, V_ID, K );
        submatrix_cc = submatrix_cc + submatrix_cond;
        submatrix_cond.reset();
      }
      //---------------------------------------------
      
      //********  End Forming Matrices ********************
      
      //---Update the result matrix for the clique---
      //K->set_sub_mat(p_clique, clique_ID, submatrix_cc);
      K( clique_ID, clique_ID ) = submatrix_cc;
      //---------------------------------------------

      //----- Clean up from updating this clique-----
      clique_ID.reset();
      if(p_clique < p & p_clique > 1){
        V_ID.reset();
      }
      submatrix_cc.reset(); 
      submatrix_D.reset();
      //----------------------------------------------

    }//End Loop through cliques
  //--We've now looped through all cliques in the graph and updated---

  return(K);  
}


#endif //_GWISH_HPP

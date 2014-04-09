#ifndef _CLIQUES_HPP
#define _CLIQUES_HPP

template< typename Tm, template <typename> class ARMA_MATRIX_TYPE, 
          typename Tv, template <typename> class ARMA_VECTOR_TYPE >
int FindCliques( ARMA_MATRIX_TYPE<Tm>& cliques, 
                 ARMA_VECTOR_TYPE<Tv>& cliquesDimens, 
                 const ARMA_MATRIX_TYPE<Tm> A )
{
  int p     = A.n_rows;
  ARMA_VECTOR_TYPE<Tv> Q( p );    //int* Q = new int[p];
  Q.zeros();
  ARMA_VECTOR_TYPE<Tv> SUBG( p ); //int* SUBG = new int[p];
  ARMA_VECTOR_TYPE<Tv> CAND( p ); //int* CAND = new int[p];
  SUBG.ones();
  CAND.ones();

  int ncliques = 0; 
  findcliquesExpand( Q, SUBG, CAND,
                      cliques, cliquesDimens, ncliques, A );
  return( ncliques );
}

template< typename Tm, template <typename> class ARMA_MATRIX_TYPE, 
          typename Tv, template <typename> class ARMA_VECTOR_TYPE >
void findcliquesExpand( ARMA_VECTOR_TYPE<Tv> Q, 
                        ARMA_VECTOR_TYPE<Tv> SUBG, 
                        ARMA_VECTOR_TYPE<Tv> CAND,
                        ARMA_MATRIX_TYPE<Tm>& cliques, 
                        ARMA_VECTOR_TYPE<Tv>& cliquesDimens, 
                        int& ncliques, 
                        const ARMA_MATRIX_TYPE<Tm>  A)
{
  int i,j;
  int p = A.n_rows;
  if( findcliquesIsEmpty( SUBG ) ){
    //found a new clique
    j = 0;
    for(i = 0; i < p; i++){
      if( Q(i) ){
        cliques( ncliques, j ) = i;
        j++;
      }
    }
    cliquesDimens( ncliques ) = j;
    ncliques++;
  } else {
    //find u
    int u;
    for( u = 0; u < p; u++ ){
      //   printf("u = %d\n",u);
      if( SUBG( u ) ){
        break;
      }
    }
    if( u == p ){
      Rprintf( "reach u = graph->nVertices\n" );
      exit( 1 );
    }
    
    int maxintersect = 0;
    for( i = 0; i < p; i++ ){
      maxintersect += CAND( i ) * A( u, i );
    }
    for( j = 0; j < p; j++ ){
      if( SUBG( j ) ){ 
        int intersect = 0;
        for( i = 0; i < p; i++ ){
          intersect += CAND( i ) * A( j, i );
        }
        if( intersect > maxintersect ){
          u = j;
          maxintersect = intersect;
        }
      }
    }
    // printf("[u = %d]\n",u+1);
    
    int q = -1;

    while( -1 != ( q = findcliquesGETQ( u, CAND, A ) ) ){
      //      printf("%d,",q+1);
      Q( q ) = 1;
      ARMA_VECTOR_TYPE<Tv> qSUBG(p); //int* qSUBG = new int[p];
      ARMA_VECTOR_TYPE<Tv> qCAND(p); //int* qCAND = new int[p];
      for( i = 0; i < p; i++ ){
          qSUBG( i ) = SUBG( i ) * A( q, i );
          qCAND( i ) = CAND( i ) * A( q, i );
      }
      findcliquesExpand( Q, qSUBG, qCAND,
                         cliques, cliquesDimens, ncliques, A );
      qCAND.reset();
      qSUBG.reset();
      CAND( q ) = 0;
      Q( q ) = 0;
      // printf("back,");
    }
  }
  return;
}

template< typename Tv, template <typename> class ARMA_VECTOR_TYPE >
int findcliquesIsEmpty( const ARMA_VECTOR_TYPE<Tv> x )
{
  for( size_t i = 0; i < x.n_elem; i++ ){
    if( x(i) ){
      return( 0 );
    }
  }
  return( 1 );
}


template< typename Tm, template <typename> class ARMA_MATRIX_TYPE, 
          typename Tv, template <typename> class ARMA_VECTOR_TYPE >
int findcliquesGETQ( const int u, 
                     const ARMA_VECTOR_TYPE<Tv>& CAND, 
                     const ARMA_MATRIX_TYPE<Tm>& A )
{
  int p = A.n_rows;
  int q;
  for( q = 0; q < p; q++ ){
    if( ( 1 == CAND( q ) ) && ( 0 == A( u, q ) ) ){
      return( q );
    }
  }
  return( -1 );
}

#endif //_CLIQUES_HPP
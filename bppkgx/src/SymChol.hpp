#ifndef _SYMCHOL_HPP
#define _SYMCHOL_HPP
/*
A Symbolic Cholesky Factorization and Completion Object.
Alex Lenkoski
 */
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

typedef class SymChol_CLASS* SymChol;//The Symbolic Cholesky Factor Class
typedef class Comp_CLASS* Comp;//The "Completion" class.

class Comp_CLASS
{
 public:
  int d_i,d_j;
  int n_dep;
  ivec dep;
  Comp_CLASS(int, int, imat);
  ~Comp_CLASS();
  void print();
};

class SymChol_CLASS
{
 public:
  int p;
  int n_diff;
  imat F;
  Comp* diff;
  SymChol_CLASS(imat);
  ~SymChol_CLASS();

  int MakeFillIn(imat);
  void print();
  template< class Matrix >
  void CompletePhi( Matrix );
};

template< class Matrix >
void SymChol_CLASS::CompletePhi( Matrix Phi )
{
  int k,l;
  double temp;
  Comp diff_k;
  int t_i,t_j;
  //------ Loop through matrix ---------------------
  for( k = 0; k < n_diff; k++ ){
    diff_k = diff[k];
    t_i = diff_k->d_i;
    t_j = diff_k->d_j;
    temp = 0;
    if(diff_k->n_dep > 0){
      for(l = 0; l < diff_k->n_dep; l++){
        temp -= Phi(diff_k->dep[l],t_i) * Phi(diff_k->dep[l],t_j);
      }
    }
    Phi( t_i,t_j) = temp / Phi(t_i, t_i);
  }
  //-------------------------------------------------------
  //  return(1);

}
#endif //_SYMCHOL_HPP
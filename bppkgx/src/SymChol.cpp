#include "SymChol.hpp"

// Constructor
Comp_CLASS::Comp_CLASS( int i, int j, imat F ){
  d_i = i;
  d_j = j;
  dep = ivec( d_i );
  n_dep = 0;

  for( int l = 0; l < d_i; l++){
    if( F( l, d_i ) && F( l, d_j ) ){
      dep( n_dep ) = l;
      n_dep++;
    }
  }
}
// Destructor
Comp_CLASS::~Comp_CLASS()
{
}
// Print
void Comp_CLASS::print()
{
  int i;
  Rprintf("Nodes %d and %d, depend on %d\n", d_i + 1, d_j + 1, n_dep );
  Rprintf("Which are: ");
  for(i = 0; i < n_dep; i++)
    {
      Rprintf("%d ",dep(i) + 1);
    }
  Rprintf("\n");

}

// Constructor
SymChol_CLASS::SymChol_CLASS( imat A )
{
  int i,j,k;
  p = A.n_rows;
  F=zeros<imat>( p, p );
  n_diff = MakeFillIn( A );
  //  F->print();

  if( n_diff > 0 ){
      diff = new Comp[n_diff];
      k = 0;
      for(i = 1; i < p - 1; i++){
        for(j = i + 1; j < p;j++){
          if( F(i,j) && !A(i,j) ){
            diff[k] = new Comp_CLASS(i,j,F);
            k++;
          }
        }
      }
    }
  }

SymChol_CLASS::~SymChol_CLASS()
{
  if( n_diff > 0){
    for(int i = 0; i < n_diff; i++ ){
      delete diff[i];
    }
    delete[] diff;
  }
}
int SymChol_CLASS::MakeFillIn( imat A )
{
  int i,j,k;
  ivec Elimination( p );
  int n_diff = 0;
  int temp;
  //--------------  Set-Up  ------------------
  for(i = 0; i < p; i++)Elimination(i) = -1;

  for(i =0;i < p - 1;i++)
    {
      for(j= i + 1; j < p; j++)
     {

       F( i, j ) =A( i, j );
     }
    }
  //-----------------------------------------

  //------ First Node is Not a Parent -------
  for(j = 1; j < p; j++){
    if( F( 0, j ) ){
      Elimination(0) = j;
      break;
    }
  }
  //-----------------------------------------

  //--- Loop Through the Remaining Nodes ---
  for(i = 1; i < p - 1; i++){
    //------- Is this a Parent? -----------
    temp = 0;
    for(j = 0; j < i; j++){
      if(Elimination(j) == i){
        for(k = i + 1; k < p; k++){
          if( F( j, k ) ){
            F( i, k ) = 1;
            if( !A( i, k ) ){
             temp = 1;
            }
          }
        }
      }
    }

    //      if(temp)n_diff++;
    //--------------------------------------

    //------ Find This Parent --------------
    for(j = i + 1; j < p; j++){
      if( F(i,j) ){
        Elimination(i) = j;
        break;
      }
    }
    //-------------------------------------
  }

  //------ End Looping Through Nodes ---------
  for(i = 0; i < p - 1; i++){
    for(j = i+1; j < p; j++){
       if(F(i,j) && !A(i,j)){
        n_diff++;
      }
    }
  }
  return( n_diff );
}

void SymChol_CLASS::print()
{
  int i;
    Rprintf("Fill In\n");
    F.print();
    Rprintf("\n");
    Rprintf("Hello %d\n",n_diff);
  Comp temp;
  if( n_diff == 0 )Rprintf("There are no differences\n");
  for( i = 0; i < n_diff;i++){
    temp = diff[i];
    Rprintf("Completion Edge %d of %d\n",i+1,n_diff);
    temp->print();
  }
}

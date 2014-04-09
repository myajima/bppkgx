/* priorpkpg.cpp
* @author Masanao Yajima
*/
#include <RcppArmadillo.h>
#include "modelpkpg.hpp"
#include "priorpkpg.hpp"

PkPgPrior::PkPgPrior()
{

}

PkPgPrior::~PkPgPrior(void)
{
}
int PkPgPrior::Init( PkPgModel& model ) {
    N  = model.N;    // No. subjects
    T  = model.T;    // No. of sampling times
    K  = model.K;    // No. of observed compartments
    Q  = model.Q;    // No. of SNPs
    Ppk= model.Ppk;  // No. of PK covariates
    Ppg= model.Ppg;  // No. of SNPs covariates
    V  = model.V;    // No. PK parameters

    // beta pk
    if( model.loadbetapk0 == true  ){
      if( !beta_pk_0.load((char*)model.betapk0file.c_str(),arma_ascii) ){
        Rf_error("Loading Beta Pk 0 failed!" );
      }
    } else {
        beta_pk_0= zeros( Ppk, V );
    }
    // B pk
    if( model.loadBpk0 == true  ){
      if( !B_pk0.load((char*)model.Bpk0file.c_str(),arma_ascii) ){
        Rf_error("Loading Beta Pk 0 failed!" );
      }
    } else {
        B_pk0 = eye( Ppk, Ppk );
    }

    // beta pg
    if( model.loadbetapg0 == true  ){
      if( !beta_pg_0.load((char*)model.betapg0file.c_str(),arma_ascii) ){
        Rf_error("Loading Beta Pg 0 failed!" );
      }
    } else {
        beta_pg_0= zeros( Ppg, Q );
    }
    // B pk
    if( model.loadBpg0 == true  ){
      if( !B_pg0.load((char*)model.Bpg0file.c_str(),arma_ascii) ){
        Rf_error("Loading Beta Pk 0 failed!" );
      }
    } else {
        B_pg0 = eye( Ppg, Ppg );
    }
    // Rho
    rho_0= zeros( Q, V );
    R_0 = eye( Q, Q );

    r1=1;
    r2=1;
    return( 0 );
}
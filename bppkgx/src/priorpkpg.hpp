#ifndef _PRIOR_PkPg_HPP
#define _PRIOR_PkPg_HPP
/* priorpkpg.hpp
 * @author Masanao Yajima
 */
#include <R.h>
#include <Rinternals.h>
#include <RcppArmadillo.h>
 
using namespace Rcpp;
using namespace arma;

class PkPgModel;
class PkPgPrior
{
public:   PkPgPrior();             // constructor
         ~PkPgPrior( void );       // destructor
    int Init( PkPgModel& model ); // 
public:
    int N;    // No. subjects
    int T;    // No. of sampling times
    int K;    // No. of observed compartments
    int Q;    // No. of SNPs
    int Ppk;  // No. of PK covariates
    int Ppg;  // No. of SNPs covariates
    int V;    // No. PK parameters
    // beta pk N( beta_pk_0, B_pk0 = inv( invB_pk0 ) )
    mat beta_pk_0;
    mat B_pk0;
    mat invB_pk0;
    // beta pg
    mat beta_pg_0;
    mat B_pg0;
    // Rho
    mat rho_0;
    mat R_0;
    // 
    int Delta_pk;
    int Delta_pg;
    mat D_prior_pk;
    mat D_prior_pg;

    // Sigma IG(r1/2, r2/2)
    int r1;
    int r2;

};
#endif //_PRIOR_PkPg_HPP

#ifndef _PkPg_MCMC_HPP
#define _PkPg_MCMC_HPP

#include <RcppArmadillo.h>
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
//#include "modelpkpg.hpp"
//#include "mcmcpkpg.hpp"

//#include "armadillo.h"
using namespace Rcpp;
using namespace arma;

RcppExport SEXP mcmc_c(     SEXP vecY_s,   SEXP W_s,     SEXP Time_s, SEXP Xpk_s,  SEXP Xpg_s,  
                            SEXP Th0_s,    SEXP PKpar,   SEXP N,      SEXP T,      SEXP G,     
                            SEXP K,        SEXP Ppk,     SEXP Ppg,    SEXP SnpFlg_s,
                            SEXP n_chain,  SEXP n_mcmc,  SEXP burn,   SEXP thin,   
                            SEXP dyn_eval, SEXP rho) ;
// RcppExport SEXP pkpg_c(     SEXP vecY_s,   SEXP W_s,     SEXP Time_s, SEXP Xpk_s,  SEXP Xpg_s,  
//                             SEXP Th0_s,    SEXP PKpar,   SEXP N,      SEXP T,      SEXP G,     
//                             SEXP K,        SEXP Ppk,     SEXP Ppg,    
//                             SEXP n_chain,  SEXP n_mcmc,  SEXP burn,   SEXP thin,   
//                             SEXP dyn_eval, SEXP rho, SEXP cluster) ;


//RcppExport void mcmc_main( SEXP dyn_eval, SEXP rho );
//RcppExport void mcmc_main( MCMCPkPg& mcmc, PkPgModel& model, PkPgResult& Result );
//RcppExport void update_fit( MCMCPkPg& mcmc, PkPgModel& model, PkPgResult& Result );
RcppExport void Save_Output();


#endif //_PkPg_MCMC_HPP

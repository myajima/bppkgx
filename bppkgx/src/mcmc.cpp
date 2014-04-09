
/***************************************************************************************/
/**** LIBRARIES ************************************************************************/
/***************************************************************************************/
#include "utilarmadillo.hpp"
#include "utilR.hpp"
#include "mcmc.hpp"
#include <R.h>
#include <Rinternals.h>
#include "mcmcpkpg.hpp"
#include "resultpkpg.hpp"
#include "modelpkpg.hpp"
//#include "gwish.hpp"

/***************************************************************************************/
/* FILES & CONSTANTS *******************************************************************/
/***************************************************************************************/
#define file_input1   "Data/PK_init.txt"
#define file_fit      "PostPkPg/fit.txt"
#define file_Gsnp     "PostPkPg/Gsnp.txt"
#define file_Gtheta   "PostPkPg/Gtheta.txt"
#define file_Grho     "PostPkPg/Grho.txt"
#define file_Rho      "PostPkPg/Rho.txt"
#define file_theta    "PostPKPg/theta.txt"
#define file_error    "PostPKPg/error.txt"
#define file_Btheta   "PostPKPg/Btheta.txt"

/***************************************************************************************/
/**** STRUCTURES & GLOBAL VARIABLES ****************************************************/
/***************************************************************************************/

/***************************************************************************************/    
/******************  R LINK  main()  ***************************************************/  
/***************************************************************************************/
SEXP  mcmc_c( SEXP vecY_s,   SEXP W_s,     SEXP Time_s, SEXP Xpk_s,  SEXP Xpg_s,  
              SEXP Th0_s,    SEXP PKpar_s, SEXP N_s,    SEXP T_s,    SEXP Q_s,     
              SEXP K_s,      SEXP Ppk_s,   SEXP Ppg_s,  SEXP SnpFlg_s,  
              SEXP n_chain,  SEXP n_mcmc,  SEXP burn,   SEXP thin,   
              SEXP dyn_eval, SEXP rho )   
{

  PkPgModel model;
  MCMCPkPg mcmc;
  PkPgResult Result;
    // Read data summaries -----------------------------------------------------------
  // model.N   = INTEGER(     N_s )[0];  // No. subjects    Fit.save( last_Fit, arma_ascii );
  // model.T   = INTEGER(     T_s )[0];  // No. of sampling times
  // model.K   = INTEGER(     K_s )[0];  // No. of observed compartments
  // model.Q   = INTEGER(     Q_s )[0];  // No. of SNPs
  // model.Ppk = INTEGER(   Ppk_s )[0];  // No. of PK covariates
  // model.Ppg = INTEGER(   Ppg_s )[0];  // No. of SNPs covariates
  // model.V   = INTEGER( PKpar_s )[0];  // No. PK parameters

  // // Read MCMC specs ---------------------------------------------------------------       
  mcmc.nCHAIN  = INTEGER( n_chain )[0];   // No. of separate mcmc chains              
  mcmc.nSAMPLE = INTEGER(  n_mcmc )[0];   // No. of mcmc simulations 
  mcmc.nBURN   = INTEGER(    burn )[0];   // No. of discarded samples
  mcmc.nTHIN   = INTEGER(    thin )[0];   // No. thinning factor 
  
  mcmc.flgSNP  = INTEGER( SnpFlg_s )[0];   // 
  // Create data structures --------------------------------------------------------
  //model.Xpg  = Rcpp::as<mat>( Xpg_s);
  //model.Xpk  = Rcpp::as<mat>( Xpk_s);
  //model.W    = Rcpp::as<imat>( W_s );
  //model.W.print();
  //model.W.save(  "./test.txt", arma_ascii );
  //model.Time = Rcpp::as<mat>(Time_s);
  //Rcpp::NumericVector vecY_r(vecY_s);           // creates Rcpp vector from SEXP
  //model.Y    = cube( vecY_r.begin(), model.T, model.K, model.N,  false );             // reuses memory and avoids extra copy
  //model.Y.print();
  model.PkModel = dyn_eval;
  model.env     = rho; 
    // model.loadSigma=false;
  model.Load( "parameters.txt" );
  mcmc.LoadData( model, Result ); 
  //model.InitializeSummary();
   //model.Y.print();
  model.NegIdxQ = field<uvec>( model.Q );

  // index for -g
  for( size_t i = 0; i < model.Q; i++ ){
    int idx = 0;
    uvec negidx(model.Q-1);
    for( size_t j = 0; j < model.Q; j++ ){
      if( i != j ){ 
        negidx( idx ) = j;
        idx += 1;
      }
    }
    model.NegIdxQ( i ) = negidx;
  } 
  model.Print();
  int resultinit_flg = Result.Init( model );
  int  priorload_flg = model.prior.Init( model );
  //int resultload_flg = Result.LoadInits( model );

  // Output files ------------------------------------------------------------------
  Rprintf("N=%d, T=%d, V=%d, K=%d, Q=%d, Pk=%d, Pg=%d \n", 
    model.N, model.T, model.V, model.K, model.Q, model.Ppk, model.Ppg);
  Rprintf("N mcmc=%d, N burn=%d, N thin=%d \n", 
    mcmc.nSAMPLE, mcmc.nBURN, mcmc.nTHIN );
 
   
 
  // Set random seed ---------------------------------------------------------------
  GetRNGstate();    
  
  // Read prior --------------------------------------------------------------------
  //Read_Prior();

  // Initialize chain parameters ---------------------------------------------------
  //Initialize_Chain(Th0);  
  
  // Gibss Sampling and Output -----------------------------------------------------
  //mcmc_main( dyn_eval, rho );
  // model.Y.print("Y");
  // model.W.print("W");
  mcmc.sample( model, Result );
  // Set random seed closure -------------------------------------------------------
  PutRNGstate();
  
  // Save data ---------------------------------------------------------------------
  //Save_Output();//Y,W);
  
  // Close C session ---------------------------------------------------------------
  // return (lag);
  return Rcpp::List::create(
          Rcpp::Named(    "ThetaAccept" ) = wrap( mcmc.acceptTheta    )
         // Rcpp::Named(  "Theta" ) = wrap( Result.Theta  ),
         // Rcpp::Named(  "Theta" ) = wrap( Result.Theta  ),
         // Rcpp::Named( "Otheta" ) = wrap( Result.Otheta ),
         // Rcpp::Named(   "Osnp" ) = wrap( Result.Osnp   ),
         // Rcpp::Named(    "Rho" ) = wrap( Result.Rho    ),
         // Rcpp::Named( "BetaPk" ) = wrap( Result.BetaPk ),
         // Rcpp::Named( "BetaPg" ) = wrap( Result.BetaPg )
         
  );
}



// SEXP  pkpg_c( SEXP vecY_s,   SEXP W_s,     SEXP Time_s, SEXP Xpk_s,  SEXP Xpg_s,  
//               SEXP Th0_s,    SEXP PKpar_s, SEXP N_s,    SEXP T_s,    SEXP Q_s,     
//               SEXP K_s,      SEXP Ppk_s,   SEXP Ppg_s,    
//               SEXP n_chain,  SEXP n_mcmc,  SEXP burn,   SEXP thin,   
//               SEXP dyn_eval, SEXP rho, SEXP cluster )   
// {

//   PkPgModel model;
//   MCMCPkPg mcmc;
//   PkPgResult Result;
//     // Read data summaries -----------------------------------------------------------
//   model.N   = INTEGER(     N_s )[0];  // No. subjects
//   model.T   = INTEGER(     T_s )[0];  // No. of sampling times
//   model.K   = INTEGER(     K_s )[0];  // No. of observed compartments
//   model.Q   = INTEGER(     Q_s )[0];  // No. of SNPs
//   model.Ppk = INTEGER(   Ppk_s )[0];  // No. of PK covariates
//   model.Ppg = INTEGER(   Ppg_s )[0];  // No. of SNPs covariates
//   model.V   = INTEGER( PKpar_s )[0];  // No. PK parameters

//   // Read MCMC specs ---------------------------------------------------------------       
//   mcmc.nCHAIN  = INTEGER( n_chain )[0];   // No. of separate mcmc chains              
//   mcmc.nSAMPLE = INTEGER(  n_mcmc )[0];   // No. of mcmc simulations 
//   mcmc.nBURN   = INTEGER(    burn )[0];   // No. of discarded samples
//   mcmc.nTHIN   = INTEGER(    thin )[0];   // No. thinning factor 
  
//   // Create data structures --------------------------------------------------------
//   model.Xpg  = Rcpp::as<mat>( Xpg_s);
//   model.Xpk  = Rcpp::as<mat>( Xpk_s);
//   model.W    = Rcpp::as<imat>(   W_s);
//   model.Time = Rcpp::as<mat>(Time_s);
//   Rcpp::NumericVector vecY_r(vecY_s);           // creates Rcpp vector from SEXP
//   model.Y    = cube( vecY_r.begin(), model.T, model.K, model.N,  false );             // reuses memory and avoids extra copy
//   model.Y.print();
//   model.PkModel = dyn_eval;
//   model.env     = rho;
//   model.env     = cluster;  
//     // std::cout << model.loadSigma <<std::endl;
//     // model.loadSigma=false;
//   model.Load("parameters.txt");
//   // if(model.loadSigma== true ) 
//   //     std::cout << "t" <<std::endl;
//   //   else
//   //     std::cout << "f" <<std::endl;
//   mcmc.LoadData( model, Result );
//    model.Y.print();
//   model.NegIdxQ = field<uvec>( model.Q );

//   // index for -g
//   for( int i = 0; i < model.Q; i++ ){
//     int idx = 0;
//     uvec negidx(model.Q-1);
//     for( int j = 0; j < model.Q; j++ ){
//       if( i != j ){ 
//         negidx( idx ) = j;
//         idx += 1;
//       }
//     }
//     model.NegIdxQ( i ) = negidx;
//   }
   
//   model.Print();
//   mcmc.LoadInits( model, Result );
//   //Result.init();

//   // Output files ------------------------------------------------------------------
//   Rprintf("N=%d, T=%d, V=%d, K=%d, Q=%d, Pk=%d, Pg=%d \n", 
//     model.N, model.T, model.V, model.K, model.Q, model.Ppk, model.Ppg);
    
//   // Initialize Output files -------------------------------------------------------
//   //Init_Output();
  
//   // Set random seed ---------------------------------------------------------------
//   GetRNGstate();    
  
//   // Read prior --------------------------------------------------------------------
//   //Read_Prior();

//   // Initialize chain parameters ---------------------------------------------------
//   //Initialize_Chain(Th0);  
  
//   // Gibss Sampling and Output -----------------------------------------------------
//   //mcmc_main( dyn_eval, rho );

//   //mcmc.sample( model, Result );
//   mcmc.update_fit_parallel( model, Result );
//   // Set random seed closure -------------------------------------------------------
//   PutRNGstate();
  
//   // Save data ---------------------------------------------------------------------
//   //Save_Output();//Y,W);
  
//   // Close C session ---------------------------------------------------------------
//   // return (lag);
//   return Rcpp::List::create(
//          Rcpp::Named("Y") = wrap(model.Y),
//          Rcpp::Named("W") = model.W
//   ) ;
// }

// /***************************************************************************************/
// /* Function: mcmc_theta1()                                                             */
// /***************************************************************************************/
// double f_theta0(double *par,  int i)
// {
//   static int     initialize = 0;
//   static double  *th0, *th0P; 
  
//   int    j, k, t;
//   SEXP   rfit;
//   double prior, like, post;
  
//   // Alloc temporary memory ----------------------------------------------------------
//   if(initialize == 0){
//     th0 = dvector(cV, 0.0);     th0P = dvector(cV, 0.0);
//     initialize = 1;
//   }
//   // Prior ---------------------------------------------------------------------------
//   for(j=0; j<cV; j++) th0[j] = par[j] - th->mth[i][j];
//   Ax(PT, th0, th0P, cV);
//   xy(th0, th0P, &prior, cV);  
//   // Likelihood ----------------------------------------------------------------------
//   for(like=0.0, k=0; k<cK; k++){
//     for(t=0; t<cT; t++){ 
//     if(is_na(conc->Y[i][k][t]) == 0)
//     like += -th->Se[k]/2*R_pow((conc->Y[i][k][t] - log(conc->fit[i][k][t] + 1.0)), 2);
//     }// t
//   }// k 
//   // Conditional Posterior -----------------------------------------------------------
//   post = -prior + like;
//   return post;
// }
// double f_theta1(double *par, SEXP dyn_eval, SEXP rho, int i)
// {
//   static int     initialize = 0;
//   static double  *th0, *th0P, *par1, **fit; 
  
//   int    v, k, t;
//   SEXP   fit1;
//   double prior, like, post;
  
//   // Alloc temporary memory ----------------------------------------------------------
//   PROTECT(fit1 = allocMatrix(REALSXP, cK, cT));
//   if(initialize == 0){
//     th0  = dvector(cV, 0.0);     th0P = dvector(cV, 0.0);
//     par1 = dvector(cV, 0.0);     fit  = dmatrix(cK, cT); 
//     initialize = 1;
//   }
//   // Prior ---------------------------------------------------------------------------
//   for(v=0; v<cV; v++) th0[v] = par[v] - th->mth[i][v];
//   Ax(PT, th0, th0P, cV);
//   xy(th0, th0P, &prior, cV);  
//   // Transform parameter vector ------------------------------------------------------
//   for(v=0; v<cV; v++) par1[v] = exp(par[v]);
//     // Likelihood ---------------------------------------------------------------------
//   fit1 = fev(dyn_eval, rho, par1, cV, i); 
//   for(k=0; k<cK; k++){
//    for(t=0; t<cT; t++){ 
//     fit[k][t] = fmax2(REAL(fit1)[k*cT + t], 0.001);
//   }}// t
//   for(like = 0.0, k=0; k<cK; k++){
//     for(t=0; t<cT; t++){ 
//       if(is_na(conc->Y[i][k][t]) == 0)
//         like += -th->Se[k]/2 * R_pow( conc->Y[i][k][t] - log(fit[k][t]+1.0), 2);
//     }// t
//   }// k
//   UNPROTECT(1);
//   // Conditional Posterior -----------------------------------------------------------
//   post = -prior + like;
//   //Rprintf("\n theta1: %f , %f \n", prior, like);
//   return post;
// }
// void mcmc_theta1(SEXP dyn_eval, SEXP rho)
// {
//   int i, v, j;
//   double fx0, fx1, ratio, u;
  
//   static int initialize = 0;
//   static double *x1, **Vx, **Vx1;
  
//   // Alloc temporaty memory ---------------------------------------------------------
//   if(initialize == 0){
//     x1  = dvector(cV, 0.0);   Vx = dmatrix(cV, cV);
//     Vx1 = dmatrix(cV, cV);    initialize = 1;
//     // Define the proposal distribution -------------------------------------------
//     V_Cov(th->theta, cN, cV, Vx);
//   }
//   // Loop subjects and parameters ---------------------------------------------------
//   for(i=0; i<cN; i++){
//     fx0   = f_theta0(th->theta[i], i);                 // post old
//     rA(mc->c_th, Vx, Vx1, cV, cV);                   // proposal Covariance
//     r_mvnorm_C(x1,  Vx1, th->theta[i], cV);            // proposal
//     fx1   = f_theta1(x1, dyn_eval, rho, i);            // post new
//     ratio = fx1 - fx0;                                 // M-H ratio
//     u     = Rf_runif( 0, 1 );
//     if( log( u ) < ratio ){
//       for(v=0; v<cV; v++) th->theta[i][v] = x1[v];   // Update theta_v
//     }    
//   }//i
// } // END mcmc_theta -------------------------------------------------------------------


// /***************************************************************************************/
// /* Function: Save_Output()                                                             */
// /***************************************************************************************/
// void Save_Output()//cube Y, mat W)
// {
//   //A.save("A1.mat");  // default save format is arma_binary
//   //A.save("A2.mat", arma_ascii);
//   //A.save("A3.mat", raw_binary);
//   //mat B;
//   // automatically detect format type
//   //B.load("A1.mat");
//   //colvec Wc =  W.col(2);
//   //rowvec Wr =  W.row(2);
//   //fit_list.save(file_fit, arma_ascii);
//   //theta_list.save(file_theta, arma_ascii);
    
  
//   //Y.save("PostPKPg/Y.txt", arma_ascii); 
//   //Y.print("Y:");
//   //W.save("PostPKPg/W.txt", arma_ascii); 
// } // END Output_Samples ---------------------------------------------------------------


// /***************************************************************************************/
// /* Function: Output_Samples()                                                          */
// /***************************************************************************************/
// void Output_Samples()
// { 
//   int i, j, t, k; 
//   // Fit ------------------------------------------------------------------------------
//   for(t=0; t<cT; t++){
//    for(i=0; i<cN; i++){
//     for(k=0; k<cK; k++) fprintf(out1, "%f ", conc->fit[i][k][t]);
//    }//i
//    fprintf(out1, "\n");
//   }//t
//   // G_snp - (Upper triang) -----------------------------------------------------------
//   for(i=0; i<(cG-1); i++){
//    for(j=(i+1); j<cG; j++)  fprintf(out2, "%d ", GG->snp[i][j]);
//   }//i
//   fprintf(out2, "\n");
//   // G_theta - (Upper triang) ---------------------------------------------------------
//   for(i=0; i<(cV-1); i++){
//    for(j=(i+1); j<cV; j++)  fprintf(out3, "%d ", GG->th[i][j]);
//   }//i
//   fprintf(out3, "\n"); 
//   // G_rho - (By row) -----------------------------------------------------------------
//   for(i=0; i<cG; i++){
//    for(j=0; j<cV; j++)  fprintf(out4, "%d ", GG->rho[i][j]);
//   }//i
//   fprintf(out4, "\n");
//   // Rho - (By row) -------------------------------------------------------------------
//   for(i=0; i<cG; i++){
//    for(j=0; j<cV; j++)  fprintf(out5, "%f ", th->Rho[i][j]);
//   }//i
//   fprintf(out5, "\n");
//   // theta - (By row) -----------------------------------------------------------------
//   for(i=0; i<cN; i++){
//   for(j=0; j<cV; j++)  fprintf(out6, "%f ", th->theta[i][j]);
//   }//i
//   fprintf(out6, "\n");
//   // Error ----------------------------------------------------------------------------
//   for(j=0; j<cK; j++)  fprintf(out7, "%f ", th->Se[j]);
//   fprintf(out7, "\n");
//   // Btheta - (By row) -----------------------------------------------------------------
//   for(i=0; i<cPpk; i++){
//     for(j=0; j<cV; j++)  fprintf(out8, "%f ", th->Bth[i][j]);
//   }//i
//   fprintf(out8, "\n");
  
// } // END Output_Samples ---------------------------------------------------------------


/****************************************************************************************/

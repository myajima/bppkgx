#include "mcmcpkpg.hpp"
//#include "resultpkpg.hpp"
//#include "modelpkpg.hpp"
#include "BayesGLassoGDP.hpp"
#include "BayesGLasso_Columnwise.hpp"
#include "BayesALasso.hpp"
#include "BayesLasso.hpp"
/*************************************************************************
 * Constructor                                                           *
 *************************************************************************/
MCMCPkPg::MCMCPkPg(){

}

/*************************************************************************
 * Destructor                                                            *
 *************************************************************************/
MCMCPkPg::~MCMCPkPg(void){

}

/*************************************************************************************/
/* Main MCMC Function: sample()                                                      */
/*************************************************************************************/
void MCMCPkPg::sample( PkPgModel& model, PkPgResult& Result )//(SEXP dyn_eval, SEXP rho)
{
  cIter = 0;
  //int tidx = 0;
  //int chain = 0;
  // Result.InitializeSummary( model );
  int total_iter = nBURN + nTHIN * nSAMPLE;
  acceptTheta=zeros<mat>( model.N, total_iter );
  for( cIter = 0; cIter < total_iter; cIter++ ){      
    //std::cout << "iteration: "<< cIter << std::endl;
    Rprintf("iteration: %d\n", cIter+1);
    //for( chain = 0; chain < NCHAIN; chain++ ){
    //Rprintf("Update Omega theta\n");
    sampleOtheta( model, Result );

    if(flgSNP == 1){
        //Rprintf("Update Omega SNP\n");
        sampleOsnp( model, Result );
    
        //Rprintf("Update Z\n");
        sampleZ( model, Result );
    
        //Rprintf("Update Beta Pg\n");
        sampleBetaPg( model, Result );
        //Rprintf("Update Rho\n");
        sampleRho( model, Result );
    }

    //Rprintf("Update Beta Pk\n");
    sampleBetaPk( model, Result );

    // update Z*Rho
    // SNPs ---------------------------------------------------------------------

    //Rprintf("Update means\n");
    Result.update_mth( model );

    // PK -----------------------------------------------------------------------
    //Rprintf("Sample Sigma\n");
    sampleSigma( model, Result );

    //Rprintf("Update Fit\n");
    update_fit( model, Result );

    //Rprintf("Sample Theta\n");
    sampleTheta( model, Result );  // make sure Result.mTheta is up to date

    //Rprintf("Update Theta Bar\n");
    update_theta_bar( model, Result );

    //Rprintf("Update Sigma Theta\n");
    update_sigma_theta( model, Result );
    // update_fit_parallel( model, Result );

    //mcmc_Sth();
    // RJ -----------------------------------------------------------------------
    // 

    // record
    if( ( cIter >= nBURN ) && ( ( cIter - nBURN ) % nTHIN == 0 ) ) {
      Result.SaveDraws( model );
    }
  }
  Result.SaveFinal( model );
 } // END MCMC_main ---------------------------------------------------------------------

// SNPs ---------------------------------------------------------------------
int MCMCPkPg::sampleBetaPg( PkPgModel& model, PkPgResult& Result ){
  mat I_Q       = eye( model.Q, model.Q ); 
  mat Sb        = inv( (1/model.N)*I_Q + Result.Osnp );
  mat betahat   = model.invXtXpg  * model.Xpg.t() * Result.Z;
  mat pcbeta    =   ( betahat * Result.Osnp  + (1/model.N) * model.prior.beta_pg_0 ) * Sb ;
  Result.BetaPg = M_norm( pcbeta,  model.invXtXpg, Sb );
  Result.XBetaPg= model.Xpg * Result.BetaPg;
  return 0;
}

int MCMCPkPg::sampleOsnp( PkPgModel& model, PkPgResult& Result ){
   cube Ores;   
   cube Sres;
   mat betapgresid = Result.Z - Result.XBetaPg;
   mat D_post_pg   = betapgresid.t() * betapgresid;

  if( 1 ){// Adaptive Bayesian Graphical Lasso
    mat Lres;
    BayesGLassoGDP( D_post_pg, model.N, Result.Ssnp, Result.Osnp, 1, 1, 1.0e-2, 1.0e-4, Sres, Ores, Lres );
    Result.LAMBDAsnp = Lres.row( 0 );
  } else {// Bayesian Graphical Lasso
    vec Lres;
    BayesGLasso_Columnwise(D_post_pg, model.N, Result.Ssnp, Result.Osnp, 1, 1, 
                            1, 0.1,
                          Sres, Ores, Lres );
    Result.LAMBDAsnp = Lres( 0 );
  }
  // G-Wishart
  // gwish_dmh_update_all( Result.Osnp, Result.Gsnp, model.prior.Delta_pg, 
  //                       model.N, model.prior.D_prior_pg, D_post_pg, 
  //                       Kres, Gres );

  Result.Osnp = Ores.slice( 0 );
  Result.Ssnp = Sres.slice( 0 );
  //Result.Gsnp = Gres;
  return 0;
}
double ZgivenW( int w, double mz, double sz)
{
  double tempZ;
  if( w == 0){ 
    tempZ = rand_truncated_normal_rej_above( mz,sz, 0 );
  } else if( w == 1 ) { 
    tempZ = rand_tn( mz, sz, 0, 1 );
  } else{ 
    tempZ = rand_truncated_normal_rej_below( mz, sz, 1 );
  }
  return(tempZ);
}

int MCMCPkPg::sampleZ( PkPgModel& model, PkPgResult& Result ){
 //zeros<mat>( model.N, model.Q ); // temporaly sample Z_ig - mu_ig
  mat RO = Result.Rho * Result.Otheta;
  mat RORt = RO * Result.Rho.t();
  mat SZi  = Result.Osnp + RORt;
  mat cZU  = chol( SZi );

  //std::cout <<"g0"<<std::endl;
  uvec aidx( 1 );
  umat bidx_temp( model.Q, 1 );
  for( size_t g = 0; g < model.Q; g++ ){
    bidx_temp(  g, 0 ) = g;
  }

  mat clogt = ( Result.logTheta- Result.XBetaPk );
  mat tempMz0 = Result.Osnp * Result.XBetaPg.t() + RO *clogt.t();
  mat v0      = solve( trimatl( trans( cZU ) ), tempMz0 );
  mat mz0     = trans(solve( trimatu( cZU ), v0 ) ); 
  //std::cout <<"g"<<std::endl;
  for( size_t i = 0; i < model.N; i++ ){
     vec z   = Result.Z.row(i).t();
      for( size_t j = 0; j < model.Q; j++ ){
      //std::cout <<"g2"<<std::endl;
      // mat tempMz = Result.Osnp * Result.XBetaPg.row(i).t() + RO *clogt.row(i).t() ;
      // mat v      = solve( trimatl( trans( cZU ) ), tempMz );
      // vec mz     = solve( trimatu( cZU ), v ) ; 
      vec mz     = mz0.row(i).t();
      aidx( 0 ) = j;
      //mz.print("mz");
      //std::cout <<"g1"<<std::endl;
      umat bidxt = bidx_temp;
      //std::cout <<"g2"<<std::endl;
      bidxt.shed_row( j );
      //std::cout <<"g3"<<std::endl;
      uvec bidx = bidxt.col( 0 );//conv_to< uvec >::from( bidxt);
      //std::cout <<"g4"<<std::endl;
      double sdz1 =pow( SZi( j, j ), - 0.5 );
      //      std::cout <<"g5"<<std::endl;
      double mzc1  = as_scalar( SZi.submat(aidx,bidx)*(z.elem(bidx)-mz.elem(bidx)));
      double mzc   = mzc1/SZi( j, j );
       ///     std::cout <<"g6"<<std::endl;
      z( j ) = ZgivenW( model.W( i, j ),mzc, sdz1 );
      // std::cout <<"mz( i, j )"<<mz( j )<<std::endl;
      // std::cout <<"z( i, j )"<<z( i, j )<<std::endl;
      //      std::cout <<"g3"<<std::endl;
      Result.Z( i, j ) = z( j ) ; 
      //  std::cout <<"g10"<<std::endl;
    }

  }


  // for( size_t i = 0; i < model.N; i++ ){
  //   int Qidx = model.Q - 1;
  //   double sdz1 = cZU( Qidx, Qidx );//1.0 / cZU( Qidx, Qidx );
  //   double mzc1 = mz(  i, Qidx);//double mzc1 = mz( i, Qidx );
  //   z( i, Qidx ) = ZgivenW(model.W( i, Qidx ),mzc1, sdz1 );
  //   for( int g = ( model.Q - 1 ); g >= 0 ; g-- ){
  //     double mzc = 0;//double mzc = mz( i, g );
  //     for( size_t j = g + 1; j < model.Q ; j++ ){
  //       //mzc -= ( z( i, j ) - mz( i, j ) ) * cZU( g, j )/cZU( g, g );
  //       mzc -= ( z( i, j ) - mz( i, j ) ) * cZU( g, j );//cZU( g, g )
  //     } // j
  //     double sdz = cZU( g, g ) ;//1.0 / cZU( g, g ) ;
  //     z( i, g ) = ZgivenW( model.W( i, g ), mzc, sdz ) + mz( i, g );
  //     Result.Z( i, g ) = z( i, g ) ;  
  //   }//g
  // }//i    

  Result.ZtZ = Result.Z.t() * Result.Z; // update ZtZ also
  Result.ZRho=  Result.Z * Result.Rho;
  return 0;
}
// int MCMCPkPg::sampleZ( PkPgModel& model, PkPgResult& Result ){
// // Based on Gaussin Markov Random Fields Theory and Applications by Rue Et.al
// //            Z_a - mZ_a | Z_b ~ Nc( -Q_ab( Z_b - mZ_b ), Q_aa )
// //            mZ : marginal mean of Z
// //            Q  : precision matrix of Z
//   mat z( model.N, model.Q ); // temporaly sample Z_ig - mu_ig
//   mat RtOR = Result.Rho.t() * Result.Otheta * Result.Rho;
//   colvec One_n( model.N ); 
//   mat hatZ = inv( RtOR )* Result.Rho.t() * Result.Otheta * ( Result.Theta - One_n * Result.Alpha.t() - Result.XBetaPk );
//   for( int i = 0; i < model.N; i++ ){
//     mat SZi  = Result.Osnp + Result.Tau(i) * RtOR;
//     // posterior mean of Z_i = SZ_i^-1 * ( Osnp %*% BetaPg %*% t( XBetaPg_i ) + tau_i*t( Rho )%*%Otheta%*% Rho%*%t(hatZ) )  
//     vec pmZi = inv( SZi ) * ( Result.Osnp * Result.XBetaPg.row( i ).t() + Result.Tau( i ) * RtOR * hatZ.row( i ).t() );
//     vec Zi = Result.Z.row( i );
//     for( int g = 0; g < model.Q; g++ ){
//     //mat mZb = Result.Z.submat( span::all, Neg ) - Result.PmZ.submat(span::all,Neg);
//       uvec gidx(1);
//       gidx(0) = g;
//       double mzc = as_scalar( -SZi( gidx, model.NegIdxQ( g ) ) * ( Zi.elem( model.NegIdxQ( g ) ) - pmZi.elem( model.NegIdxQ( g ) ) ) );
//       double sdz = pow( SZi( g, g ), 0.5 );
//       if( model.W( i, g ) == 0){ 
//         //z = norm_trunc(mzc, hg, -Inf0, alpha);
//         z( i, g ) = rand_truncated_normal_rej_above( mzc, sdz, 0 );
//       } else if(model.W(i,g) == 1) { 
//         //z = norm_trunc(mzc, hg, alpha, beta);
//         z( i, g ) = rand_tn( mzc, sdz, 0, 1 );
//       } else{ 
//         //z = norm_trunc(mzc, hg, beta, Inf0);
//         z( i, g ) = rand_truncated_normal_rej_below( mzc, sdz, 1 );
//       } 
//       // if(is_na(snp->W[i][g])) snp->Z[i][g] = rnorm(0, 2.33);
//       // if(snp->Z[i][g] >  200) snp->Z[i][g] = rnorm(200, 5);
//       // if(snp->Z[i][g] < -200) snp->Z[i][g] = rnorm(-200, 5);
//       //Rprintf("(%d,%d) w=%d, z=%f\n", i, g, snp->W[i][g], snp->Z[i][g]);  
//       Result.Z = z( i, g ) + pmZi( g );  
//     }//g
//   }//i    
 
//   Result.ZtZ = Result.Z.t() * Result.Z; // update ZtZ also
//   return 0;
// }
/***************************************************************************************/
/* Function: sampleRho()                                                                */
/***************************************************************************************/
int MCMCPkPg::sampleRho( PkPgModel& model, PkPgResult& Result ){

  mat ymat   = Result.logTheta - Result.XBetaPk;
  vec yvec   = reshape( ymat, ymat.n_elem, 1 );
  mat Xmat   = kron( eye<mat>( model.V, model.V ), Result.Z );
  vec rhovec = reshape( Result.Rho, Result.Rho.n_elem, 1 );
  if( 1 ){
    // adaptive lasso
    mat resrho;
    vec resdelta, resiota;
    Bayes_AL( yvec, Xmat, 
              0, 1, 1,
              rhovec, Result.delta, Result.iota,
              resrho, resdelta, resiota,
              false, true, true ); 
    Result.Rho   = reshape( resrho.row( 0 ), Result.Rho.n_rows, Result.Rho.n_cols );
    Result.delta = resdelta( 0 );
    Result.iota  = resiota( 0 );
    Result.ZRho  = Result.Z * Result.Rho;
  } else{
    // lasso
    mat resrho;
    vec sigma2Samples;  
    mat invTau2Samples; 
    vec lambdaSamples;  
  
    int initialize = 0;
    if( Result.Rhoinit == 1 ){ // one time only
      initialize = 1;
      Result.Rhoinit = 0;
    }

    Bayes_Lasso( yvec, Xmat, 0, 1, 1, 
                 initialize,
                 rhovec, Result.Rhos2, Result.Rhot2, Result.Rholb,
                 resrho, sigma2Samples, invTau2Samples, lambdaSamples  );
  
    Result.Rho   = reshape( resrho.row( 0 ), Result.Rho.n_rows, Result.Rho.n_cols );
    Result.Rhos2 = sigma2Samples( 0 );
    Result.Rhot2 = invTau2Samples.row( 0 ).t();
    Result.Rholb = lambdaSamples( 0 ); 
    Result.ZRho  =  Result.Z * Result.Rho;
  }
  return 0;
} // END sampleRho --------------------------------------------------------------------

// PK -----------------------------------------------------------------------
int MCMCPkPg::sampleSigma( PkPgModel& model, PkPgResult& Result ){
  for( size_t i = 0; i < model.N; i++ ){
    mat Y   = model.Y.slice( i );      // log( Y )
    mat Fit = Result.Fit.slice( i ); // need to check if it's logged
    for( size_t k = 0; k < model.K; k++ ){
      double pr1 = 0.5 * ( model.prior.r1 + model.prior.T );
      vec resid  = log( Y.col( k ) ) - log( Fit.col( k ) );
      double ssr = dot(resid, resid );
      double pr2 = 0.5 * ( model.prior.r2 + ssr );
      Result.Sigma( i, k ) = rinvgamma( pr1, pr2 );
    }
  }
  return 0;
}
/***************************************************************************************/
/* sampleTheta()                                                                       */  
/*     Metropolis Hasting step for Theta                                               */
/***************************************************************************************/
// Posterior Density Function to sample Theta 
double f_theta( colvec logTheta, colvec mTheta, mat Otheta, double Tau, mat Y, mat Fit, rowvec Sigma, double logEHRTime ) 
{
  double prior, like, post;
  colvec cTheta = ( logTheta - mTheta ) ;
  prior = - 0.5 * Tau * as_scalar( cTheta.t() * Otheta * cTheta );
  mat D =  diagmat( 1 / Sigma );
  like  = - 0.5 * accu( pow ( log( Y ) - log( Fit ), 2 ) * D ); // Need to figure out what to do with Sigma
  // Conditional Posterior -----------------------------------------------------------
  post = prior + like + Rf_dnorm4( logEHRTime, 0, 1, 1 );
  return post;
}

int MCMCPkPg::sampleTheta( PkPgModel& model, PkPgResult& Result ){
  //int i, v, j;
  double fx0, fx1, ratio, u;
  mat proplogTheta( model.N, model.V ); // Proposal for Theta.  In matrix with vision of parallelizing in the future
  mat propTheta( model.N, model.V ); 
  // Loop subjects ---------------------------------------------------
  for( size_t i = 0; i < model.N; i++ ){
    // current posterior log likelihood
    colvec CurrentlogTheta = Result.logTheta.row( i ).t();
    fx0 = f_theta( CurrentlogTheta,     Result.mTheta.row( i ).t(), 
                   Result.Otheta,    1.0, //Result.Tau( i ), 
                   model.Y.slice( i ), Result.Fit.slice( i ), 
                   Result.Sigma.row( i ) , log(Result.EHRTime(i,0)));              
    // propoal ----------- start -----------
    colvec theta_bar_i = Result.logTheta_bar.row( i ).t();
    mat sigma_theta_i  = Result.Sigma_theta.slice( i );
    //colvec ProplogThetaCenter = CurrentlogTheta + 2 * ( Result.logTheta_MLE.row(i).t() - CurrentlogTheta );
    //colvec logtemp = MVNORM( 0, ProplogThetaCenter , sigma_theta_i ) ;      // proposal

    colvec logtemp = MVNORM( 0, CurrentlogTheta , sigma_theta_i ) ;      // proposal
    proplogTheta.row( i ) = logtemp.t();
    //logtemp.print("logtemp");
    propTheta.row( i )   = exp( proplogTheta.row( i ) );
    rowvec tempPThetaRow(15);
    for( int j = 0; j< 10; j++)tempPThetaRow(j) = propTheta( i, j );
    double proplogEHRTime = Rf_rnorm(log(Result.EHRTime(i,0)),0.1);
    double propEHRTime = exp( proplogEHRTime );
    tempPThetaRow(10) =  propEHRTime; //Result.EHRTime(i,0);
    for( int j = 10; j< 14; j++) tempPThetaRow(j+1) = propTheta( i, j );
    rowvec param        = tempPThetaRow;
    mat TempFit         = fev( model.PkModel, model.env, model.K, model.T, param, i ); // update the likelihood
    if( TempFit.n_cols==model.K){
      uvec idx            = find( TempFit < 0.00001 );
      TempFit.elem( idx ) = ones<vec>( idx.n_elem ) * 0.00001;
      // propoal ----------- end -----------
      // proposed posterior log likelihood  
      fx1 = f_theta( proplogTheta.row( i ).t(), Result.mTheta.row( i ).t(), 
                     Result.Otheta,    1.0, //Result.Tau(i),
                     model.Y.slice( i ), TempFit, 
                     Result.Sigma.row( i ), proplogEHRTime ); 
      // Accept or Reject
      ratio = fx1 - fx0;     // M-H ratio                    
      u = Rf_runif( 0, 1 ); 
      if( log( u ) < ratio ){
        // for( v = 0; v < model.V; v++ ) {
        //   Result.Theta( i, v ) = propTheta( i, v );   // Update theta_v
        // }
        acceptTheta(i,cIter)=1;
        Result.logTheta_old.row( i ) = Result.logTheta.row( i );
        Result.logTheta.row( i )     = proplogTheta.row( i );
        Result.Theta.row( i )        = propTheta.row( i );
        Result.Fit.slice( i )        = TempFit;
        Result.EHRTime(i,0)          = propEHRTime;
      }  else {
        Result.logTheta_old.row( i ) = Result.logTheta.row( i );
        Result.logTheta.row( i )     = Result.logTheta.row( i );
        Result.Theta.row( i )        = Result.Theta.row( i );
        Result.Fit.slice( i )        = Result.Fit.slice( i );
      }
    } else{
      acceptTheta(i,cIter)=-1;
      Result.logTheta_old.row( i ) = Result.logTheta.row( i );
      Result.logTheta.row( i )     = Result.logTheta.row( i );
      Result.Theta.row( i )        = Result.Theta.row( i );
      Result.Fit.slice( i )        = Result.Fit.slice( i );
    }

  } // end loop i

  return 0;
}


// int MCMCPkPg::sampleTheta( PkPgModel& model, PkPgResult& Result ){
//   //int i, v, j;
//   double fx0, fx1, ratio, u;
//   mat proplogTheta( model.N, model.V ); // Proposal for Theta.  In matrix with vision of parallelizing in the future
//   mat propTheta( model.N, model.V ); 
//   // Loop subjects ---------------------------------------------------
//   for( int i = 0; i < model.N; i++ ){
//     // current posterior log likelihood
//     colvec CurrentlogTheta = Result.logTheta.row( i ).t();
//     fx0 = f_theta( CurrentlogTheta,     Result.mTheta.row( i ).t(), 
//                    Result.Otheta,    1.0, //Result.Tau( i ), 
//                    model.Y.slice( i ), Result.Fit.slice( i ), 
//                    Result.Sigma.row( i ) , log(Result.EHRTime(i,0)));              
//     // propoal ----------- start -----------
//     colvec theta_bar_i = Result.logTheta_bar.row( i ).t();
//     mat sigma_theta_i  = Result.Sigma_theta.slice( i );
//     //std::cout << i << " : " ;
//     //theta_bar_i.print("thetabar");
//     //colvec ProplogThetaCenter = CurrentlogTheta + 2 * ( Result.logTheta_MLE.row(i).t() - CurrentlogTheta );
//     //colvec logtemp = MVNORM( 0, ProplogThetaCenter , sigma_theta_i ) ;      // proposal
//     colvec logtemp = MVNORM( 0, CurrentlogTheta , sigma_theta_i ) ;      // proposal
//     proplogTheta.row( i ) = logtemp.t();
//     //logtemp.print("logtemp");
//     propTheta.row( i )   = exp( proplogTheta.row( i ) );
//     rowvec tempPThetaRow(15);
//     for( int j = 0; j< 10; j++)tempPThetaRow(j) = propTheta( i, j );
//     double proplogEHRTime = Rf_rnorm(log(Result.EHRTime(i,0)),0.1);
//     double propEHRTime = exp( proplogEHRTime );
//     tempPThetaRow(10) =  propEHRTime; //Result.EHRTime(i,0);
//     for( int j = 10; j< 14; j++) tempPThetaRow(j+1) = propTheta( i, j );
//     rowvec param        = tempPThetaRow;
//     //param.print("param");
//     //param.print("param");
//     mat TempFit         = fev( model.PkModel, model.env, model.K, model.T, param, i ); // update the likelihood
//     if(TempFit.n_cols==model.K){
//       uvec idx            = find( TempFit < 0.00001 );
//       TempFit.elem( idx ) = ones<vec>( idx.n_elem ) * 0.00001;
//       // propoal ----------- end -----------
//       // proposed posterior log likelihood  
//       fx1 = f_theta( proplogTheta.row( i ).t(), Result.mTheta.row( i ).t(), 
//                      Result.Otheta,    1.0, //Result.Tau(i),
//                      model.Y.slice( i ), TempFit, 
//                      Result.Sigma.row( i ), proplogEHRTime ); 
//       // Accept or Reject
//       ratio = fx1 - fx0;     // M-H ratio                    
//       u = Rf_runif( 0, 1 ); 
//       if( log( u ) < ratio ){
//         // for( v = 0; v < model.V; v++ ) {
//         //   Result.Theta( i, v ) = propTheta( i, v );   // Update theta_v
//         // }
//         //std::cout << "accept"<< std::endl;
//         acceptTheta(i,cIter)=1;
//         Result.logTheta_old.row( i ) = Result.logTheta.row( i );
//         Result.logTheta.row( i )     = proplogTheta.row( i );
//         Result.Theta.row( i )        = propTheta.row( i );
//         Result.Fit.slice( i )        = TempFit;
//         Result.EHRTime(i,0)          = propEHRTime;
//       }  else {
//         //std::cout << "fail"<< std::endl;
//         Result.logTheta_old.row( i ) = Result.logTheta.row( i );
//         Result.logTheta.row( i )     = Result.logTheta.row( i );
//         Result.Theta.row( i )        = Result.Theta.row( i );
//         Result.Fit.slice( i )        = Result.Fit.slice( i );
//       }
//     } else{
//       //std::cout << i <<std::endl;
//       acceptTheta(i,cIter)=-1;
//       Result.logTheta_old.row( i ) = Result.logTheta.row( i );
//       Result.logTheta.row( i )     = Result.logTheta.row( i );
//       Result.Theta.row( i )        = Result.Theta.row( i );
//       Result.Fit.slice( i )        = Result.Fit.slice( i );
//     }

//   } // end loop i

//   return 0;
// }
/***************************************************************************************/
/* update_fit()                                                                        */
/***************************************************************************************/
void MCMCPkPg::update_fit( PkPgModel& model, PkPgResult& Result )
{
  // Compute lsoda solution ---------------------------------------------------------
  mat ThetaTemp = Result.Theta;
  ThetaTemp.insert_cols( 10, Result.EHRTime );
  for( size_t i = 0; i < model.N; i++ ){
    rowvec param          = ThetaTemp.row( i ); 
    mat tempMat           = fev( model.PkModel, model.env, model.K, model.T, param, i ); 
    uvec idx              = find( tempMat < 0.00001 );
    tempMat.elem( idx )   = ones<vec>( idx.n_elem ) * 0.00001;
    Result.Fit.slice( i ) = tempMat;
  } //i
} // END update_fit -------------------------------------------------------------------
/***************************************************************************************/
/* update_fit()                                                                        */
/***************************************************************************************/
void MCMCPkPg::update_fit_parallel( PkPgModel& model, PkPgResult& Result )
{
  int           i ;
  Rprintf("starting update fit");
  // // Compute lsoda solution ---------------------------------------------------------
  // for( i = 0; i < model.N; i++ ){
  //   rowvec param = exp( Result.Theta.row(i) );
  //     //param.print();
  //     //mat TEST= fev( model.PkModel, model.env, model.K, model.T, param, i );
  //     //TEST.print();
  //     //Result.Fit.slice(i) = fev( model.PkModel, model.env, model.K, model.T, param, i ); 
  //     Result.Fit(span( i*model.K, i*model.K+model.K-1 ),span::all) = fev( model.PkModel, model.env, model.K, model.T, param, i ); 
  // }//i
  Result.Fit = fev_parallel( model.PkModel, model.env, model.cluster, model.K, model.T, model.N, Result.Theta, i ); 
} // END update_fit -------------------------------------------------------------------

void MCMCPkPg::update_theta_bar( PkPgModel& model, PkPgResult& Result ){
  Result.logTheta_bar_old = Result.logTheta_bar;// Theta_bar(l-1)
  Result.logTheta_bar     = ( Result.logTheta_bar * ( cIter + 1 ) + Result.logTheta )/( cIter + 2 );// Theta_bar(l)
}
void MCMCPkPg::update_sigma_theta( PkPgModel& model, PkPgResult& Result ){
  mat eIv( model.V, model.V ); 
  eIv.eye();
  eIv = eIv * model.epsilon;
  double sd =( 5.6644 / model.V ); 
  for(size_t i = 0; i< model.N; i++ ){
      rowvec theta_bar_old = Result.logTheta_bar_old.row( i );
      rowvec theta_bar_now = Result.logTheta_bar.row( i );
      rowvec theta_now     = Result.logTheta.row( i );
      rowvec theta_old     = Result.logTheta_old.row( i );
      if( cIter < 5000 ){ 
         Result.SS_theta.slice( i ) = Result.SS_theta.slice( i ) + theta_now.t() * theta_now;
      } else if( cIter >= 5000 ){
         Result.SS_theta.slice( i ) = Result.SS_theta.slice( i ) + theta_now.t() * theta_now;
         if( Rf_rbinom(1,0.01) == 1 ){ 
                 Result.Sigma_theta.slice(i)  =  sd * ( (  Result.SS_theta.slice( i ) 
                                                         - ( cIter + 2 ) * theta_bar_now.t()*theta_bar_now )/( cIter + 1 ) 
                                                          + eIv ) ; //sd *  
          } else {
                  Result.Sigma_theta.slice(i) = Result.SigmaTheta0;
          }                                  // /(cIter);
                                                              //+ sd * eIv;
        //Result.Sigma_theta.slice(i).print("Result.Sigma_theta.slice(i)");                                                    
      } //else{ 
        // if( cIter> 20 ){ //T0=10
        //   mat to = cIter*theta_old.t() * theta_bar_old;
        //   mat tn=( cIter + 1 ) * theta_bar_now.t()*theta_bar_now;
        //   mat tt=theta_now.t() * theta_now;
          // Result.Sigma_theta.slice( i ) = (  (cIter) / cIter+1 ) * Result.Sigma_theta.slice( i ) 
          //                               + ( sd / ( cIter + 1 ) ) 
          //                               * (  ( cIter + 1 ) * theta_bar_old.t() * theta_bar_old
          //                                  - ( cIter + 2 ) * theta_bar_now.t() * theta_bar_now
          //                                  + theta_now.t() * theta_now
          //                                  + eIv );
        // // } else {
          
        //} 
    //}
  }
}
int MCMCPkPg::sampleBetaPk( PkPgModel& model, PkPgResult& Result  ){
  mat Sb        = inv( model.prior.B_pk0 + model.XtXpk );
  mat betahat   = model.invXtXpk * model.Xpk.t() * ( Result.logTheta - Result.ZRho );
  mat pcbeta    = Sb * ( model.prior.B_pk0 * model.prior.beta_pk_0 + model.XtXpk * betahat );
  Result.BetaPk = M_norm( pcbeta, Sb, Result.Stheta );
  Result.XBetaPk= model.Xpk * Result.BetaPk; // update summary
  return 0;
}

// int MCMCPkPg::sampleAlpha( PkPgModel& model, PkPgResult& Result  ){
//     double alphahat = as_scalar( Result.Tau.t() * ( Result.Theta - Result.XBetaPk - Result.ZRho ));
//     double sumtau = accu( Result.Tau );
//     alphahat = alphahat / sumtau;
//     mat sumtOth = sumtau * Result.Otheta;
//     mat Pa      = model.prior.invA_0 + sumtOth;  // precision 
//     mat Sa      = inv( Pa );                    // covariance
//     vec pcalpha  = Sa * ( model.prior.invA_0 * model.prior.alpha_0 *  + sumtOth * alphahat );
//     vec propAlpha = MVNORM( 1, pcalpha, Pa ); //1 for precision
//     return 0;
// }


int MCMCPkPg::sampleOtheta( PkPgModel& model,PkPgResult& Result ){
  cube Sres, Ores ;
  imat Gres;
  mat betapkresid = Result.BetaPk - model.prior.beta_pk_0;
  //betapkresid.print("betapkresid ");

  Result.logTheta_Resid = Result.logTheta - (Result.XBetaPk + Result.ZRho);//Result.mTheta;
  mat D_post_pk   = Result.logTheta_Resid.t() * Result.logTheta_Resid
                    + betapkresid.t() * model.prior.B_pk0 * betapkresid; //model.prior.Delta_pk 

  // Adaptive
  if( 1 ){
    mat Lres;
    BayesGLassoGDP( D_post_pk, model.N + model.Ppk, Result.Stheta, Result.Otheta, 1, 1, 
      1.0e-2, 1.0e-4, Sres, Ores, Lres );
    Result.LAMBDAtheta = Lres.row( 0 );
  } else {// GLasso
    vec Lres;
    // Result.Stheta.print("S"); 
    // Result.Otheta.print("O");
    // D_post_pk.print("D_post_pk");
    BayesGLasso_Columnwise( D_post_pk,  model.N + model.Ppk, Result.Stheta, Result.Otheta, 1, 1, 
                            1, 0.1,
                            Sres, Ores, Lres );
    // Result.Stheta.print("S"); 
    // Result.Otheta.print("O");
    Result.LAMBDAtheta = Lres( 0 );
  }
  // G-Wishart
  // gwish_dmh_update_all( Result.Otheta, Result.Gtheta, model.prior.Delta_pk, 
  //                       model.N + model.Ppk, model.prior.D_prior_pk, D_post_pk, 
  //                       Kres, Gres );
  Result.Otheta      = Ores.slice( 0 );
  Result.Stheta      = Sres.slice( 0 );
  return 0;
}


bool MCMCPkPg::LoadData( PkPgModel& model, PkPgResult& Result )
{
  Rprintf("Loading Data \n");
  if( !model.Y.load((char*)model.mstrYDataFile.c_str(),arma_ascii) ){
      Rf_error("Loading Y failed!" );
  }
  if( !model.Time.load((char*)model.mstrTimeDataFile.c_str(),arma_ascii) ){
      Rf_error("Loading Time failed!" );
  }
  if( !model.W.load((char*)model.mstrWDataFile.c_str(),arma_ascii) ){
      Rf_error("Loading W failed!" );
  }
  if( !model.Xpk.load((char*)model.mstrXpkDataFile.c_str(),arma_ascii) ){
      Rf_error("Loading Xpk failed!" );
  }
  if( !model.Xpg.load((char*)model.mstrXpgDataFile.c_str(),arma_ascii) ){
       Rf_error("Loading Xpg failed!" );
  }
  // Preprocess
  model.XtXpk    = model.Xpk.t() * model.Xpk;
  model.XtXpg    = model.Xpg.t() * model.Xpg;
  model.invXtXpg = inv( model.XtXpg );
  model.invXtXpk = inv( model.XtXpk );
  Rprintf( "Loading Data Finished\n" );
  return true;
}


#include "resultpkpg.hpp"
#include "modelpkpg.hpp"
/***************************************************************************************/
/* FILES & CONSTANTS *******************************************************************/
/***************************************************************************************/
//#define file_input1   "Data/PK_init.txt"
// #define file_Fit      "./PostPkPg/Fit.txt"
// #define file_Gsnp     "./PostPkPg/Gsnp.txt"
// #define file_Gtheta   "./PostPkPg/Gtheta.txt"
// #define file_Grho     "./PostPkPg/Grho.txt"
// #define file_Rho      "./PostPkPg/Rho.txt"
// #define file_Theta    "./PostPKPg/Theta.txt"
// #define file_Sigma    "./PostPKPg/Error.txt"
// #define file_BetaPk   "./PostPKPg/BetaPk.txt"
// #define file_BetaPg   "./PostPKPg/BetaPg.txt"

// #define last_Fit      "./SavePkPg/Fit.txt"
// #define last_Gsnp     "./SavePkPg/Gsnp.txt"
// #define last_Gtheta   "./SavePkPg/Gtheta.txt"
// #define last_Grho     "./SavePkPg/Grho.txt"
// #define last_Rho      "./SavePkPg/Rho.txt"
// #define last_Theta    "./SavePKPg/Theta.txt"
// #define last_Sigma    "./SavePKPg/Sigma.txt"
// #define last_BetaPk   "./SavePKPg/BetaPk.txt"
// #define last_BetaPg   "./SavePKPg/BetaPg.txt"

/*************************************************************************
 * Constructor                                                           *
 *************************************************************************/
PkPgResult::PkPgResult(){

    // postFitfile.open( file_Fit.c_str() );
    // if(postFitfile.fail()){
    //     //std::cout << "Failed to create " << file_Fit << endl;
    //     Rf_error("Failed to create %s", file_Fit);
    // }
    // postGsnpfile.open(file_Gsnp.c_str());
    // if(postGsnpfile.fail()){
    //     //std::cout << "Failed to create " << file_Gsnp << endl;
    //     Rf_error("Failed to create %s", file_Gsnp);
    // }
    // postGthetafile.open(file_Gtheta.c_str());
    // if(postGthetafile.fail()){
    //     //std::cout << "Failed to create " << file_Gtheta << endl;
    //     Rf_error("Failed to create %s", file_Gtheta );
    // }
    // postGrhofile.open(file_Grho.c_str());
    // if(postGrhofile.fail()){
    //     //std::cout << "Failed to create " << file_Grho << endl;
    //     Rf_error("Failed to create %s", file_Grho );
    // }
    // postRhofile.open(file_Rho.c_str());
    // if(postRhofile.fail()){
    //     //std::cout << "Failed to create " << file_Rho << endl;
    //     Rf_error("Failed to create %s", file_Rho );
    // }
    // postThetafile.open(file_Theta.c_str());
    // if(postThetafile.fail()){
    //     //std::cout << "Failed to create " << file_Theta << endl;
    //     Rf_error("Failed to create %s", file_Theta );
    // }
    // postSigmafile.open(file_Sigma.c_str());
    // if(postSigmafile.fail()){
    //     //std::cout << "Failed to create " << file_Sigma << endl;
    //     Rf_error("Failed to create %s", file_Sigma );
    // }
    // postBetaPkfile.open(file_BetaPk.c_str());
    // if(postBetaPkfile.fail()){
    //     //std::cout << "Failed to create " << file_BetaPk << endl;
    //     Rf_error("Failed to create %s", file_BetaPk );
    // }
    // postBetaPgfile.open(file_BetaPg.c_str());
    // if(postBetaPgfile.fail()){
    //     //std::cout << "Failed to create " << file_BetaPg << endl;
    //     Rf_error("Failed to create %s", file_BetaPg);
    // }
}

/*************************************************************************
 * Destructor                                                            *
 *************************************************************************/
PkPgResult::~PkPgResult(void)
{
    //delete [] W;
    //delete [] K;
    postFitfile.close();
    // postGsnpfile.close();
    // postGthetafile.close();
    // postGrhofile.close();
    postZfile.close();
    postOsnpfile.close();
    postOthetafile.close();
    postRhofile.close();
    postThetafile.close();
    postSigmafile.close();
    postBetaPkfile.close();
    postEHRTimefile.close();
}

int PkPgResult::Init( PkPgModel& model )
{

    postFitfile.open( model.file_Fit.c_str() );
    if(postFitfile.fail()){
        Rf_error( "Failed to create %s", model.file_Fit.c_str() );
    }
    postOsnpfile.open( model.file_Osnp.c_str() );
    if(postOsnpfile.fail()){
        Rf_error( "Failed to create %s", model.file_Osnp.c_str() );
    }
    postOthetafile.open( model.file_Otheta.c_str() );
    if(postOthetafile.fail()){
        Rf_error( "Failed to create %s", model.file_Otheta.c_str() );
    }
    // postGrhofile.open( model.file_Grho.c_str() );
    // if(postGrhofile.fail()){
    //     Rf_error( "Failed to create %s", model.file_Grho.c_str() );
    // }
    postZfile.open( model.file_Z.c_str() );
    if(postZfile.fail()){
        Rf_error( "Failed to create %s", model.file_Z.c_str() );
    }
    postRhofile.open( model.file_Rho.c_str() );
    if(postRhofile.fail()){
        Rf_error( "Failed to create %s", model.file_Rho.c_str() );
    }
    postThetafile.open( model.file_Theta.c_str() );
    if(postThetafile.fail()){
        Rf_error( "Failed to create %s", model.file_Theta.c_str() );
    }
    postSigmafile.open( model.file_Sigma.c_str() );
    if(postSigmafile.fail()){
        Rf_error( "Failed to create %s", model.file_Sigma.c_str() );
    }
    postBetaPkfile.open( model.file_BetaPk.c_str() );
    if(postBetaPkfile.fail()){
        Rf_error( "Failed to create %s", model.file_BetaPk.c_str() );
    }
    postBetaPgfile.open( model.file_BetaPg.c_str() );
    if(postBetaPgfile.fail()){
        Rf_error( "Failed to create %s", model.file_BetaPg.c_str() );
    }
    postEHRTimefile.open( model.file_EHRTime.c_str() );
    if(postEHRTimefile.fail()){
        Rf_error( "Failed to create %s", model.file_EHRTime.c_str() );
    }

//     return 0;
// }

// int PkPgResult::LoadInits( PkPgModel& model )
// {
  Rprintf("Loading Initial Values\n");
  if( model.loadFit == true  ){
    if( !Fit.load((char*)model.initFitfile.c_str(),arma_ascii) ){
      Rf_error("Loading Fit failed!" );
    }
  }
  if( model.loadTheta == true  ){
    if( !Theta.load((char*)model.initThetafile.c_str(),arma_ascii) ){
      Rf_error("Loading Theta failed!" );
    }
  }
  if( model.loadZ == true  ){
    if( !Z.load((char*)model.initZfile.c_str(),arma_ascii) ){
      Rf_error("Loading Z failed!" );
    }
  }
  if( model.loadBetaPk == true  ){
    if( !BetaPk.load((char*)model.initBetaPkfile.c_str(),arma_ascii) ){
      Rf_error("Loading BetaPk failed!" );
    }
  }
  if( model.loadBetaPg == true  ){
    if( !BetaPg.load((char*)model.initBetaPgfile.c_str(),arma_ascii) ){
      Rf_error("Loading BetaPg failed!" );
    }
  }

  if( model.loadRho == true  ){
    if( !Rho.load((char*)model.initRhofile.c_str(),arma_ascii) ){
      Rf_error("Loading Rho failed!" );
    }
  }
    //std::cout << model.loadSigma<<std::endl;
  if( model.loadSigma == true ){
    if( !Sigma.load((char*)model.initSigmafile.c_str(),arma_ascii) ){
      Rf_error("Loading Sigma failed!" );
    }
  }
  if( model.loadEHRTime == true ){
    if( !EHRTime.load((char*)model.initEHRTimefile.c_str(),arma_ascii) ){
      Rf_error("Loading EHRT failed!" );
    }
  }
  EHRTime.print("EHRT");
  if( model.loadOsnp == true  ){
    if( !Osnp.load((char*)model.initOsnpfile.c_str(),arma_ascii) ){
      Rf_error("Loading Osnp failed!" );
    }
  }
  Ssnp = inv( Osnp );

  if( model.loadOtheta == true  ){
    if( !Otheta.load((char*)model.initOthetafile.c_str(),arma_ascii) ){
      Rf_error("Loading Otheta failed!" );
    }
  }
  Stheta = inv( Otheta );

  if( model.loadSigmaTheta0 == true  ){
    if( !SigmaTheta0.load((char*)model.initSigmaTheta0file.c_str(),arma_ascii) ){
      Rf_error("Loading SigmaTheta0 failed!" );
    }
  }

 
  //Other parameters 
  delta = 1;
  iota  = 1;
  Rhoinit = 1;
  //Rprintf( "Initializing Summary\n" );
  logTheta_bar_old = zeros<mat>( model.N,  model.V );// Theta_bar(l-1)
  logTheta         = log( Theta + 0.01 );
  logTheta_bar     = logTheta; //zeros<mat>( model.N,  model.V ); // Theta_bar(0)
  logTheta_old     = zeros<mat>( model.N,  model.V );
  logTheta_MLE     = logTheta;
  Sigma_theta   = zeros<cube>( model.V, model.V, model.N );// Sigma_theta(0)
  SS_theta      = zeros<cube>( model.V, model.V, model.N );// Sigma_theta(0)
  for ( size_t i = 0; i < model.N; i++ ){ 
      SS_theta.slice(i)    = logTheta.row(i).t() * logTheta.row(i);//zeros<mat>( model.V, model.V);
       if( model.loadSigmaTheta0 == true  ){ 
         Sigma_theta.slice(i) = SigmaTheta0;
       } else { 
         Sigma_theta.slice(i) = 0.01* eye<mat>( model.V, model.V ); 
      }
      //Sigma_theta.slice(i).print("sigmatheta");
  }
  XBetaPg     = model.Xpg * BetaPg;
  XBetaPk     = model.Xpk * BetaPk;
  ZRho        =         Z * Rho;
  mTheta      = XBetaPk + ZRho;   // mean of prior Theta Xpk*BetaPk + Z*Rho
  logTheta_Resid = logTheta - mTheta;
  //Rprintf("Finished Initializing Summary\n");



  Rprintf("Finished Loading Initial Values\n");
  return( 0 );
}

// int PkPgResult::InitializeSummary( PkPgModel& model ){
//     Rprintf("Initializing Summary\n");
//     Theta_bar_old = zeros<mat>( model.N,  model.V );// Theta_bar(l-1)
//     Theta_bar     = Theta; // Theta_bar(0)
//     Theta_old     = Theta;
//     Sigma_theta   = zeros<cube>( model.V, model.V, model.N );// Sigma_theta(0)
//     SS_theta      = zeros<cube>( model.V, model.V, model.N );// Sigma_theta(0)
//     for ( int i = 0; i < model.N; i++ ){ 
//         SS_theta.slice(i) = Theta.row(i).t() * Theta.row(i);
//     }
//     for ( int i =0; i<model.N; i++){ 
//         Sigma_theta.slice(i) = eye<mat>( model.V, model.V ); 
//     }
//     XBetaPg     = model.Xpg * BetaPg;
//     XBetaPk     = model.Xpk * BetaPk;
//     ZRho        =         Z * Rho;
//     mTheta      = XBetaPk + ZRho;   // mean of prior Theta Xpk*BetaPk + Z*Rho
//     Theta_Resid = Theta - mTheta;
//     Rprintf("Finished Initializing Summary\n");
//     return(0);
// }
void PkPgResult::update_mth( PkPgModel& model ){
    //XBetaPg     = model.Xpg * BetaPg;
    mTheta      = XBetaPk + ZRho;  
}
/************************************************************************
 * Functions for logging the predictive density during each step of mcmc *
 ***********************************************************************/
bool PkPgResult::SaveDraws( PkPgModel& model ){
  //cube FitMat = Fit;
  //FitMat.reshape( model.N*model.T, model.K, 1 );
  for( size_t i=0; i< Fit.n_slices; i++ ){
    Fit.slice(i).raw_print( postFitfile );
  }   
   
   //   Gsnp.raw_print( postGsnpfile   );
   // Gtheta.raw_print( postGthetafile );
   //   Grho.raw_print( postGrhofile   );
        Z.raw_print( postZfile      );
      Rho.raw_print( postRhofile    );
    Theta.raw_print( postThetafile  );
    Sigma.raw_print( postSigmafile  );
   BetaPk.raw_print( postBetaPkfile );
   BetaPg.raw_print( postBetaPgfile );
     Osnp.raw_print( postOsnpfile   );
   Otheta.raw_print( postOthetafile );
   trans(EHRTime).raw_print( postEHRTimefile );
   
 //    Fit.print( postFitfile    );
 //   Gsnp.print( postGsnpfile   );
 // Gtheta.print( postGthetafile );
 //   Grho.print( postGrhofile   );
 //    Rho.print( postRhofile    );
 //  Theta.print( postThetafile  );
 //  Sigma.print( postSigmafile  );
 // BetaPk.print( postBetaPkfile );
 // BetaPg.print( postBetaPgfile );
  return( true );
}
/*********************************************************************
 * functions for saving all parameter values needed to continue mcmc
 ********************************************************************/
bool PkPgResult::SaveFinal( PkPgModel& model ) {
       Fit.save( model.last_Fit   , arma_ascii );
    //   Gsnp.save( model.last_Gsnp  , arma_ascii );
    // Gtheta.save( model.last_Gtheta, arma_ascii );
    //   Grho.save( model.last_Grho  , arma_ascii );
       Osnp.save( model.last_Osnp  , arma_ascii );
     Otheta.save( model.last_Otheta, arma_ascii );
          Z.save( model.last_Z      , arma_ascii );
        Rho.save( model.last_Rho   , arma_ascii );
      Theta.save( model.last_Theta , arma_ascii );
      Sigma.save( model.last_Sigma , arma_ascii );
     BetaPk.save( model.last_BetaPk, arma_ascii );
     BetaPg.save( model.last_BetaPg, arma_ascii );
     EHRTime.save( model.last_EHRTime, arma_ascii );
    return( true );
}



int PkPgResult::Reload( PkPgModel& model )
{
    postFitfile.open( model.file_Fit.c_str() );
    if(postFitfile.fail()){
        Rf_error( "Failed to create %s", model.file_Fit.c_str() );
    }
    postOsnpfile.open( model.file_Osnp.c_str() );
    if(postOsnpfile.fail()){
        Rf_error( "Failed to create %s", model.file_Osnp.c_str() );
    }
    postOthetafile.open( model.file_Otheta.c_str() );
    if(postOthetafile.fail()){
        Rf_error( "Failed to create %s", model.file_Otheta.c_str() );
    }
    postZfile.open( model.file_Z.c_str() );
    if(postZfile.fail()){
        Rf_error( "Failed to create %s", model.file_Z.c_str() );
    }
    postRhofile.open( model.file_Rho.c_str() );
    if(postRhofile.fail()){
        Rf_error( "Failed to create %s", model.file_Rho.c_str() );
    }
    postThetafile.open( model.file_Theta.c_str() );
    if(postThetafile.fail()){
        Rf_error( "Failed to create %s", model.file_Theta.c_str() );
    }
    postSigmafile.open( model.file_Sigma.c_str() );
    if(postSigmafile.fail()){
        Rf_error( "Failed to create %s", model.file_Sigma.c_str() );
    }
    postBetaPkfile.open( model.file_BetaPk.c_str() );
    if(postBetaPkfile.fail()){
        Rf_error( "Failed to create %s", model.file_BetaPk.c_str() );
    }
    postBetaPgfile.open( model.file_BetaPg.c_str() );
    if(postBetaPgfile.fail()){
        Rf_error( "Failed to create %s", model.file_BetaPg.c_str() );
    }
//     return 0;
// }

// int PkPgResult::LoadInits( PkPgModel& model )
// {
  Rprintf("Loading Saved Values\n");
  //if( model.loadFit == true  ){
  if( !Fit.load((char*)model.last_Fit.c_str(),arma_ascii) ){
      Rf_error("Loading Fit failed!" );
  }
  //}
  //if( model.loadTheta == true  ){
  if( !Theta.load((char*)model.last_Theta.c_str(),arma_ascii) ){
      Rf_error("Loading Theta failed!" );
  }
  //}
  //if( model.loadZ == true  ){
    if( !Z.load((char*)model.last_Z.c_str(),arma_ascii) ){
      Rf_error("Loading Z failed!" );
    }
  // }
  // if( model.loadBetaPk == true  ){
    if( !BetaPk.load((char*)model.last_BetaPk.c_str(),arma_ascii) ){
      Rf_error("Loading BetaPk failed!" );
    }
  // }
  // if( model.loadBetaPg == true  ){
    if( !BetaPg.load((char*)model.last_BetaPg.c_str(),arma_ascii) ){
      Rf_error("Loading BetaPg failed!" );
    }
  // }

  // if( model.loadRho == true  ){
    if( !Rho.load((char*)model.last_Rho.c_str(),arma_ascii) ){
      Rf_error("Loading Rho failed!" );
    }
  // }
    //std::cout << model.loadSigma<<std::endl;
  // if( model.loadSigma == true ){
    if( !Sigma.load((char*)model.last_Sigma.c_str(),arma_ascii) ){
      Rf_error("Loading Sigma failed!" );
    }
  // }
  // if( model.loadOsnp == true  ){
    if( !Osnp.load((char*)model.last_Osnp.c_str(),arma_ascii) ){
      Rf_error("Loading Osnp failed!" );
    }
  // }
  Ssnp = inv( Osnp );

  // if( model.loadOtheta == true  ){
    if( !Otheta.load((char*)model.last_Otheta.c_str(),arma_ascii) ){
      Rf_error("Loading Otheta failed!" );
    }
  // }
  Stheta = inv( Otheta );

  //if( model.loadSigmaTheta0 == true  ){
  // if( !SigmaTheta0.load((char*)model.last_SigmaTheta0.c_str(),arma_ascii) ){
  //   Rf_error("Loading SigmaTheta0 failed!" );
  // }
  //}

 
  //Other parameters 
  delta = 1;
  iota  = 1;
  Rhoinit = 1;
  //Rprintf( "Initializing Summary\n" );
  logTheta_bar_old = zeros<mat>( model.N,  model.V );// Theta_bar(l-1)
  logTheta         = log( Theta + 0.00001 );
  logTheta_bar     = logTheta; //zeros<mat>( model.N,  model.V ); // Theta_bar(0)
  logTheta_old     = zeros<mat>( model.N,  model.V );
  logTheta_MLE     = logTheta;
  Sigma_theta   = zeros<cube>( model.V, model.V, model.N );// Sigma_theta(0)
  SS_theta      = zeros<cube>( model.V, model.V, model.N );// Sigma_theta(0)
  for ( size_t i = 0; i < model.N; i++ ){ 
      SS_theta.slice(i)    = logTheta.row(i).t() * logTheta.row(i);//zeros<mat>( model.V, model.V);
       if( model.loadSigmaTheta0 == true  ){ 
         Sigma_theta.slice(i) = SigmaTheta0;
       } else { 
         Sigma_theta.slice(i) = 0.01* eye<mat>( model.V, model.V ); 
      }
      //Sigma_theta.slice(i).print("sigmatheta");
  }
  XBetaPg     = model.Xpg * BetaPg;
  XBetaPk     = model.Xpk * BetaPk;
  ZRho        =         Z * Rho;
  mTheta      = XBetaPk + ZRho;   // mean of prior Theta Xpk*BetaPk + Z*Rho
  logTheta_Resid = logTheta - mTheta;
  //Rprintf("Finished Initializing Summary\n");



  Rprintf("Finished Loading Saved Values\n");
  return( 0 );
}


// bool PkPgResult::SaveFitDraw(string FileName){  
//     return true;
// }
// bool PkPgResult::SaveGsnpDraw(string FileName){
//     return true;
// }
// bool PkPgResult::SaveGthetaDraw(string FileName){
//     return true;
// }
// bool PkPgResult::SaveGrhoDraw(string FileName){
//     return true;
// }
// bool PkPgResult::SaveRhoDraw(string FileName){
//     return true;
// }
// bool PkPgResult::SaveThetaDraw(string FileName){
//     return true;
// }
// bool PkPgResult::SaveErrorDraw(string FileName){   
//     return true;
// }
// bool PkPgResult::SaveBthetaDraw(string FileName){
//     return true;
// }

// bool PkPgResult::SavePDraw(){
//   for(int j=0;j<J;j++)
//     {
//       for(int t=0;t<T;t++)
//     {
//       postpfile << q[j]*p[j][t] << "\t";
//     }
//       postpfile<<endl;
//     }
//   return true;
// }

// bool PkPgResult::SaveMuDraw(){
//   for (int i = 0; i < J*T; i++) {
//     for (int j = 0; j < D; j++) {
//       postmufile << mu[i][j] << "\t";
//     }
//     postmufile << endl;
//   }
//   return true;
// }

// bool PkPgResult::SaveSigmaDraw(){
//   for (int i = 0; i < J*T; i++) {
//     for (int j = 0; j < D; j++) {
//       for (int k = 0; k <= j; k++) {
//     postSigmafile << Sigma[i][j][k] << "\t";
//       }
//       for (int k = j+1; k < D; k++) {
//     postSigmafile << Sigma[i][k][j] << "\t";
//       }
//     }
//     postSigmafile << endl;
//   }
//   return true;
// }

// bool PkPgResult::SaveQDraw(){
//   for(int j=0;j<J;j++)
//     {
//       postqfile << q[j] << "\t";
//     }
//   postqfile << endl;
//   return true;
// }

// bool PkPgResult::SaveMDraw(){
//   for (int i = 0; i < J; i++) {
//     for (int j = 0; j < D; j++) {
//       postmfile << m[i][j] << "\t";
//     }
//     postmfile << endl;
//   }
//   return true;
// }

// bool PkPgResult::SavePhiDraw()
// {
//   for (int i = 0; i < J; i++) {
//     for (int j = 0; j < D; j++) {
//       for (int k = 0; k <= j; k++) {
//     postPhifile << Phi[i][j][k] << "\t";
//       }
//       for (int k = j+1; k < D; k++) {
//     postPhifile << Phi[i][k][j] << "\t";
//       }
//     }
//     postPhifile << endl;
//   }
//   return true;
// }



// bool PkPgResult::SaveFit( string FileName ){
//   ofstream theFile( FileName.c_str() );
//   if( theFile.fail() ){
//         std::cout << "Failed to create file " << FileName.c_str()  << endl;
//         exit(1);
//   }
//   //theFile << Fit << endl;
//   theFile.close();
//   return true;
// }
// bool PkPgResult::SaveGsnp(string FileName){
//   ofstream theFile( FileName.c_str() );
//   if( theFile.fail() ){
//         std::cout << "Failed to create file " << FileName.c_str()  << endl;
//         exit(1);
//   }
//   //theFile << Fit << endl;
//   theFile.close();
//   return true;
// }
// bool PkPgResult::SaveGtheta(string FileName){
//   ofstream theFile( FileName.c_str() );
//   if( theFile.fail() ){
//         std::cout << "Failed to create file " << FileName.c_str()  << endl;
//         exit(1);
//   }
//   //theFile << Fit << endl;
//   theFile.close();
//   return true;
// }
// bool PkPgResult::SaveGrho(string FileName){
//   ofstream theFile( FileName.c_str() );
//   if( theFile.fail() ){
//         std::cout << "Failed to create file " << FileName.c_str()  << endl;
//         exit(1);
//   }
//   //theFile << Fit << endl;
//   theFile.close();
//   return true;
// }
// bool PkPgResult::SaveRho(string FileName){
//   ofstream theFile( FileName.c_str() );
//   if( theFile.fail() ){
//         std::cout << "Failed to create file " << FileName.c_str()  << endl;
//         exit(1);
//   }
//   //theFile << Fit << endl;
//   theFile.close();
//   return true;
// }
// bool PkPgResult::SaveTheta(string FileName){
//   ofstream theFile( FileName.c_str() );
//   if( theFile.fail() ){
//         std::cout << "Failed to create file " << FileName.c_str()  << endl;
//         exit(1);
//   }
//   //theFile << Fit << endl;
//   theFile.close();
//   return true;
// }
// bool PkPgResult::SaveError(string FileName){
//   ofstream theFile( FileName.c_str() );
//   if( theFile.fail() ){
//         std::cout << "Failed to create file " << FileName.c_str()  << endl;
//         exit(1);
//   }
//   //theFile << Fit << endl;
//   theFile.close();
//   return true;
// }
// bool PkPgResult::SaveBtheta(string FileName){
//   ofstream theFile( FileName.c_str() );
//   if( theFile.fail() ){
//         std::cout << "Failed to create file " << FileName.c_str()  << endl;
//         exit(1);
//   }
//   //theFile << Fit << endl;
//   theFile.close();
//   return true;
// }

// bool PkPgResult::SaveAlpha0(string FileName) {
//   ofstream theFile(FileName.c_str());
//   if (theFile.fail()) {
//         std::cout << "Failed to create file " << FileName.c_str()  << endl;
//         exit(1);
//     }
//   theFile << alpha0 << endl;
//   theFile.close();
//   return true;
// }

// bool PkPgResult::SaveAlpha(string FileName) {
//   ofstream theFile(FileName.c_str());
//   if (theFile.fail()) {
//         std::cout << "Failed to create file " << FileName.c_str()  << endl;
//         exit(1);
//     }
//   for(int j=0;j<J;j++)
//     {
//       theFile << alpha[j] << "\t";
//     }
//   theFile << endl;
//   theFile.close();
//   return true;
// }

// bool PkPgResult::SaveW(string FileName) {
//     ofstream theFile(FileName.c_str());
//     if (theFile.fail()) {
//         std::cout << "Failed to create file " << FileName.c_str()  << endl;
//         exit(1);
//     }
//     for (int i = 0; i < N; i++) {
//       theFile << W[i] +1<< endl;

//     }
//     //theFile << endl;
//     theFile.close();
//     return true;
// }

// bool PkPgResult::SaveM(string FileName) {
//     ofstream theFile(FileName.c_str());
//     if (theFile.fail()) {
//         std::cout << "Failed to create file " << FileName.c_str()  << endl;
//         exit(1);
//     }
//     //    int D = mu[0].Ncols();
//     for (int i = 0; i < J; i++) {
//         for (int j = 0; j < D; j++) {
//             theFile << m[i][j] << "\t";
//         }
//         theFile << endl;
//     }
//     theFile.close();
//     return true;
// }

// bool PkPgResult::SavePhi(string FileName) {
//     ofstream theFile(FileName.c_str());
//     if (theFile.fail()) {
//         std::cout << "Failed to create file " << FileName.c_str()  << endl;
//         exit(1);
//     }
//     //    int D = mu[0].Ncols();
//     for (int i = 0; i < J; i++) {
//         for (int j = 0; j < D; j++) {
//             for (int k = 0; k <= j; k++) {
//                 theFile << Phi[i][j][k] << "\t";
//             }
//             for (int k = j+1; k < D; k++) {
//                 theFile << Phi[i][k][j] << "\t";
//             }
//         }
//         theFile << endl;
//     }
//     theFile.close();
//     return true;
// }

// bool PkPgResult::SaveQ(string FileName) {
//   ofstream theFile(FileName.c_str());
//   if (theFile.fail()) {
//         std::cout << "Failed to create file " << FileName.c_str()  << endl;
//         exit(1);
//     }
//   for(int j=0;j<J;j++)
//     {
//       theFile << q[j] << "\t";
//     }
//   theFile << endl;
//   theFile.close();
//   return true;
// }

// bool PkPgResult::SaveqV(string FileName) {
//   ofstream theFile(FileName.c_str());
//   if (theFile.fail()) {
//         std::cout << "Failed to create file " << FileName.c_str()  << endl;
//         exit(1);
//     }
//   for(int j=0;j<J;j++)
//     {
//       theFile << qV[j] << "\t";
//     }
//   theFile << endl;
//   theFile.close();
//   return true;
// }

// bool PkPgResult::SaveK(string FileName) {
//     ofstream theFile(FileName.c_str());
//     if (theFile.fail()) {
//         std::cout << "Failed to create file " << FileName.c_str()  << endl;
//         exit(1);
//     }
//     for (int i = 0; i < N; i++) {
//       theFile << K[i]+1 << endl;
//     }
//     //theFile << endl;
//     theFile.close();
//     return true;
// }


// bool PkPgResult::SaveMu(string FileName) {
//     ofstream theFile(FileName.c_str());
//     if (theFile.fail()) {
//         std::cout << "Failed to create file " << FileName.c_str()  << endl;
//         exit(1);
//     }
//     //    int D = mu[0].Ncols();
//     for (int i = 0; i < J*T; i++) {
//         for (int j = 0; j < D; j++) {
//             theFile << mu[i][j] << "\t";
//         }
//         theFile << endl;
//     }
//     theFile.close();
//     return true;
// }

// bool PkPgResult::SaveSigma(string FileName) {
//     ofstream theFile(FileName.c_str());
//     if (theFile.fail()) {
//         std::cout << "Failed to create file " << FileName.c_str()  << endl;
//         exit(1);
//     }
//     //    int D = mu[0].Ncols();
//     for (int i = 0; i < J*T; i++) {
//         for (int j = 0; j < D; j++) {
//             for (int k = 0; k <= j; k++) {
//                 theFile << Sigma[i][j][k] << "\t";
//             }
//             for (int k = j+1; k < D; k++) {
//                 theFile << Sigma[i][k][j] << "\t";
//             }
//         }
//         theFile << endl;
//     }
//     theFile.close();
//     return true;
// }

// bool PkPgResult::SaveP(string FileName) {
//   ofstream theFile(FileName.c_str());
//   if (theFile.fail()) {
//         std::cout << "Failed to create file " << FileName.c_str()  << endl;
//         exit(1);
//   }
//   //  int blah=1;
//   for(int j=0;j<J;j++)
//     {
//       for(int t=0;t<T;t++)
//     {
//       theFile << p[j][t] << "\t";
//     }
//       theFile<<endl;
//     }
//   theFile.close();
//   return true;
// }

// bool PkPgResult::SavepV(string FileName) {
//   ofstream theFile(FileName.c_str());
//   if (theFile.fail()) {
//         std::cout << "Failed to create file " << FileName.c_str()  << endl;
//         exit(1);
//   }
//   //  int blah=1;
//   for(int j=0;j<J;j++)
//     {
//       for(int t=0;t<T;t++)
//     {
//       theFile << pV[j][t] << "\t";
//     }
//       theFile<<endl;
//     }
//   theFile.close();
//   return true;
// }


/*******************************************************
 * functions for maintaining and logging running means
 *******************************************************/
// bool PkPgResult::SaveBar() {
//     SaveXMbar("postxmbar.txt");
//     SaveXMubar("postxmubar.txt");
//     return true;
// }

// bool PkPgResult::SaveXMbar(string FileName)
// {
//     ofstream theFile(FileName.c_str());
//     if (theFile.fail()) {
//         std::cout << "Failed to create file " << FileName.c_str()  << endl;
//         exit(1);
//     }
//     //    int D = mu[0].Ncols();
//     for (int i = 0; i < N; i++) {
//         for (int j = 0; j < D; j++) {
//             theFile << xmbar[i][j]/nmcits << "\t";
//         }
//         theFile << endl;
//     }
//     theFile.close();
//     return true;

// }
// bool PkPgResult::SaveXMubar(string FileName)
// {
//     ofstream theFile(FileName.c_str());
//     if (theFile.fail()) {
//         std::cout << "Failed to create file " << FileName.c_str()  << endl;
//         exit(1);
//     }
//     //    int D = mu[0].Ncols();
//     for (int i = 0; i < N; i++) {
//         for (int j = 0; j < D; j++) {
//             theFile << xmubar[i][j]/nmcits << "\t";
//         }
//         theFile << endl;
//     }
//     theFile.close();
//     return true;

// }

// void PkPgResult::UpdateMeans()
// {
//   int j;
//   RowVector foo(D);foo=0;
//   SymmetricMatrix foo2(D);foo2=0;
//   RowVector foo3(T);foo3=0;
//   // if this is the first update, zero out the parameters
//   if(nmcits==0)
//     {
//       for(j=0;j<N;j++)
//     {
//       xmbar.push_back(foo);
//       xmubar.push_back(foo);
//     }
//     }
  
//   for(j=0;j<N;j++)
//     {
//       xmbar[j]+=m[W[j]];
//       xmubar[j]+=mu[GetIndex(W[j],K[j])];
//     }
//   nmcits++;
// }


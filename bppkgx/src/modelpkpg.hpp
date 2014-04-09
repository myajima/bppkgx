#ifndef _MODEL_PkPg_HPP
#define _MODEL_PkPg_HPP
/* 
 * model.hpp 
 */
#include <RcppArmadillo.h>
// #include <map>
#include <cstring>
#include <string>
//#include <cstdio>
// #include <cstdlib>
// #include <fstream>
// #include <iostream>
// extern "C" {
// #include <RcppArmadillo.h>
#include <R.h>
#include <Rinternals.h>
#include "priorpkpg.hpp"
//#include "resultpkpg.hpp"

using namespace Rcpp;
using namespace arma;
using namespace std;

//class MCMCPkPg;
//class PkPgResult;
//class PkPgPrior;
class PkPgModel
{
public:
    PkPgModel();
    ~PkPgModel(void);
    bool Load(string FileName);
    bool Save(string FileName);
    bool InitializeSummary();
    bool Print();
public:
    size_t N;    // No. subjects
    size_t T;    // No. of sampling times
    size_t K;    // No. of observed compartments
    size_t Q;    // No. of SNPs
    size_t Ppk;  // No. of PK covariates
    size_t Ppg;  // No. of SNPs covariates
    size_t V;    // No. PK parameters
    double epsilon;
    // Data 
    mat Time;
    imat W;

    // Covariates
    mat Xpg;        // U: matrix( Q, Ppg )
    mat XtXpg;      // t(U)%*%U
    //mat XXtpg;      // t(U)%*%U
    mat invXtXpg;   // solve( t(U)%*%U )
    //mat invXXtpg;   // solve( t(U)%*%U )

    mat Xpk;        // X: matrix( Q, Ppg )
    mat XtXpk;      // t(X)%*%X
    mat invXtXpk;
    //mat XXtpk;      // t(X)%*%X
    cube Y;
    field<uvec> NegIdxQ;
    field<uvec> NegIdxV;

    // 
    SEXP PkModel;
    SEXP env;
    SEXP cluster;    

    PkPgPrior prior; // prior parameters
//  string mstralgorithm;
//  string mstrDataFile;
    
//  int mnN;
//  int mnD;
//  int mnJ;
//  int mnT;
    
//  RowVector mdm0;
//  double mdlambda0;
//  double mdphi0;
//  double mdnu0;
//  double mdgamma;
//  double mdnu;
//  double mde0;
//  double mdf0;
//  double mdee;
//  double mdff;
//  double mdaa;

//  int mnCompSwaps;
//  int mnBurnin;
//  int mnIter;
    
//  int mnPrintout;
//  double mnErrorPerTol;         // EM: Percent Difference stopping rule for diffence in posterior.

//  int mnSeed;

    // vec ToRowVector(string str);
    static double ToDouble(string str);
    static int ToInt(string str);
    static string ToLower(string str);
    static string ToString(int value);
    static string ToString(double value);
    static string trim(string s);

//  // load sets of parameters
    string mstrYDataFile;
    string mstrTimeDataFile;
    string mstrWDataFile;
    string mstrXpgDataFile;
    string mstrXpkDataFile;

    // Flags to determine if variable is loaded from the file
    // priors
    bool loadA0;
    bool loadalpha0;
    bool loadbetapg0;
    bool loadbetapk0;
    bool loadBpg0;
    bool loadBpk0;

    // initial values
    bool loadFit;
    bool loadBetaPk;
    bool loadBetaPg;
    bool loadAlpha;
    bool loadOsnp;
    bool loadOtheta;
    bool loadSigmaTheta0;
    bool loadRho;
    bool loadSigma;
    bool loadTau;
    bool loadTheta;
    bool loadZ;
    bool loadEHRTime;
    // bool loadGtheta;
    // bool loadGsnp;
    // bool loadGrho;


    bool saveFit   ; 
    bool saveOsnp  ; 
    bool saveOtheta; 
    //bool saveGrho  ; 
    bool saveZ   ; 
    bool saveRho   ; 
    bool saveTheta ; 
    bool saveSigma ; 
    bool saveBetaPk; 
    bool saveBetaPg; 
    bool saveEHRTime;
    // file names 
    // priors
    string A0file;
    string alpha0file;
    string betapg0file;
    string betapk0file;
    string Bpg0file;
    string Bpk0file;
    // initial values
    string initFitfile;
    string initAlphafile;
    string initBetaPgfile;
    string initBetaPkfile;
    string initOsnpfile;
    string initOthetafile;
    string initSigmaTheta0file;
    string initRhofile;
    string initSigmafile;
    string initTaufile;
    string initThetafile;
    string initZfile;
    string initEHRTimefile;
    // string initGsnpfile;
    // string initGthetafile;
    // string initGrhofile;

    // files to save the current values
    string file_Fit   ; 
    // string file_Gsnp  ; 
    // string file_Gtheta; 
    // string file_Grho  ; 
    string file_Z     ; 
    string file_Osnp  ; 
    string file_Otheta; 
    string file_Rho   ; 
    string file_Theta ; 
    string file_Sigma ; 
    string file_BetaPk; 
    string file_BetaPg; 
    string file_EHRTime;
    // files to save the last value
    string last_Fit   ; 
    string last_Z     ; 
    // string last_Gsnp  ; 
    // string last_Gtheta; 
    // string last_Grho  ;
    string last_Osnp  ; 
    string last_Otheta;  
    string last_Rho   ; 
    string last_Theta ; 
    string last_Sigma ; 
    string last_BetaPk; 
    string last_BetaPg; 
    string last_EHRTime;

};
#endif //_MODEL_PkPg_HPP


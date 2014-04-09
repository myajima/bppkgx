#ifndef _RESULT_PkPg_HPP
#define _RESULT_PkPg_HPP

#include <RcppArmadillo.h>
#include <vector>
#include <R.h>
#include <Rinternals.h>
//#include "modelpkpg.hpp"
using namespace std;
using namespace Rcpp;
using namespace arma;
class PkPgModel;
class MCMCPkPg;


class PkPgResult
{
    friend class MCMCPkPg;
    public:
        //PkPgResult(int nclusters, int ncomponents, int npoints,int dimension);
        PkPgResult();
        ~PkPgResult(void); 
        int Init( PkPgModel& model );
        int Reload( PkPgModel& model );
        //int LoadInits( PkPgModel& model );
    public:
        cube Fit;
        mat logFit;
        // imat Gtheta;   //= mat(   V, V );
        // imat Gsnp;     //= mat(   Q, Q );
        // imat Grho;     //= mat(   Q, V );
        // Likelihood parameters
        mat Sigma;     //= mat(   N, K )

        // Pk parameters
        mat logTheta;
        mat Theta;    //= mat(   N, V );      
        mat EHRTime;
        mat logTheta_old;    //= mat(   N, V ); 
        mat Alpha;    //= mat( Ppk, V );  
        mat BetaPk;   //= mat( Ppk, V );  
        mat Rho;      //= mat( Q,   V);
        double delta;
        double iota; 
        int Rhoinit; // initialization for the 1st round
        double Rhos2;
        vec Rhot2;
        double Rholb; 
        mat Stheta;   // Sigma_theta: mat( V,   V);
        mat Otheta;   // Ometa_theta= inv( Stheta ): mat( V,   V);
        mat LAMBDAtheta;
        //vec Tau; // vec( N,   0.01 );

        // Pg parameters
        mat Z;        //
        mat mZ;       //
        mat BetaPg;     //= mat( Ppg, Q ); 
        mat Ssnp;       // Sigma_snp: mat( Q,   Q);
        mat Osnp;       // Omega_snp = inv( Ssnp ): mat( Q,   Q);

        mat LAMBDAsnp;
        mat PsZ;   // Posterior inverse Covariance Matrix 
                   // Osnp + Tau(i) Rho.t()*Otheta* Rho
                   // need updating before updating Z
        mat PmZ;   // Posterior mean of Z
                   // PsZ^-1[Osnp BetaPg * Xpg(i,) + Tau(i) Rho.t()*Otheta* Rho hat Z]
        mat Zhat;  // Tau(i) Rho.t()*Otheta* Rho


        mat Pth;      //= mat( V,   V);
        mat th_Pz;       //= mat( Q,   Q);
        mat th_h_th;     //= vec( V,   0.01 );
        mat th_h_z;      //= vec( Q,   0.01 );


        // Summary Statst (check if some of them can be removed)
        // is this a result or the model... 
        mat logTheta_bar;
        mat logTheta_bar_old;
        cube Sigma_theta;
        cube SS_theta;
        mat SigmaTheta0;
        mat ZtZ;      // t(Z) %*% Z      ( update when Z is updated )
        mat XBetaPg;  // t(Xpg)%*%Betapg ( update when BetaPg is updated )
        mat XBetaPk;  // t(Xpk)%*%Betapk ( update when BetaPk is updated )
        mat ZRho;     // t(Z)%*%Rho      ( update when Rho and Z are updated )
        mat mTheta;   // mean of prior Theta Xpk*BetaPk + Z*Rho
        mat logTheta_Resid; // logTheta - XbetaPk - ZRho
        mat D_post_pk;
        mat D_post_pg;
        mat logTheta_MLE;
    //int InitializeSummary(PkPgModel& model );

    // log the current state of the mcmc
    bool SaveDraws( PkPgModel& model );
    bool SaveFinal( PkPgModel& model );

    /* functions to log the *bar values above */
    void update_mth(PkPgModel& model);
    void UpdateMeans();
    bool SaveXMbar(string FileName);
    bool SaveXMubar(string FileName);
    bool SaveBar();
    void OpenPostFiles();


    /* output files for logging predictive distribution*/
    ofstream postFitfile   ;
    // ofstream postGsnpfile  ;
    // ofstream postGthetafile;
    // ofstream postGrhofile  ;
    ofstream postOsnpfile  ;
    ofstream postOthetafile;
    ofstream postRhofile   ;
    ofstream postZfile     ;
    ofstream postThetafile ;
    ofstream postSigmafile ;
    ofstream postBetaPkfile;
    ofstream postBetaPgfile;
    ofstream postEHRTimefile;
};

#endif //_RESULT_PkPg_HPP

    // int isEM;
    // double alpha0; 
    // int r;
    
    // /* mean values of component/cluster means/covariances */
    // int nmcits; // number of mcmc iterations upon which the means are calculated

    // vector<RowVector> xmbar;
    // vector<RowVector> xmubar;
    // RowVector alpha; //J by 1

    // #if defined(CDP_TBB)
    //     concurrent_vector<RowVector> m; //J by D        
    //     concurrent_vector<SymmetricMatrix> Phi;    //J by D * D
        
    //     concurrent_vector<RowVector> mu; //(J by T) by D
    //     concurrent_vector<SymmetricMatrix> Sigma; //(J by T) by (D by D)
        
    //     concurrent_vector<RowVector> p; //J by T
    //     concurrent_vector<RowVector> pV; // J by T
    
    //     concurrent_vector<RowVector> eta; // J by T 
        
    //     //work variable
    //     concurrent_vector<UpperTriangularMatrix> Phi_T_i;   //J by D * D
    //     concurrent_vector<double> Phi_log_det;
    //     concurrent_vector<LowerTriangularMatrix> L_i; //(J by T) by (D by D)
    //     concurrent_vector<double> Sigma_log_det;
    // #else
    //     vector<RowVector> m; //J by D        
    //     vector<SymmetricMatrix> Phi;    //J by D * D
        
    //     vector<RowVector> mu; //(J by T) by D
    //     vector<SymmetricMatrix> Sigma; //(J by T) by (D by D)
        
    //     vector<RowVector> p; //J by T
    //     vector<RowVector> pV; // J by T
    
    //     vector<RowVector> eta; // J by T 
        
    //     //work variable
    //     vector<UpperTriangularMatrix> Phi_T_i;   //J by D * D */
    //      vector<double> Phi_log_det; 
    //     vector<LowerTriangularMatrix> L_i; //(J by T) by (D by D)
    //     vector<double> Sigma_log_det;
    // #endif
    
    // bool loadZ; // are we tracking reference classification?
    // int* Z; //N by 1 -- Classification Vector
    // int* refZ; //N by 1 -- Reference Classification Vector
    // int* refZobs; //T by 1 -- Number of obs in each Classification Component
        
    // RowVector q;    // J by 1
    // RowVector qV; //J by 1

    // int* W; //N by 1
    // int* K; //N by 1
    
    // inline int GetIndex(int j, int t) {
    //     return (j * T + t);
    // };
    // int J;
    // int T;
    // int N;
    // int D;


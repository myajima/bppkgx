#include <RcppArmadillo.h>
#include <R.h>
#include <Rinternals.h>
#include <cstring>
#include <string>

#include "utilarmadillo.hpp"
#include "utilR.hpp"

// #include "SymChol.Rcpp"
// #include "cliques.Rcpp"
#include "resultpkpg.hpp"
#include "modelpkpg.hpp"
#include "gwish.hpp"
//class PkPgResult;
//class PkPgModel;
class MCMCPkPg
{
  public:
    MCMCPkPg();
    ~MCMCPkPg(void);
    int nCHAIN;  // No. of separate mcmc chains              
    int nSAMPLE; // No. of mcmc simulations 
    int nBURN;   // No. of discarded samples
    int nTHIN;   // No. thinning factor 
    int flgSNP;
    int cIter;
    mat acceptTheta;
    //int initialize( int n_chain, int n_sample, int n_burn, int n_thin );
    bool LoadData( PkPgModel& model, PkPgResult& Result ) ;

    //bool Initialize( PkPgModel& model,PkPgResult& Result )
    void sample( PkPgModel& model, PkPgResult& Result );
    // SNPs ---------------------------------------------------------------------
    int sampleBetaPg( PkPgModel& model, PkPgResult& Result );
    int sampleOsnp( PkPgModel& model, PkPgResult& Result );
    int sampleZ( PkPgModel& model, PkPgResult& Result ); 
    // PK -----------------------------------------------------------------------
    int sampleSigma( PkPgModel& model, PkPgResult& Result );
    int sampleTheta( PkPgModel& model, PkPgResult& Result ); //dyn_eval, rho
    void update_fit( PkPgModel& model, PkPgResult& Result ); //dyn_eval, rho
    void update_fit_parallel( PkPgModel& model, PkPgResult& Result ); //dyn_eval, rho
    void update_theta_bar( PkPgModel& model, PkPgResult& Result );
    void update_sigma_theta( PkPgModel& model, PkPgResult& Result );

    //int sampleAlpha( PkPgModel& model,PkPgResult& Result);
    int sampleBetaPk( PkPgModel& model, PkPgResult& Result );
    int sampleRho( PkPgModel& model, PkPgResult& Result );
    //update_mth();
    int sampleOtheta( PkPgModel& model, PkPgResult& Result );
//  double sampleAlpha(double* V, double e, double f, MTRand& mt);
//  void sampleMuSigma(Matrix& x, int n, double nu, double gamma,RowVector& m, SymmetricMatrix& Phi, RowVector& PostMu, SymmetricMatrix& PostSigma,UpperTriangularMatrix& uti, double& logdet, MTRand& mt);
//  void sampleMuSigma(vector<int>& indexes, double nu, double gamma,RowVector& m, SymmetricMatrix& Phi, RowVector& PostMu,SymmetricMatrix& PostSigma,UpperTriangularMatrix& uti, double& logdet, MTRand& mt);
//  void sampleP(vector<int>& p, int n, double gamma, int T, RowVector& postp, RowVector& postV, MTRand& mt);
//  void sampleP(int* p, int n, double gamma, int T, RowVector& postp, R
//  void cov(Matrix& x, RowVector& mu,int mul, SymmetricMatrix& result);
// owVector& postV, MTRand& mt);
//  void sampleM(int n, double gamma, vector<RowVector>& mu, vector<SymmetricMatrix>& Sigma, RowVector& m0,
//      SymmetricMatrix& Phi0, RowVector& newm, MTRand& mt);
//  void samplePhi(int n, double nu, vector<SymmetricMatrix>& Sigma, double nu0, SymmetricMatrix& Lambda0,
//      SymmetricMatrix& newphi, MTRand& mt);
//  //MTRand mt;
//  SpecialFunctions2 msf;
//  //Matrix mX;
//  double** mX;
//  CDPPrior prior;
//  double precalculate;


};

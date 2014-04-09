#include <RcppArmadillo.h>
#include <R.h>
#include <Rinternals.h>
#include <cstring>
#include <string>
#include "resultpkpg.hpp"
#include "modelpkpg.hpp"
#include "utilarmadillo.hpp"
#include "utilR.hpp"

class MCMCPkPg
{
public:
    MCMCPkPg();
    ~MCMCPkPg(void);
    int nCHAIN;    // No. of separate mcmc chains              
    int nSAMPLE;   // No. of mcmc simulations 
    int nBURN;   // No. of discarded samples
    int nTHIN;   // No. thinning factor 
    //int initialize( int n_chain, int n_sample, int n_burn, int n_thin );
    void sample( PkPgModel& model, PkPgResult& Result );   
    // SNPs ---------------------------------------------------------------------
    int sampleBz();
    int sampleSz();
    int sampleZ(); 
    // PK -----------------------------------------------------------------------
    int sampleSe();
    int sampletheta1(); //dyn_eval, rho
    void update_fit( PkPgModel& model, PkPgResult& Result ); //dyn_eval, rho
    int sampleBth();
    int sampleRho();
    //update_mth();
    int sampleSth();
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

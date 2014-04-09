/**************************************************
 * model.cpp: implementation of the Model class.
 **************************************************/
#include "modelpkpg.hpp"
//#include "resultpkpg.hpp"
//#include "priorpkpg.hpp"

 //Construction/Destruction
PkPgModel::PkPgModel(){
   // Flags to determine if variable is loaded from the file
    // priors
   loadA0     = false;
   loadalpha0 = false;
   loadbetapg0= false;
   loadbetapk0= false;
   loadBpg0   = false;
   loadBpk0   = false;

   //itial values
   loadFit   = false;
   loadBetaPk= false;
   loadBetaPg= false;
   loadAlpha = false;
   loadOsnp  = false;
   loadOtheta= false;
   loadRho   = false;
   loadSigma = false;
   loadTau   = false;
   loadTheta = false;
   loadZ     = false;
   loadOtheta= false;
   loadOsnp  = false;
   loadSigmaTheta0=false;
   loadEHRTime = false;
   //loadGrho  = false;   

   //constant
   epsilon = 0.0001;
}
// Model::Model(int n, int t, int k, int q, int ppk, int ppg, int v){
//     N   = n;    // No. subjects
//     T   = t;    // No. of sampling times
//     K   = k;    // No. of observed compartments
//     Q   = q;    // No. of SNPs
//     Ppk = ppk;  // No. of PK covariates
//     Ppg = ppg;  // No. of SNPs covariates
//     V   = v;    // No. PK parameters
// }

PkPgModel::~PkPgModel(){

}
bool PkPgModel::InitializeSummary(){
  Rprintf("in init model");
  XtXpk = Xpk.t() * Xpk;
  invXtXpk = inv( XtXpk );
  XtXpg = Xpg.t() * Xpg;
  invXtXpg = inv( XtXpg );
  // XXtpk = Xpk.t() * Xpk;
  // XXtpg = Xpg.t() * Xpg;
  // invXXtpg = inv( XXtpg );
    Rprintf("out init model");
  return false;
}
string PkPgModel::ToLower(string str)
{
    //new implementation for GNU
    char *newstr = strdup(str.c_str());
    int i = 0;
    while (newstr[i] != '\0') {
        newstr[i] = tolower(newstr[i]);
        i++;
    }
    return newstr;
}

int PkPgModel::ToInt(string str)
{
    return atoi(str.c_str());
}

double PkPgModel::ToDouble(string str)
{
    return atof(str.c_str());
}

string PkPgModel::ToString(int value) {
    char  buffer[10]; 
    sprintf(buffer, "%d", value);
    string result(buffer);
    return result;
    //cout << result.c_str() << endl;
}


string PkPgModel::ToString(double value) {
    char  buffer[10]; 
    sprintf(buffer, "%4.1f", value);
    string result(buffer);
    return result;
}

string PkPgModel::trim(string s)
{
    if (s.empty()) return s;
    string ret;
    for (int i = 0; i < (int)s.length();  i++) {
        if (!isspace(s.at(i)) && !iscntrl(s.at(i)))
            ret.append(s.substr(i,1));
    }
    return ret;
}
bool PkPgModel::Load(string FileName){
  int BufferSize = 4096;
  char* theLine = new char[BufferSize];
  ifstream theFile(FileName.c_str());
  if (theFile.fail()) {
    Rf_error( "Failed to open the model file!" );
    return false;
  }
  
  int nLineCount = 0;
  while (!theFile.eof()) {
    theFile.getline(theLine, BufferSize);
    nLineCount++;
    string theline(theLine);
    string Name(""), Value("");
    theline = trim(theline);            //so space is not allowed, should be improved later
    if (theline.length() && (theline.c_str()[0] != '#'))
      {
    int pos = 0;
    if ((pos = theline.find("=")) != -1) {
      Name = theline.substr(0, pos);
      Value = theline.substr(pos + 1);
    }
    if (Name == "" && Value == "") {
    } else if (Name == "" || Value == "") {
      Rf_error( "Invalid Parameter: %s", theLine );
      return false;
    } else {
      string name = ToLower(Name);
      string value = ToLower(Value);
      if( name == "yfile" ) {
        mstrYDataFile = Value;
      } else if( name == "timefile" ) {
        mstrTimeDataFile = Value;
      } else if( name == "wfile" ) {
        mstrWDataFile = Value;
      } else if( name == "xpkfile" ) {
        mstrXpkDataFile = Value;
      } else if( name == "xpgfile" ) {
        mstrXpgDataFile = Value;
      } else if( name=="alpha0file" ){
        loadalpha0=true;
        alpha0file = Value;
      } else if( name=="a0file" ){
        loadA0=true;
        A0file = Value;
      } else if(name=="betapk0file"){
        loadbetapk0=true;
        betapk0file=Value;
      } else if(name=="bpk0file"){
        loadBpk0=true;
        Bpk0file=Value;
      } else if(name=="betapg0file"){
        loadbetapg0=true;
        betapg0file=Value;
      } else if(name=="bpg0file"){
        loadBpg0=true;
        Bpg0file=Value;
      } else if(name=="initfitfile"){
        loadFit=true;
        initFitfile=Value;    
      } else if(name=="initthetafile"){
        loadTheta=true;
        initThetafile=Value;    
      } else if(name=="initzfile"){
        loadZ=true;
        initZfile=Value;
      } else if(name=="initalphafile"){
        loadAlpha=true;
        initAlphafile=Value;
      } else if(name=="initbetapkfile"){
        loadBetaPk=true;
        initBetaPkfile=Value;
      } else if(name=="initbetapgfile"){
        loadBetaPg=true;
        initBetaPgfile=Value;
      } else if(name=="initrhofile"){
        loadRho=true;
        initRhofile=Value;
      } else if(name=="inittaufile"){
        loadTau=true;
        initTaufile=Value;
      } else if(name=="initsigmafile"){
        loadSigma=true;
        initSigmafile=Value;
      } else if(name=="initosnpfile"){
        loadOsnp=true;
        initOsnpfile=Value;
      } else if(name=="initothetafile"){
        loadOtheta=true;
        initOthetafile=Value;
      // } else if(name=="initgsnpfile"){
      //   loadGsnp=true;
      //   initGsnpfile=Value;
      // } else if(name=="initgthetafile"){
      //   loadGtheta=true;
      //   initGthetafile=Value;
      // } else if(name=="initgrhofile"){
      //   loadGrho=true;
      //   initGrhofile=Value;
      } else if(name=="initsigmatheta0file"){
        loadSigmaTheta0=true;
        initSigmaTheta0file=Value;
      } else if(name=="initehrtimefile"){
        loadEHRTime=true;
        initEHRTimefile=Value;
      } else if( name == "savefitfile"){
        saveFit     = true; 
        file_Fit    = Value; 
      } else if( name == "saveosnpfile"){
        saveOsnp    = true; 
        file_Osnp   = Value; 
      } else if( name == "saveothetafile"){
        saveOtheta  = true; 
        file_Otheta = Value; 
      // } else if( name == "saveGrhofile"){
      //   saveGrho    = true; 
      //   file_Grho   = Value; 
      } else if( name == "savezfile"){
        saveZ     = true; 
        file_Z    = Value; 
      } else if( name == "saverhofile"){
        saveRho     = true; 
        file_Rho    = Value; 
      } else if( name == "savethetafile"){
        saveTheta   = true; 
        file_Theta  = Value; 
      } else if( name == "savesigmafile"){
        saveSigma   = true; 
        file_Sigma  = Value; 
      } else if( name == "savebetapkfile"){
        saveBetaPk  = true; 
        file_BetaPk = Value; 
      } else if( name == "savebetapgfile"){
        saveBetaPg  = true; 
        file_BetaPg = Value; 
      } else if( name == "saveehrtimefile"){
        saveEHRTime  = true; 
        file_EHRTime = Value; 
      } else if( name == "lastfitfile"){
        last_Fit    = Value; 
      } else if( name == "lastosnpfile"){
        last_Osnp   = Value; 
      } else if( name == "lastothetafile"){
        last_Otheta = Value; 
      // } else if( name == "lastGsnpfile"){
      //   last_Gsnp   = Value; 
      // } else if( name == "lastGthetafile"){
      //   last_Gtheta = Value; 
      // } else if( name == "lastGrhofile"){
      //   last_Grho   = Value; 
      } else if( name == "lastzfile"){
        last_Z      = Value; 
      } else if( name == "lastrhofile"){
        last_Rho    = Value; 
      } else if( name == "lastthetafile"){
        last_Theta  = Value; 
      } else if( name == "lastsigmafile"){
        last_Sigma  = Value; 
      } else if( name == "lastbetapkfile"){
        last_BetaPk = Value; 
      } else if( name == "lastbetapgfile"){
        last_BetaPg = Value; 
      } else if( name == "lastehrtimefile"){
        last_EHRTime = Value; 
      }  else if( name == "n" ){
        N = ToInt( value );
      } else if( name == "t" ){
        T = ToInt( value );
      } else if( name == "q" ){
        Q = ToInt( value );
      } else if( name == "k" ){
        K = ToInt(value);
      } else if( name == "v" ){
        V = ToInt( value );
      } else if( name == "ppk" ){
        Ppk = ToInt(value);
      } else if( name == "ppg" ){
        Ppg = ToInt(value);
      } 
      else {
        Rf_warning( "Unknown Parameter: %s", theLine );
        //return false;
      }
      
    }
      }
  }
  //if(setm0){mdm0 = RowVector(mnD);mdm0=0;}
  delete[] theLine;
  theLine = NULL;
  // if (mstrDataFile == "") {
  //   std::cout << endl << "Warning: Text file "  << FileName.c_str() << " might come from a different platform" << endl;
  //   std::cout << "Suggestion: Use dos2unix/unis2dos or similar tools to convert the file first" << endl;
  //   return false;
  // }
  //    Print();
  return true;
}

bool PkPgModel::Print(){
Rprintf("## Current Settings ##\n");
Rprintf("");
Rprintf("#data section");
Rprintf("No. subjects N = %d \n", N);
Rprintf("No. of sampling times T = %d \n", T);
Rprintf("No. of observed compartments K = %d \n", K);
Rprintf("No. of SNPs Q = %d \n", Q);
Rprintf("No. of PK covariates Ppk = %d \n", Ppk);
Rprintf("No. of SNPs covariates Ppg= %d \n", Ppg);
Rprintf("No. PK parameters V = %d \n", V);

//  cout << "N = " << mnN << endl;
//  cout << "D = " << mnD << endl;
//  if (mstrDataFile != "") {
//      cout << "datafile = " << mstrDataFile.c_str() << endl;
//  } else {
//      cout << "#datafile = " << endl;
//  }   
//  cout << endl;

//  cout << "#prior section" << endl;
//  cout << "J = " << mnJ << endl;
//  cout << "T = " << mnT << endl;
//  cout << "m0 = ";
//  for(int i=0;i<mnD;i++)
//    cout << mdm0[i] << " ";
//  cout << endl;
//  cout << "nu0 = " << mdnu0 << endl;
//  cout << "gamma = " << mdgamma << endl;
//  cout << "phi0 = " << mdphi0 << endl;
//  cout << "lambda0 = " << mdlambda0 << endl;
//  cout << "e0 = " << mde0 << endl;
//  cout << "f0 = " << mdf0 << endl;
//  cout << "nu = " << mdnu << endl;
//  cout << "ee = " << mdee << endl;
//  cout << "ff = " << mdff << endl;
//  cout << endl;
    
//  cout << "#MCMC section" << endl;
//  cout << "burnin = " << mnBurnin << endl;
//  cout << "iter = " << mnIter << endl;
//  cout << "seed = " << mnSeed << endl;
//  cout << "swaps = " << mnCompSwaps << endl << endl;

//  cout << "#Turn on/off individual sampling steps"<<endl;
//  cout << "samplem = " << samplem << endl;
//  cout << "samplePhi = " << samplePhi << endl;
//  cout << "samplew = " << samplew << endl;
//  cout << "sampleq = " << sampleq << endl;
//  cout << "samplealpha0 = " << samplealpha0 << endl;
//  cout << "samplemu = " << samplemu << endl;
//  cout << "sampleSigma = " << sampleSigma << endl;
//  cout << "samplek = " << samplek << endl;
//  cout << "samplep = " << samplep << endl;
//  cout << "samplealpha = " << samplealpha << endl << endl;
        
//  cout << "#Load initial MCMC values from files"<<endl;
//  if (loadAlpha0){
//    cout << "Alpha0file = " << Alpha0file.c_str() << endl;
//  } else {
//    cout << "#Alpha0file = " << endl;
//  }   

//  if (loadM){
//    cout << "Mfile = " << Mfile.c_str() << endl;
//  } else {
//    cout << "#Mfile = " << endl;
//  }
  if( loadTheta == true  ){
    Rprintf("Thetafile = %s \n", initThetafile.c_str());
  } else{
    Rprintf("#Thetafile =  \n");  
  }
  if( loadZ == true  ){
    Rprintf("Zfile = %s \n", initZfile.c_str());
  } else{
    Rprintf("#Zfile =  \n");  
  }
  if( loadAlpha == true  ){
    Rprintf("Alphafile = %s \n", initAlphafile.c_str());
  } else{
    Rprintf("#Alphafile =  \n");  
  }

  if( loadBetaPk == true  ){
    Rprintf("BetaPkfile = %s \n", initBetaPkfile.c_str());
  } else{
    Rprintf("#BetaPkfile =  \n");  
  }
  if( loadBetaPg == true  ){
    Rprintf("BetaPgfile = %s \n", initBetaPgfile.c_str());
  } else{
    Rprintf("#BetaPgfile =  \n");  
  }
  if( loadRho == true  ){
    Rprintf("Rhofile = %s \n", initRhofile.c_str());
  } else{
    Rprintf("#Rhofile =  \n");  
  }

  if( loadSigma == true  ){
    Rprintf("Sigmafile = %s \n", initSigmafile.c_str());
  } else{
    Rprintf("#Sigmafile =  \n");   
  }
  if( loadTau == true  ){
    Rprintf("Taufile = %s \n", initTaufile.c_str());
  } else{
    Rprintf("#Taufile =  \n");   
  }

  if( loadOsnp == true  ){
    Rprintf("Osnpfile = %s \n", initOsnpfile.c_str());
  } else{
    Rprintf("#Osnpfile =  \n");   
  }
    if( loadOtheta == true  ){
    Rprintf("Othetafile = %s \n", initOthetafile.c_str());
  } else {
    Rprintf("#Othetafile =  \n");   
  }

  // if( loadGtheta == true  ){
  //   Rprintf("Gthetafile = %s \n", initGthetafile.c_str());
  // } else{
  //   Rprintf("#Gthetafile =  \n");   
  // }
  // if( loadGsnp == true  ){
  //   Rprintf("Gsnpfile = %s \n", initGsnpfile.c_str());
  // } else{
  //   Rprintf("#Gsnpfile =  \n");   
  // }
  // if( loadGrho == true  ){
  //   Rprintf("Grhofile = %s \n", initGrhofile.c_str());
  // } else{
  //   Rprintf("#Grhofile =  \n");   
  // }
//  if (loadPhi){
//    cout << "Phifile = " << Phifile.c_str() << endl;
//  } else {
//    cout << "#Phifile = " << endl;
//  }   
//  if (loadW){
//    cout << "Wfile = " << Wfile.c_str() << endl;
//  } else {
//    cout << "#Wfile = " << endl;
//  }   
//  if (loadQ){
//    cout << "Qfile = " << Qfile.c_str() << endl;
//  } else {
//    cout << "#Qfile = " << endl;
//  }   
//  if (loadqV){
//    cout << "qVfile = " << qVfile.c_str() << endl;
//  } else {
//    cout << "#qVfile = " << endl;
//  }   
//  if (loadAlpha){
//    cout << "Alphafile = " << Alphafile.c_str() << endl;
//  } else {
//    cout << "#Alphafile = " << endl;
//  }   
//  if (loadMu){
//    cout << "Mufile = " << Mufile.c_str() << endl;
//  } else {
//    cout << "#Mufile = " << endl;
//  }   
//  if (loadSigma){
//    cout << "Sigmafile = " << Sigmafile.c_str() << endl;
//  } else {
//    cout << "#Sigmafile = " << endl;
//  }   
//  if (loadK){
//    cout << "Kfile = " << Kfile.c_str() << endl;
//  } else {
//    cout << "#Kfile = " << endl;
//  }   
//  if (loadP){
//    cout << "Pfile = " << Pfile.c_str() << endl;
//  } else {
//    cout << "#Pfile = " << endl;
//  }   
//  if (loadpV){
//    cout << "pVfile = " << pVfile.c_str() << endl;
//  } else {
//    cout << "#pVfile = " << endl;
//  }   
    
    
//          cout << endl;

//  cout.flush();
  return true;
}
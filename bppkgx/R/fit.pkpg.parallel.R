# --------------------------------------------------------------------------------------------- #
# Masanao Yajima                                                                           #
# --------------------------------------------------------------------------------------------- #
# FUNCTION :  fit.pkpg()                                                                        #
# Purpose  :  MCMC sampler for the PK/PGx model of Telesca/Rosner/Muller                        # 
# Version  :  0.0.1                                                                             #
# Date     :  April 2012                                                                        # 
# --------------------------------------------------------------------------------------------- #

fit.pkpgpall <- function( Conc,  
                      SNP,  
                      Time, 
                      Xpk     = NULL, 
                      Xpg     = NULL,  
                      Model,  
                      Jacobian= NULL, 
                      Theta0,
                      PK_V, 
                      Dose0,   
                      gen.init= TRUE,
                      n_chain = 1,
                      n_mcmc  = 1000,  
                      n_burn  = 1000,  
                      n_lag   = 1,
                      workdir = ".", 
                      prior   = 'T')
{
  # CREATE PRIOR FILE -----------------------------------------------------------------------
  # ---- To add an interactive window for prior specification

  # CREATE OUTPUT DIRECTORY -----------------------------------------------------------------
  files <- list.files();
  dir   <- files == "PostPkPg";
  if( length( files[dir] ) == 0 ){ 
    cat("Directory PostPkPg created: \n")
    dir.create("PostPkPg"); 
  }
  # CREATE SAVE DIRECTORY -----------------------------------------------------------------
  if( length( files[files == "SavePkPg"] ) == 0 ){ 
    cat("Directory SavePkPg created: \n")
    dir.create("SavePkPg"); 
  }
  # GET DATA SUMMARIES ----------------------------------------------------------------------
  # SNPs --------------------- #
  Nsnp <- nrow(SNP);           # Number of subject in SNPs matrix
  Q    <- ncol(SNP) - 1;       # Number of plymorphisms
  # Concentrations ----------- #
  dConc<- dim( Conc )
  T    <- dConc[1] #unique( table( Conc$ID ) )
  K    <- dConc[2] #ncol(Conc) - 2;      # Number of compartments
  Nconc<- dConc[3] # length( unique( Conc$ID ) ) # Number of subjects in Concentration file #max(Conc$ID)  
  #TOT  <- nrow(Conc);             # Total number of concentration values
  # Sampling Times ----------- #
  #Time    <- TOT/ Nconc;    # Number of sampling times

  # Covariates --------------- # 
  # need to tidy up a little more using R standard functions.
  if( is.null( Xpk ) ){
    Xpk  <- rep(1, Nconc);  
    Ppk  <- 1;  
  } else {
    if( !all( Xpk[,1] == 1 ) ){
      Xpk <- cbind( rep( 1, nrow( Xpk ) ), Xpk );
    } 
    Ppk <- ncol(Xpk);  
  }
  if( is.null( Xpg ) ){
    Xpg  <- rep(1, Nconc);
    Ppg  <- 1;
  } else {
    if( !all( Xpg[,1] == 1 ) ){
      Xpg <- cbind( rep( 1, nrow( Xpg ) ), Xpg );
    }
    Ppg  <- ncol(Xpg);  
  }
  # Initial Dosage and Theta Values ---------------------------------------------------------
  # PKparam  <- names(Theta0)[3:ncol(Theta0)]; 
  Th0    <- log( Theta0  + 0.0001);
  if( PK_V != ncol(Theta0) ) stop("PK_V and the column of Theta0 should match");

  # Priors
  #alpha0  = rep( 0, PK_V )
  #A0      = diag( PK_V )
  betapk0 = matrix(0, Ppk, PK_V )
  Bpk0    = diag(Ppk)
  betapg0 = matrix(0, Ppg, Q )
  Bpg0    = diag(Ppk)
  # Initial values ( think about how to let users specify this )

  if( gen.init == TRUE ){
    library( mvtnorm )
    iFit   <- array( rnorm( Nconc * K * T ), c( T, K, Nconc) )
    iGsnp  <- diag( Q )
    iGtheta<- diag( PK_V )
    iGrho  <- matrix( 0, Q, PK_V )
    iZ     <- matrix( rnorm( Nsnp * Q ), Nsnp, Q )
    #iAlpha <- matrix( rmvnorm( 1, alpha0, A0 ), PK_V, 1 )
    iBetaPk<- matrix( rmvnorm( PK_V, rep(0,Ppk), Bpk0 ), Ppk, PK_V )
    iBetaPg<- matrix( rmvnorm(    Q, rep(0,Ppg), Bpg0 ), Ppg, Q  )
    iRho   <- matrix( rnorm(Q * PK_V),   Q, PK_V )
    iTau   <- matrix( rmvnorm( 1, alpha0, A0 ), PK_V, 1 )
    iSigma <- NULL
    iGsnp  <- diag( Q )
    iGtheta<- diag( PK_V )
    iGrho  <- matrix( 0, Q, PK_V )
    iOsnp  <- diag( Q )
    iOtheta<- diag( PK_V )

  }
  # VERIFY DIMENTIONAL COMPATIBILTY ---------------------------------------------------------
  #if(T%%round(T) != 0) stop("Check that all subject have the same number of sampling times");    
  #if( length(Time) != 1 ) stop( "Check that all subject have the same number of sampling times");  
  if( Nconc != Nsnp )  stop( gettextf("N in Conc is %d, while N in SNP is %d", Nconc, Nsnp), domain = NA);    
  # DEFINE data structures ------------------------------------------------------------------
  #c1 <- log( as.matrix( Conc )[,3:ncol( Conc )] + 0.1);
  c1 <- log( Conc + 0.1)
  s1 <- as.matrix( SNP )[,2:ncol( SNP )];    
  # DEFINE LSODA C ENVELOPE PARALLEL---------------------------------------------------------
  library( deSolve );
  library( parallel );
  D0   <- Dose0;                                                    # D0 dose
  time <- Time;#time <- t( matrix( Conc[,2], ncol = Nconc ) );                               # time matrix 
  .runlsoda<- function( i ){ #, param ) {
     param <- c( param[i,], D0[i] ); 
        t  <- c( 0.0, time[i,] );
      fit  <- lsoda( rep( 0, 7 ), t, Model, parms=param, rtol = 0.5, atol = 0.5, jacfunc= Jacobian )
      return( fit / param[length(param) - 1 ] )
  }
  n.cluster = 2
  cl <- makeCluster( n.cluster, methods = FALSE )
  dyn.eval <- function( param, cl ){
      #parms  <- c( param, D0[i] ); 
      #t      <- c( 0.0, time[i,] );
      print(param)
      clusterExport( cl, c("param") )
      #out    <- lsoda( rep( 0, 7 ), t, Model, parms=parms, rtol = 0.5, atol = 0.5, jacfunc= Jacobian );
      #print(out)
      res    <- parSapply( cl, seq_len(86), .runlsoda )
      out1   <- matrix( res[-c( 1:16, 17:32, 49:64, 113:128 ), ], ,4 )

      #as.matrix(out)[2:nrow(out),c(2,4,5,6)];
      #as.double(out1/param[length(param)-1]);

      #need to fix this so that array of ( T, K, N ) gets passed back.
      as.double( t( out1 ) );
  }
  dyn.env<- new.env( dyn.eval )
  assign("time", time, dyn.env)
  assign("D0", D0, dyn.env)
  assign("Model", Model, dyn.env)
  assign("Jacobian", Jacobian, dyn.env)
  assign(".runlsoda", .runlsoda, dyn.env)
  assign("lsoda", lsoda, dyn.env)
  clusterExport( cl, c("time", "D0", "Model", "Jacobian", "lsoda", ".runlsoda"), envir = dyn.env )
  ###########################################################################################
  info = Sys.info()
  
  if(info[1]=="Windows"){
      dirsep="\\"
  } else {
      dirsep="/"
  }
  ##############################################################################
  # Parameters
  ##############################################################################
  filename = "parameters.txt";8
  cat("#data section \n", file = filename )
  cat("N   = ",Nconc, "\n", file = filename, append = TRUE)
  cat("T   = ",    T, "\n", file = filename, append = TRUE)
  cat("K   = ",    K, "\n", file = filename, append = TRUE)
  cat("Q   = ",    Q, "\n", file = filename, append = TRUE)
  cat("V   = ", PK_V, "\n", file = filename, append = TRUE)
  cat("Ppk = ",  Ppk, "\n", file = filename, append = TRUE)
  cat("Ppg = ",  Ppg, "\n", file = filename, append = TRUE)
  ##############################################################################
  # data
  ##############################################################################
  if( !( is.null( c1 ) ) ){
      Yfilename = paste( c( workdir, "SavePkPg", "Y.txt" ), collapse = dirsep )
      #write.table( t(c1), file = Yfilename,
      #             sep = "\t", row.names = FALSE, col.names = FALSE )
      # array( as.double(t(c1)), c(T, K, Nconc ) )
      save_arma_cube( c1, filename = Yfilename )
      #save_arma_cube( t(c1), filename = Yfilename )
      cat("Yfile = ", Yfilename, "\n",file = filename, append = TRUE )
  }
  if( !( is.null( Time ) ) ){
      Timefilename = paste( c( workdir, "SavePkPg", "Time.txt" ), collapse = dirsep )
      #write.table( t(c1), file = Yfilename,
      #             sep = "\t", row.names = FALSE, col.names = FALSE )
      # array( as.double(t(c1)), c(T, K, Nconc ) )
      save_arma_matrix( Time, filename = Timefilename )
      #save_arma_cube( t(c1), filename = Yfilename )
      cat("Timefile = ", Timefilename, "\n",file = filename, append = TRUE )
  }
  if( !( is.null( s1 ) ) ){
      Wfilename = paste( c( workdir, "SavePkPg", "W.txt" ), collapse = dirsep )
      #write.table( s1, file = Wfilename,
      #             sep = "\t", row.names = FALSE, col.names = FALSE )
      save_arma_matrix( s1, filename = Wfilename )
      cat("Wfile = ", Wfilename, "\n",file = filename, append = TRUE )
  }
  if( !( is.null( Xpk ) ) ){
      Xpkfilename = paste( c( workdir, "SavePkPg", "Xpk.txt" ), collapse = dirsep )
      #write.table( Xpk, file = Xpkfilename,
      #             sep = "\t", row.names = FALSE, col.names = FALSE )
      save_arma_matrix( Xpk, filename = Xpkfilename )
      cat("Xpkfile = ", Xpkfilename, "\n",file = filename, append = TRUE )
  }
  if( !( is.null( Xpg ) ) ){
      Xpgfilename = paste( c( workdir, "SavePkPg", "Xpg.txt" ), collapse = dirsep )
      #write.table( Xpg, file = Xpgfilename,
      #             sep = "\t", row.names = FALSE, col.names = FALSE )
      save_arma_matrix( Xpg, filename = Xpgfilename )
      cat("Xpgfile = ", Xpgfilename, "\n",file = filename, append = TRUE )
  }
  ##############################################################################
  # prior
  ##############################################################################
  cat("#prior section \n",  file = filename, append = TRUE)
  # if( !( is.null( alpha0 ) ) ){
  #     alpha0filename = paste( c( workdir, "SavePkPg", "prioralpha0.txt" ), collapse = dirsep )
  #     #write.table( alpha0, file = alpha0filename,
  #     #             sep = "\t", row.names = FALSE, col.names = FALSE )
  #     save_arma_matrix( alpha0, filename = alpha0filename )
  #     cat("alpha0file = ", alpha0filename, "\n",file = filename, append = TRUE )
  # }
  # if( !is.null( A0 ) ){
  #     A0filename = paste( c( workdir, "SavePkPg", "priorA0.txt" ), collapse = dirsep )
  #     #write.table( A0, file = A0filename,
  #     #             sep = "\t", row.names = FALSE, col.names = FALSE )
  #     save_arma_matrix( A0, filename = alpha0filename )
  #     cat("A0file = ", A0filename, "\n",file = filename, append = TRUE )
  # }
  if( !is.null( betapk0 ) ){
      betapk0filename = paste( c( workdir, "SavePkPg", "priorbetapk0.txt" ), collapse = dirsep )
      #write.table( betapk0, file = betapk0filename,
      #             sep = "\t", row.names = FALSE, col.names = FALSE )
      save_arma_matrix( betapk0, filename = betapk0filename )
      cat("betapk0file = ", betapk0filename, "\n",file = filename, append = TRUE )
  }
  if( !is.null( Bpk0 ) ) {
      Bpk0filename = paste( c( workdir, "SavePkPg", "priorBpk0.txt" ), collapse = dirsep )
      #write.table( Bpk0, file = Bpk0filename,
      #             sep = "\t", row.names = FALSE, col.names = FALSE )
      save_arma_matrix( Bpk0, filename = Bpk0filename )
      cat("Bpk0file = ", Bpk0filename, "\n",file = filename, append = TRUE )
  }
  if( !is.null( betapg0 ) ){
      betapg0filename = paste( c( workdir, "SavePkPg", "priorbetapg0.txt" ), collapse = dirsep )
      #write.table( betapg0, file = betapg0filename,
      #             sep = "\t", row.names = FALSE, col.names = FALSE )
      save_arma_matrix( betapg0, filename = betapg0filename )
      cat("betapg0file = ", betapg0filename, "\n",file = filename, append = TRUE )
  }
  if( !is.null( Bpg0 ) ){
      Bpg0filename = paste( c( workdir, "SavePkPg", "priorBpg0.txt" ), collapse = dirsep )
      write.table( Bpg0, file = Bpg0filename,
                   sep = "\t", row.names = FALSE, col.names = FALSE )
      save_arma_matrix( Bpg0, filename = Bpg0filename )
      cat("Bpg0file = ", Bpg0filename, "\n",file = filename, append = TRUE )
  }
  ##############################################################################
  # initial values 
  ##############################################################################
  cat("#initial value section \n",  file = filename, append = TRUE)
  if( !is.null( Th0 ) ){
      Thetafilename = paste( c( workdir, "SavePkPg", "initTheta.txt" ), collapse = dirsep )
      save_arma_matrix( Th0, filename = Thetafilename )
      cat("initThetafile = ", Thetafilename, "\n",file = filename, append = TRUE )
  }
  if( !is.null( iFit ) ){
      initFitfilename = paste( c( workdir, "SavePkPg", "initFit.txt" ), collapse = dirsep )
      save_arma_cube( iFit, filename = initFitfilename )
      cat("initFitfile = ", initFitfilename, "\n",file = filename, append = TRUE )
  }
  if( !is.null( iZ ) ){
      initZfilename = paste( c( workdir, "SavePkPg", "initZ.txt" ), collapse = dirsep )
      save_arma_matrix( iZ, filename = initZfilename )
      cat("initZfile = ", initZfilename, "\n",file = filename, append = TRUE )
  }
  # if( !is.null( iAlpha ) ){
  #     initAlphafilename = paste( c( workdir, "SavePkPg", "initAlpha.txt" ), collapse = dirsep )
  #     save_arma_matrix( iAlpha, filename = initAlphafilename )
  #     cat("initAlphafile = ", initAlphafilename, "\n",file = filename, append = TRUE )
  # }
  if( !is.null( iBetaPk ) ){
      initBetaPkfilename = paste( c( workdir, "SavePkPg", "initBetaPk.txt" ), collapse = dirsep )
      save_arma_matrix( iBetaPk, filename = initBetaPkfilename )
      cat("initBetaPkfile = ", initBetaPkfilename, "\n",file = filename, append = TRUE )
  }
  if( !is.null( iBetaPg ) ){
      initBetaPgfilename = paste( c( workdir, "SavePkPg", "initBetaPg.txt" ), collapse = dirsep )
      save_arma_matrix( iBetaPg, filename = initBetaPgfilename )
      cat("initBetaPgfile = ", initBetaPgfilename, "\n",file = filename, append = TRUE )
  }
  if( !is.null( iRho ) ){
      initRhofilename = paste( c( workdir, "SavePkPg", "initRho.txt" ), collapse = dirsep )
      save_arma_matrix( iRho, filename = initRhofilename )
      cat("initRhofile = ", initRhofilename, "\n",file = filename, append = TRUE )
  }
  if( !is.null( iTau ) ){
      initTaufilename = paste( c( workdir, "SavePkPg", "initTau.txt" ), collapse = dirsep )
      save_arma_matrix( iTau, filename = initTaufilename )
      cat("initTaufile = ", initTaufilename, "\n",file = filename, append = TRUE )
  }
  if( !is.null( iSigma ) ){
      initSigmafilename = paste( c( workdir, "SavePkPg", "initSigma.txt" ), collapse = dirsep )
      save_arma_matrix( iSigma, filename = initSigmafilename )
      cat("initSigmafile = ", initSigmafilename, "\n",file = filename, append = TRUE )
  }
  if( !is.null( iOsnp ) ){
      initOsnpfilename = paste( c( workdir, "SavePkPg", "initOsnp.txt" ), collapse = dirsep )
      save_arma_matrix( iOsnp, filename = initOsnpfilename )
      cat("initOsnpfile = ", initOsnpfilename, "\n",file = filename, append = TRUE )
  }
  if( !is.null( iOtheta ) ){
      initOthetafilename = paste( c( workdir, "SavePkPg", "initOtheta.txt" ), collapse = dirsep )
      save_arma_matrix( iOtheta, filename = initOthetafilename )
      cat("initOthetafile = ", initOthetafilename, "\n",file = filename, append = TRUE )
  }
  if( !is.null( iGsnp ) ){
      initGsnpfilename = paste( c( workdir, "SavePkPg", "initGsnp.txt" ), collapse = dirsep )
      save_arma_int_matrix( iGsnp, filename = initGsnpfilename )
      cat("initGsnpfile = ", initGsnpfilename, "\n",file = filename, append = TRUE )
  }
  if( !is.null( iGtheta ) ){
      initGthetafilename = paste( c( workdir, "SavePkPg", "initGtheta.txt" ), collapse = dirsep )
      save_arma_int_matrix( iGtheta, filename = initGthetafilename )
      cat("initGthetafile = ", initGthetafilename, "\n",file = filename, append = TRUE )
  }
  if( !is.null( iGrho ) ){
      initGrhofilename = paste( c( workdir, "SavePkPg", "initGrho.txt" ), collapse = dirsep )
      save_arma_int_matrix( iGrho, filename = initGrhofilename )
      cat("initGrhofile = ", initGrhofilename, "\n",file = filename, append = TRUE )
  }

    print("main")
  ###########################################################################################
  # CALL MCMC SAMPLER (C function) ----------------------------------------------------------                    
  obj<-.Call("pkpg_c",
             as.double(      as.double(c1) ), # Drug concentration    
             as.matrix(                 s1 ), # Polymorph                    
             as.matrix(               Time ), # Time matrix
             as.matrix(                Xpk ), # PK covariates
             as.matrix(                Xpg ), # Polymorph covariates
             as.matrix(                Th0 ), # Initial values PK
             as.integer(              PK_V ), # No. PK parameters
             as.integer(             Nconc ), # No. subjects
             as.integer(                 T ), # No. of sampling times
             as.integer(                 Q ), # No. of SNPs
             as.integer(                 K ), # No. of observed compartments
             as.integer(               Ppk ), # No. of PK covariates
             as.integer(               Ppg ), # No. of SNPs covariates
             as.integer(           n_chain ), # No. of mcmc simulations  
             as.integer(            n_mcmc ), # No. of mcmc simulations
             as.integer(            n_burn ), # No. of discarded samples
             as.integer(             n_lag ), # No. thinning factor
             body(                dyn.eval ), # ODE Model
             dyn.env, #new.env(             dyn.eval ), # R environment
                                        cl  ,
             PACKAGE = "bppkgx" );
  stopCluster(cl)  
  # RETURN PkPg OBJECT ---------------------------------------------------------------------

}           
# END fit.pkpg ------------------------------------------------------------------------------- #  
# save_arma_matrix <- function( mat, filename, sep = "\t" ){
#   write( "ARMA_MAT_TXT_FN008",file = filename )
#   write( paste( dim( mat ), sep=" ", collapse=" "), file = filename, append = TRUE )
#   write.table( mat, file = filename,
#                sep = sep, row.names = FALSE, col.names = FALSE, append = TRUE )
# }
# save_arma_int_matrix <- function( mat, filename, sep = "\t" ){
#   write( "ARMA_MAT_TXT_IU004",file = filename )
#   write( paste( dim( mat ), sep=" ", collapse=" "), file = filename, append = TRUE )
#   write.table( mat, file = filename,
#                sep = sep, row.names = FALSE, col.names = FALSE, append = TRUE )
# }


# save_arma_cube <- function( arr, filename, sep = "\t" ){
#   write( "ARMA_CUB_TXT_FN008",file = filename )
#   write( paste( dim( arr ), sep=" ", collapse=" "), file = filename, append = TRUE )
#   write.table( matrix(unlist(aperm(arr,c(2,1,3))),,ncol=dim(arr)[2],byrow=T), file = filename,
#                sep = sep, row.names = FALSE, col.names = FALSE, append = TRUE )
# }
# SEXP vecY_s, 
# SEXP W_s,    
# SEXP Time_s, 
# SEXP Xpk_s,  
# SEXP Xpg_s,  
# SEXP Th0,   
# SEXP PKpar,    
# SEXP N,      
# SEXP T,      
# SEXP G,     
# SEXP K,
# SEXP Ppk,      
# SEXP Ppg,
# SEXP n_chain,     
# SEXP n_mcmc, 
# SEXP burn,  
# SEXP thin,   
# SEXP dyn_eval, 
# SEXP rho  


# # --------------------------------------------------------------------------------------------- #
# # Masanao Yajima                                                                           #
# # --------------------------------------------------------------------------------------------- #
# # FUNCTION :  fit.pkpg()                                                                        #
# # Purpose  :  MCMC sampler for the PK/PGx model of Telesca/Rosner/Muller                        # 
# # Version  :  0.0.1                                                                             #
# # Date     :  April 2012                                                                        # 
# # --------------------------------------------------------------------------------------------- #

# fit.pkpg.parallel <- function( Conc,  
#                       SNP,  
#                       Xpk     = NULL, 
#                       Xpg     = NULL,  
#                       Model,  
#                       Jacobian= NULL, 
#                       Theta0,
#                       gen.init= TRUE,
#                       n_chain = 1,
#                       n_mcmc  = 1000,  
#                       n_burn  = 1000,  
#                       n_lag   = 1,
#                       workdir = ".", 
#                       prior   = 'T')
# {
#   # CREATE PRIOR FILE -----------------------------------------------------------------------
#   # ---- To add an interactive window for prior specification

#   # CREATE OUTPUT DIRECTORY -----------------------------------------------------------------
#   files <- list.files();
#   dir   <- files == "PostPkPg";
#   if( length( files[dir] ) == 0 ){ 
#     cat("Directory PostPkPg created: \n")
#     dir.create("PostPkPg"); 
#   }
#   # CREATE SAVE DIRECTORY -----------------------------------------------------------------
#   if( length( files[files == "SavePkPg"] ) == 0 ){ 
#     cat("Directory SavePkPg created: \n")
#     dir.create("SavePkPg"); 
#   }
#   # GET DATA SUMMARIES ----------------------------------------------------------------------
#   # SNPs --------------------- #
#   Nsnp <- nrow(SNP);           # Number of subject in SNPs matrix
#   Q    <- ncol(SNP) - 1;       # Number of plymorphisms
#   # Concentrations ----------- #
#   K    <- ncol(Conc) - 2;      # Number of compartments
#   Nconc<- length( unique( Conc$ID ) ) # Number of subjects in Concentration file #max(Conc$ID)  
#   TOT  <- nrow(Conc);             # Total number of concentration values
#   # Sampling Times ----------- #
#   #Time    <- TOT/ Nconc;    # Number of sampling times
#   Time  <- unique( table( Conc$ID ) )
#   # Covariates --------------- # 
#   if( is.null( Xpk ) && is.null( Xpg ) ){
#     Xpk  <- rep(1, Nconc);  
#     Xpg  <- rep(1, Nconc);
#     Ppk  <- 1;              
#     Ppg  <- 1;
#   }
#   else {
#     Ppk  <- ncol(Xpk);  
#     Ppg  <- ncol(Xpg);  
#   }
#   # Initial Dosage and Theta Values ---------------------------------------------------------
#   # PKparam  <- names(Theta0)[3:ncol(Theta0)]; 
#   PK_V   <- ncol(Theta0) - 2;

#   # Priors
#   alpha0  = rep( 0, PK_V )
#   A0      = diag( PK_V )
#   betapk0 = matrix(0, Ppk, PK_V )
#   Bpk0    = diag(Ppk)
#   betapg0 = matrix(0, Ppg, Q )
#   Bpg0    = diag(Ppk)
#   # Initial values ( think about how to let users specify this )
#   Th0    <- log(as.matrix(Theta0[,3:ncol(Theta0)]) + 0.0001);

#   if( gen.init == TRUE ){
#     library( mvtnorm )
#     iFit   <- matrix( rnorm( Nconc * K * Time ), K * Nconc, Time )
#     iGsnp  <- diag( Q )
#     iGtheta<- diag( PK_V )
#     iGrho  <- matrix( 0, Q, PK_V )
#     iZ     <- matrix( rnorm( Nsnp * Q ), Nsnp, Q )
#     iAlpha <- matrix( rmvnorm( 1, alpha0, A0 ), PK_V, 1 )
#     iBetaPk<- matrix( rmvnorm( PK_V, rep(0,Ppk), Bpk0 ), Ppk, PK_V )
#     iBetaPg<- matrix( rmvnorm(    Q, rep(0,Ppg), Bpg0 ), Ppg, Q  )
#     iRho   <- matrix( rnorm(Q * PK_V),   Q, PK_V )
#     iTau   <- matrix( rmvnorm( 1, alpha0, A0 ), PK_V, 1 )
#     iSigma <- NULL
#     iGsnp  <- diag( Q )
#     iGtheta<- diag( PK_V )
#     iGrho  <- matrix( 0, Q, PK_V )
#     iOsnp  <- diag( Q )
#     iOtheta<- diag( PK_V )

#   }
#   # VERIFY DIMENTIONAL COMPATIBILTY ---------------------------------------------------------
#   #if(T%%round(T) != 0) stop("Check that all subject have the same number of sampling times");    
#   if( length(Time) != 1 ) stop( "Check that all subject have the same number of sampling times");  
#   if( Nconc != Nsnp )  stop( gettextf("N in Conc is %d, while N in SNP is %d", Nconc, Nsnp), domain = NA);    
#   # DEFINE data structures ------------------------------------------------------------------
#   c1 <- log( as.matrix( Conc )[,3:ncol( Conc )] + 0.1);
#   s1 <- as.matrix( SNP )[,2:ncol( SNP )];    

#   # DEFINE LSODA C ENVELOPE -----------------------------------------------------------------
#   library( deSolve );
#   D0   <- Theta0[,2];                                                    # D0 dose
#   time <- t( matrix( Conc[,2], ncol = Nconc ) );                               # time matrix 
#   dyn.eval <- function( i, param ){
#       parms  <- c( param, D0[i] ); 
#       t      <- c( 0.0, time[i,] );
#       out    <- lsoda( rep( 0, 7 ), t, Model, parms=parms, rtol = 0.5, atol = 0.5, jacfunc= Jacobian );
#       #print(out)
#       out1   <- as.matrix(out)[2:nrow(out),c(2,4,5,6)];
#       #as.double(out1/param[length(param)-1]);
#       as.matrix(t( out1/param[length(param)-1]));
#   }

#   # DEFINE LSODA C ENVELOPE PARALLEL---------------------------------------------------------
#   library( deSolve );
#   library( parallel );
#   D0   <- Theta0[,2];                                                    # D0 dose
#   time <- t( matrix( Conc[,2], ncol = Nconc ) );                               # time matrix 

#   .runlsoda<- function( i, param ) {
#      param <- c( param[i,], D0[i] ); 
#         t  <- c( 0.0, time[i,] );
#       fit  <- lsoda( rep( 0, 7 ), t, Model, parms=param, rtol = 0.5, atol = 0.5, jacfunc= Jacobian )
#       return( fit / param[length(param) - 1 ] )
#   }

#   n.cluster = 2
#   cl <- makeCluster( n.cluster, methods = FALSE )
#   clusterExport( cl, c("time", "D0", "Model", "Jacobian", "lsoda", ".runlsoda") )
#   dyn.eval <- function( param ){
#       #parms  <- c( param, D0[i] ); 
#       #t      <- c( 0.0, time[i,] );
#       clusterExport( cl, c("param") )
#       #out    <- lsoda( rep( 0, 7 ), t, Model, parms=parms, rtol = 0.5, atol = 0.5, jacfunc= Jacobian );
#       #print(out)
#       res    <- parSapply( cl, seq_len(86), .runlsoda )
#       out1   <- matrix( res[-c( 1:16, 17:32, 49:64, 113:128 ), ], ,4 )

#       #as.matrix(out)[2:nrow(out),c(2,4,5,6)];
#       #as.double(out1/param[length(param)-1]);

#       #need to fix this so that array of ( T, K, N ) gets passed back.
#       as.double( t( out1 ) );
#   }
#   dyn.env<- new.env( dyn.eval )
#   assign("time", time, dyn.env)
#   assign("D0", D0, dyn.env)
#   assign("Model", Model, dyn.env)
#   assign("Jacobian", Jacobian, dyn.env)
#   assign(".runlsoda", .runlsoda, dyn.env)
#   assign("lsoda", lsoda, dyn.env)
#   ###########################################################################################
#   info = Sys.info()
  
#   if(info[1]=="Windows"){
#       dirsep="\\"
#   } else {
#       dirsep="/"
#   }
#   ##############################################################################
#   # Parameters
#   ##############################################################################
#   filename = "parameters.txt";8
#   cat("#data section \n", file = filename )
#   cat("N   = ",Nconc, "\n", file = filename, append = TRUE)
#   cat("T   = ", Time, "\n", file = filename, append = TRUE)
#   cat("K   = ",    K, "\n", file = filename, append = TRUE)
#   cat("Q   = ",    Q, "\n", file = filename, append = TRUE)
#   cat("V   = ", PK_V, "\n", file = filename, append = TRUE)
#   cat("Ppk = ",  Ppk, "\n", file = filename, append = TRUE)
#   cat("Ppg = ",  Ppg, "\n", file = filename, append = TRUE)
#   ##############################################################################
#   # data
#   ##############################################################################
#   if( !( is.null( c1 ) ) ){
#       Yfilename = paste( c( workdir, "SavePkPg", "Y.txt" ), collapse = dirsep )
#       write.table( t(c1), file = Yfilename,
#                    sep = "\t", row.names = FALSE, col.names = FALSE )
#       cat("Yfile = ", Yfilename, "\n",file = filename, append = TRUE )
#   }
#   if( !( is.null( s1 ) ) ){
#       Wfilename = paste( c( workdir, "SavePkPg", "W.txt" ), collapse = dirsep )
#       write.table( s1, file = Wfilename,
#                    sep = "\t", row.names = FALSE, col.names = FALSE )
#       cat("Wfile = ", Wfilename, "\n",file = filename, append = TRUE )
#   }
#   if( !( is.null( Xpk ) ) ){
#       Xpkfilename = paste( c( workdir, "SavePkPg", "Xpk.txt" ), collapse = dirsep )
#       write.table( Xpk, file = Xpkfilename,
#                    sep = "\t", row.names = FALSE, col.names = FALSE )
#       cat("Xpkfile = ", Xpkfilename, "\n",file = filename, append = TRUE )
#   }
#   if( !( is.null( Xpg ) ) ){
#       Xpgfilename = paste( c( workdir, "SavePkPg", "Xpg.txt" ), collapse = dirsep )
#       write.table( Xpg, file = Xpgfilename,
#                    sep = "\t", row.names = FALSE, col.names = FALSE )
#       cat("Xpgfile = ", Xpgfilename, "\n",file = filename, append = TRUE )
#   }
#   ##############################################################################
#   # prior
#   ##############################################################################
#   cat("#prior section \n",  file = filename, append = TRUE)
#   if( !( is.null( alpha0 ) ) ){
#       alpha0filename = paste( c( workdir, "SavePkPg", "prioralpha0.txt" ), collapse = dirsep )
#       write.table( alpha0, file = alpha0filename,
#                    sep = "\t", row.names = FALSE, col.names = FALSE )
#       cat("alpha0file = ", alpha0filename, "\n",file = filename, append = TRUE )
#   }
#   if( !is.null( A0 ) ){
#       A0filename = paste( c( workdir, "SavePkPg", "priorA0.txt" ), collapse = dirsep )
#       write.table( A0, file = A0filename,
#                    sep = "\t", row.names = FALSE, col.names = FALSE )
#       cat("A0file = ", A0filename, "\n",file = filename, append = TRUE )
#   }
#   if( !is.null( betapk0 ) ){
#       betapk0filename = paste( c( workdir, "SavePkPg", "priorbetapk0.txt" ), collapse = dirsep )
#       write.table( betapk0, file = betapk0filename,
#                    sep = "\t", row.names = FALSE, col.names = FALSE )
#       cat("betapk0file = ", betapk0filename, "\n",file = filename, append = TRUE )
#   }
#   if( !is.null( Bpk0 ) ) {
#       Bpk0filename = paste( c( workdir, "SavePkPg", "priorBpk0.txt" ), collapse = dirsep )
#       write.table( Bpk0, file = Bpk0filename,
#                    sep = "\t", row.names = FALSE, col.names = FALSE )
#       cat("Bpk0file = ", Bpk0filename, "\n",file = filename, append = TRUE )
#   }
#   if( !is.null( betapg0 ) ){
#       betapg0filename = paste( c( workdir, "SavePkPg", "priorbetapg0.txt" ), collapse = dirsep )
#       write.table( betapg0, file = betapg0filename,
#                    sep = "\t", row.names = FALSE, col.names = FALSE )
#       cat("betapg0file = ", betapg0filename, "\n",file = filename, append = TRUE )
#   }
#   if( !is.null( Bpg0 ) ){
#       Bpg0filename = paste( c( workdir, "SavePkPg", "priorBpg0.txt" ), collapse = dirsep )
#       write.table( Bpg0, file = Bpg0filename,
#                    sep = "\t", row.names = FALSE, col.names = FALSE )
#       cat("Bpg0file = ", Bpg0filename, "\n",file = filename, append = TRUE )
#   }
#   ##############################################################################
#   # initial values 
#   ##############################################################################
#   cat("#initial value section \n",  file = filename, append = TRUE)
#   if( !is.null( Th0 ) ){
#       Thetafilename = paste( c( workdir, "SavePkPg", "initTheta.txt" ), collapse = dirsep )
#       save_arma_matrix( Th0, filename = Thetafilename )
#       cat("initThetafile = ", Thetafilename, "\n",file = filename, append = TRUE )
#   }
#   if( !is.null( iFit ) ){
#       initFitfilename = paste( c( workdir, "SavePkPg", "initFit.txt" ), collapse = dirsep )
#       save_arma_matrix( iFit, filename = initFitfilename )
#       cat("initFitfile = ", initFitfilename, "\n",file = filename, append = TRUE )
#   }
#   if( !is.null( iZ ) ){
#       initZfilename = paste( c( workdir, "SavePkPg", "initZ.txt" ), collapse = dirsep )
#       save_arma_matrix( iZ, filename = initZfilename )
#       cat("initZfile = ", initZfilename, "\n",file = filename, append = TRUE )
#   }
#   if( !is.null( iAlpha ) ){
#       initAlphafilename = paste( c( workdir, "SavePkPg", "initAlpha.txt" ), collapse = dirsep )
#       save_arma_matrix( iAlpha, filename = initAlphafilename )
#       cat("initAlphafile = ", initAlphafilename, "\n",file = filename, append = TRUE )
#   }
#   if( !is.null( iBetaPk ) ){
#       initBetaPkfilename = paste( c( workdir, "SavePkPg", "initBetaPk.txt" ), collapse = dirsep )
#       save_arma_matrix( iBetaPk, filename = initBetaPkfilename )
#       cat("initBetaPkfile = ", initBetaPkfilename, "\n",file = filename, append = TRUE )
#   }
#   if( !is.null( iBetaPg ) ){
#       initBetaPgfilename = paste( c( workdir, "SavePkPg", "initBetaPg.txt" ), collapse = dirsep )
#       save_arma_matrix( iBetaPg, filename = initBetaPgfilename )
#       cat("initBetaPgfile = ", initBetaPgfilename, "\n",file = filename, append = TRUE )
#   }
#   if( !is.null( iRho ) ){
#       initRhofilename = paste( c( workdir, "SavePkPg", "initRho.txt" ), collapse = dirsep )
#       save_arma_matrix( iRho, filename = initRhofilename )
#       cat("initRhofile = ", initRhofilename, "\n",file = filename, append = TRUE )
#   }
#   if( !is.null( iTau ) ){
#       initTaufilename = paste( c( workdir, "SavePkPg", "initTau.txt" ), collapse = dirsep )
#       save_arma_matrix( iTau, filename = initTaufilename )
#       cat("initTaufile = ", initTaufilename, "\n",file = filename, append = TRUE )
#   }
#   if( !is.null( iSigma ) ){
#       initSigmafilename = paste( c( workdir, "SavePkPg", "initSigma.txt" ), collapse = dirsep )
#       save_arma_matrix( iSigma, filename = initSigmafilename )
#       cat("initSigmafile = ", initSigmafilename, "\n",file = filename, append = TRUE )
#   }
#   if( !is.null( iOsnp ) ){
#       initOsnpfilename = paste( c( workdir, "SavePkPg", "initOsnp.txt" ), collapse = dirsep )
#       save_arma_matrix( iOsnp, filename = initOsnpfilename )
#       cat("initOsnpfile = ", initOsnpfilename, "\n",file = filename, append = TRUE )
#   }
#   if( !is.null( iOtheta ) ){
#       initOthetafilename = paste( c( workdir, "SavePkPg", "initOtheta.txt" ), collapse = dirsep )
#       save_arma_matrix( iOtheta, filename = initOthetafilename )
#       cat("initOthetafile = ", initOthetafilename, "\n",file = filename, append = TRUE )
#   }
#   if( !is.null( iGsnp ) ){
#       initGsnpfilename = paste( c( workdir, "SavePkPg", "initGsnp.txt" ), collapse = dirsep )
#       save_arma_int_matrix( iGsnp, filename = initGsnpfilename )
#       cat("initGsnpfile = ", initGsnpfilename, "\n",file = filename, append = TRUE )
#   }
#   if( !is.null( iGtheta ) ){
#       initGthetafilename = paste( c( workdir, "SavePkPg", "initGtheta.txt" ), collapse = dirsep )
#       save_arma_int_matrix( iGtheta, filename = initGthetafilename )
#       cat("initGthetafile = ", initGthetafilename, "\n",file = filename, append = TRUE )
#   }
#   if( !is.null( iGrho ) ){
#       initGrhofilename = paste( c( workdir, "SavePkPg", "initGrho.txt" ), collapse = dirsep )
#       save_arma_int_matrix( iGrho, filename = initGrhofilename )
#       cat("initGrhofile = ", initGrhofilename, "\n",file = filename, append = TRUE )
#   }
#   ###########################################################################################
#   # CALL MCMC SAMPLER (C function) ----------------------------------------------------------                    
#   obj<-.Call("mcmc_c",
#              as.double(              t(c1) ), # Drug concentration    
#              as.matrix(                 s1 ), # Polymorph                    
#              as.matrix(           Conc[,2] ), # Time matrix
#              as.matrix(                Xpk ), # PK covariates
#              as.matrix(                Xpg ), # Polymorph covariates
#              as.matrix(                Th0 ), # Initial values PK
#              as.integer(              PK_V ), # No. PK parameters
#              as.integer(             Nconc ), # No. subjects
#              as.integer(              Time ), # No. of sampling times
#              as.integer(                 Q ), # No. of SNPs
#              as.integer(                 K ), # No. of observed compartments
#              as.integer(               Ppk ), # No. of PK covariates
#              as.integer(               Ppg ), # No. of SNPs covariates
#              as.integer(           n_chain ), # No. of mcmc simulations  
#              as.integer(            n_mcmc ), # No. of mcmc simulations
#              as.integer(            n_burn ), # No. of discarded samples
#              as.integer(             n_lag ), # No. thinning factor
#              body(                dyn.eval ), # ODE Model
#                                      dyn.env, # R environment
#              PACKAGE = "testpackage" );

#   stopCluster(cl)  
#   # RETURN PkPg OBJECT ---------------------------------------------------------------------
#   return( obj )
# }           
# # END fit.pkpg ------------------------------------------------------------------------------- #  
# save_arma_matrix <- function( mat, filename, sep = "\t" ){
#   write( "ARMA_MAT_TXT_FN008",file = filename )
#   write( paste( dim( mat ), sep=" ", collapse=" "), file = filename, append = TRUE )
#   write.table( mat, file = filename,
#                sep = sep, row.names = FALSE, col.names = FALSE, append = TRUE )
# }
# save_arma_int_matrix <- function( mat, filename, sep = "\t" ){
#   write( "ARMA_MAT_TXT_IU004",file = filename )
#   write( paste( dim( mat ), sep=" ", collapse=" "), file = filename, append = TRUE )
#   write.table( mat, file = filename,
#                sep = sep, row.names = FALSE, col.names = FALSE, append = TRUE )
# }

# save_arma_cube <- function( mat, filename, sep = "\t" ){
#   write( "ARMA_CUB_TXT_FN008",file = filename )
#   write( paste( dim( mat ), sep=" ", collapse=" "), file = filename, append = TRUE )
#   write.table( mat, file = filename,
#                sep = sep, row.names = FALSE, col.names = FALSE, append = TRUE )
# }
# # SEXP vecY_s, 
# # SEXP W_s,    
# # SEXP Time_s, 
# # SEXP Xpk_s,  
# # SEXP Xpg_s,  
# # SEXP Th0,   
# # SEXP PKpar,    
# # SEXP N,      
# # SEXP T,      
# # SEXP G,     
# # SEXP K,
# # SEXP Ppk,      
# # SEXP Ppg,
# # SEXP n_chain,     
# # SEXP n_mcmc, 
# # SEXP burn,  
# # SEXP thin,   
# # SEXP dyn_eval, 
# # SEXP rho  

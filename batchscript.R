starttime <- Sys.time()
Sys.getenv(c("SLURM_SUBMIT_DIR"))
Sys.getenv(c("HOST", "SLURM_JOB_ID", "SLURM_NODELIST", "SLURM_NNODES", "SLURM_NTASKS", "SLURM_CPUS_PER_TASK", "SLURM_CPUS_ON_NODE", "SLURM_NTASKS_PER_NODE",  "SLURM_TASK_PID",  "SLURM_PARTITION"))

library(mvtnorm)
library(doParallel)
library(foreach)
library(plyr)
cl <- makeCluster(10)  # Use 3 cores
registerDoParallel(cl) # register these 3
dirs=NULL
for(i in 1:10){
	dirtosave=file.path(".", i)
	dir.create(dirtosave)
	dirs=c(dirs,dirtosave)
}
remindex <- 1:10

foo<-function(tp,dirs,remindex){
	library(bppkgx)
	setwd(dirs[tp])
	dyn.load("/shared/silo_researcher/Gottardo_R/myajima_working/RExp/bppkgx/testdirectory/mymod.so")
	deleteIdx <-remindex[tp]
	Nconc     <- Nconc - length(deleteIdx)
	ConcArray <- ConcArray[,,-deleteIdx]
	Time      <- Time[-deleteIdx, ]
	T0        <- T0[-deleteIdx, ]
	SNP       <- SNP[-deleteIdx,]
	X1        <- X1[-deleteIdx,]
	EHRTime   <- matrix(Theta0[-deleteIdx,11],length(Theta0[-deleteIdx,11]),1)
	Theta0    <- Theta0[-deleteIdx,-c(11)]
	Stheta0   <- cov(log(Theta0+0.01))*0.001/dim(Theta0)[2]
	PK_V      <- ncol(Theta0);
	fit <- fit.pkpg( Conc   = ConcArray,  SNP = SNP,  Time = Time,    Xpk = X1,        Xpg = X1,  
                 Model  = EHRpk, Jacobian = jac1,  
                 Theta0 = Theta0,  PK_V= PK_V, Dose0 = T0[,2], EHRT0 =EHRTime, #Volume = V0,
                 iSigmaTheta0 = Stheta0,
                 n_chain  = Nchain, n_mcmc = Nmcmc, n_burn   = Nburn,  
                 n_lag    = Nlag,   prior  = 'T')
	return(NULL)
}

foreach(i=1:10) %dopar% foo(i,dirs,remindex)

q(save = "yes")

Sys.time() - starttime


sbatch -J "ppg" --cpus-per-task=10 --time=10-4 --wrap="R --save --restore < /shared/silo_researcher/Gottardo_R/myajima_working/RExp/bppkgx/testdirectory/crossvalidation/fullfold10/pkrunCV10.R"
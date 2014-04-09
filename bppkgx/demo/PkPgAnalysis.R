# ----------------------------------------------------------------------------------- #
#   Pk/Pgx Model (Telesca / Rosner / Muller, 2010)                                    #
# ----------------------------------------------------------------------------------- #
#setwd("~/Dropbox/PKPG")
# Load Data ------------------------------------------------------------------------- #
Conc.raw <- read.table("~/Dropbox/PKPG/03_Data/Concentration_fixed.txt", header=T);
#Conc <- subset(Conc, Conc$ID <= 30);

Conc <- Conc.raw
Conc[is.na(Conc)]<- -999999
for(i in 3:6){
  for(j in nrow(Conc):1){
      if( Conc[j, i] <= 0.15 ){ 
            # if(!is.na(Conc[j-1, i]) && !is.na(Conc[j+1, i]) && Conc[j, 1]==Conc[j-1, 1]&&Conc[j, 1]==Conc[j+1, 1]){
            #     Conc[j, i]  <- (Conc[(j-1),i] + Conc[(j+1),i])/2; # average
            # } else 
            if(Conc[j-1, i] > 0.15 && Conc[j, 1]==Conc[j-1, 1]){
                Conc[j, i] <- Conc[(j-1),i] - 1; # last value carry forward
            } else if(Conc[j+1, i] > 0.15 && Conc[j, 1]==Conc[j+1, 1]){
                Conc[j, i] <- Conc[(j+1),i] + 1; # next value carried backward
            }
	    }
	}
}

for(j in 1:nrow(Conc)){
  if( Conc[j, 2] < 0 ) Conc[j, 2] =Conc[j-1, 2]+1
}

#Conc[is.na(Conc)] <- 919.919;
Nconc<- length( unique( Conc$ID ) )
ConcArray <- array(NA, c(15, 4, Nconc) )

for(i in 1:Nconc){
  ConcArray[,,i]=as.array( unlist( Conc[(15*(i-1)+1):(15*(i)),c(3:6)]))
}

Time <- matrix( Conc[,2], 86, 15, byrow = TRUE )
#
SNP  <- read.table("~/Dropbox/PKPG/03_Data/snp.txt", header = T) - 1;
SNP  <- SNP[,names(SNP)!="ID"]
#SNP[is.na(SNP)] <- 919.919;
SNP[is.na(SNP)] <- 0
# 
X    <- read.csv("~/Dropbox/PKPG/03_Data/covariate.csv", header=T); 
race <- rep(1, nrow(X));
race[X$race != "African American"] = 0;
sex  <- X$sex - 1; 
age  <- X$age - mean(X$age, na.rm=T);
X1   <- cbind(rep(1, length(age)), race, sex, age); 
#X1[is.na(X1)] <- 919.919;
#
T0     <- as.matrix(read.table("~/Dropbox/PKPG/03_Data/T01.txt", header=T));
Theta0 <- as.matrix(T0[,3:ncol(T0)]);


# delete observation 47
deleteIdx <- c(47,62,81)
Nconc     <- Nconc - length(deleteIdx)
ConcArray <- ConcArray[,,-deleteIdx]
Time      <- Time[-deleteIdx, ]
T0        <- T0[-deleteIdx, ]
SNP       <- SNP[-deleteIdx,]
X1        <- X1[-deleteIdx,]
#V0       <- Theta0[-deleteIdx,15]
EHRTime   <- matrix(Theta0[-deleteIdx,11],length(Theta0[-deleteIdx,11]),1)
Theta0    <- Theta0[-deleteIdx,-c(11)]
PK_V      <- ncol(Theta0);
# Param.raw <- read.csv("~/Dropbox/PKPG/03_Data/PkPar.csv", header=T);
# Param.names<- c( 'Ke',  'Kcp', 'Kpc', 'Kesn', 'K30', 'K35', 'K50', 'Keapc', 'K70', 'K3B', 'EHRT', 'KBG', 'KBG1', 'KG3', 'V')
# Param.init <- Param.raw[-deleteIdx,Param.names]
# Define PK model ------------------------------------------------------------------- #
#"1"  "Ke"   #"2"  "Kcp"  
#"3"  "Kpc"  #"4"  "Kesn" 
#"5"  "K30"  #"6"  "K35"  
#"7"  "K50"  #"8"  "Keapc"
#"9"  "K70"  #"10" "K3B"  
#"11" "EHRT" #"12" "KBG"  
#"13" "KBG1" #"14" "KG3"  
# ---------------------------------------------------------------------------------- #
# EHRpk <- function(t, y, p){
#     if(t <= 1.5) yd1 <- p[ 3]*y[2] - (p[ 2]+p[ 4]+p[ 8]+p[1])*y[1] + p[16];
#     if(t >  1.5) yd1 <- p[ 3]*y[2] - (p[ 2]+p[ 4]+p[ 8]+p[1])*y[1];
#                  yd2 <- p[ 3]*y[1] -              p[ 2]      *y[2];
#                  yd3 <- p[ 4]*y[1] - (p[ 5]+p[ 6]+p[10]     )*y[3] + p[14]*y[7]; 
#                  yd4 <- p[ 6]*y[3] -              p[ 7]      *y[4];
#                  yd5 <- p[ 8]*y[1] -              p[ 9]      *y[5];
#     if( t > p[11] && t < ( p[11] + 1 ) ){ 
#                  yd6 <- p[10]*y[3] - (p[12]+p[13]           )*y[6];
#                  yd7 <-              (p[12]+p[13]           )*y[6] - p[14]*y[7];
#     } else {
#                  yd6 <- p[10]*y[3] -  p[12]*y[6];
#                  yd7 <-               p[12]*y[6]                   - p[14]*y[7];
#     }
#     list( c( yd1, yd2, yd3, yd4, yd5, yd6, yd7 ) );
# }
EHRpk  <- function  (  t, y, p ) {
    ydot <- vector(len = 7)
    if(t <= 1.5) ydot[1] <- p[ 3]*y[2] - (p[ 2]+p[ 4]+p[ 8]+p[1])*y[1] + p[16];
    if(t >  1.5) ydot[1] <- p[ 3]*y[2] - (p[ 2]+p[ 4]+p[ 8]+p[1])*y[1];
                 ydot[2] <- p[ 2]*y[1] -  p[ 3]*y[2] ;
                 ydot[3] <- p[ 4]*y[1] - (p[ 5]+p[ 6]+p[10]     )*y[3] + p[14]*y[7]; 
                 ydot[4] <- p[ 6]*y[3] -              p[ 7]      *y[4];
                 ydot[5] <- p[ 8]*y[1] -              p[ 9]      *y[5];
    if( t > p[11] && t < ( p[11] + 1 ) ){ 
                 ydot[6] <- p[10]*y[3] - (p[12]+p[13]           )*y[6];
                 ydot[7] <-              (p[12]+p[13]           )*y[6] - p[14]*y[7];
    } else {
                 ydot[6] <- p[10]*y[3] -  p[12]*y[6];
                 ydot[7] <-               p[12]*y[6]                   - p[14]*y[7];
    }
       return(list(ydot))
    }

# EHRpk  <- function  (  t, y, p ) {
#     ydot <- vector(len = 7)
#     if(t <= 1.5) ydot[1] <- p[ 3]*y[2] - (p[ 2]+p[ 4]+p[ 8]+p[1])*y[1] + p[16];
#     if(t >  1.5) ydot[1] <- p[ 3]*y[2] - (p[ 2]+p[ 4]+p[ 8]+p[1])*y[1];
#                  ydot[2] <- p[ 3]*y[1] -              p[ 2]      *y[2];
#                  ydot[3] <- p[ 4]*y[1] - (p[ 5]+p[ 6]+p[10]     )*y[3] + p[14]*y[7]; 
#                  ydot[4] <- p[ 6]*y[3] -              p[ 7]      *y[4];
#                  ydot[5] <- p[ 8]*y[1] -              p[ 9]      *y[5];
#     if( t > p[11] && t < ( p[11] + 1 ) ){ 
#                  ydot[6] <- p[10]*y[3] - (p[12]+p[13]           )*y[6];
#                  ydot[7] <-              (p[12]+p[13]           )*y[6] - p[14]*y[7];
#     } else {
#                  ydot[6] <- p[10]*y[3] -  p[12]*y[6];
#                  ydot[7] <-               p[12]*y[6]                   - p[14]*y[7];
#     }
#        return(list(ydot))
#     }

# Define Jacobian (optional) -------------------------------------------------------- #
jac1 <- function(t, y, p){

    if(t > p[11] && t < (p[11]+1)){ 
        jstar = matrix(nrow = 7, ncol = 7, byrow = FALSE,
                       data=
                       c(-(p[2]+p[4]+p[8]+p[1]),  p[2],               p[4],     0,  p[8],                0,               0,
                                           p[3], -p[3],                  0,     0,     0,                0,               0,
                                              0,     0, -(p[5]+p[6]+p[10]),  p[6],     0,            p[10],               0,   
                                              0,     0,                  0, -p[7],     0,                0,               0, 
                                              0,     0,                  0,     0, -p[9],                0,               0,
                                              0,     0,                  0,     0,     0, -(p[12] + p[13]), (p[12] + p[13]),
                                              0,     0,              p[14],     0,     0,                0,         -p[14] )
                       )
  }else{
    jstar = matrix(nrow = 7, ncol = 7, byrow = FALSE,
                       data= 
                       c(-(p[2]+p[4]+p[8]+p[1]),  p[2],               p[4],     0,  p[8],                0,              0,
                                           p[3],  -p[3],                  0,     0,     0,                0,              0,
                                              0,     0, -(p[5]+p[6]+p[10]),  p[6],     0,            p[10],              0,   
                                              0,     0,                  0, -p[7],     0,                0,              0, 
                                              0,     0,                  0,     0, -p[9],                0,              0,
                                              0,     0,                  0,     0,     0,  -p[12]         ,  p[12]        ,
                                              0,     0,              p[14],     0,     0,                0,         -p[14])
                       )
  }
    return(jstar);
}
# jac1 <- function(t, y, p){

#     if(t > p[11] && t < (p[11]+1)){ 
#         jstar = matrix(nrow = 7, ncol = 7, byrow = FALSE,
#                        data=
#                        c(-(p[2]+p[4]+p[8]+p[1]),  p[3],               p[4],     0,  p[8],                0,               0,
#                                            p[3], -p[2],                  0,     0,     0,                0,               0,
#                                               0,     0, -(p[5]+p[6]+p[10]),  p[6],     0,            p[10],               0,   
#                                               0,     0,                  0, -p[7],     0,                0,               0, 
#                                               0,     0,                  0,     0, -p[9],                0,               0,
#                                               0,     0,                  0,     0,     0, -(p[12] + p[13]), (p[12] + p[13]),
#                                               0,     0,              p[14],     0,     0,                0,         -p[14] )
#                        )
#   }else{
#     jstar = matrix(nrow = 7, ncol = 7, byrow = FALSE,
#                        data= 
#                        c(-(p[2]+p[4]+p[8]+p[1]),  p[3],               p[4],     0,  p[8],                0,              0,
#                                            p[3], -p[2],                  0,     0,     0,                0,              0,
#                                               0,     0, -(p[5]+p[6]+p[10]),  p[6],     0,            p[10],              0,   
#                                               0,     0,                  0, -p[7],     0,                0,              0, 
#                                               0,     0,                  0,     0, -p[9],                0,              0,
#                                               0,     0,                  0,     0,     0,  -p[12]         ,  p[12]        ,
#                                               0,     0,              p[14],     0,     0,                0,         -p[14])
#                        )
#   }
#     return(jstar);
# }

#Stheta0 <- round( cov(log(Theta0[-c(15,28, 33,44,57,60,68,73,83),]+0.01))*5.6644/dim(Theta0)[2],2)
Stheta0 <- cov(log(Theta0+0.01))*0.001/dim(Theta0)[2]
#Stheta0 <- diag(diag(cov(log(Theta0+0.000001))))*0.1^2/dim(Theta0)[2]
#Stheta0 <- diag(dim(Theta0)[2])*0.1^2/dim(Theta0)[2]
Nmcmc =10000
Nburn =3000
Nchain=15
Nlag  =200
# Call library and fit model -------------------------------------------------------- #
#library(PkPg);
library(testpackage);
tpm <- proc.time();
workdir="/Users/masanaoyajima/Documents/Test/Rarmadillo/testdirectory/test201303016"
dir.create(workdir); 
setwd(workdir)
dyn.load("/Users/masanaoyajima/Documents/Test/Rarmadillo/testdirectory/mymod.so")
if(0){
fit <- fit.pkpg( Conc   = ConcArray,  SNP = SNP,  Time = Time,    Xpk = X1,        Xpg = X1,  
         Model  = EHRpk, Jacobian = jac1,  Theta0 = Theta0,  PK_V= PK_V, Dose0 = T0[,2], EHRT0 =EHRTime, #Volume = V0,
         iSigmaTheta0 = Stheta0,
           n_chain= Nchain, n_mcmc = Nmcmc, n_burn   = Nburn,  
           n_lag    = Nlag,   prior  = 'T')
}
setwd(paste(workdir,"/PostPkPg",sep=""))

#Nmcmc<-1596

N<-Nconc
K<-dim(ConcArray)[2]
Nsnp <- dim(SNP)[2]
Ntime <- dim(Time)[2]
# read result
Sigma=read.table("Sigma.txt")
Sigma<-array(unlist(Sigma), c(N,Nmcmc,K))
Sigma <- aperm(Sigma,c(1,3,2))

Fit = read.table("Fit.txt")
Fit <- array(unlist(Fit), c(Ntime, N,Nmcmc,K))
Fit <- aperm(Fit,c(2, 1,4,3))

Theta = read.table("Theta.txt")  # plot(Theta[seq(15,dim(Theta)[1],by=83),10])
Theta<- array(unlist(Theta), c( N,Nmcmc,PK_V))
Theta <- aperm(Theta,c(1,3,2))

Otheta = read.table("Otheta.txt")# plot(Otheta[seq(12,dim(Otheta)[1],by=14),10])
Otheta <- array(unlist(Otheta), c(PK_V, Nmcmc, PK_V))
Otheta <-aperm(Otheta,c(1, 3, 2))


Osnp = read.table("Osnp.txt") #plot(Osnp[seq(1,dim(Osnp)[1],by=41),10])
Osnp <- array(unlist(Osnp), c(Nsnp, Nmcmc, Nsnp))
Osnp <-aperm(Osnp,c(1, 3, 2))


Rho = read.table("Rho.txt")#plot(Rho[seq(1,dim(Rho)[1],by=41),10])
Rho <- array(unlist(Rho), c(Nsnp, Nmcmc, PK_V))
Rho <-aperm(Rho,c(1, 3, 2))

Z = read.table("Z.txt")#plot(Z[seq(1,dim(Z)[1],by=83),2])
Z <- array(unlist(Z), c(N, Nmcmc, Nsnp))
Z <-aperm(Z,c(1, 3, 2))

BetaPk = read.table("BetaPk.txt")# plot(BetaPk[seq(1,dim(BetaPk)[1],by=4),2])
BetaPk <- array(unlist(BetaPk), c(K, Nmcmc, PK_V))
BetaPk <-aperm(BetaPk,c(1, 3, 2))

BetaPg = read.table("BetaPg.txt")# plot(BetaPg[seq(12,dim(BetaPg)[1],by=4),2])
BetaPg <- array(unlist(BetaPg), c(K , Nmcmc,Nsnp ))
BetaPg <-aperm(BetaPg,c(1, 3, 2))

EHRTime = read.table("EHRTime.txt")

# Summary
Sigma.mean <- apply( Sigma, c( 1, 2 ), mean )
Otheta.mean<- apply(Otheta, c( 1, 2 ), mean )
Osnp.mean  <- apply(  Osnp, c( 1, 2 ), mean )
Otheta.median<- apply(Otheta, c( 1, 2 ), median)
Osnp.median  <- apply(  Osnp, c( 1, 2 ), median )
Theta.mean <- apply( Theta, c( 1, 2 ), mean )




Rho.mean  <-apply(Rho,c(1,2),mean)
Rho.median<-apply(Rho,c(1,2),median)
dimnames(Rho.mean)<-list(names(SNP),dimnames(Theta0)[[2]])
dimnames(Rho.median)<-list(names(SNP),dimnames(Theta0)[[2]])
Z.mean  <-apply(Z,c(1,2),mean)
  
BetaPk.mean<-apply(BetaPk,c(1,2),mean)
BetaPg.mean<-apply(BetaPg,c(1,2),mean)
dimnames(BetaPk.mean)<-list(dimnames(X1)[[2]],dimnames(Theta0)[[2]])
dimnames(BetaPg.mean)<-list(dimnames(X1)[[2]],names(SNP))
Mean.Fit.Patient <- array(NA,c(N,Ntime,K))
for(i in 1:N){
  for(j in 1:K){
    Mean.Fit.Patient[i,,j]<- rowMeans( Fit[i,,j,])
  }
}

ci.Fit.Patient <- array(NA,c(N,Ntime,K,7))
for(i in 1:N){
  for(j in 1:K){
    for(t in 1:Ntime){
      ci.Fit.Patient [i,t,j,]<- quantile( Fit[i,t,j,],c(0.01,0.025,0.05,0.5,0.95,0.975,0.99))
    }
  }
}

MSE.Fit <- array(NA,c(N,Nmcmc,K))
for(i in 1:N){
  for(j in 1:K){
    for(k in 1:Nmcmc){
      MSE.Fit[i,k,j] <- sum( (Fit[i,,j,k]-Mean.Fit.Patient[i,,j])^2 ) / Ntime
    }
  }
}

rank.MSE.Fit <- array(NA,c(N,Nmcmc,K))
for(i in 1:N){
  for(j in 1:K){
    rank.MSE.Fit[i,,j] <-rank(MSE.Fit[i,,j],ties.method="first")
  }
}
# plot 
dir.create("plot"); 
# plot(c(1,K), range(Sigma.mean),type="n")
# for(i in 1:N){
#   lines(1:K,Sigma.mean[i,] ,col=rgb(0,0,1,alpha=0.3))
# }
library(RColorBrewer)
library(MASS)
compname=c("irinotecan", "SN-38", "SN-38G (glucuronide)", "APC")
linecolors = colorRampPalette( brewer.pal( 9, "PuBuGn" ) )(Nmcmc)
linecolors = colorRampPalette( brewer.pal( 9, "YlOrBr" ) )(Nmcmc)
linecolors = colorRampPalette( brewer.pal( 9, "Oranges" ) )(Nmcmc)
linecolors = colorRampPalette( brewer.pal( 9, "Greys" ) )(Nmcmc)
linecolors.trans=rev(linecolors)
for(i in 1:Nmcmc) {
 rgbcol=col2rgb(linecolors.trans[i])
 linecolors.trans[i] = rgb(rgbcol[1],rgbcol[2],rgbcol[3],alpha=50, maxColorValue =255)
}
  for(id in 1:N){
    pdf(paste("./plot/ConcFitG_",id,".pdf",sep=""),width=12,height=2.7)
    par(oma=c(0,0,0,0),mgp=c(2,.5,0), tck=-.04)
    par(mar=c(3.5,3.5,3,1))
    par(mfrow=c(1,K))
    for(k in 1:K){
      plot(range(Time[id,]),range(c((Fit[id,,k,]),(ConcArray[,k,id]))),type="n",main = compname[k],ylab="concentration",xlab="hours after infusion started")
      rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = "gray85")
      oMSE <- order( MSE.Fit[id,,k] ,decreasing=TRUE)
      for(i in seq(1,Nmcmc,by=10)){ lines(Time[id,],(Fit[id,,k,oMSE[i]]) ,col=linecolors.trans[rank.MSE.Fit[id,oMSE[i],k]]) }
            #for(i in 1:Nmcmc){ lines(1:Ntime,(Fit[id,,k,i]) ,col=rgb(0,0,1,alpha=0.3)) }
      #lines(Time[id,], (Mean.Fit.Patient[id,,k]),col="red",lty=3,lwd=3)
      lines(Time[id,], (ci.Fit.Patient[id,,k,4]),col="red",lty=3,lwd=3)
      lines(Time[id,], (ci.Fit.Patient[id,,k,2]),col="orange",lty=2,lwd=1)
      lines(Time[id,], (ci.Fit.Patient[id,,k,6]),col="orange",lty=2,lwd=1)
      # lines(Time[id,], (ci.Fit.Patient[id,,k,1]),col="yellow",lty=3,lwd=0.5)
      # lines(Time[id,], (ci.Fit.Patient[id,,k,7]),col="yellow",lty=3,lwd=0.5)
      lines(Time[id,], (ConcArray[,k,id]),col=rgb(0, 191, 255,alpha=255,maxColorValue = 255),lty=1,lwd=2.5)#

    }
    dev.off()
  }


   for(id in 1:N){
      pdf(paste("./plot/ParcoordTheta_",id,".pdf",sep=""),width=16,height=4)
      parcoord(rbind(t(Theta[id,,]), Theta0[id,]),col=c(linecolors,rgb(0, 191, 255,alpha=255,maxColorValue = 255)),var.label=TRUE)
      dev.off()
    }








ci.snp.Theta <- array(NA,c(3,PK_V,7))
for(i in 0:2){
  for(j in 1:PK_V){
      ci.snp.Theta[i+1,j,]<- quantile( log(Theta[SNP[,1]==i,j,]),c(0.01,0.025,0.05,0.5,0.95,0.975,0.99))
  }
}
dimnames(ci.snp.Theta)<-list(c(0,1,2),Pk.names,c("1%","2.5%","5%","50%","95%","97.5%","99%"))

ccsnp <-brewer.pal(9,"Oranges")
ccsnp3 <-ccsnp[c(2,5,8)]
      pdf(paste("./plot/logParcoordThetaSNP.pdf",sep=""),width=14,height=3.5)
          par(oma=c(0,0,0,0),mgp=c(2,.5,0), tck=-.04)
    par(mar=c(3,0.2,1,0.2))
      parcoord(rbind( ci.snp.Theta[,,4], ci.snp.Theta[,,2], ci.snp.Theta[,,6] ),col=c(ccsnp3,ccsnp3,ccsnp3) ,var.label=TRUE
        ,lty=c(rep(1,3),rep(2,3),rep(2,3)),lwd=c(rep(2,3),rep(0.5,3),rep(0.5,3)))
                  rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = "gray85") 
      parcoord(rbind( ci.snp.Theta[,,4], ci.snp.Theta[,,2], ci.snp.Theta[,,6] ),col=c(ccsnp3,ccsnp3,ccsnp3) ,var.label=TRUE
        ,lty=c(rep(1,3),rep(2,3),rep(2,3)),lwd=c(rep(2,3),rep(0.5,3),rep(0.5,3)),add=T)
      dev.off()

ci.snp.Theta2 <- array(NA,c(9,PK_V,7))
for(i in 1:9){
  for(j in 1:PK_V){
      ci.snp.Theta2[i,j,]<- quantile( log(Theta[uuscore==i,j,]),c(0.01,0.025,0.05,0.5,0.95,0.975,0.99))
  }
}
dimnames(ci.snp.Theta2)<-list(1:9,Pk.names,c("1%","2.5%","5%","50%","95%","97.5%","99%"))

      pdf(paste("./plot/logParcoordThetaSNP2.pdf",sep=""),width=14,height=3.5)

          par(oma=c(0,0,0,0),mgp=c(2,.5,0), tck=-.04)
    par(mar=c(3,0.2,1,0.2))
      parcoord(rbind( ci.snp.Theta2[,,4], ci.snp.Theta2[,,2], ci.snp.Theta2[,,6] ),col=c(ccsnp,ccsnp,ccsnp) ,var.label=TRUE
        ,lty=c(rep(1,9),rep(2,9),rep(2,9)),lwd=c(rep(2,9),rep(0.5,9),rep(0.5,9)))
                  rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = "gray85") 
      parcoord(rbind( ci.snp.Theta2[,,4], ci.snp.Theta2[,,2], ci.snp.Theta2[,,6] ),col=c(ccsnp,ccsnp,ccsnp) ,var.label=TRUE
        ,lty=c(rep(1,9),rep(2,9),rep(2,9)),lwd=c(rep(2,9),rep(0.5,9),rep(0.5,9)),add=T)
      dev.off()
  
# id=84
# plot(c(1,PK_V),range(c(log(Fit[id,,k,]),log(ConcArray[,k,id]))),type="n")
#     for(i in 1:Nmcmc){ lines(1:PK_V,log(Fit[id,,k,i]) ,col=rgb(0,0,1,alpha=0.3)) }
#       lines(1:PK_V, log(ConcArray[,k,id]),col="red")

plot(c(1,PK_V), range(Theta),type="n")
for(i in 1:N){
  lines(1:PK_V,Theta.mean[i,] ,col=rgb(0,0,1,alpha=0.3))
}
# Theta
for(i in 1:N){
  pdf(paste("./plot/ThetaHist_",i,"_",j,".pdf",sep=""),width=20,height=12)
  par(mfrow=c(3,5))
  for(j in 1:PK_V){
    hist((Theta[i,j,]),xlim=range(c((Theta[i,j,]),Theta0[i,j])))
    abline(v=Theta0[i,j],col="red")
  } 
  dev.off()
}
 
#Rho
for(i in 1:Nsnp){
  for(j in 1:PK_V){
  pdf(paste("./plot/RhoHist_",i,"_",j,".pdf",sep=""),width=4,height=4)
  hist(Rho[i,j,])
  dev.off()
}
}

# Z
for(i in 1:N){
  pdf(paste("./plot/ZHist_",i,"_",j,".pdf",sep=""),width=20,height=20)
  par(mfrow=c(6,7))
  for(j in 1:Nsnp){
    hist(Z[i,j,])
  } 
  dev.off()
}
Pk.names<- c( 'Ke',  'Kcp', 'Kpc', 'Ksn', 'K30', 'K35', 'K50', 'Kapc', 'K70', 'K3B', 'KBG', 'KBG1', 'KG3', 'V')
pdf("./plot/density_Otheta.pdf",width=30,heigh=30)
par(mfrow=c(PK_V,PK_V),oma=c(0,0,0,0),mgp=c(2,.5,0), tck=-.04)
for( i in 1:PK_V) for(j in 1:PK_V){
if(i!=j){
par(mar=c(1,1,0.2,0))
plot(density(Otheta[i,j,]),xlab="",ylab="",main="",xlim=range(Otheta[i,j,]),axes=FALSE, cex.main=3) #main=bms[i,j],
#hist(Otheta[i,j,],xlab="",ylab="",main="",xlim=range(Otheta),axes=FALSE, cex.main=3) #main=bms[i,j],
axis( side=1, at=0, label="", tick=TRUE,)
#if( adj_select[i,j]!=0 ){ box(lwd=3,col="red") } else{box(col=rgb(0,0,0,alpha=0.2))}
#text(-1.1,max(density( beta[ifelse(i>=j,i-1,i),j,])$y)*0.95,labels=bms[i,j],cex=1.2)
box(col=rgb(0,0,0,alpha=0.3))
}else{
  plot.new()
  nnamchar <- nchar(Pk.names[i])
  if(nnamchar>7){
    text(0.5,0.7,substring(Pk.names[i],1,7),cex=3)
    text(0.5,0.3,substring(Pk.names[i],8,nnamchar),cex=3)
  } else {
    text(0.5,0.5,Pk.names[i],cex=3)
  }
}
}
#mtext( expression(paste("posterior mean of ", beta) ),
#      NORTH<-3, line=1, adj=0.5, cex=4 , outer=TRUE)
dev.off()
#png("density_beta_rppa.png",width=5000,heigh=5000)
pdf("./plot/density_Osnp.pdf",width=80,heigh=80)
par(mfrow=c(Nsnp,Nsnp),oma=c(0,0,0,0),mgp=c(2,.5,0), tck=-.04)
for( i in 1:Nsnp) for(j in 1:Nsnp){
if(i!=j){
par(mar=c(1,1,0.2,0))
plot(density(Osnp[i,j,]),xlab="",ylab="",main="",xlim=range(Osnp[i,j,]),axes=FALSE, cex.main=3) #main=bms[i,j],
axis( side=1, at=0, label="", tick=TRUE,)
#if( adj_select[i,j]!=0 ){ box(lwd=3,col="red") } else{box(col=rgb(0,0,0,alpha=0.2))}
#text(-1.1,max(density( beta[ifelse(i>=j,i-1,i),j,])$y)*0.95,labels=bms[i,j],cex=1.2)
#box(col=rgb(0,0,0,alpha=0.3))
}else{
  plot.new()
  nnamchar <- nchar(names(SNP)[i])
  if(nnamchar>7){
    text(0.5,0.7,substring(names(SNP)[i],1,7),cex=3)
    text(0.5,0.3,substring(names(SNP)[i],8,nnamchar),cex=3)
  } else {
    text(0.5,0.5,names(SNP)[i],cex=3)
  }
}
}
#mtext( expression(paste("posterior mean of ", beta) ),
#      NORTH<-3, line=1, adj=0.5, cex=4 , outer=TRUE)
dev.off()


pdf("./plot/density_Rho.pdf",width=80,heigh=80)
par(mfrow=c(Nsnp,PK_V),oma=c(0,0,0,0),mgp=c(2,.5,0), tck=-.04)
for( i in 1:Nsnp) for(j in 1:PK_V){
#if(i!=j){
par(mar=c(1,1,0.2,0))
plot(density(Rho[i,j,]),xlab="",ylab="",main="",xlim=range(Rho[i,j,]),axes=FALSE, cex.main=3) #main=bms[i,j],
axis( side=1, at=0, label="", tick=TRUE,)
}
dev.off()

Rho.ci <- array( NA, c( Nsnp,PK_V, 5) )
for( i in 1:Nsnp){ 
  for(j in 1:PK_V){
    Rho.ci[i,j,]<- quantile(Rho[i,j,],c(0.01,0.025,0.5,0.975,0.99))
  }
}

pdf("./plot/lineplot_Rho.pdf",width=21,heigh=6)
par(mfrow=c(2,7), oma=c(0,0,0,0),mgp=c(2,.5,0), tck=-.04)
for(j in 1:PK_V){
#par(mar=c(1,1,0.2,0))
plot(c(1,Nsnp),range(Rho.ci[,j,]),type="n",xlab="",ylab="", cex.main=3, main=Pk.names[j])
for( i in 1:Nsnp) {
  cc = ifelse(sign(Rho.ci[i,j,1])!=sign(Rho.ci[i,j,5]), 1, 2)
  points( i, Rho.ci[i,j,3],col= cc , pch=19)
  lines( c(i,i), c(Rho.ci[i,j,2],Rho.ci[i,j,4]),lwd=3,col= cc )
  lines( c(i,i), c(Rho.ci[i,j,1],Rho.ci[i,j,5]),lwd=1,col= cc )
  abline(h=0,lty=2)
  #axis( side=1, at=0, label="", tick=TRUE,)
}
}
dev.off()

pcor.Otheta=Otheta
for(i in 1:Nmcmc){
  D=diag( 1/sqrt(diag(Otheta[,,i])))
  pcor.Otheta[,,i] = D%*%Otheta[,,i]%*%D
}
pcor.Otheta.mean = apply(pcor.Otheta,c(1,2),mean)
pcor.Otheta.median = apply(pcor.Otheta,c(1,2),median)
diag(pcor.Otheta.mean)=0
diag(pcor.Otheta.median)=0
pcor.Otheta.ci <- array( NA, c( PK_V,PK_V, 5) )
for( i in 1:PK_V){ 
  for(j in 1:PK_V){
    pcor.Otheta.ci[i,j,]<- quantile(pcor.Otheta[i,j,],c(0.01,0.025,0.5,0.975,0.99))
  }
}
sum(sign(pcor.Otheta.ci[,,1])==sign(pcor.Otheta.ci[,,5]))
sum(sign(pcor.Otheta.ci[,,2])==sign(pcor.Otheta.ci[,,4]))

pcor.Osnp=Osnp
for(i in 1:Nmcmc){
  D=diag( 1/sqrt(diag(Osnp[,,i])))
  pcor.Osnp[,,i] = D%*%Osnp[,,i]%*%D
}
pcor.Osnp.mean = apply(pcor.Osnp,c(1,2),mean)
diag(pcor.Osnp.mean)=0
pcor.Osnp.median = apply(pcor.Osnp,c(1,2),median)
diag(pcor.Osnp.median)=0
pcor.Osnp.ci <- array( NA, c( Nsnp,Nsnp, 5) )
for( i in 1:Nsnp){ 
  for(j in 1:Nsnp){
    pcor.Osnp.ci[i,j,]<- quantile(pcor.Osnp[i,j,],c(0.01,0.025,0.5,0.975,0.99))
  }
}
sum(sign(pcor.Osnp.ci[,,1])==sign(pcor.Osnp.ci[,,5]))
sum(sign(pcor.Osnp.ci[,,2])==sign(pcor.Osnp.ci[,,4]))

pdf("./plot/lineplot_pcorOsnp.pdf",width=21,heigh=20)
par(mfrow=c(6,7), oma=c(0,0,0,0),mgp=c(2,.5,0), tck=-.04)
for(j in 1:Nsnp){
#par(mar=c(1,1,0.2,0))
plot(c(1,Nsnp),range(pcor.Osnp.ci[,j,]),type="n",xlab="",ylab="", cex.main=3, main=names(SNP)[j])
for( i in 1:Nsnp) {
  if(i!=j){
  cc = ifelse(sign(pcor.Osnp.ci[i,j,1])!=sign(pcor.Osnp.ci[i,j,5]), 1, 2)
  points( i, pcor.Osnp.ci[i,j,3],col= cc , pch=19)
  lines( c(i,i), c(pcor.Osnp.ci[i,j,2],pcor.Osnp.ci[i,j,4]),lwd=3,col= cc )
  lines( c(i,i), c(pcor.Osnp.ci[i,j,1],pcor.Osnp.ci[i,j,5]),lwd=1,col= cc )
  abline(h=0,lty=2)
}
  #axis( side=1, at=0, label="", tick=TRUE,)
}
}
dev.off()
# rescaled rho
sRho=Rho*0
for(i in 1:Nmcmc){
  sRho[,,i] = Rho[,,i]*matrix(apply(Z[,,i],2,sd),Nsnp,PK_V)
}
sRho.ci <- array( NA, c( Nsnp,PK_V, 5) )
for( i in 1:Nsnp){ 
  for(j in 1:PK_V){
    sRho.ci[i,j,]<- quantile(sRho[i,j,],c(0.01,0.025,0.5,0.975,0.99))
  }
}
pdf("./plot/lineplot_sRho.pdf",width=21,heigh=6)
par(mfrow=c(2,7), oma=c(0,0,0,0),mgp=c(2,.5,0), tck=-.04)
for(j in 1:PK_V){
#par(mar=c(1,1,0.2,0))
plot(c(1,Nsnp),range(sRho.ci[,j,]),type="n",xlab="",ylab="", cex.main=3, main=Pk.names[j])
for( i in 1:Nsnp) {
  cc = ifelse(sign(sRho.ci[i,j,1])==sign(sRho.ci[i,j,5]), "red",ifelse(sign(sRho.ci[i,j,2])==sign(sRho.ci[i,j,4]), "orange","black"))
  points( i, sRho.ci[i,j,3],col= cc , pch=19)
  lines( c(i,i), c(sRho.ci[i,j,2],sRho.ci[i,j,4]),lwd=3,col= cc )
  lines( c(i,i), c(sRho.ci[i,j,1],sRho.ci[i,j,5]),lwd=1,col= cc )
  abline(h=0,lty=2)
  #axis( side=1, at=0, label="", tick=TRUE,)
}
}
dev.off()
sRho.mean <- apply(sRho,c(1,2),mean)
sRho.median <- apply(sRho,c(1,2),median)


sRho2=Rho*0
for(i in 1:Nmcmc){
  sRho2[,,i] = matrix(1/apply(log(Theta[,,i]),2,sd),Nsnp,PK_V,byrow=T)*Rho[,,i]*matrix(apply(Z[,,i],2,sd),Nsnp,PK_V)
}
sRho2.ci <- array( NA, c( Nsnp,PK_V, 5) )
for( i in 1:Nsnp){ 
  for(j in 1:PK_V){
    sRho2.ci[i,j,]<- quantile(sRho2[i,j,],c(0.01,0.025,0.5,0.975,0.99))
  }
}
sRho2.mean <- apply(sRho2,c(1,2),mean)
sRho2.median <- apply(sRho2,c(1,2),median)
pdf("lineplot_sRho2.pdf",width=21,heigh=6)
par(mfrow=c(2,7), oma=c(0,0,0,0),mgp=c(2,.5,0), tck=-.04)
for(j in 1:PK_V){
#par(mar=c(1,1,0.2,0))
plot(c(1,Nsnp),range(sRho2.ci[,j,]),type="n",xlab="",ylab="", cex.main=3, main=Pk.names[j])
for( i in 1:Nsnp) {
  cc = ifelse(sign(sRho2.ci[i,j,1])==sign(sRho2.ci[i,j,5]), "red",ifelse(sign(sRho2.ci[i,j,2])==sign(sRho2.ci[i,j,4]), "orange","black"))
  points( i, sRho2.ci[i,j,3],col= cc , pch=19)
  lines( c(i,i), c(sRho2.ci[i,j,2],sRho2.ci[i,j,4]),lwd=3,col= cc )
  lines( c(i,i), c(sRho2.ci[i,j,1],sRho2.ci[i,j,5]),lwd=1,col= cc )
  abline(h=0,lty=2)
  #axis( side=1, at=0, label="", tick=TRUE,)
}
}
dev.off()
library(sna)
library(network)
Otheta_select <- ifelse(sign(pcor.Otheta.ci[,,2])==sign(pcor.Otheta.ci[,,4]), 1, 0)#abs(pcor.Otheta.mean)>0.2
Osnp_select   <- ifelse(sign(pcor.Osnp.ci[,,2])==sign(pcor.Osnp.ci[,,4]), 1, 0)#abs(pcor.Osnp.mean)>0.5
Rho_select    <- ifelse(sign(sRho.ci[,,2])==sign(sRho.ci[,,4]), 1, 0)  #abs(Rho.mean)>0.1
dimnames(Rho_select)<-list(names(SNP), Pk.names)
CombGraph <- matrix(0, PK_V+Nsnp, PK_V+Nsnp)
CombGraph[1:PK_V,1:PK_V]=Otheta_select
CombGraph[(PK_V+1):(PK_V+Nsnp),(PK_V+1):(PK_V+Nsnp)]=Osnp_select
CombGraph[(PK_V+1):(PK_V+Nsnp), 1:PK_V]=Rho_select
RhoGraph <- matrix(0, PK_V+Nsnp, PK_V+Nsnp)
RhoGraph[(PK_V+1):(PK_V+Nsnp), 1:PK_V]=Rho_select
dimnames(CombGraph)<-list(c(Pk.names,names(SNP)),c(Pk.names,names(SNP)))
dimnames(RhoGraph)<-list(c(Pk.names,names(SNP)),c(Pk.names,names(SNP)))
network_Otheta <- network(Otheta_select, directed=FALSE)
network_Osnp   <- network(Osnp_select, directed=FALSE)
  network_Comb <- network(CombGraph, directed=TRUE)
  network_Rho <- network(RhoGraph, directed=TRUE)
  col_Comb <- c(rep("white",PK_V),rep("black",Nsnp))
  pdf("graph.pdf",width=10,height=7)
   cord_Comb <-gplot(network_Comb, gmode = "graph", label=dimnames(CombGraph)[[1]],displaylabels=T,main="",vertex.col=col_Comb ,cex.main=3,displayisolates = FALSE,vertex.cex = 2,edge.col =color[za])
        rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = "gray85") 
  gplot(network_Comb, gmode = "graph", label=dimnames(CombGraph)[[1]],displaylabels=T,main="",vertex.col=col_Comb ,cex.main=3,coord=cord_Comb, displayisolates = FALSE,vertex.cex = 2,edge.col =color[za], new=FALSE)

  cord_Rho  <-gplot(network_Rho,  gmode = "digraph",displaylabels=F,main="",vertex.col=rgb(0,0,0,alpha=0),vertex.border = 0,coord=cord_Comb, new=FALSE,displayisolates = FALSE,vertex.cex = 2,arrowhead.cex = 2,edge.col =color[za])
  dev.off()

CombGraphCol <- matrix(0, PK_V+Nsnp, PK_V+Nsnp)
CombGraphCol[1:PK_V,1:PK_V]=pcor.Otheta.median
CombGraphCol[(PK_V+1):(PK_V+Nsnp),(PK_V+1):(PK_V+Nsnp)]=pcor.Osnp.median
CombGraphCol[(PK_V+1):(PK_V+Nsnp), 1:PK_V]=sRho.median
color = rev( brewer.pal( 5, "RdBu" ))
color = color[-3]
z.breaks = seq( -0.6, 0.6, length.out = 4 + 1 )
za = array( as.double( cut( CombGraphCol, breaks = z.breaks, labels = 1:4) ), dim( CombGraphCol ) )
  
Rho_selected_only = Rho_select[rowSums(Rho_select)>0,]
a=hclust(dist( Rho_select[rowSums(Rho_select)>0,],method="manhattan"))
Rho_selected_only[a$order,]
Rho_select_mean <- Rho.mean*Rho_select
z.breaks = seq( -0.5, 0.5, length.out = 11 + 1 )
za = array( as.double( cut( Rho_select_mean, breaks = z.breaks, labels = 1:11 ) ), dim( Rho.mean ) )
color = rev( brewer.pal( 11, "RdBu" ))
adj_col <- array(color[za],dim( Rho.mean ))

pdf("network.pdf", width=14, height=7)
#par(mfrow=c(1,3))
layout( matrix( c( 1, 2, 3 ), 1, 3, byrow = FALSE ), c( 10.5, 10.5, 1.5 ) )
cord <- gplot(nbase,label=names(Y),displaylabels=T,main="refractory patients",vertex.col="white",edge.col=(adj_col),edge.lty = adj_lty,cex.main=3)                # A less useful layout....
gplot(ndiff,label=names(Y),displaylabels=T,main="relapsed patients", coord = cord,vertex.col="white",edge.col=(adjz_col),edge.lty =adjz_lty,cex.main=3)                # A less useful layout....

par.label = list( mar = c( 0.5, 0.1, 2, 0.1 ), pty = "m" )
    op <- par( names( par.label ) )
    par( par.label )
    plot( c( 0, 1 ), c( min( z.breaks ), max( z.breaks ) ), type = "n",
        bty = "n", xlab = "", ylab = "", xaxt = "n", yaxt = "n" )
    for( i in 2:( length( z.breaks ) ) ) {
        rect( xleft = 0.5, ybottom = z.breaks[i - 1], 
              xright = 1, ytop = z.breaks[i], col = color[i - 1] )
        text( x = 0.45, y = z.breaks[i - 1], 
              labels = format(round(z.breaks[i - 1], 2)), 
              cex =  1, adj = 1, xpd = TRUE)
    }
    rect( xleft = 0.5, ybottom = z.breaks[length( z.breaks )], 
          xright = 1, ytop = z.breaks[length( z.breaks )], col = color[length( color )])
    text( x = 0.45, y = z.breaks[length( z.breaks )], 
          labels = format(round(z.breaks[length( z.breaks )], 2 ) ), 
          cex =  1, adj = 1, xpd = TRUE )
    par( op )
dev.off()
# pdf("network.pdf", width=14, height=7)
# #par(mfrow=c(1,3))
# layout( matrix( c( 1, 2, 3 ), 1, 3, byrow = FALSE ), c( 10.5, 10.5, 1.5 ) )
# cord_theta <- gplot(network_Otheta ,gmode = "graph", label=Pk.names,displaylabels=T,main="",vertex.col="white",cex.main=3) #,edge.col=(adj_col),edge.lty = adj_lty)                # A less useful layout....
# cord_snp <- gplot( network_Osnp ,gmode = "graph", label=names(SNP),displaylabels=T,main="",vertex.col="white",cex.main=3) #,edge.col=(adj_col),edge.lty = adj_lty)                # A less useful layout....
# cord_snp_shift=cord_snp
# cord_snp_shift[,1]= cord_snp_shift[,1]+20
# xrange=range(c(cord_theta[,1],cord_snp_shift[,1]))
# yrange=range(c(cord_theta[,2],cord_snp_shift[,2]))
# gplot(network_Otheta ,gmode = "graph", label=Pk.names,coord = cord_theta, displaylabels=T,main="",vertex.col="white",cex.main=3,xlim=xrange,ylim=yrange) #,edge.col=(adj_col),edge.lty = adj_lty)                # A less useful layout....
# gplot( network_Osnp ,gmode = "graph", label=names(SNP),coord = cord_snp_shift,displaylabels=T,main="",vertex.col="white",cex.main=3,new=FALSE,xlim=xrange,ylim=yrange) #,edge.col=(adj_col),edge.lty = adj_lty)                # A less useful layout....

# gplot(ndiff,label=names(Y),displaylabels=T,main="", coord = cord,vertex.col="white",cex.main=3) #,edge.col=(adjz_col),edge.lty =adjz_lty,cex.main=3)                # A less useful layout....

# par.label = list( mar = c( 0.5, 0.1, 2, 0.1 ), pty = "m" )
#     op <- par( names( par.label ) )
#     par( par.label )
#     plot( c( 0, 1 ), c( min( z.breaks ), max( z.breaks ) ), type = "n",
#         bty = "n", xlab = "", ylab = "", xaxt = "n", yaxt = "n" )
#     for( i in 2:( length( z.breaks ) ) ) {
#         rect( xleft = 0.5, ybottom = z.breaks[i - 1], 
#               xright = 1, ytop = z.breaks[i], col = color[i - 1] )
#         text( x = 0.45, y = z.breaks[i - 1], 
#               labels = format(round(z.breaks[i - 1], 2)), 
#               cex =  1, adj = 1, xpd = TRUE)
#     }
#     rect( xleft = 0.5, ybottom = z.breaks[length( z.breaks )], 
#           xright = 1, ytop = z.breaks[length( z.breaks )], col = color[length( color )])
#     text( x = 0.45, y = z.breaks[length( z.breaks )], 
#           labels = format(round(z.breaks[length( z.breaks )], 2 ) ), 
#           cex =  1, adj = 1, xpd = TRUE )
#     par( op )

Stheta0 <- round( cov(log(Theta0[-c(15,33,57,60,68,73,83),]+0.01))*5.6644/dim(Theta0[,])[2],2)

Stheta0 <- round( cov(log(Theta0+0.01))*5.6644/10,2)
rsample<-exp( rmvnorm(1000,colMeans(log(Theta0[,]+0.001)),Stheta0))
  pdf(paste("./plot/ParcoordTheta_001.pdf",sep=""),width=16,height=4)
      parcoord(rbind(rsample, Theta0[,]),col=c(rep( rgb(0,0,1,alpha=0.3),1000),rep(rgb(1,0,0,alpha=0.2),nrow(Theta0))),var.label=TRUE)
      dev.off()

# Residual plot
Resid=Fit*0
for(i in 1:N){
  for(j in 1:Nmcmc){
    Resid[i,,,j]=Fit[i,,,j]-ConcArray[,,i]
  }
}


for(i in 1:N){
  pdf(paste("resid_plot_",i,".pdf"),width=16,height=4)
  par(mfrow=c(1,4))
  for(k in 1:K){
    plot( range(Time[i,]), range(Resid[i,,k,]),type="n")
    for(j in 1:Nmcmc){
      points(Time[i,],Resid[i,,k,j],col=rgb(0,0,1,alpha=0.2))
    }
  }
  dev.off()
}

# eps <- 1e-2
# pvalue.pcor.Otheta <- array( NA, c( PK_V,PK_V) )
# for( i in 1:PK_V){ 
#   for(j in 1:PK_V){
#     pvalue.pcor.Otheta[i,j] <- sum(abs(pcor.Otheta[i,j,])<eps)/Nmcmc
#   }
# }
# pOtheta=pvalue.pcor.Otheta[upper.tri(pvalue.pcor.Otheta)]
# mOtheta=pcor.Otheta.median[upper.tri(pcor.Otheta.median)]
# pOtheta.ord <- pOtheta[order(mOtheta)]
# mOtheta.ord <- mOtheta[order(mOtheta)]

# ell.Otheta<- max((1:length(pOtheta.ord ))[cumsum( pOtheta.ord ) <=  1:length(pOtheta.ord ) *0.05])
# mell.Otheta <- mOtheta.ord[ell.Otheta]
# pcor.Otheta.mean>mell.Otheta
# pcor.Otheta.select <- 1*(abs(pcor.Otheta.mean)>abs(mell.Otheta))


# eps <- 1e-2
# pvalue.pcor.Osnp <- array( NA, c( Nsnp,Nsnp) )
# for( i in 1:Nsnp){ 
#   for(j in 1:Nsnp){
#     if(i<j) { pvalue.pcor.Osnp[i,j] <- sum(abs(pcor.Osnp[i,j,])<eps)/Nmcmc
#       pvalue.pcor.Osnp[j,i] <-pvalue.pcor.Osnp[i,j]
#     }
#   }
# }
# pOsnp=pvalue.pcor.Osnp[upper.tri(pvalue.pcor.Osnp)]
# mOsnp=pcor.Osnp.median[upper.tri(pcor.Osnp.median)]
# pOsnp.ord <- pOsnp[order(mOsnp)]
# mOsnp.ord <- mOsnp[order(mOsnp)]

# ell.Osnp<- max((1:length(pOsnp.ord ))[cumsum( pOsnp.ord ) <=  1:length(pOsnp.ord ) *0.05])
# mell.Osnp <- mOsnp.ord[ell.Osnp]
# pcor.Osnp.mean>mell.Osnp
# pcor.Osnp.select <- 1*(abs(pcor.Osnp.mean)>abs(mell.Osnp))


#   eps.Rho   <- c(1e-6,1e-5,1e-4,1e-3, 1e-2, 0.025, 0.05, 0.1, 0.2,0.5 ) 
#   alpha.Rho <- c(1e-6,1e-5,1e-4,1e-3, 1e-2, 0.025, 0.05, 0.1, 0.2,0.5 )
#   pvalue.sRho <- array( NA, c( Nsnp,PK_V, length(eps.Rho) ) )
#   sRho.select <- array( NA, c( Nsnp,PK_V, length(eps.Rho), length(alpha.Rho) ) )
#   for( k in 1:length(eps.Rho) ){
#     for( i in 1:Nsnp){ 
#       for(j in 1:PK_V){
#         pvalue.sRho[i,j,k] <- sum( abs(sRho[i,j,])<eps.Rho[k])/Nmcmc
#       }
#     }
#   }
#   # dimnames(pvalue.sRho)<-list(names(SNP), Pk.names)
#   # dimnames(sRho.mean)<-list(names(SNP), Pk.names)
#   # dimnames(sRho.median)<-list(names(SNP), Pk.names)
#   for( k in 1:length(eps.Rho) ){
#     for( u in 1:length(alpha.Rho) ){
#   psRho=as.double(pvalue.sRho[,,k])
#   msRho=as.double(sRho.median)
#   psRho.ord <- psRho[order(abs(msRho),decreasing = TRUE)]
#   msRho.ord <- msRho[order(abs(msRho),decreasing = TRUE)]

#   ell.sRho<- max((1:length(psRho.ord ))[cumsum( psRho.ord ) <=  (1:length( psRho.ord )) * alpha.Rho[u]])
#   mell.sRho <- msRho.ord[ell.sRho]
#   pell.sRho <- psRho.ord[ell.sRho]
#   if(is.na(mell.sRho)) mell.sRho <-Inf
#   sRho.select[,,k,u] <- 1*(abs(sRho.median)>abs(mell.sRho))
# }
# }
# select.grid.Rho <- apply(sRho.select,c(3,4),sum)
# dimnames(select.grid.Rho)<-list(eps.Rho,alpha.Rho)

apply(sRho.select,c(1,2),sum)
  #sRho.select <- 1*(pvalue.sRho<pell.sRho)
  dimnames(sRho.select)<-list(names(SNP), Pk.names)
CombGraph <- matrix(0, PK_V+Nsnp, PK_V+Nsnp)
CombGraph[1:PK_V,1:PK_V]=pcor.Otheta.select
CombGraph[(PK_V+1):(PK_V+Nsnp),(PK_V+1):(PK_V+Nsnp)]=pcor.Osnp.select
CombGraph[(PK_V+1):(PK_V+Nsnp), 1:PK_V]=sRho.select
RhoGraph <- matrix(0, PK_V+Nsnp, PK_V+Nsnp)
RhoGraph[(PK_V+1):(PK_V+Nsnp), 1:PK_V]=sRho.select
dimnames(CombGraph)<-list(c(Pk.names,names(SNP)),c(Pk.names,names(SNP)))
dimnames(RhoGraph)<-list(c(Pk.names,names(SNP)),c(Pk.names,names(SNP)))
  network_Comb <- network(CombGraph, directed=TRUE)
  network_Rho <- network(RhoGraph, directed=TRUE)
  col_Comb <- c(rep("white",PK_V),rep("black",Nsnp))
  pdf("graph.pdf",width=10,height=9)
   cord_Comb <-gplot(network_Comb, gmode = "graph", label=dimnames(CombGraph)[[1]],displaylabels=T,main="",vertex.col=col_Comb ,cex.main=3,displayisolates = FALSE,vertex.cex = 2,edge.col =color[za])
        rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = "gray85") 
  gplot(network_Comb, gmode = "graph", label=dimnames(CombGraph)[[1]],displaylabels=T,main="",vertex.col=col_Comb ,cex.main=3,coord=cord_Comb, displayisolates = FALSE,vertex.cex = 2,edge.col =color[za], new=FALSE)

  cord_Rho  <-gplot(network_Rho,  gmode = "digraph",displaylabels=F,main="",vertex.col=rgb(0,0,0,alpha=0),vertex.border = 0,coord=cord_Comb, new=FALSE,displayisolates = FALSE,vertex.cex = 2,arrowhead.cex = 2,edge.col =color[za])
  dev.off()

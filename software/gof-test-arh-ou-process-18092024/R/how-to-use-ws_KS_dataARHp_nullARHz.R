
# load auxiliary functions and packages
rm(list = ls())
library(dplyr)
library(fda.usc)
library(goffda)
library(sde)
library(lubridate)
source("./R/OU_estimation_test.R")

# ----- ARH GENERATION -----

std_bb <- function(t = seq(0, 1, l = 101)) {
  
  bb_path <-
    as.numeric(rproc2fdata(1, t = t, sigma = "brownian")$data)
  return(bb_path)
}
BrowBridge <- function(t = seq(0, 1, l = 101)) {
  nt <- length(t) 
  I  <- rnorm(nt - 1) / sqrt(nt - 1) 
  B  <- c(0, cumsum(I)) 
  BB <- B - t * B[length(t)]
  return(BB)
}
ARH <- function(z = 1, burn = 200, n = 50, t = seq(0, 1, l = 101), c = NULL, X0) {
  
  X <- matrix(NA, n + burn, length(t))
  X[1, ] <- X0
  
  
  for (i in 2:(n + burn)) {
    
    if (z == 2 & i == 2) {
      
      X[i, ] <- BrowBridge(t)
      
    } else {
      
      bb <- BrowBridge(t)
      X[i, ] <- bb # noise
      
      # phi1 <- t |> map_dfc(function(x) { phi_ARH1(x, s = t, c = 0.500568) })
      # X[i, ] <-  X[i, ] + colSums(phi1 * X[i - 1, ]) / (length(t) - 1)
      
      for (j in 1:length(t)) {
        
        if (z == 1) { # ARH1
          
          X[i, j] <-
            X[i, j] + (sum(phi_ARH1(t = t[j], s = t, c = c) * X[i - 1, ]) /
                         (length(t) - 1))
          
          
          
        } else if (z == 2) { # ARH2
          
          X[i,j] <- X[i,j] + 
            (sum(phi_ARH2(t = t[j], s = t, c = c)[1, ] * X[i - 1, ]) / (length(t) - 1)) +
            (sum(phi_ARH2(t = t[j], s = t, c = c)[2, ] * X[i - 2, ]) / (length(t) - 1))
          
        } 
      }
    }
  }
  X <- fdata(X, argvals = t)
  return(X)
  
}
phi_ARH1 <- function(t, s, c = 0.5005669) {
  
  beta <- c*(2 - (2*t - 1)^2 - (2*s - 1)^2)
  return(beta)
  
}
phi_ARH2 <- function(t, s, c = c(0.669502, 0.40170)) {
  
  beta <- matrix(NA, 2, max(length(t), length(s)))
  beta[1, ] <- c[1]*exp(-((t^2 + s^2) / 2))
  beta[2, ] <- c[2]*exp(-((t^2 + s^2) / 2))
  return(beta)
  
}


# ----- FAR TEST -----

Order_FAR<-function(X, p_all = 4, npc = 3){
  
  nharm<-30
  if(any(grepl("fd", class(X)))) {
    M<-101
    t_eval<-seq(0,1,length=M)
    Xg = t(eval.fd(t_eval,X))
    nb<-dim(X$coefs)[1]
    fbf<-X$basis
  }else if(any(grepl("matrix", class(X)))){
    M<-dim(X)[2]
    N<-dim(X)[1]
    t_eval<-seq(0,1,length=M)
    nb<-100
    fbf = create.bspline.basis(rangeval=c(0,1), nbasis=nb) 
    Xg = X
  }else{
    stop("X needs to be a matrix or functional")
  }
  Test<-numeric(0)
  pval<-numeric(0)
  if(length(p_all) == 1){
    p_vec = 1:p_all
  }else{p_vec = p_all}
  for(p in p_vec){
    Xgp<-Xg[p:(N-1),]
    if( p > 1){
      for(j in 2:p){Xgp<-cbind(Xgp,Xg[(p-j+1):(N-j),])}
    }
    t_eval2 = seq(0,1,length=p*M)
    Xgfo = Data2fd(t_eval2, t(Xgp),  fbf)
    Ygfo = Data2fd(t_eval, t(Xg[(p+1):N,]),  fbf)
    Xgfo<-center.fd(Xgfo)
    Ygfo<-center.fd(Ygfo)
    
    qy<-npc
    qx<-qy*p
    
    
    X.fd<-pca.fd(Xgfo,qx)
    Y.fd<-pca.fd(Ygfo,qy)
    
    X<-X.fd$scores
    Y<-Y.fd$scores
    lm1<-lm(Y~X)
    Ceps<-(1/(N-p-qx-1))*(t(lm1$res)%*%lm1$res)
    
    M0<-400 #points to approx integral of harmonics over [(p-1)/(p),1] interval
    t.tmp<-seq(from=(p-1)/(p), to = 1, length.out = (M0+1))
    t.tmp<-t.tmp[2:(M0+1)]
    v.tmp<-eval.fd(t.tmp,X.fd$harmonics)
    V<-t(t(v.tmp)%*%v.tmp)/(M*p)
    eg.v<-eigen(V)
    
    qv<-length(which(eg.v$values>0.9))
    alpha<-eg.v$vectors[,(1:qv)]
    lambda<-diag(t(X)%*%X)/N
    Psi<-diag((N*lambda)^(-1))%*%(t(X)%*%Y)
    I_qy<-diag(rep(1,times=qy))
    
    av_tmp<-as.vector(t(alpha)%*%Psi)
    Tst.Stat<-(N-p)*av_tmp%*%solve( (I_qy%x%t(alpha))%*%(Ceps%x%diag(1/lambda))%*%(I_qy%x%alpha) )%*%av_tmp
    
    Test<-c(Test,Tst.Stat)
    pval<-c(pval,1-pchisq(Tst.Stat,df=(qv*qy)))
  }	
  Sig<-rep('',times=length(p_vec))
  Sig[pval < 0.1] = '*'
  Sig[pval < 0.05] = '**'
  Sig[pval < 0.01] = '***'
  Table<-data.frame(p_vec,Test,pval,Sig)
  names(Table) <-c("Order under HA","Test Stat","P-Value","Sig")
  return(Table)
}

# ----- SIMULATION PARAMETERS -----

burn <- 200
n <- c(50, 150, 250, 350, 500, 750)
B <- 1500
MC <- 3
p.value <- array(NA, c(length(n), MC, 3))
t <- seq(0, 1, l = 101)
iterations <- MC * length(n)
pb <- txtProgressBar(max = iterations, style = 3)
progress <- function(n) { setTxtProgressBar(pb, n) }
opts <- list(progress = progress)
t_init <- Sys.time()

# ----- NULL PARAMETERS -----
z <- 1:3
p_data <- 0

# ----- SIMULATIONS -----
for (j in 1:length(n)) {
  for (m in 1:MC) {
    
    # ARHp generation
    X_ARHp <- ARH(z = p_data, burn = burn, n = n[j],
                  t = t, X0 = BrowBridge(t), c = 0.5005669)
    X_ARHp <- X_ARHp[(burn + 1):(burn + n[j])]
    
    cat(paste0("\nSample size ", n[j], ", iteration ", m))
    FAR_test <- Order_FAR(X_ARHp$data, z, npc = min(30, n[j]/5))
    
    p.value[j, m, ] <- FAR_test$`P-Value`
    
  }
}


# collecting and measuring time
t_end <- Sys.time()
t_lapse <- t_end - t_init


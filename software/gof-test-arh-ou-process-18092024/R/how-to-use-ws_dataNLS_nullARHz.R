
# ----------------------
# DATA: NLS
# NULL: ARH(z)
# ----------------------

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

ARH_NLS <- function(burn = 200, n = 50, t = seq(0, 1, l = 101), c = 0.5, X0) {
  
    phi <- function(t, s) 0.669502 * exp(-((t^2 + s^2) / 2)) 

    X <- matrix(NA, n + burn, length(t))
    X[1, ] <- X0
    for (i in 2:(n + burn)) {
      
      bb <- BrowBridge(t) # noise

      for (j in 1:length(t)) {
        
        X[i, j] <- (sum(phi(t[j], t) * sqrt(abs(X[i - 1, ]))) / (length(t) - 1)) + bb[j]
      }
    }
    
    X <- fdata(X, argvals = t)
    return(X)
  
}


## ----- SIMULATION PARAMETERS -----

burn <- 200
n <- c(50, 150, 250, 350, 500, 750)
B <- 1500
MC <- 3
p.value <- matrix(NA, length(n), MC)
t <- seq(0, 1, l = 101)
iterations <- MC * length(n)
pb <- txtProgressBar(max = iterations, style = 3)
progress <- function(n) { setTxtProgressBar(pb, n) }
opts <- list(progress = progress)
t_init <- Sys.time()

# ----- NULL PARAMETERS -----
z_null <- 0


# ----- SIMULATIONS -----
for (j in 1:length(n)) {
  for (m in 1:MC) {
    
    # ARH NLS generation
    X_NLS <- ARH_NLS(burn = burn, n = n[j], t = t, X0 = BrowBridge(t), c = 0.5)
    X_NLS <- list("fstoch_proc" = X_NLS[(burn + 1):(burn + n[j])])
    
    # STAGE 1: Is X_n an ARH(z) process?
    cat(paste0("\nSample size ", n[j], ", iteration ", m))
    testing_ARHz <- ARHz_test(X_NLS, z = z_null, B = B, hyp_simp = FALSE,
                              hyp_comp = TRUE, est_method = "fpcr_l1s",
                              thre_p = 0.995, thre_q = 0.995, lambda = NULL,
                              boot_scores = TRUE, plot_dens = FALSE,
                              plot_proc = FALSE, save_fit_flm = FALSE,
                              save_boot_stats = FALSE, cv_1se = TRUE)
    p.value[j, m] <- testing_ARHz$testing_comp$p.value
    
  }
}


# collecting and measuring time
t_end <- Sys.time()
t_lapse <- t_end - t_init


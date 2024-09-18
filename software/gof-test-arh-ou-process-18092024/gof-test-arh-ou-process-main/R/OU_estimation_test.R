# ###############################################################
# Description: companion software of the article
#              «A goodness-of-fit test for functional time series
#              with applications to diffusion processes»
# Authors: A. López-Pérez and J. Álvarez-Liébana
# ###############################################################

# ###########################
# PACKAGES
# ###########################
rm(list = ls())
# setwd(dirname(rstudioapi::getSourceEditorContext()$path))
repos <- "http://cran.us.r-project.org"
if(!require(goffda)) install.packages("goffda", repos = repos)
if(!require(ggplot2)) install.packages("ggplot2", repos = repos)
if(!require(fda.usc)) install.packages("fda.usc", repos = repos)
# if(!require(tidyverse)) install.packages("tidyverse", repos = repos)
if(!require(sde)) install.packages("sde", repos = repos)


# ###########################
# AUXILIARY FUNCTIONS
# ###########################


# Description: generate functional stochastic processes
#
# Example:
# OU <- r_stoch_proc(500, seq(0, 1, l = 101),
#                    par.sde = list("alpha" = 0, "beta" = 0.5, "sigma" = 0.1),
#                    warm_up = -1, model = "OU")
r_stoch_proc <- function(n, t, par.sde = list("alpha" = 0, "beta" = 0.5,
                                              "sigma" = 1),
                         mu = par.sde[["alpha"]] / par.sde[["beta"]], X0 = mu,
                         warm_up = 50, type = "CKLS", model = "OU",
                         delta = 1/(length(t) - 1), verbose = FALSE,
                         plot = FALSE, drf = NULL, sig = NULL) {
  
  # Generating from a CKLS model
  # CKLS model: dX_t = (alpha - beta*X_t) dt + sigma * X_t^gamma dW_t
  #
  # OU model: dX_t = (alpha - beta*X_t) dt + sigma dW_t
  #           gamma = 0; beta = theta; alpha = theta * mu
  #
  if (type == "CKLS") {
    
    # X0 as initial condition at t0, N the number of points to be generated
    # between [t0, n * t1], delta denotes the discretization step and warm_up
    # the number of trajectorias to be discarded (note that N + 1 are generated)
    suppressMessages(
      sde_proc <- sde.sim(X0 = X0, theta = c(par.sde$alpha, par.sde$beta,
                                             par.sde$sigma),
                          N = (n * length(t)) + warm_up, t0 = t[1],
                          T = n * t[length(t)], delta = delta, model = model)
    )
    
    if (warm_up >= 0) {
      
      sde_proc <- sde_proc[-(1:(warm_up + 1))]
      
    }
  } else if (type == "other") {
    
    suppressMessages(
      sde_proc <- sde.sim(X0 = X0, N = (n * length(t)) + warm_up, t0 = t[1],
                          T = n * t[length(t)], delta = delta, drift = drf, sigma = sig)
    )
    
  }
  
  # Centered SDE
  mu_est    <- mean(sde_proc)
  sde_proc  <- (sde_proc - mu_est)
  
  # SDE as a functional object with n = 1
  stoch_proc <- fdata(t(as.matrix(sde_proc)),
                      argvals = seq(t[1], t[length(t)] * n, l = length(t) * n))
  
  # SDE as a set of functional trajectories (assumed to be centered)
  fstoch_proc <- fdata(t(matrix(sde_proc, nrow = length(t))), argvals = t)
  
  # 
  # Plot the functional trajectories
  if (plot) {
    
    par(mfrow = c(1, 2))
    plot(stoch_proc, main = "SDE")
    plot(fstoch_proc, main = "SDE as functional data")
    par(mfrow = c(1, 1))
    
  }
  
  # Output
  return(list("fstoch_proc" = fstoch_proc, "stoch_proc" = stoch_proc))
  
}


# Description: function to be optimized for estimating the
#              volatility (sigma) of OU process, given an OU process
#              as a vectorial one-dimensional stochastic process
#
# Example:
# OU <- r_stoch_proc(500, seq(0, 1, l = 101),
#                    par.sde = list("alpha" = 0, "beta" = 0.5, "sigma" = 0.1),
#                    warm_up = -1, model = "OU")
# optimize(mce, r = as.vector(OU$stoch_proc$data),
#          Delta = 1/(length(OU$fstoch_proc$argvals) - 1),
#          interval = c(0, 2))$minimum
mce <- function(sigma, r, Delta){
  y    <- diff(r)
  n    <- length(y)
  suma <- sum(log(sigma^2) + (y^2) / (Delta * sigma^2))
  return(suma / n)
}


# Description: function for computing the MLE estimator of theta,
# given an OU process into a fdata (fda.usc package) object with n = 1
#
# Example:
# OU <- r_stoch_proc(500, seq(0, 1, l = 101),
#                    par.sde = list("alpha" = 0, "beta" = 0.5, "sigma" = 1),
#                    warm_up = -1, model = "OU")
# MLE_theta(OU$stoch_proc)
MLE_theta <- function(sde_OU) {
  
  # Grid points where the whole OU takes values: a set of subintervals 
  t_OU <- sde_OU[["argvals"]] 
  
  # Maximum Likelihood Estimator (MLE) of theta
  MLE_theta <- sum(as.vector(sde_OU[["data"]])[2:length(t_OU)] *
                     diff(as.vector(sde_OU[["data"]]))) /
    goffda::integral1D(fx = as.vector(sde_OU[["data"]])^2, t = t_OU)
  
  # Output: MLE of theta
  return(MLE_theta)
  
}


# Description: function for converting functional stochastic process
# into a FLMFR given by Y = X_n vs X = X_{n-p}
# 
# Example:
# OU <- r_stoch_proc(500, seq(0, 1, l = 101),
#                    par.sde = list("alpha" = 0, "beta" = 0.5,
#                                   "sigma" = 0.1),
#                    warm_up = -1, model = "OU")
# OU_FLMFR <- ARH_to_FLMFR(ARHz = OU$fstoch_proc, z = 1)
# plot(OU_FLMFR$X_fdata) # Plotting functional X variable
# plot(OU_FLMFR$Y_fdata) # Plotting functional Y variable
ARH_to_FLMFR <- function(ARHz, z, centered = FALSE) {
  
  # Common settings: sample size and number of grids where it is evaluated
  n <- dim(ARHz)[1]
  n_grids <- length(ARHz[["argvals"]])
  
  # Data should be centered to be transformed
  if (!centered) {
    
    ARHz <- ARHz - func_mean(ARHz)
    
  }
  
  # We will remove the first z trajectories since Y_n = X_{n + 1}
  Y <- ARHz[(z + 1):n, ]
  
  if (z == 0) {
    
    X <- Y
    
  } else if (z == 1) {
    
    X <- ARHz[1:(n - 1), ]
    
  } else if (z > 1) {
    
    X <- fda.usc::fdata(matrix(0, n - z, n_grids), argvals = ARHz[["argvals"]])
    for (j in 1:z) {
      
      # s <- (1:n_grids)[ARHz[["argvals"]] >= (j - 1)/z & ARHz[["argvals"]] <= j/z]
      s <- ARHz[["argvals"]][ARHz[["argvals"]] >= (j - 1)/z & ARHz[["argvals"]] <= j/z]
      scaled_s <- round(ARHz[["argvals"]], z) %in% round((s * z - (j - 1)), z)
      # X[["data"]] <- X[["data"]] + ARHz[["data"]][1:(n - j), s * z - (j - 1)]
      new_X <- matrix(0, n - z, n_grids)
      new_X[, ARHz[["argvals"]] %in% s] <- ARHz[["data"]][1:(n - z),  scaled_s]
      X[["data"]] <- X[["data"]] + new_X# ARHz[["data"]][1:(n - j),  scaled_s]
      
    }
  }
  
  # Output
  return(list("X_fdata" = X, "Y_fdata" = Y))
}


# Description: function for predicting an OU as an ARH(1) process
#
# Example
# OU <- r_stoch_proc(1500, seq(0, 1, l = 101),
#                    par.sde = list("alpha" = 0, "beta" = 1.5,
#                                   "sigma" = sqrt(1e-2)),
#                    warm_up = -1, model = "OU")
# pred_OU <- ARH_pred_OU(OU, thre_p = 0.995, fpc = TRUE)
ARH_pred_OU <- function(OU, thre_p = 0.95, fpc = TRUE,
                        centered = FALSE) {
  
  # Sample size and grid points where ARH(1) paths take values
  t <- OU$fstoch_proc[["argvals"]] # Interval [a, b] (L2[a, b] as the Hilbert)
  n <- dim(OU$fstoch_proc[["data"]])[1] # Number of trajectories
  
  # Functional mean
  f_mean <- func_mean(OU$fstoch_proc)
  
  # If fpc = TRUE, the OU is firstly decomposed into the first FPC
  if (fpc) {
    
    fpc_OU <- fpc(OU$fstoch_proc, n_fpc = n) # OU$f_stoch_proc is internally centered
    s <- cumsum(fpc_OU[["d"]]^2) # Cumulative empirical explained variance
    p_thre <- 1:which(s/s[length(s)] > thre_p)[1] # The first FPC > thre_p
    
    # Rebuilding the OU process (as fdata). Note that OU_fpc is centered
    OU_fpc <- fpc_to_fdata(X_fpc = fpc_OU, coefs = fpc_OU$scores[, p_thre])
    
    # Recomputing the whole OU (as sde) as a fdata object with n = 1
    OU$stoch_proc <-
      fdata(as.vector(t((OU_fpc + f_mean)$data)),
            argvals = seq(t[1], t[length(t)] * n, l = length(t) * n))
    
  } else {
    
    OU_fpc <- OU$f_stoch_proc - f_mean
    
  }
  
  # Maximum Likelihood Estimator (MLE) of theta
  theta_est <- MLE_theta(OU$stoch_proc)
  
  # Computing the plug-in predictor \widehat{X}_n for each n, according to the
  # ARH1 framework proposed in Alvarez-Liebana et al. (2016)
  predictor_OU <- OU_fpc # Already centered!
  
  # Autocorrelation operator given by rho(X)(t) = exp(-theta_hat * t) * X(b)
  for (nn in 2:n) {
    
    predictor_OU[["data"]][nn, ] <- exp(-theta_est * t) *
      predictor_OU[["data"]][nn - 1, length(t)]
    
  }
  
  # Since f_mean was removed (predictor_OU = OU_fpc), we add it
  if (!centered) {
    
    predictor_OU <- predictor_OU + f_mean
    
  }
  
  # Functional residuals (just 2:n since eps_1 "does not exist")
  residuals <- OU$fstoch_proc[2:n, ] - predictor_OU[2:n, ]
  
  # Output list
  return(list("theta_est" = theta_est, "OU" = OU, "predictor_OU" = predictor_OU,
              "residuals" = residuals))
  
}

# ##################################
# TWO-STAGES SPECIFICATION TEST
# ##################################

# STAGE 1: GOF TEST FOR ARH(z) processes
#
# Description: function for testing if a stochastic process {X_n}
# is an ARH(z) process or not
#
# H0: ARH(z) process
# H1: unspecified alternative
#
# Example:
# OU <- r_stoch_proc(300, seq(0, 1, l = 81),
#                    par.sde = list("alpha" = 0, "beta" = 0.5, "sigma" = 0.1),
#                    warm_up = -1, model = "OU")
# testing_ARHz <- ARHz_test(OU, z = 1, B = 1000, hyp_simp = FALSE,
#                           hyp_comp = TRUE, est_method = "fpcr_l1s",
#                           thre_p = 0.995, thre_q = 0.995)
ARHz_test <- function(X, z = 1, B = 1000, hyp_simp = FALSE, hyp_comp = TRUE,
                      est_method = "fpcr", thre_p = 0.995, thre_q = 0.995,
                      lambda = NULL, boot_scores = TRUE, plot_dens = FALSE,
                      plot_proc = FALSE, save_fit_flm = TRUE,
                      save_boot_stats = TRUE, cv_1se = TRUE) {
  
  # If hyp_comp = hyp_simp = 0 ==> hyp_comp = 1
  hyp_comp <- ifelse(hyp_comp + hyp_simp == 0, 1, hyp_comp)
  
  # Convert X into an ARH(p)  process
  cat("\nConverting X into an ARH(z) process...\n")
  X_flmfr <- ARH_to_FLMFR(X[["fstoch_proc"]], z)
  
  testing_simple <- testing_comp <- NULL
  if (hyp_simp) {
    
    cat("Testing if X is an ARH(0) (rho = null operator)...\n")
    # Testing if X is an ARH(0); that is, rho = null operator (beta0 = 0)
    testing_simple <- flm_test(X_flmfr$X_fdata, X_flmfr$Y_fdata, beta0 = 0,
                               B = B, est_method = est_method,
                               thre_p = thre_p,
                               thre_q = thre_q, lambda = lambda,
                               boot_scores = boot_scores,
                               plot_dens = plot_dens, plot_proc = plot_proc,
                               save_fit_flm = save_fit_flm,
                               save_boot_stats = save_boot_stats,
                               cv_1se = cv_1se, verbose = FALSE)
    
  }
  
  if (hyp_comp) {
    
    cat("Testing if X is an ARH(z) process\n")
    # Testing if X is an ARH(z)
    testing_comp <- flm_test(X_flmfr$X_fdata, X_flmfr$Y_fdata, B = B,
                             est_method = est_method,
                             thre_p = thre_p,
                             thre_q = thre_q, lambda = lambda,
                             boot_scores = boot_scores, plot_dens = plot_dens,
                             plot_proc = plot_proc,
                             save_fit_flm = save_fit_flm,
                             save_boot_stats = save_boot_stats,
                             cv_1se = cv_1se, verbose = FALSE)
    
  }
  
  # Output
  return(list("X" = X, "X_flmfr" = X_flmfr, "testing_simple" = testing_simple,
              "testing_comp" = testing_comp))
}

# STAGE 2: F-TEST FOR OU PROCESSES
#
# Description: given a stochastic process {X_n}_{n in Z},
# this function, firstly, builds the functional model Y = rho(X), 
# with Y = X_n and X = X_{n-1}, and secondly, computes the F-statistic
# between FLMFR (unrestricted) vs OU (that is, a particular ARH(1),
# with a specific operator rho, understood as a particular case of FLMFR).
# 
# Example:
# H0: OU process (as particular ARH(1) process)
# H1: any FLMFR different to H0
#
# Example:
# OU <- r_stoch_proc(300, seq(0, 1, l = 81),
#                    par.sde = list("alpha" = 0, "beta" = 0.5, "sigma" = 0.1),
#                    warm_up = -1, model = "OU")
# testing_ARHz <- ARHz_test(OU, z = 1, B = 1000, hyp_simp = FALSE,
#                           hyp_comp = TRUE, est_method = "fpcr_l1s",
#                           thre_p = 0.995, thre_q = 0.995)
# testing_F_OU <- F_stat_OU(testing_ARHz$X_flmfr, OU, est_method = "fpcr_l1s",
#                           thre_p = 0.995, thre_q = 0.995,
#                           cv_1se = FALSE)
# F_stat <- testing_F_OU$F_stat
F_stat_OU <- function(X_flmfr, X, est_method = "fpcr_l1s", thre_p = 0.995,
                      thre_q = 0.995, lambda = NULL, cv_1se = FALSE,
                      tol = 1e-2, verbose = TRUE) {
  
  ## RESTRICTED OU MODEL: X as an ARH(1) (as a particular case of FLMFR)
  ## between X_n vs X_{n-1}
  
  # Computing predictor as an OU process, decomposed into FPC
  if (verbose) {
    
    cat("\nComputing predictor as an OU process...\n")
    
  }
  pred_OU <- ARH_pred_OU(X, thre_p = thre_p, fpc = TRUE)
  theta_est <- pred_OU[["theta_est"]]
  
  # Residual Sum of Squared L2-norms
  RSS_OU <- mean(apply(pred_OU[["residuals"]][["data"]]^2, 1,
                       FUN = "integral1D",
                       pred_OU[["residuals"]][["argvals"]]))
  
  ## UNRESTRICTED MODEL: X as a FLMFR X_n vs X_{n-1}
  if (verbose) {
    
    cat("Computing FLMFR estimator...\n")
    
  }
  FLMFR_est <- flm_est(X = X_flmfr$X_fdata, Y = X_flmfr$Y_fdata,
                       est_method = est_method, thre_p = thre_p,
                       thre_q = thre_q, lambda = lambda,
                       cv_1se = cv_1se)
  
  # Residual Sum of Squared L2-norms
  RSS_FLMFR <- mean(apply(FLMFR_est[["residuals"]][["data"]]^2, 1,
                          FUN = "integral1D",
                          X_flmfr$Y_fdata[["argvals"]]))
  
  # Computing F-statistic
  F_stat <- (RSS_OU - RSS_FLMFR) / RSS_FLMFR
  
  # Output
  return(list("F_stat" = F_stat, "RSS_FLMFR" = RSS_FLMFR,
              "RSS_OU" = RSS_OU, "pred_OU" = pred_OU,
              "FLMFR_est" = FLMFR_est))
}

# Description: function for implementing a two-stages
#              specification test for OU processes
# STAGE 1:
#     H0: X_n is an ARH(1) process (X_n vs X_{n-1} FLMFR)
#     H1: X_n is not an ARH(1) process
# STAGE 2 (F-test): 
#     H0: X_n = rho(X_{n-1}) from an OU process
#     H1: X_n = rho(X_{n-1}) an alternative FLMFR model
#
# Example: testing Ait-Sahalia (AS) model:
# drf <- expression(0.00107/x - 0.0517 + 0.877*x - 4.604*x^2)  # SDE drift function
# sig <- expression(0.8 * x^1.5)                               # SDE diffusion function
# AS_test <- test_gof_OU(150, t = seq(0, 1, l = 101), warm_up = -1,
#                        type = "other", z = 1, B = 50, X0 = 0.08,
#                        est_method = "fpcr_l1s", thre_p = 0.995, thre_q = 0.995,
#                        verbose = TRUE, drf = drf, sig = sig)
# AS_test$p.value  # p-values Stage 1 and 2
#
# Example: testing OU model
# OU_test <-
#   test_gof_OU(150, t = seq(0, 1, l = 101), warm_up = -1,
#               par.sde = list("alpha" = 0, "beta" = 0.5, "sigma" = 0.05),
#               type = "CKLS", z = 1, B = 50, est_method = "fpcr_l1s",
#               thre_p = 0.995, thre_q = 0.995, verbose = TRUE)
# OU_test$p.value  # p-values Stage 1 and 2
test_gof_OU <- function(n, t = seq(0, 1, l = 71),
                        par.sde = list("alpha" = 0, "beta" = 0.5, "sigma" = 1),
                        mu = par.sde[["alpha"]] / par.sde[["beta"]], X0 = mu,
                        warm_up = 50, type = "CKLS", model = "OU", z = 1, B = 1000,
                        hyp_simp = FALSE, hyp_comp = TRUE, est_method = "fpcr_l1s", 
                        thre_p = 0.995, thre_q = 0.995, lambda = NULL, 
                        boot_scores = TRUE, plot_dens = FALSE, plot_proc = FALSE, 
                        save_fit_flm = TRUE, save_boot_stats = TRUE, cv_1se = TRUE, 
                        verbose = FALSE, drf = NULL, sig = NULL) {
  
  # Generating sde as a functional stochastic process
  X <- r_stoch_proc(n, t = t, par.sde = par.sde, mu = mu, X0 = X0,
                    warm_up = warm_up, type = type, model = model,
                    verbose = FALSE, plot = FALSE, drf = drf, sig = sig)
  
  # STAGE 1: Is X_n an ARH(1) process (that is, X_n vs X_{n-1} FLMFR)?
  cat("\nSTAGE 1\n")
  testing_ARHz <- ARHz_test(X, z = z, B = B, hyp_simp = hyp_simp,
                            hyp_comp = hyp_comp, est_method = est_method,
                            thre_p = thre_p, thre_q = thre_q, lambda = lambda,
                            boot_scores = boot_scores, plot_dens = plot_dens,
                            plot_proc = plot_proc, save_fit_flm = save_fit_flm,
                            save_boot_stats = save_boot_stats, cv_1se = cv_1se)
  
  # STAGE 2: F-test OU vs ARH(1) (unrestricted)
  cat("\nSTAGE 2\n")
  testing_F_OU <- F_stat_OU(testing_ARHz$X_flmfr, X, est_method = est_method,
                            thre_p = thre_p, thre_q = thre_q, lambda = lambda,
                            cv_1se = FALSE, verbose = FALSE)
  F_stat <- testing_F_OU$F_stat
  F_stat_ast <- theta_est_ast <- rep(0, B)
  
  if (verbose) {
    
    pb <- txtProgressBar(max = B, style = 3)
    
  }
  
  beta_est  <- testing_F_OU$pred_OU$theta_est
  sigma_est <- abs(optimize(mce, r = as.vector(X$stoch_proc$data),
                            Delta = 1/(length(t) - 1), interval = c(0,2))$minimum)
  
  par.sde_ast       <- par.sde
  par.sde_ast$beta  <- beta_est
  par.sde_ast$alpha <- 0          # centered process, mu = 0
  par.sde_ast$sigma <- sigma_est
  
  for (b in 1:B) {
    
    # Simulating bootstrap replicates of X, with theta_est as the true parameter
    X_ast <- r_stoch_proc(n, t = t, par.sde = par.sde_ast, mu = 0, X0 = 0,
                          warm_up = warm_up, type = "CKLS", model = "OU",
                          verbose = FALSE, plot = FALSE)
    
    # Converting X_ast into a FLMFR X_n vs X_{n-1}
    X_flmfr_ast <- ARH_to_FLMFR(X_ast[["fstoch_proc"]], z)
    
    # Computing bootstrap F-statistics
    testing_F_OU_ast <- F_stat_OU(X_flmfr_ast, X_ast, est_method = est_method,
                                  thre_p = thre_p, thre_q = thre_q,
                                  lambda = lambda, cv_1se = FALSE,
                                  verbose = FALSE)
    theta_est_ast[b] <- testing_F_OU_ast$pred_OU$theta_est
    F_stat_ast[b] <- testing_F_OU_ast$F_stat
    
    # Display progress
    if (verbose) {
      
      setTxtProgressBar(pb, b)
      
    }
  }
  
  # Approximation of the p-value by MC
  p_value_F_test <- mean(F_stat < F_stat_ast)
  
  # Return htest object
  meth <- "Two-stages procedure for testing if X is as an OU"
  data_name <- paste(deparse(substitute(Y)), deparse(substitute(X)),
                     sep = " ~ ")
  result <-
    structure(list(statistic =
                     c("ARH(0)-stat" = testing_ARHz$testing_simple$statistic,
                       "ARH(1)-stat" = testing_ARHz$testing_comp$statistic,
                       "F-stat OU" = F_stat),
                   p.value = c("GOF ARH(0)" = testing_ARHz$testing_simple$p.value,
                               "GOF ARH(1)" = testing_ARHz$testing_comp$p.value,
                               "F-test OU" = p_value_F_test),
                   boot_statistics =
                     c("GOF ARH(0)" = testing_ARHz$testing_simple$boot_statistics,
                       "GOF ARH(1)" = testing_ARHz$testing_comp$boot_statistics,
                       "F-test" = F_stat_ast),
                   method = meth,
                   parameter = testing_ARHz$testing_comp$parameter,
                   p = testing_ARHz$testing_comp$p,
                   q = testing_ARHz$testing_comp$q,
                   fit_flm = testing_ARHz$testing_comp$fit_flm,
                   theta_est = testing_F_OU$pred_OU$theta_est,
                   boot_theta_est = theta_est_ast, data.name = data_name))
  class(result) <- "htest"
  return(result)
  
}


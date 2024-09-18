<!--
Goodness-of-fit test for autoregressive functional processes and specification test for OU process from a functional perspective. Software companion for "A goodness-of-fit test for functional time series with applications to diffusion processes"
Authors: Alejandra López-Pérez and Javier Álvarez-Liébana (@DadosDeLaplace)
-->

Goodness-of-fit test for autoregressive functional processes and specification test for diffusion processes
======

[![License:
GPLv3](https://img.shields.io/badge/license-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

<!-- <img src="" alt="goffda  hexlogo" align="right" width="200" style="padding: 0 15px; float: right;"/> -->

Overview
-------

Software companion for the article [A goodness-of-fit test for functional time series with applications to diffusion processes (Álvarez-Liébana, López-Perez, González-Manteiga and Febrero Bande, 2021)](https://arxiv.org/), currently submitted. It implements the proposed estimators and **goodness-of-fit tests for 
autoregressive functional processes of order one in Hilbert spaces (ARH(1) processes)**, as well as a **specification test for diffusion processes**, focused on testing [Ornstein-Uhlenbeck processes](https://www.sciencedirect.com/science/article/abs/pii/S016771521630044X).

It also allows to replicate simulations and the data application presented.



```r
# Install and load packages
if(!require(sde)) install.packages("sde", repos = repos)
if(!require(ggplot2)) install.packages("ggplot2", repos = repos)
if(!require(ggthemes)) install.packages("ggthemes", repos = repos)
if(!require(latex2exp)) install.packages("latex2exp", repos = repos)
if(!require(goffda)) install.packages("goffda", repos = repos)
if(!require(fda.usc)) install.packages("fda.usc", repos = repos)
if(!require(tidyverse)) install.packages("tidyverse", repos = repos)
if(!require(lubridate)) install.packages("lubridate", repos = repos)
if(!require(glue)) install.packages("glue", repos = repos)
if(!require(reshape2)) install.packages("reshape2", repos = repos)
```

![](https://github.com/dadosdelaplace/gof-test-arh-ou-process/blob/main/plots/fig2a_plotOU.png)
*Continuous path of an OU process splitted into subintervals [0,h] (h=1), such that `X_i(t) = \xi_{ih + t}`, with i=1, ..., 5*

![](https://github.com/dadosdelaplace/gof-test-arh-ou-process/blob/main/plots/fig2b_plotOU.png)
*OU process characterized as a set of ARH(1) trajectories `{X_n(t): t in [0,1]}`, with i=1, ...,5*

<details><summary>Code</summary>

```r
 # ####################################
# Description: plotting trajectories of OU process as an ARH(1) process
# Authors: A. López-Pérez and J. Álvarez-Liébana
# Article: «A goodness-of-fit test for functional time series: applications
#          to specification test for stochastic diffusion models» (submitted)
# ####################################

# Packages and working space
rm(list = ls())
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
repos <- "http://cran.us.r-project.org"
if(!require(sde)) install.packages("sde", repos = repos)
if(!require(ggplot2)) install.packages("ggplot2", repos = repos)
if(!require(latex2exp)) install.packages("latex2exp", repos = repos)
if(!require(tidyverse)) install.packages("tidyverse", repos = repos)
if(!require(ggthemes)) install.packages("ggthemes", repos = repos)

# Fixing the seed for the randomness
set.seed(123)

# We simulate paths of continuous-time zero-mean stochastic
# OU process as CKLS process with kappa = 2, sigma = 0.001,
# with 250 grid points, trajectories splitted into subintervals
# [0, h], with h = 50.
x <- sde.sim(X0 = 0, theta = c(0, 2, sqrt(0.001)), N = 250, T = 5,
             t0 = 0, delta = 1/50, model = "OU")
t <- seq(0, 5, l = length(x))

# We interpolate the trajectories into 10 001 points
interp <- spline(t, x, n = 1e4 + 1)
OU <- as.tibble(data.frame("t" = interp$x, "x" = interp$y,
                           "interv" = as.factor(pmin(floor(interp$x), 4))))

# Plot the trajectories
theme_set(theme_bw(base_family = "Poppins"))
theme_update(text = element_text(family = "Poppins", size = 15,
                                 color = "black"))
             
fig2a <- 
  ggplot(OU, aes(x = t, y = x, color = interv)) +
  geom_line(size = 1.3) + # width of line
  # Vertical lines
  geom_vline(aes(xintercept = 1), linetype = "dashed", size = 1.2) +
  geom_vline(aes(xintercept = 2), linetype = "dashed", size = 1.2) +
  geom_vline(aes(xintercept = 3), linetype = "dashed", size = 1.2) +
  geom_vline(aes(xintercept = 4), linetype = "dashed", size = 1.2) +
  # Color pattern from tableau appareance according to intervals
  scale_color_tableau() +
  # Annotations: curve + text
  annotate(geom = "curve", x = 0.33, y = 0.02, xend = 0.8, yend = 0.01, 
           color = "firebrick", curvature = -0.5, size = 1.1,
           arrow = arrow(length = unit(3.5, "mm"))) +
  annotate(geom = "text", x = 0.27, y = 0.02,
           label = TeX("$X_{1}(t)$", bold = TRUE),
           hjust = "right", family = "Poppins", size = 6) +
  annotate(geom = "curve", x = 1.55, y = 0.026, xend = 1.4, yend = 0.017, 
           color = "firebrick", curvature = 0.4, size = 1.1,
           arrow = arrow(length = unit(3.5, "mm"))) +
  annotate(geom = "text", x = 1.47, y = 0.029,
           label = TeX("$X_{2}(t)$", bold = TRUE),
           hjust = "left", family = "Poppins", size = 6) +
  annotate(geom = "curve", x = 2.65, y = 0.01, xend = 2.83, yend = 0, 
           color = "firebrick", curvature = -0.7, size = 1.1,
           arrow = arrow(length = unit(3.5, "mm"))) +
  annotate(geom = "text", x = 2.6, y = 0.01,
           label = TeX("$X_{3}(t)$", bold = TRUE),
           hjust = "right", family = "Poppins", size = 6) +
  annotate(geom = "curve", x = 3.25, y = 0.025, xend = 3.45, yend = 0.01, 
           color = "firebrick", curvature = 0.3, size = 1.1,
           arrow = arrow(length = unit(3.5, "mm"))) +
  annotate(geom = "text", x = 3.2, y = 0.027,
           label = TeX("$X_{4}(t)$", bold = TRUE),
           hjust = "left", family = "Poppins", size = 6)  +
  annotate(geom = "curve", x = 4.7, y = 0.025, xend = 4.55, yend = 0.005, 
           color = "firebrick", curvature = 0.5, size = 1.1,
           arrow = arrow(length = unit(3.5, "mm"))) +
  annotate(geom = "text", x = 5, y = 0.028,
           label = TeX("$X_{5}(t)$", bold = TRUE),
           hjust = "right", family = "Poppins", size = 6) +
  # Axis
  labs(x = "t",  y = "") +
  # Title (splitted in two lines)
  ggtitle("Splitting trajectories\nof an OU process") +
  # Font settings
  theme(axis.title.x = element_text(family = "Poppins", hjust = .5,
                                    size = 13, face = "bold"),
        axis.title.y = element_text(family = "Poppins", hjust = .5,
                                    size = 13, face = "bold"),
        plot.title = element_text(face = "bold", family = "Poppins",
                                  size = 13),
        legend.position = "none") # without legend

# Save as 8x4 png 
ggsave("fig2a_plotOU.png", width = 8, height = 5)

# OU plotted as ARH(1) process
OU_ARH1 <- data.frame("t" = rep(seq(0, 1, l = 2000), 5),
                      "x" = rev(rev(OU$x)[-1]),
                      "interv" = rev(rev(OU$interv)[-1]))
fig2b <- 
  ggplot(OU_ARH1, aes(x = t, y = x, color = interv)) +
  geom_line(size = 1.3) + # width of line
  # Color pattern from tableau appareance according to intervals
  scale_color_tableau() +
  annotate(geom = "text", x = 0.34, y = 0.02,
           label = TeX("$X_{1}(t)$", bold = TRUE),
           hjust = "right", family = "Poppins", size = 6) +
  annotate(geom = "text", x = 0.94, y = 0.03,
           label = TeX("$X_{2}(t)$", bold = TRUE),
           hjust = "right", family = "Poppins", size = 6) +
  annotate(geom = "text", x = 0.83, y = -0.03,
           label = TeX("$X_{3}(t)$", bold = TRUE),
           hjust = "left", family = "Poppins", size = 6) +
  annotate(geom = "text", x = 0.515, y = 0.0155,
           label = TeX("$X_{4}(t)$", bold = TRUE),
           hjust = "right", family = "Poppins", size = 6) +
  annotate(geom = "text", x = 0.045, y = -0.011,
           label = TeX("$X_{5}(t)$", bold = TRUE),
           hjust = "right", family = "Poppins", size = 6) +
  # Axis
  labs(x = "t", y = "") +
  # Title (splitted in two lines)
  ggtitle("OU stochastic model characterized\nas an ARH(1) process") +
  # Font settings
  theme(axis.title.x = element_text(family = "Poppins", hjust = .5,
                                    size = 13, face = "bold"),
        axis.title.y = element_text(family = "Poppins", hjust = .5,
                                    size = 13, face = "bold"),
        plot.title = element_text(face = "bold", family = "Poppins",
                                  size = 13),
        legend.position = "none") # without legend

# Save as 8x4 png 
ggsave("fig2b_plotOU.png", width = 8, height = 5)




```
</details>

Installation
------------

Get the released version from [GitHub (Download ZIP button)](https://github.com/dadosdelaplace/gof-test-arh-ou-process/archive/refs/heads/main.zip), including [R codes](https://github.com/dadosdelaplace/gof-test-arh-ou-process/tree/main/R), [plots](https://github.com/dadosdelaplace/gof-test-arh-ou-process/tree/main/plots) and [data](https://github.com/dadosdelaplace/gof-test-arh-ou-process/tree/main/data). The installation of the following packages will be required:

* `sde`: package to simulate Stochastic Differential Equations (SDEs).
* `ggplot2`, `latex2exp`, `ggthemes`: packages related with visualizing data.
* `tidyverse`: a collection of R packages designed for data science.
* [`goffda`](https://github.com/egarpor/goffda): package for performing a goodness-of-fit test for the functional linear model with functional response
* `fda.usc`: package for managing functional data objects in R.


Auxiliary functions
------------

#### `r_stoch_proc()`

* Description: function to simulate stochastic processes from a functional perspective, splitted the whole trajectories into subintervals
* Main inputs:
  * `n`: sample size (numer of functional trajectories valued in `t`).
  * `t`: grid points where functional trajectories are valued.
  * `par.sde`: a list of parameters to introduce if the model is generated nested in the CKLS parametric family (`type = CKLS`). Otherwise, wnevaluated expressions for drift (`drf`) and sigma (`sig`) functions can be also provided.
  * `warm_up`: number of trajectories to be discarded.
  * `model`: type of CKLS model to be simulate.
* Outputs: a list with the following elements
  * `fstoch_proc`: stochastic process as a `fdata` object splitted `n` trajectories.
  * `stoch_proc`: stochastic process as a `fdata` object but with a single trajectories.

```r
# Example
OU <- r_stoch_proc(500, seq(0, 1, l = 101), par.sde = list("alpha" = 0, "beta" = 0.5, "sigma" = 0.1),
                   warm_up = -1, model = "OU")
```

<details><summary>Code</summary>

```r
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

```
</details>

#### `mce()`

* Description: function to be optimized for estimating the volatility (sigma) of OU process, given an OU process as a vectorial one-dimensional stochastic process
* Main inputs:
  * `sigma`: grid of values to be chosen.
  * `r`: OU process as a vectorial one-dimensional stochastic process
  * `Delta`: discretization step.
* Output: the root-n consistent estimator of sigma, based on the proposal by [Corradi and White, 1999](https://onlinelibrary.wiley.com/doi/abs/10.1111/1467-9892.00136)


```r
# Example:
OU <- r_stoch_proc(500, seq(0, 1, l = 101), par.sde = list("alpha" = 0, "beta" = 0.5, "sigma" = 0.1),
                   warm_up = -1, model = "OU")
optimize(mce, r = as.vector(OU$stoch_proc$data), Delta = 1/(length(OU$fstoch_proc$argvals) - 1),
         interval = c(0, 2))$minimum
```

<details><summary>Code</summary>

```r
mce <- function(sigma, r, Delta){
  y    <- diff(r)
  n    <- length(y)
  suma <- sum(log(sigma^2) + (y^2) / (Delta * sigma^2))
  return(suma / n)
}
```
</details>
  

#### `MLE_theta()`

* Description: function for computing the MLE estimator of theta of an OU process
* Main inputs:
  * `sde_OU`: OU process given as a fdata (fda.usc package) object with n = 1
* Output: the MLE strongly-consistent estimator of theta, based on the proposal by [Álvarez-Liébana et al., 2016](https://www.sciencedirect.com/science/article/pii/S016771521630044X)


```r
# Example:
OU <- r_stoch_proc(500, seq(0, 1, l = 101), par.sde = list("alpha" = 0, "beta" = 0.5, "sigma" = 1),
                   warm_up = -1, model = "OU")
MLE_theta(OU$stoch_proc)

```

<details><summary>Code</summary>

```r
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
               
```
               
</details>
  

##### `ARH_to_FLMFR()`
  
* Description: function for converting functional stochastic process into a FLMFR
* Main inputs:
  * `ARHz`: ARH process as a `fdata` object with `n` trajectories.
  * `z`: order of the ARH process (ARH(z) process).
  * `centered`: flag to indicate if ARH process is centered; default to `FALSE`.
* Output: `X_fdata` and `Y_fdata`, as functional covariates and responses, according to the characterization proposed by [López-Perez et al., 2021](arxiv.org/...)
  
```r
# Example:
OU <- r_stoch_proc(500, seq(0, 1, l = 101),  par.sde = list("alpha" = 0, "beta" = 0.5, "sigma" = 0.1),
                   warm_up = -1, model = "OU")
OU_FLMFR <- ARH_to_FLMFR(OU$fstoch_proc, 1)
plot(OU_FLMFR$X_fdata) # Plotting functional X variable
plot(OU_FLMFR$Y_fdata) # Plotting functional Y variable
```

<details><summary>Code</summary>

```r
  
ARH_to_FLMFR <- function(ARHz, z, centered = FALSE) {
  
  # Common settings: sample size and number of grids where it is evaluated
  n <- dim(ARHz)[1]
  n_grids <- length(ARHz[["argvals"]])
  
  # Data should be centered to be transformed
  if (!centered) {
    
    ARHz <- ARHz - func_mean(ARHz)
    
  }
  
  # We will remove the first p trajectories since Y_n = X_{n + 1}
  Y <- ARHz[(z + 1):n, ]
  
  if (z == 1) {
    
    X <- ARHz[1:(n - 1), ]
    
  } else {
    
    X <- fda.usc::fdata(matrix(0, n - z, n_grids), argvals = ARHz[["argvals"]])
    for (j in 1:z) {
      
      s <- (1:n_grids)[ARHz[["argvals"]] >= (j - 1)/z & ARHz[["argvals"]] <= j/z]
      X[["data"]] <- X[["data"]] + ARHz[["data"]][1:(n - j), s * z - (j - 1)]
      
    }
  }
  
  # Output
  return(list("X_fdata" = X, "Y_fdata" = Y))
}
```
</details>
  
##### `ARH_pred_OU()`
  
* Description: function for predicting an OU as an ARH(1) process
* Main inputs:
  * `OU`: OU process as the list provided by `r_stoch_proc()`.
  * `thre_p`: order of the ARH process (ARH(z) process).
  * `fpc`: flag to indicate if the functional process should be firstly decomposed into a truncated set of Functional Principal Components (FPC).
  * `centered = FALSE`: flag to indicate if ARH process is centered; default to `FALSE`.
* Output: a list with the following elements
  * `theta_est`: MLE estimator of theta, performed by function `MLE_theta()`.
  * `OU`: input OU process.
  * `predictor_OU`: OU plug-in predictor according to [Álvarez-Liébana et al., 2016](https://www.sciencedirect.com/science/article/pii/S016771521630044X)
  * `residuals`: functional errors.

  
```r
# Example
OU <- r_stoch_proc(1500, seq(0, 1, l = 101), par.sde = list("alpha" = 0, "beta" = 1.5, "sigma" = sqrt(1e-2)),
                   warm_up = -1, model = "OU")
pred_OU <- ARH_pred_OU(OU, thre_p = 0.995, fpc = TRUE)
```
  
<details><summary>Code</summary>

```r
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
```
</details>



Main functions: two-stages specification test
------------
 

##### `ARHz_test()` (stage 1: GoF test for ARH(z) processes)
  

* Description: function for testing if a stochastic process is an ARH(z) process or not
  * H0: ARH(z) process
  * H1: unspecified alternative
* Main inputs:
  * `X`: stochastic process as the list provided by `r_stoch_proc()`.
  * `z`: order of the ARH process (ARH(z) process).
  * `B`: bootstrap replicates.
  * `hyp_simp`/ `hyp_comp`: flags to indicate if simple and/or composite hypothesis are tested.
  * `est_method`: FLMFR estimation method.
  * `thre_p`, `thre_q`: amount of variance to be captured by the initial FPC selection.
* Output: a list with the following elements
  * `X`: input stochastic process.
  * `X_flmfr`: X splitted and characterized as FLMFR
  * `testing_simple`/ `testing_comp`: test outputs for simple and/or composite hypothesis test.
 
  
```r
# Example
OU <- r_stoch_proc(300, seq(0, 1, l = 81), par.sde = list("alpha" = 0, "beta" = 0.5, "sigma" = 0.1),
      warm_up = -1, model = "OU")
testing_ARHz <- ARHz_test(OU, z = 1, B = 1000, hyp_simp = FALSE, hyp_comp = TRUE, est_method = "fpcr_l1s",
                          thre_p = 0.995, thre_q = 0.995)
```
  
<details><summary>Code</summary>

```r
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
```
</details>
 
#### `F_stat_OU()` (stage 2: F-test for OU processes):
 
* Description: given a stochastic process, this function, firstly, builds the functional model Y = rho(X), with Y = X_n and X = X_{n-1}, and secondly, computes the F-statistic between FLMFR (unrestricted) vs OU (that is, a particular ARH(1), with a specific operator rho, understood as a particular case of FLMFR).
  * H0: OU process (as particular ARH(1) process)
  * H1: any FLMFR different to H0
* Main inputs:
  * `X_flmfr`: stochastic process characterized as a FLMFR
  * `X`: stochastic process as the list provided by `r_stoch_proc()`.
  * `est_method`: FLMFR estimation method.
  * `thre_p`, `thre_q`: amount of variance to be captured by the initial FPC selection.
* Output: a list with the following elements
  * `F_stat`: F-statistic of the test
  * `RSS_FLMFR`: residual sum of squared norms under the unrestricted FLMFR hypothesis
  * `RSS_OU`: residual sum of squared norms under the OU model
  * `pred_OU`: prediction of the stochastic process as an OU process (via ARH(1) model)
  * `FLMFR_est`: prediction of the stochastic process as a FLMFR process
 

```r
# Example
OU <- r_stoch_proc(300, seq(0, 1, l = 81),
                   par.sde = list("alpha" = 0, "beta" = 0.5, "sigma" = 0.1),
                   warm_up = -1, model = "OU")
testing_ARHz <- ARHz_test(OU, z = 1, B = 1000, hyp_simp = FALSE,
                          hyp_comp = TRUE, est_method = "fpcr_l1s",
                          thre_p = 0.995, thre_q = 0.995)
testing_F_OU <- F_stat_OU(testing_ARHz$X_flmfr, OU, est_method = "fpcr_l1s",
                          thre_p = 0.995, thre_q = 0.995,
                          cv_1se = FALSE)
F_stat <- testing_F_OU$F_stat
```
 
 
 <details><summary>Code</summary>

```r
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
```
</details>
  
  
#### `test_gof_OU()`:
  
* Description: function for implementing a two-stages OU specification test
  * STAGE 1:
    * H0: X_n is an ARH(1) process (X_n vs X_{n-1} FLMFR)
    * H1: X_n is not an ARH(1) process
  * STAGE 2 (F-test): 
    * H0: X_n = rho(X_{n-1}) from an OU process
    * H1: X_n = rho(X_{n-1}) an alternative FLMFR model
* Main inputs:
  * `n`: sample size (numer of functional trajectories valued in `t`).
  * `t`: grid points where functional trajectories are valued.
  * `par.sde`: a list of parameters to introduce if the model is generated nested in the CKLS parametric family (`type = CKLS`). Otherwise, unevaluated expressions for drift (`drf`) and sigma (`sig`) functions can be also provided.
  * `warm_up`: number of trajectories to be discarded.
  * `model`: type of CKLS model to be simulated.
  * `hyp_simp`/ `hyp_comp`: flags to indicate if simple and/or composite hypothesis are tested.
  * `est_method`: FLMFR estimation method.
  * `thre_p`, `thre_q`: amount of variance to be captured by the initial FPC selection.
  * `z`: order to be test (ARH(z) model)
* Outputs: an htest object with (among others elements)
  * `statistics` (`ARH(0)-stat`, `ARH(1)-stat`, `F-stat OU`): statistics for the simple/composite hypothesis of the stage 1, and F-statistic of the stage 2.
  * `p.value` (`GOF ARH(0)`, `GOF ARH(1)`, `F-test OU`): p-values of the corresponding tests
  * `boot_statistics` (`GOF ARH(0)`, `GOF ARH(1)`, `F-test`): bootstrapped statistics.
  class(result) <- "htest"


```r
# Example: testing Ait-Sahalia (AS) model:
drf <- expression(0.00107/x - 0.0517 + 0.877*x - 4.604*x^2)  # SDE drift function
sig <- expression(0.8 * x^1.5)                               # SDE diffusion function
AS_test <- test_gof_OU(150, t = seq(0, 1, l = 101), warm_up = -1,  type = "other", z = 1, B = 50,
                       X0 = 0.08, est_method = "fpcr_l1s", thre_p = 0.995, thre_q = 0.995,
                       verbose = TRUE, drf = drf, sig = sig)
AS_test$p.value  # p-values Stage 1 and 2

# Example: testing OU model
OU_test <- test_gof_OU(150, t = seq(0, 1, l = 101), warm_up = -1, par.sde = list("alpha" = 0, "beta" = 0.5, "sigma" = 0.05),
                       type = "CKLS", z = 1, B = 50, est_method = "fpcr_l1s",
                       thre_p = 0.995, thre_q = 0.995, verbose = TRUE)
OU_test$p.value  # p-values Stage 1 and 2
```
  
<details><summary>Code</summary>

```r
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
```
</details>

  

Real-data application
------------

Motivated by the extensive use of diffusion processes, the three datasets considered consist on currency pair rates Euro-British pound (EURGBP), Euro-US Dollar (EURUSD), British pound-US Dollar (GBPUSD), determined in the [foreign exchange market (the largest and most liquid financial market)](www.histdata.com), such that data was recorded every 5 minutes from January 1, 2019 to December 31, 2019.

* [`EURGBP_2019_5min.csv` file](https://github.com/dadosdelaplace/gof-test-arh-ou-process/blob/main/data/EURGBP_2019_5min.csv): Euro-British pound (EURGBP) high-frequency (5-minutes) exchange rates.
* [`EURUSD_2019_5min.csv` file](https://github.com/dadosdelaplace/gof-test-arh-ou-process/blob/main/data/EURUSD_2019_5min.csv): Euro-US Dollar (EURUSD) high-frequency (5-minutes) exchange rates.
* [`GBPUSD_2019_5min.csv` file](https://github.com/dadosdelaplace/gof-test-arh-ou-process/blob/main/data/GBPUSD_2019_5min.csv): British pound-US Dollar (GBPUSD) high-frequency (5-minutes) exchange rates.

Motivated by the extensive use of diffusion processes in finance, commonly used to model currency exchange rates [Ball and Roma, 1994] and intra-day patterns in the foreign exchange market [Andersen and Bollerslev, 1997], we apply our OU specification test, in the context of **high-frequency financial data**. Unlike univariate frameworks, where a vast time window is required for getting properly sample sizes, and thus, certain properties essential for long horizon asymptotics (e.g., ergodicity) would not be achieved, we consider a finite time observation window, where the **FDA scheme allows to capture the dynamic of the process**. It is also noteworthy that, due to the high-frequency, observations are subject to market microstructure noise [Aït-Sahalia et al., 2005], interacting with the sampling frequency, that could lead to reject the null hypothesis.
    	
The daily curves are valued in the interval [0,1], accounts for a 1-day window, discretized in 288 equispaced grid points.
 
<details><summary>Code: loading and preprocessing the data</summary>
 
```r
# ###########################
# PACKAGES
# ###########################
rm(list = ls())
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
repos <- "http://cran.us.r-project.org"
if(!require(ggplot2)) install.packages("ggplot2", repos = repos)
if(!require(fda.usc)) install.packages("fda.usc", repos = repos)
if(!require(tidyverse)) install.packages("tidyverse", repos = repos)
if(!require(lubridate)) install.packages("lubridate", repos = repos)
if(!require(glue)) install.packages("glue", repos = repos)
if(!require(reshape2)) install.packages("reshape2", repos = repos)
if(!require(latex2exp)) install.packages("latex2exp", repos = repos)

# Companion software
source("./OU_estimation_test.R")

# ###########################
# DATA
# ###########################

# URL of data freely available at GitHub repository
# https://github.com/dadosdelaplace/gof-test-arh-ou-process
repo <- paste0("https://raw.githubusercontent.com/dadosdelaplace/",
               "gof-test-arh-ou-process/main/data")
url <- c(glue("{repo}/EURGBP_2019_5min.csv"),
         glue("{repo}/EURUSD_2019_5min.csv"),
         glue("{repo}/GBPUSD_2019_5min.csv"))

# Data was extracted from histdata.com/download-free-forex-data/
EURGBP <- read_delim(file = url[1], delim = ";")
EURUSD <- read_delim(file = url[2], delim = ";")
GBPUSD <- read_delim(file = url[3], delim = ";")

# Preparing dates
dates <- as.Date(as.character(EURGBP$date), format = "%Y%m%d")
hours <- floor(as.numeric(EURGBP$hour) / 1e4)
minutes <- floor((as.numeric(EURGBP$hour) - hours * 1e4) / 1e2)
hourtime <- hms::as_hms(hours * 60 * 60 + minutes*60)

# Build data
data <- data.frame("dates" = dates, "hourtime" = hourtime,
                   "day" = yday(dates), "EURGBP" = EURGBP[, 3],
                   "EURUSD" = EURUSD[, 3], "GBPUSD" = GBPUSD[, 3])
names(data) <- c("dates", "hourtime", "day",
                 "EURGBP", "EURUSD", "GBPUSD")

# Settings of working space: 288 daily curves recorded each 5 minutes
t <- seq(0, 1, by = 1/(60*24/5 - 1))
fda_data <- list()
fda_data$EURGBP <-
  fdata(t(matrix(data$EURGBP, nrow = length(t))), argvals = t)
fda_data$EURUSD <-
  fdata(t(matrix(data$EURUSD, nrow = length(t))), argvals = t)
fda_data$GBPUSD <-
  fdata(t(matrix(data$GBPUSD, nrow = length(t))), argvals = t)

```
</details>
 
 
![](https://github.com/dadosdelaplace/gof-test-arh-ou-process/blob/main/plots/fig2a_real_data.png)
 
![](https://github.com/dadosdelaplace/gof-test-arh-ou-process/blob/main/plots/fig2b_real_data.png) 
 
![](https://github.com/dadosdelaplace/gof-test-arh-ou-process/blob/main/plots/fig2c_real_data.png)
 
*EURGBP, EURUSD and GBPUSD exchange rates: whole trajectories of the observed stochastic process.*

![](https://github.com/dadosdelaplace/gof-test-arh-ou-process/blob/main/plots/fig2d_real_data.png) 
 
![](https://github.com/dadosdelaplace/gof-test-arh-ou-process/blob/main/plots/fig2e_real_data.png) 
 
![](https://github.com/dadosdelaplace/gof-test-arh-ou-process/blob/main/plots/fig2f_real_data.png)
 
*EURGBP, EURUSD and GBPUSD exchange rates: trajectories splitted in daily curves as an ARH(1) process*
 
<details><summary>Code: plotting the data</summary>
 
```r
# ###########################
# PLOTS
# ###########################

# Theme
theme_set(theme_bw(base_family = "Times New Roman"))
theme_update(text = element_text(family = "Times New Roman", size = 15),
             axis.title.x =
               element_text(family = "Times New Roman",
                            hjust = .5, size = 13),
             axis.title.y = element_text(family = "Times New Roman",
                                         hjust = .5, size = 13),
             plot.title =
               element_text(face = "bold", family = "Times New Roman",
                            size = 15),
             legend.position = "none") # without legend

# Plotting the whole trajectories as one-dimensional sde
fig2a <- 
  ggplot(data, aes(x = dates, y = EURGBP, color = day)) +
  geom_line(size = 1.3) + # width of line
  labs(x = "dates",  y = "EURGBP")
# Save as 8x4 png 
ggsave("fig2a_real_data.png", width = 8, height = 5)

fig2b <- 
  ggplot(data, aes(x = dates, y = EURUSD, color = day)) +
  geom_line(size = 1.3) + # width of line
  labs(x = "dates",  y = "EURUSD")
# Save as 8x4 png 
ggsave("fig2b_real_data.png", width = 8, height = 5)

fig2c <- 
  ggplot(data, aes(x = dates, y = GBPUSD, color = day)) +
  geom_line(size = 1.3) + # width of line
  labs(x = "dates",  y = "GBPUSD")
# Save as 8x4 png 
ggsave("fig2c_real_data.png", width = 8, height = 5)

# Plotting as FD
EURGBP_fd <- data.frame("t" = rep(fda_data$EURGBP$argvals,
                                  dim(fda_data$EURGBP)[1]),
                        "x" = as.numeric(t(fda_data$EURGBP$data)),
                        "n" = rep(1:dim(fda_data$EURGBP)[1],
                                  each = dim(fda_data$EURGBP)[2]))
EURUSD_fd <- data.frame("t" = rep(fda_data$EURUSD$argvals,
                                  dim(fda_data$EURUSD)[1]),
                        "x" = as.numeric(t(fda_data$EURUSD$data)),
                        "n" = rep(1:dim(fda_data$EURUSD)[1],
                                  each = dim(fda_data$EURUSD)[2]))
GBPUSD_fd <- data.frame("t" = rep(fda_data$GBPUSD$argvals,
                                  dim(fda_data$GBPUSD)[1]),
                        "x" = as.numeric(t(fda_data$GBPUSD$data)),
                        "n" = rep(1:dim(fda_data$GBPUSD)[1],
                                  each = dim(fda_data$GBPUSD)[2]))

fig2d <- 
  ggplot(EURGBP_fd, aes(x = t, y = x, group = n, color = n)) +
  geom_line(size = 0.7) + # width of line
  labs(x = "t",  y = TeX("$X_{n}(t)$ (EURGBP)"))
# Save as 8x4 png 
ggsave("fig2d_real_data.png", width = 8, height = 5)

fig2e <- 
  ggplot(EURUSD_fd, aes(x = t, y = x, group = n, color = n)) +
  geom_line(size = 0.7) + # width of line
  labs(x = "t",  y = TeX("$X_{n}(t)$ (EURUSD)"))
# Save as 8x4 png 
ggsave("fig2e_real_data.png", width = 8, height = 5)

fig2f <- 
  ggplot(GBPUSD_fd, aes(x = t, y = x, group = n, color = n)) +
  geom_line(size = 0.7) + # width of line
  labs(x = "t",  y = TeX("$X_{n}(t)$ (GBPUSD)"))
# Save as 8x4 png 
ggsave("fig2f_real_data.png", width = 8, height = 5)
```
</details>

The following steps are implemented:
 * **Stage 0**: one-dimensional SDE trajectories are converted into ARH(1) processes.
 * **Stage 1**: testing if the set of trajectories follows an ARH(1) model, calibrating the test by a wild boostrap.
 * **Stage 2**: testing if the parametric FLMFR model is, in particular, an OU model (characterized as an ARH(1) processes), calibrating the test by a parametric bootstrap.
 
<details><summary>Code: OU specification test</summary>
 
```r
 
# ###########################
# TESTING EURGBP DATA
# ###########################

# Parameters
set.seed(123)
z <- 1 # order to test
hyp_simp <- FALSE # we test composite hypothesis
hyp_comp <- TRUE
est_method  <- "fpcr_l1s" # FLMFR estimation method
thre_p <- thre_q <- 0.995 # Initial (p, q) for capturing 99.5% of variance
cv_1se <- TRUE # optimal lambda be the lambda.1se, as returned by cv.glmnet?
B <- 1000 # Bootstrap replicates

# STAGE 1: Is X_n an ARH(1) process (that is, X_n vs X_{n-1} FLMFR)?
X <- list("fstoch_proc" = fda_data$EURGBP)
testing_ARHz <- ARHz_test(X, z = z, B = B, hyp_simp = hyp_simp,
                          hyp_comp = hyp_comp, est_method = est_method,
                          thre_p = thre_p, thre_q = thre_q, cv_1se = cv_1se)
testing_ARHz$testing_comp$p.value

# STAGE 2: F-test OU vs ARH(1) (unrestrictedd)
testing_F_OU <- F_stat_OU(testing_ARHz$X_flmfr, X, est_method = est_method,
                          thre_p = thre_p, thre_q = thre_q, cv_1se = FALSE)
F_stat <- testing_F_OU$F_stat
F_stat_ast <- theta_est_ast <- rep(0, B) 

# Parameter estimates
beta_est <- testing_F_OU$pred_OU$theta_est
sigma_est <- abs(optimize(mce, r = as.vector(data$EURGBP),
                          Delta = 1/(length(t) - 1),
                          interval = c(0, 2))$minimum)
par.sde_ast <- list("alpha" = 0, "beta" = beta_est,
                    "sigma" = sigma_est) # centered process, alpha = mu = 0

# Parametric boostrap
pb <- txtProgressBar(max = B, style = 3)
for (b in 1:B) {
  
  # Simulating bootstrap replicates of X
  X_ast <- r_stoch_proc(dim(fda_data$EURGBP$data)[1], t = t,
                        par.sde = par.sde_ast, mu = 0, X0 = 0,
                        warm_up = -1, type = "CKLS", model = "OU",
                        verbose = FALSE, plot = FALSE)
  
  # Converting X_ast into a FLMFR X_n vs X_{n-1}
  X_flmfr_ast <- ARH_to_FLMFR(X_ast[["fstoch_proc"]], z)
  
  # Computing bootstrap F-statistics
  testing_F_OU_ast <- F_stat_OU(X_flmfr_ast, X_ast, est_method = est_method,
                                thre_p = thre_p, thre_q = thre_q,
                                cv_1se = FALSE, verbose = FALSE)
  theta_est_ast[b] <- testing_F_OU_ast$pred_OU$theta_est
  F_stat_ast[b]    <- testing_F_OU_ast$F_stat
  
  setTxtProgressBar(pb, b)
}

# Approximation of the p-value by MC
p_value_F_test <- mean(F_stat < F_stat_ast)

# p-values
c(testing_ARHz$testing_comp$p.value, p_value_F_test)


# ###########################
# TESTING EURUSD DATA
# ###########################

# Parameters
set.seed(123)
z <- 1 # order to test
hyp_simp <- FALSE # we test composite hypothesis
hyp_comp <- TRUE
est_method  <- "fpcr_l1s" # FLMFR estimation method
thre_p <- thre_q <- 0.995 # Initial (p, q) for capturing 99.5% of variance
cv_1se <- TRUE # optimal lambda be the lambda.1se, as returned by cv.glmnet?
B <- 1000 # Bootstrap replicates

# STAGE 1: Is X_n an ARH(1) process (that is, X_n vs X_{n-1} FLMFR)?
X <- list("fstoch_proc" = fda_data$EURUSD)
testing_ARHz <- ARHz_test(X, z = z, B = B, hyp_simp = hyp_simp,
                          hyp_comp = hyp_comp, est_method = est_method,
                          thre_p = thre_p, thre_q = thre_q, cv_1se = cv_1se)
testing_ARHz$testing_comp$p.value

# STAGE 2: F-test OU vs ARH(1) (unrestrictedd)
testing_F_OU <- F_stat_OU(testing_ARHz$X_flmfr, X, est_method = est_method,
                          thre_p = thre_p, thre_q = thre_q, cv_1se = FALSE)
F_stat <- testing_F_OU$F_stat
F_stat_ast <- theta_est_ast <- rep(0, B) 

# Parameter estimates
beta_est <- testing_F_OU$pred_OU$theta_est
sigma_est <- abs(optimize(mce, r = as.vector(data$EURUSD),
                          Delta = 1/(length(t) - 1),
                          interval = c(0, 2))$minimum)
par.sde_ast <- list("alpha" = 0, "beta" = beta_est,
                    "sigma" = sigma_est) # centered process, alpha = mu = 0

# Parametric boostrap
pb <- txtProgressBar(max = B, style = 3)
for (b in 1:B) {
  
  # Simulating bootstrap replicates of X
  X_ast <- r_stoch_proc(dim(fda_data$EURUSD$data)[1], t = t,
                        par.sde = par.sde_ast, mu = 0, X0 = 0,
                        warm_up = -1, type = "CKLS", model = "OU",
                        verbose = FALSE, plot = FALSE)
  
  # Converting X_ast into a FLMFR X_n vs X_{n-1}
  X_flmfr_ast <- ARH_to_FLMFR(X_ast[["fstoch_proc"]], z)
  
  # Computing bootstrap F-statistics
  testing_F_OU_ast <- F_stat_OU(X_flmfr_ast, X_ast, est_method = est_method,
                                thre_p = thre_p, thre_q = thre_q,
                                cv_1se = FALSE, verbose = FALSE)
  theta_est_ast[b] <- testing_F_OU_ast$pred_OU$theta_est
  F_stat_ast[b]    <- testing_F_OU_ast$F_stat
  
  setTxtProgressBar(pb, b)
}

# Approximation of the p-value by MC
p_value_F_test <- mean(F_stat < F_stat_ast)

# p-values
c(testing_ARHz$testing_comp$p.value, p_value_F_test)


# ###########################
# TESTING GBPUSD DATA
# ###########################

# Parameters
set.seed(123)
z <- 1 # order to test
hyp_simp <- FALSE # we test composite hypothesis
hyp_comp <- TRUE
est_method  <- "fpcr_l1s" # FLMFR estimation method
thre_p <- thre_q <- 0.995 # Initial (p, q) for capturing 99.5% of variance
cv_1se <- TRUE # optimal lambda be the lambda.1se, as returned by cv.glmnet?
B <- 1000 # Bootstrap replicates

# STAGE 1: Is X_n an ARH(1) process (that is, X_n vs X_{n-1} FLMFR)?
X <- list("fstoch_proc" = fda_data$GBPUSD)
testing_ARHz <- ARHz_test(X, z = z, B = B, hyp_simp = hyp_simp,
                          hyp_comp = hyp_comp, est_method = est_method,
                          thre_p = thre_p, thre_q = thre_q, cv_1se = cv_1se)
testing_ARHz$testing_comp$p.value

# STAGE 2: F-test OU vs ARH(1) (unrestrictedd)
testing_F_OU <- F_stat_OU(testing_ARHz$X_flmfr, X, est_method = est_method,
                          thre_p = thre_p, thre_q = thre_q, cv_1se = FALSE)
F_stat <- testing_F_OU$F_stat
F_stat_ast <- theta_est_ast <- rep(0, B) 

# Parameter estimates
beta_est <- testing_F_OU$pred_OU$theta_est
sigma_est <- abs(optimize(mce, r = as.vector(data$GBPUSD),
                          Delta = 1/(length(t) - 1),
                          interval = c(0, 2))$minimum)
par.sde_ast <- list("alpha" = 0, "beta" = beta_est,
                    "sigma" = sigma_est) # centered process, alpha = mu = 0

# Parametric boostrap
pb <- txtProgressBar(max = B, style = 3)
for (b in 1:B) {
  
  # Simulating bootstrap replicates of X
  X_ast <- r_stoch_proc(dim(fda_data$GBPUSD$data)[1], t = t,
                        par.sde = par.sde_ast, mu = 0, X0 = 0,
                        warm_up = -1, type = "CKLS", model = "OU",
                        verbose = FALSE, plot = FALSE)
  
  # Converting X_ast into a FLMFR X_n vs X_{n-1}
  X_flmfr_ast <- ARH_to_FLMFR(X_ast[["fstoch_proc"]], z)
  
  # Computing bootstrap F-statistics
  testing_F_OU_ast <- F_stat_OU(X_flmfr_ast, X_ast, est_method = est_method,
                                thre_p = thre_p, thre_q = thre_q,
                                cv_1se = FALSE, verbose = FALSE)
  theta_est_ast[b] <- testing_F_OU_ast$pred_OU$theta_est
  F_stat_ast[b]    <- testing_F_OU_ast$F_stat
  
  setTxtProgressBar(pb, b)
}

# Approximation of the p-value by MC
p_value_F_test <- mean(F_stat < F_stat_ast)

# p-values
c(testing_ARHz$testing_comp$p.value, p_value_F_test)

```
</details>
 

License
----------

[![License:
GPLv3](https://img.shields.io/badge/license-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

<!-- <img src="" alt="goffda  hexlogo" align="right" width="200" style="padding: 0 15px; float: right;"/> -->

**Feel free to use** any of the contents but don't forget to cite them.

References
----------

Álvarez-Liébana, J., López-Perez, A., González-Manteiga, W. and Febrero-Bande, M. (2021). 
A goodness-of-fit test for functional time series with applications to diffusion processes. *arXiv: *
<a href="https://arxiv.org" class="uri">https://arxiv.org</a>

García-Portugués, E., Álvarez-Liébana, J., Álvarez-Pérez, G. and
González-Manteiga, W. (2021). A goodness-of-fit test for the functional
linear model with functional response. *Scand J Statist.* 2021; 48: 502-528.
<a href="https://doi.org/10.1111/sjos.12486" class="uri">https://doi.org/10.1111/sjos.12486</a>

Álvarez-Liébana, J., Bosq, D., and Ruiz-Medina, M. D. (2016). Consistency of the plug-in functional predictor of the Ornstein–Uhlenbeck process in Hilbert and Banach spaces. *Stat. Probab. Letters* 2016; 117: 12-22. <a href="https://doi.org/10.1016/j.spl.2016.04.023" class="uri">https://doi.org/10.1016/j.spl.2016.04.023</a>

Corradi, V. and White, H. (1999). Specification Tests for the Variance of a Diffusion. *J. Time Ser. Anal. 1999; 20: 253-270. <a href="https://doi.org/10.1111/1467-9892.00136" class="uri">https://doi.org/10.1111/1467-9892.00136</a>


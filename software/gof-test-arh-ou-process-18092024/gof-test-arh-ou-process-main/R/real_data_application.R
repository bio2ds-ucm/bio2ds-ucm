# ###############################################################
# Description: companion software to replicate the real
#              data application of the article
#              «A goodness-of-fit test for functional time series
#              with applications to diffusion processes»
# Authors: A. López-Pérez and J. Álvarez-Liébana
# ###############################################################

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

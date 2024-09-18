
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
             
fig1a <- 
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
  # annotate(geom = "curve", x = 0.33, y = 0.02, xend = 0.8, yend = 0.01, 
  #          color = "firebrick", curvature = -0.5, size = 1.1,
  #          arrow = arrow(length = unit(3.5, "mm"))) +
  annotate(geom = "text", x = 0.5, y = 0.04,
           label = TeX("$X_{1}(t)$"), hjust = "center",
           family = "Times New Roman", size = 7) +
  # annotate(geom = "curve", x = 1.55, y = 0.026, xend = 1.4, yend = 0.017, 
  #          color = "firebrick", curvature = 0.4, size = 1.1,
  #          arrow = arrow(length = unit(3.5, "mm"))) +
  annotate(geom = "text", x = 1.5, y = 0.04,
           label = TeX("$X_{2}(t)$"), hjust = "center",
           family = "Times New Roman", size = 7) +
  # annotate(geom = "curve", x = 2.65, y = 0.01, xend = 2.83, yend = 0, 
  #          color = "firebrick", curvature = -0.7, size = 1.1,
  #          arrow = arrow(length = unit(3.5, "mm"))) +
  annotate(geom = "text", x = 2.5, y = 0.04,
           label = TeX("$X_{3}(t)$"), hjust = "center",
           family = "Times New Roman", size = 7) +
  # annotate(geom = "curve", x = 3.25, y = 0.025, xend = 3.45, yend = 0.01, 
  #          color = "firebrick", curvature = 0.3, size = 1.1,
  #          arrow = arrow(length = unit(3.5, "mm"))) +
  annotate(geom = "text", x = 3.5, y = 0.04,
           label = TeX("$X_{4}(t)$"), hjust = "center",
           family = "Times New Roman", size = 7)  +
  # annotate(geom = "curve", x = 4.7, y = 0.025, xend = 4.55, yend = 0.005, 
  #          color = "firebrick", curvature = 0.5, size = 1.1,
  #          arrow = arrow(length = unit(3.5, "mm"))) +
  annotate(geom = "text", x = 4.5, y = 0.04,
           label = TeX("$X_{5}(t)$"), hjust = "center",
           family = "Times New Roman", size = 7) +
  # Axis
  labs(x = "t",  y = "") #+
  # Title (splitted in two lines)
  #ggtitle("Splitting trajectories\nof an OU process")

# Save as 8x4 png 
ggsave("fig1a_plotOU.png", width = 8, height = 5)

# OU plotted as ARH(1) process
OU_ARH1 <- data.frame("t" = rep(seq(0, 1, l = 2000), 5),
                      "x" = rev(rev(OU$x)[-1]),
                      "interv" = rev(rev(OU$interv)[-1]))
fig1b <- 
  ggplot(OU_ARH1, aes(x = t, y = x, color = interv)) +
  geom_line(size = 1.3) + # width of line
  # Color pattern from tableau appareance according to intervals
  scale_color_tableau() +
  annotate(geom = "text", x = 0.34, y = 0.02,
           label = TeX("$X_{1}(t)$"),
           hjust = "right", family = "Times New Roman", size = 7) +
  annotate(geom = "text", x = 0.94, y = 0.03,
           label = TeX("$X_{2}(t)$"),
           hjust = "right", family = "Times New Roman", size = 7) +
  annotate(geom = "text", x = 0.85, y = -0.03,
           label = TeX("$X_{3}(t)$"),
           hjust = "left", family = "Times New Roman", size = 7) +
  annotate(geom = "text", x = 0.515, y = 0.0155,
           label = TeX("$X_{4}(t)$"),
           hjust = "right", family = "Times New Roman", size = 7) +
  annotate(geom = "text", x = 0.037, y = -0.011,
           label = TeX("$X_{5}(t)$"),
           hjust = "right", family = "Times New Roman", size = 7) +
  # Axis
  labs(x = "t", y = "") #+
  # # Title (splitted in two lines)
  # ggtitle("OU stochastic model characterized\nas an ARH(1) process") 
# Save as 8x4 png 
ggsave("fig1b_plotOU.png", width = 8, height = 5)




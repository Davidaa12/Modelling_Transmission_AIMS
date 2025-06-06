# ====================================================================
# Concave-Down Quadratic Model Fit for Vector Competence (bc.succ)
# Species: Anopheles stephensi
# Author: David A. Aderinkola
# Email: davidaa@aims.ac.za 
# Description: Fits a concave-down quadratic model to the 'bc.succ' trait
# using Bayesian inference (Stan) and visualises posterior estimates.
# ====================================================================

# ---- Load Required Libraries ----
library(tidyverse)   #Load a collection of R packages for data manipulation (e.g. dplyr, tidyr) and visualisation (ggplot2)
library(rstan)       # For Bayesian modelling via Stan

# ---- Load and Filter Trait Data ----
data <- read.csv("traits.csv")

# Filter for 'bc.succ' trait of An. stephensi with non-zero values
data_bc.succ <- filter(data, trait.name == "bc.succ" & trait != 0 & specie == "An. stephensi")

# ---- Prepare Data for Stan ----
x_t <- seq(5, 40, 0.5)  # Sequence of temperatures for posterior predictions

data_T <- list(
  x = data_bc.succ$T,                            # Observed temperatures
  y = data_bc.succ$trait / 100,                  # Trait values scaled to [0,1]
  N = nrow(data_bc.succ),                        # Number of observations
  x_t = x_t,                                     # Prediction temperatures
  N_t = length(x_t)                              # Number of prediction points
)

# ---- Fit the Quadratic Model via Stan ----
Quad_Model <- stan(
  file = "David_bc_succ.stan",                  # Stan file (concave-down quadratic function)
  data = data_T,
  iter = 2000,
  warmup = 1000,
  chains = 4,
  control = list(adapt_delta = 0.8)             # Improve sampling stability
)

# ---- Summarise Posterior Output ----
print(Quad_Model)

# ---- Extract Posterior Predictions ----
posterior_samples <- rstan::extract(Quad_Model)
mu_t_samples <- posterior_samples$mu_t           # Matrix of posterior draws for predicted trait values (rows = iterations, columns = temperatures)

# Save x_t and mu_t_samples as an RDS file
saveRDS(list(
  x_t = x_t,
  mu_t_samples = mu_t_samples
), file = "bc_model.rds")

# ---- Compute Posterior Summaries ----
mu_t_mean  <- apply(mu_t_samples, 2, mean)        # Posterior mean of the predicted trait values across temperatures
mu_t_lower <- apply(mu_t_samples, 2, quantile, probs = 0.025)        # 2.5% CI
mu_t_upper <- apply(mu_t_samples, 2, quantile, probs = 0.975)        # 97.5% CI

# Combine results into a tidy data frame for plotting
df_mu <- data.frame(
  x_t = x_t,
  mu_mean = mu_t_mean,
  mu_lower = mu_t_lower,
  mu_upper = mu_t_upper
)


# ---- Visualise Posterior Predictions with Credible Intervals ----
ggplot(df_mu, aes(x = x_t, y = mu_mean)) +
  geom_line(color = "red", linewidth = 0.7) +
  geom_ribbon(aes(ymin = mu_lower, ymax = mu_upper), fill = "blue", alpha = 0.3) +
  geom_point(data = data_bc.succ, aes(x = T, y = trait / 100), colour = "black") +
  scale_x_continuous(breaks = seq(5, 40, by = 5)) +  # Set x-axis labels at 5-unit intervals
  labs(
    x = "Temperature (°C)",
    y = "bc trait value",
    #y = "Predicted Trait Value (mu_t, scaled)",
    title = "Vector Competence (bc) with 95% Credible Interval"
  ) +
  theme_minimal()

#ggsave("bc.succ_plot.png", width = 7, height = 5, dpi = 300)

ggsave("bc.succ_plot.pdf", width = 7, height = 5)

##########################################################################################################################

# ---- Extract Posterior Samples for Parameters ----
posterior_samples <- rstan::extract(Quad_Model)

# Extract Tmin, Tmax, and compute Topt for each posterior draw
Tmin_samples <- posterior_samples$Tmin
Tmax_samples <- posterior_samples$Tmax
Topt_samples <- (Tmin_samples + Tmax_samples) / 2

# Summarise Topt
Topt_mean <- mean(Topt_samples)
Topt_ci <- quantile(Topt_samples, probs = c(0.025, 0.975))

# Summarise Tmin
Tmin_mean <- mean(Tmin_samples)
Tmin_ci <- quantile(Tmin_samples, probs = c(0.025, 0.975))

# Summarise Tmax
Tmax_mean <- mean(Tmax_samples)
Tmax_ci <- quantile(Tmax_samples, probs = c(0.025, 0.975))

# Print results
cat("Tmin: ", Tmin_mean, " [", Tmin_ci[1], ", ", Tmin_ci[2], "]\n")
cat("Tmax: ", Tmax_mean, " [", Tmax_ci[1], ", ", Tmax_ci[2], "]\n")
cat("Topt: ", Topt_mean, " [", Topt_ci[1], ", ", Topt_ci[2], "]\n")



###########################################################################################################################

# ====================================================================
# Posterior Diagnostics and Visualisation for Concave-Down Quadratic Model
# Species: Anopheles stephensi
# Description: Generates trace and density plots for key parameters
# from the quadratic thermal response model using bayesplot and patchwork.
# ====================================================================

# ---- Load Required Packages ----
library(bayesplot)
#library(ggplot2)
library(patchwork)  # for combining plots

posterior_array <- as.array(Quad_Model)  # Convert stanfit object to 3D array [iterations, chains, parameters]

# --------------------------------------------------------------------
#                         TRACEPLOTS
# --------------------------------------------------------------------
# Combined traceplot for Tmin, Tmax, and sigma
mcmc_trace(posterior_array, pars = c("Tmin", "Tmax", "f_max", "sigma")) +
  ggtitle("bc Traceplots for Tmin, Tmax, f_max and \u03C3")

ggsave("bc_tp.pdf", width = 7, height = 5)


# --------------------------------------------------------------------
#                        POSTERIOR DENSITY PLOTS
# --------------------------------------------------------------------
# Create individual plots
p1 <- mcmc_areas(posterior_array, pars = "Tmin", prob = 0.95) +
  ggtitle("Posterior Distribution of Tmin")

p2 <- mcmc_areas(posterior_array, pars = "Tmax", prob = 0.95) +
  ggtitle("Posterior Distribution of Tmax")

p3 <- mcmc_areas(posterior_array, pars = "f_max", prob = 0.95) +
  ggtitle("Posterior Distribution of f_max")

p4 <- mcmc_areas(posterior_array, pars = "sigma", prob = 0.95) +
  ggtitle("Posterior Distribution of \u03C3")  # Unicode for σ

# Combine plots in a 2x2 grid
combined_plot <- (p1 | p2) / (p3 | p4)+
  plot_annotation(title = "bc Posterior Distributions of Key Parameters")


# Display the combined plot
print(combined_plot)

# Save to PDF using cairo_pdf to support Unicode
ggsave("bc_dp.pdf", plot = combined_plot, width = 7, height = 5, device = cairo_pdf)

# ====================================================================
# Concave-Up Quadratic Model Fit for Mortality Rate (mu)
# Species: Anopheles stephensi
# Author: David A. Aderinkola
# Email: davidaa@aims.ac.za 
# Description: Fits a U-shaped quadratic model to the 'mu' trait
# using Bayesian inference (Stan) and visualises posterior output.
# ====================================================================

# ---- Load Required Libraries ----
library(tidyverse)   # Load a collection of R packages for data manipulation (e.g. dplyr, tidyr) and visualisation (ggplot2)
library(rstan)       # For Bayesian modelling in Stan

# ---- Set working directory ----
## setwd("...")

# ---- Load and Filter Trait Data ----
data <- read.csv("traits.csv")

# Filter for mortality rate (mu) for An. stephensi with non-zero values
data_mu <- filter(data, trait.name == "mu" & trait != 0 & specie == "An. stephensi")

# ---- Prepare Data for Stan Model ----
x_t <- seq(5, 40, 0.5)  # Sequence of temperatures for prediction

data_T <- list(
  x = data_mu$T,                        # Observed temperatures
  y = data_mu$trait,                    # Observed trait values
  N = nrow(data_mu),                    # Number of observations
  x_t = x_t,                            # Prediction temperatures
  N_t = length(x_t)                     # Number of prediction points
)

# ---- Fit Quadratic-Up Stan Model ----
QuadUp_Model <- stan(
  file = "David_mu.stan",              # Stan file implementing concave-up quadratic function
  data = data_T,
  iter = 2000,
  warmup = 1000,
  chains = 4,
  control = list(adapt_delta = 0.99, max_treedepth = 12)  # Enhance sampling stability
)

# ---- Summarise Posterior Output ----
print(QuadUp_Model)

# ---- Extract Posterior Predictions ----
posterior_samples <- rstan::extract(QuadUp_Model)
mu_t_samples <- posterior_samples$mu_t                     # Posterior draws of predicted values

# ---- Posterior Summary Statistics ----
mu_t_mean  <- apply(mu_t_samples, 2, mean)                           # Compute the posterior mean of the trait at each temperature
mu_t_lower <- apply(mu_t_samples, 2, quantile, probs = 0.025)        # 2.5% CI
mu_t_upper <- apply(mu_t_samples, 2, quantile, probs = 0.975)        # 97.5% CI

# Save x_t and mu_t_samples as an RDS file
saveRDS(list(
  x_t = x_t,
  mu_t_samples = mu_t_samples
), file = "mu_model.rds")

# Combine into a tidy data frame for plotting
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
  geom_point(data = data_mu, aes(x = T, y = trait), colour = "black") +
  scale_x_continuous(breaks = seq(5, 40, by = 5)) +  # Set x-axis labels at 5-unit intervals
  labs(
    x = "Temperature (°C)",
    y = "mu trait value",
    title = "Mortality Rate (mu) with 95% Credible Interval"
  ) +
  theme_minimal()

#ggsave("mu_plot.png", width = 7, height = 5, dpi = 300)

ggsave("mu_plot.pdf", width = 7, height = 5)

##########################################################################################################################
##########################################################################################################################

# ---- Extract Posterior Samples for Concave-Up Quadratic Model ----
posterior_samples <- rstan::extract(QuadUp_Model)

# Extract Topt directly from the posterior
Topt_samples <- posterior_samples$Topt

# Extract q and Topt to compute beta and gamma
q_samples <- posterior_samples$q
f_min_samples <- posterior_samples$f_min

# Compute beta and gamma from q and Topt
beta_samples <- -2 * q_samples * Topt_samples
gamma_samples <- f_min_samples + q_samples * Topt_samples^2

# Topt is already sampled; compute Tmin and Tmax using quadratic roots near f_min + threshold
compute_bounds <- function(q, beta, gamma, f_min, threshold = 0.01) {
  roots <- polyroot(c(gamma - f_min - threshold, beta, q))
  real_roots <- Re(roots[abs(Im(roots)) < 1e-8])
  sort(real_roots)
}

bounds <- t(mapply(compute_bounds, q_samples, beta_samples, gamma_samples, f_min_samples))
Tmin_samples <- bounds[, 1]
Tmax_samples <- bounds[, 2]

# Summary statistics function
summary_stats <- function(samples) {
  c(mean = mean(samples), lower = quantile(samples, 0.025), upper = quantile(samples, 0.975))
}

# Apply summary function
Topt_summary <- summary_stats(Topt_samples)
Tmin_summary <- summary_stats(Tmin_samples)
Tmax_summary <- summary_stats(Tmax_samples)

# Output all
list(Topt = Topt_summary, Tmin = Tmin_summary, Tmax = Tmax_summary)


###########################################################################################################################

# ====================================================================
# Posterior Diagnostics and Visualisation for Concave-Up Quadratic Model
# Description: Generates trace and density plots for key parameters
# in the Bayesian quadratic using the bayesplot & patchwork packages.
# ====================================================================

# ---- Load Required Packages ----
library(bayesplot)
#library(ggplot2)
library(patchwork)  # for combining plots

posterior_array <- as.array(QuadUp_Model)  # Convert stanfit object to 3D array [iterations, chains, parameters]

# --------------------------------------------------------------------
#                         TRACEPLOTS
# --------------------------------------------------------------------
# Combined traceplot for Tmin, Tmax, and sigma
mcmc_trace(posterior_array, pars = c("beta", "gamma", "f_min", "sigma")) +
  ggtitle("mu Traceplots for \u03B2, \u03B3, f_max and \u03C3")

ggsave("mu_tp.pdf", width = 7, height = 5)



# --------------------------------------------------------------------
#                        POSTERIOR DENSITY PLOTS
# --------------------------------------------------------------------
# Create individual plots
p1 <- mcmc_areas(posterior_array, pars = "beta", prob = 0.95) +
  ggtitle("Posterior Distribution of \u03B2")

p2 <- mcmc_areas(posterior_array, pars = "gamma", prob = 0.95) +
  ggtitle("Posterior Distribution of \u03B3")

p3 <- mcmc_areas(posterior_array, pars = "f_min", prob = 0.95) +
  ggtitle("Posterior Distribution of f_min")

p4 <- mcmc_areas(posterior_array, pars = "sigma", prob = 0.95) +
  ggtitle("Posterior Distribution of \u03C3")  # Unicode for σ

# Combine plots in a 2x2 grid
combined_plot <- (p1 | p2) / (p3 | p4)+
  plot_annotation(title = "mu Posterior Distributions of Key Parameters")


# Display the combined plot
print(combined_plot)

# Save to PDF using cairo_pdf to support Unicode
ggsave("mu_dp.pdf", plot = combined_plot, width = 7, height = 5, device = cairo_pdf)

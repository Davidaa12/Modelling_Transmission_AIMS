# ====================================================================
# Briere Model Fit for Parasite Development Rate (pdr)
# Species: Anopheles stephensi
# Author: David A. Aderinkola
# Email: davidaa@aims.ac.za
# Description: Fits a Briere thermal response model to trait 'pdr'
# using Bayesian inference (Stan) and visualises posterior estimates.
# ====================================================================

# ---- Load Required Libraries ----
library(tidyverse)   # Load a collection of R packages for data manipulation (e.g. dplyr, tidyr) and visualisation (ggplot2)
library(rstan)       # Bayesian modelling in Stan

# ---- Load and Filter Trait Data ----
data <- read.csv("traits.csv")

# Filter for 'pdr' trait of An. stephensi with non-zero values
data_pdr <- filter(data, trait.name == "pdr" & trait != 0 & specie == "An. stephensi")

# ---- Prepare Data for Stan Model ----
x_t <- seq(5, 40, 0.5)  # Temperature range for predictions

data_T <- list(
  x = data_pdr$T,                       # Observed temperatures
  y = data_pdr$trait,                   # Observed trait values (pdr)
  N = nrow(data_pdr),                   # Number of observations
  x_t = x_t,                            # Temperatures for prediction
  N_t = length(x_t)                     # Number of prediction points
)

# ---- Fit the Briere Model ----
Briere_Model <- stan(
  file = "David_pdr.stan",             # Stan model file (Briere function)
  data = data_T,
  iter = 2000,
  warmup = 1000,
  chains = 4,
  control = list(adapt_delta = 0.8)    # Improve sampling stability
)

# ---- Summarise Model Output ----
print(Briere_Model)

# ---- Extract Posterior Predictions ----
posterior_samples <- rstan::extract(Briere_Model)
mu_t_samples <- posterior_samples$mu_t         # Matrix of posterior draws for predicted trait values (rows = iterations, columns = temperatures)
#mu_t_samples
# ---- Compute Posterior Summaries ----
mu_t_mean  <- apply(mu_t_samples, 2, mean)                           # Compute the posterior mean of the trait at each temperature
mu_t_lower <- apply(mu_t_samples, 2, quantile, probs = 0.025)        # Lower 95% CI
mu_t_upper <- apply(mu_t_samples, 2, quantile, probs = 0.975)        # Upper 95% CI

# Save x_t and mu_t_samples as an RDS file
saveRDS(list(
  x_t = x_t,
  mu_t_samples = mu_t_samples
), file = "pdr_model.rds")

# Combine into a tidy data frame for plotting
df_mu <- data.frame(
  x_t = x_t,
  mu_mean = mu_t_mean,
  mu_lower = mu_t_lower,
  mu_upper = mu_t_upper
)

# ---- Visualise Posterior Predictions ----
ggplot(df_mu, aes(x = x_t, y = mu_mean)) +
  geom_line(color = "red", size = 0.7) +
  geom_ribbon(aes(ymin = mu_lower, ymax = mu_upper), fill = "blue", alpha = 0.3) +
  geom_point(data = data_pdr, aes(x = T, y = trait), colour = "black") +
  scale_x_continuous(breaks = seq(5, 40, by = 5)) +  # Set x-axis labels at 5-unit intervals
  labs(
    x = "Temperature (°C)",
    y = "pdr trait value",
    title = "Parasite Development Rate (pdr) with 95% Credible Interval"
  ) +
  theme_minimal()

#ggsave("pdr_plot.png", width = 7, height = 5, dpi = 300)

ggsave("pdr_plot.pdf", width = 7, height = 5)

##########################################################################################################################

# ---- Extract Posterior Samples for Parameters ----
posterior_samples <- rstan::extract(Briere_Model)

# Extract Tmin, Tmax, and compute Topt for each posterior draw
Tmin_samples <- posterior_samples$Tmin
Tmax_samples <- posterior_samples$Tmax

# Compute Topt using the analytical formula from the Stan function
Topt_samples <- (4*Tmax_samples + 3*Tmin_samples + sqrt(16*Tmax_samples^2 + 9*Tmin_samples^2 - 16*Tmax_samples*Tmin_samples)) / 10


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


##########################################################################################################################
# ====================================================================
# Posterior Diagnostics and Visualisation for Briere Model
# Description: Generates trace plots and density plots for key parameters
# in the Bayesian Briere model using the bayesplot & patchwork packages.
# ====================================================================

# ---- Load Required Packages ----
library(bayesplot)
#library(ggplot2)
library(patchwork)  # for combining plots

posterior_array <- as.array(Briere_Model)  # Convert stanfit object to 3D array [iterations, chains, parameters]

# --------------------------------------------------------------------
#                         TRACEPLOTS
# --------------------------------------------------------------------
# Combined traceplot for Tmin, Tmax, and sigma
mcmc_trace(posterior_array, pars = c("Tmin", "Tmax", "f_max", "sigma")) +
  ggtitle("pdr Traceplots for Tmin, Tmax, f_max and \u03C3")

ggsave("pdr_tp.pdf", width = 7, height = 5)



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
  plot_annotation(title = "pdr Posterior Distributions of Key Parameters")


# Display the combined plot
print(combined_plot)

# Save to PDF using cairo_pdf to support Unicode
ggsave("pdr_dp.pdf", plot = combined_plot, width = 7, height = 5, device = cairo_pdf)

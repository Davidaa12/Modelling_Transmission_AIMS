# ====================================================================
# Computation of Temperature-Dependent Relative R_0(T)
# Based on vector and parasite traits modelled in Stan
# Author: David A. Aderinkola
# Email: davidaa@aims.ac.za
# ====================================================================

# ---- Load Required Libraries ----
library(tidyverse)

# ---- Load Posterior samples of Traits from RDS Files ----
mu_samps  <- readRDS("mu_model.rds")$mu_t_samples
a_samps   <- readRDS("a_model.rds")$mu_t_samples
bc_samps  <- readRDS("bc_model.rds")$mu_t_samples
pdr_samps <- readRDS("pdr_model.rds")$mu_t_samples
efd_samps <- readRDS("efd_model.rds")$mu_t_samples
e2a_samps <- readRDS("e2a_model.rds")$mu_t_samples
mdr_samps <- readRDS("mdr_model.rds")$mu_t_samples

# ---- Load temperature vector (assumed consistent across models) ----
x_t <- readRDS("mu_model.rds")$x_t    # Consistency across all traits
x_t

# ---- Define Constants ----
r    <- 0.03        # Human recovery rate
mu_p <- 0.083       # Pre-adult mortality

# ---- Compute R_0(T) using Posterior Samples ----
R0_samples <- (
  (a_samps^2 * bc_samps * exp(-mu_samps / pdr_samps)) / (mu_samps * r)
) *
  ((mdr_samps * e2a_samps) / mu_samps) *
  ((efd_samps * mdr_samps * e2a_samps) / mu_samps - (mu_p + mdr_samps))

# ---- Normalise R_0(T) Samples to [0, 1] ----
R0_samples_norm <- sweep(R0_samples, 1, apply(R0_samples, 1, max), FUN = "/")

# ---- Summarise Posterior for Normalised R_0(T) ----
R0_mean  <- apply(R0_samples_norm, 2, mean)
R0_lower <- apply(R0_samples_norm, 2, quantile, probs = 0.025)
R0_upper <- apply(R0_samples_norm, 2, quantile, probs = 0.975)

# ---- Create Data Frame for Plotting ----
df_R0 <- data.frame(
  temperature = x_t,
  R0_mean = R0_mean,
  R0_lower = R0_lower,
  R0_upper = R0_upper
)

# ---- Plot Normalised R_0(T) with 95% Credible Interval ----
ggplot(df_R0, aes(x = temperature, y = R0_mean)) +
  geom_line(color = "red", size = 0.7) +
  geom_ribbon(aes(ymin = R0_lower, ymax = R0_upper), fill = "blue", alpha = 0.3) +
  scale_x_continuous(breaks = seq(5, 40, by = 5)) +
  labs(
    x = "Temperature (°C)",
    y = expression(R[0]^rel*(T)),
    title = expression("Basic Reproduction Number " ~ R[0]^rel*(T) ~ "with 95% Credible Interval")
  ) +
  theme_minimal()

ggsave("R0_Relative_Normalised.pdf", width = 7, height = 5)

# ---- Optimal Temperature and CI ----
opt_index <- apply(R0_samples_norm, 1, which.max)
opt_temps <- x_t[opt_index]

opt_temp_mean <- mean(opt_temps)
opt_temp_ci <- quantile(opt_temps, probs = c(0.025, 0.975))

print(paste("Optimal Temperature:", round(opt_temp_mean, 2)))
print(paste("95% CI:", round(opt_temp_ci[1], 2), "-", round(opt_temp_ci[2], 2)))

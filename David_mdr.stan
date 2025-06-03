functions {
  // Modified Briere function to model thermal performance curve
  // Inputs: 
  //   x      - temperature
  //   f_max  - peak trait value
  //   Tmin   - lower thermal threshold
  //   Tmax   - upper thermal threshold
  // Computes an empirical thermal response curve with peak at approximated Topt
  real briere(real x, real f_max, real Tmin, real Tmax) {
    real out;
    real Topt;
    // Analytical approximation of temperature at peak trait value
    Topt = (4*Tmax + 3*Tmin + sqrt( 16*Tmax^2 + 9*Tmin^2  - 16*Tmax*Tmin ) ) / 10;
    // Scaling coefficient for the function shape
    real q = f_max / (Topt * (Topt - Tmin) * sqrt(Tmax - Topt));
    
    // Evaluate function only within thermal limits
    if (x > Tmin && x < Tmax && x > 0) {
      out = q * x * (x - Tmin) * sqrt(Tmax - x);
    } else {
      out = 0;
    }
    return out;
  }
}

//f_max=q*Topt*(Topt-Tmin)*sqrt(Tmax-Topt)

data {
  int N;               // Number of observations
  vector[N] x;         // Observed temperatures
  vector[N] y;         // Observed trait values (e.g. mosquito development rate)
  int N_t;             // Number of temperatures for prediction
  vector[N_t] x_t;     // Temperatures for prediction
}

parameters {
  real<lower=0> f_max;     // Maximum trait value (peak height)
  real Tmin;               // Minimum temperature threshold
  real Tmax;               // Maximum temperature threshold (constrained > Tmin)
  real<lower=0> sigma;     // Observation error standard deviation
}

transformed parameters {
  vector[N] mu;  // Expected trait values
  for (i in 1:N) {
    mu[i] = briere(x[i], f_max, Tmin, Tmax);  // Apply Brière function to each observation
  }
}

model {
  // Priors for thermal limits and peak trait value
  Tmin ~ normal(5, 10);
  Tmax ~ normal(40, 10);
  f_max ~ normal(1, 1);
  
  // Gaussian likelihood for observed data
  y ~ normal(mu, sigma);
}

generated quantities {
  vector[N_t] mu_t;  // Predicted trait values at new temperatures
  for (i in 1:N_t) {
    mu_t[i] = briere(x_t[i], f_max, Tmin, Tmax);
  }
}

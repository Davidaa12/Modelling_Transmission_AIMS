functions {
  // Concave-down quadratic thermal response function
  // Models trait performance with a symmetric peak at Topt
  // Inputs:
  //   x      - temperature
  //   f_max  - peak trait value
  //   Tmin   - minimum temperature threshold
  //   Tmax   - maximum temperature threshold
  real quadratic(real x, real f_max, real Tmin, real Tmax) {
    real out;
    real Topt;
    Topt = (Tmax + Tmin) / 2;  // Peak temperature (symmetric midpoint)
    real q = f_max / ((Topt - Tmin) * (Tmax - Topt));  // Scaling coefficient
    
    // Evaluate response only within thermal limits
    if (x > Tmin && x < Tmax && x > 0) {
      out = q * (x - Tmin) * (Tmax - x);
    } else {
      out = 0;
    }
    
    // Optional upper bound on output
    // if (out > 1) { out = 1; }
    
    return out;
  }
}

// f_max = q * (Topt - Tmin) * (Tmax - Topt)

data {
  int N;               // Number of observations
  vector[N] x;         // Observed temperatures
  vector[N] y;         // Observed trait values
  int N_t;             // Number of temperatures for prediction
  vector[N_t] x_t;     // Temperatures for prediction
}

parameters {
  real<lower=0> f_max;     // Maximum trait value (peak height)
  real Tmin;               // Lower thermal threshold
  real Tmax;               // Upper thermal threshold
  real<lower=0> sigma;     // Observation error standard deviation
}

transformed parameters {
  vector[N] mu;  // Expected trait values at observed temperatures
  for (i in 1:N) {
    mu[i] = quadratic(x[i], f_max, Tmin, Tmax);
  }
}

model {
  // Priors for parameters
  Tmin ~ normal(5, 10);
  Tmax ~ normal(40, 10);
  f_max ~ normal(1, 1);

  // Likelihood: Gaussian error model
  y ~ normal(mu, sigma);
}

generated quantities {
  vector[N_t] mu_t;  // Predicted trait values at new temperatures
  for (i in 1:N_t) {
    mu_t[i] = quadratic(x_t[i], f_max, Tmin, Tmax);
  }
}

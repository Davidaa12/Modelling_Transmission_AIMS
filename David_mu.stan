functions {
  // Concave-up quadratic function: f(T) = q*T^2 + beta*T + gamma
  // Models U-shaped thermal performance curves
  real quad_up(real x, real q, real beta, real gamma) {
    return q * square(x) + beta * x + gamma;
  }
}

data {
  int N;                // Number of observations
  vector[N] x;          // Observed temperatures
  vector[N] y;          // Observed trait values
  int N_t;              // Number of prediction temperatures
  vector[N_t] x_t;      // Temperatures for prediction
}

parameters {
  real q;                      // Quadratic coefficient (must be > 0 for concave-up shape)
  real Topt;                   // Temperature where trait reaches minimum (vertex of the parabola)
  real<lower=0> f_min;         // Minimum trait value at Topt
  real<lower=0> sigma;         // Observation noise standard deviation
}

transformed parameters {
  real beta;                   // Linear term in quadratic function
  real gamma;                  // Intercept term (function value at T = 0)
  
  // Derived from setting f′(Topt) = 0: Topt = -beta / (2q)
  beta = -2 * q * Topt;

  // Compute gamma to ensure f(Topt) = f_min
  gamma = f_min + q * square(Topt);

  vector[N] mu;                // Expected trait values at observed temperatures
  for (i in 1:N) {
    mu[i] = quad_up(x[i], q, beta, gamma);
  }
}

model {
  // Priors for parameters
  gamma ~ normal(0.5, 0.5);     
  Topt ~ normal(25, 10);      
  f_min ~ normal(0.5, 0.5);    
 
  // Gaussian likelihood for trait observations
  y ~ normal(mu, sigma);
}

generated quantities {
  vector[N_t] mu_t;            // Predicted trait values at new temperatures
  for (i in 1:N_t) {
    mu_t[i] = quad_up(x_t[i], q, beta, gamma);
  }
}

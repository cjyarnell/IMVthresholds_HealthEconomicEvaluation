// ARIVe project
// Stan code for Bayesian multistate model 

functions {
}

data {
  
  int<lower=0> P; // number of beta parameters including intercept

  int<lower=0> N_cen; // total number of censored observations
  int<lower=0> N_obs; // total number of observed transitions
  
  // covariates 
  matrix[N_cen, P] X_cen; 
  matrix[N_obs, P] X_obs; 
  
  // observed times
  array[N_cen] real y_cen;
  array[N_obs] real y_obs;
}


parameters {
    vector[P] beta;
    real logalpha;
}

transformed parameters{
}

model {
  // priors
  beta ~ normal(0,0.5);
  logalpha ~ normal(0,0.4);
  {
      // model Weibull rate as function of covariates
    vector[N_cen] lambda_cen;
    vector[N_obs] lambda_obs;
    
    real alpha = exp(logalpha);
    
    // standard weibull AFT re-parameterization
    lambda_cen = exp((X_cen*beta)*alpha);
    lambda_obs = exp((X_obs*beta)*alpha);
  
    // likelihood function
    target += weibull_lcdf( y_obs | alpha, lambda_obs) + 
              weibull_lccdf(y_cen | alpha, lambda_cen);

  }
}

generated quantities {
// alpha
  real alpha = exp(logalpha);
// intercept
  real intercept = exp(beta[1] * alpha);
  vector[P-1] hazard_ratio;
  for (i in 2:P){hazard_ratio[i-1] = exp(beta[i] * alpha);}
}

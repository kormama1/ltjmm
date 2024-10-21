// Stacked mixed model with independent normal distributions on random effects

data{
	// Define variables in data
	int<lower=0> N_obs;              // Number of observations
	int<lower=0, upper=N_obs> N_sub; // Number of subjects
	int<lower=0, upper=N_obs> N_out; // Number of outcomes

	int<lower=0> N_X;                
	matrix[N_obs, N_X] X;            // fixed effects

	vector[N_obs] y;                 // outcome
	vector[N_obs] obs_time;          // observation times
  
	int<lower=1, upper=N_out> outcome[N_obs];      // outcome index
	int<lower=1, upper=N_sub> subject[N_obs];      // subject index
}

parameters{
	matrix[N_out, N_X] beta;        // parameter vector of observed covariates
	vector[N_out] gamma;
	vector<lower=0>[N_out] sigma_y; // population standard deviation 
 
	matrix[N_sub, N_out] alpha1;        // random slopes
	matrix[N_sub, N_out] alpha0;        // random intercepts 
 
	vector<lower=0>[N_out] sigma_alpha1;    // standard deviation of random slopes
	vector<lower=0>[N_out] sigma_alpha0;    // standard deviation of random intercepts
}

transformed parameters{
	matrix[N_obs, N_X] beta_aug;
	vector[N_obs] gamma_aug;
	vector[N_obs] alpha0_aug;
	vector[N_obs] alpha1_aug;
	vector[N_obs] sigma_y_aug;
	vector[N_obs] mu;

	for(n_obs in 1:N_obs){ 
		beta_aug[n_obs] = beta[outcome[n_obs]];
		gamma_aug[n_obs] = gamma[outcome[n_obs]];
		alpha1_aug[n_obs] = alpha1[subject[n_obs], outcome[n_obs]];
		alpha0_aug[n_obs] = alpha0[subject[n_obs], outcome[n_obs]];    
		sigma_y_aug[n_obs] = sigma_y[outcome[n_obs]];
		mu[n_obs] = dot_product(X[n_obs], beta_aug[n_obs]) + gamma_aug[n_obs] * obs_time[n_obs] + alpha0_aug[n_obs] + alpha1_aug[n_obs] * obs_time[n_obs];
	}
}

model{	
	// Priors
	gamma ~ cauchy(0, 2.5);
	sigma_y ~ cauchy(0, 2.5);
	sigma_alpha0 ~ cauchy(0, 2.5);
	sigma_alpha1 ~ cauchy(0, 2.5);
	  
	for(n_out in 1:N_out){
		beta[n_out] ~ normal(0, 10.0);
		alpha0[,n_out] ~ normal(0, sigma_alpha0[n_out]);
		alpha1[,n_out] ~ normal(0, sigma_alpha1[n_out]);
	}

	// Likelihood
	for(n_obs in 1:N_obs){
		y[n_obs] ~ normal(mu[n_obs], sigma_y_aug[n_obs]);
	} 
}

generated quantities{
	vector[N_obs] log_lik;           
	// log_likelihood (observations)
	for(n_obs in 1:N_obs){
		log_lik[n_obs] = normal_lpdf(y[n_obs] | mu[n_obs], sigma_y_aug[n_obs]);
	}
}

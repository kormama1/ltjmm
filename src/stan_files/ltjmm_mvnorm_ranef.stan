// LTJMM with multinormal distribution on random effects

data{

	// Define variables in data
	int<lower=0> N_obs;              // Number of observations
	int<lower=0, upper=N_obs> N_sub; // Number of subjects
	int<lower=0, upper=N_obs> N_out; // Number of outcomes

	int<lower=0> N_X;                // Dimension of fixed effects
	matrix[N_obs, N_X] X;            // Fixed effects

	vector[N_obs] y;                 // Outcome
	vector[N_obs] obs_time;          // Observation times
  
	int<lower=1, upper=N_out> outcome[N_obs];      // Outcome index
	int<lower=1, upper=N_sub> subject[N_obs];      // Subject index
}

parameters{
	matrix[N_out, N_X] beta;
	vector<lower=0>[N_out] gamma;    // Restrict gamma to be positive for identifiability
	vector<lower=0>[N_out] sigma_y;
	real<lower=0> sigma_delta;
	vector[N_sub] delta;
	cholesky_factor_corr[2*N_out-1] Lcorr;
	vector<lower=0>[2*N_out-1] sigma_diag;
	matrix[(2*N_out-1), N_sub] z_alpha;
}

transformed parameters{ 
	matrix[2*N_out-1, 2*N_out-1] sigma_L;
	matrix[N_sub, 2*N_out-1] alpha;
	matrix[N_sub, N_out-1] alpha0_raw;  // Raw random intercepts
	matrix[N_sub, N_out] alpha0;        
	matrix[N_sub, N_out] alpha1;        // Random slopes

	vector[N_obs] delta_aug;            
	vector[N_obs] alpha0_aug;
	vector[N_obs] alpha1_aug;
	matrix[N_obs, N_X] beta_aug;
	vector[N_obs] gamma_aug;
	vector[N_obs] sigma_y_aug;
	vector[N_obs] mu;
  
	sigma_L = diag_pre_multiply(sigma_diag, Lcorr);  // sigma_L is lower triangular matrix with positive diagonal elements
	alpha = (sigma_L * z_alpha)';                    // Random effects'
	alpha0_raw = alpha[, 1:(N_out-1)];
	alpha1 = alpha[, N_out:(2*N_out-1)];

	// Sum-to-zero constraint for subject-specific random intercepts over outcomes
	for(n_sub in 1:N_sub){
		for(n_out in 1:(N_out-1)){
			alpha0[n_sub, n_out] = alpha0_raw[n_sub, n_out];
		}
		alpha0[n_sub, N_out] = -sum(alpha0_raw[n_sub,]);
	}

	for(n_obs in 1:N_obs){ 
		beta_aug[n_obs] = beta[outcome[n_obs]];
		gamma_aug[n_obs] = gamma[outcome[n_obs]];
		delta_aug[n_obs] = delta[subject[n_obs]];
		alpha0_aug[n_obs] = alpha0[subject[n_obs], outcome[n_obs]]; 
		alpha1_aug[n_obs] = alpha1[subject[n_obs], outcome[n_obs]];    
		sigma_y_aug[n_obs] = sigma_y[outcome[n_obs]];
		mu[n_obs] = dot_product(X[n_obs], beta_aug[n_obs]) + gamma_aug[n_obs] * (obs_time[n_obs] + delta_aug[n_obs]) + alpha0_aug[n_obs] + alpha1_aug[n_obs] * obs_time[n_obs];
	}
}

model{
	// Priors
	gamma ~ cauchy(0, 2.5);
	for(n_out in 1:N_out){
		beta[n_out] ~ normal(0, 10);
	}
  
	to_vector(z_alpha) ~ normal(0, 1);  // implies delta_alpha[n_sub] ~ multi_normal_cholesky(0, sigma_L);
	Lcorr ~ lkj_corr_cholesky(2);
	sigma_diag ~ cauchy(0, 2.5);
	sigma_y ~ cauchy(0, 2.5);
	sigma_delta ~ cauchy(0, 2.5);
	delta ~ normal(0, sigma_delta); 

	// Likelihood
	for(n_obs in 1:N_obs){
		y[n_obs] ~ normal(mu[n_obs], sigma_y_aug[n_obs]);
	}
}

generated quantities{
	matrix[2*N_out-1, 2*N_out-1] Omega;  // correlation matrix 
	matrix[2*N_out-1, 2*N_out-1] Sigma;  // covariance matrix
	vector[N_obs] log_lik;               // log_likelihood (observations)

	Omega = multiply_lower_tri_self_transpose(Lcorr);
	Sigma = quad_form_diag(Omega, sigma_diag);
	
	// Computing log_likelihood for each subject
	for(n_obs in 1:N_obs){
		log_lik[n_obs] = normal_lpdf(y[n_obs] | mu[n_obs], sigma_y_aug[n_obs]);
	}
}

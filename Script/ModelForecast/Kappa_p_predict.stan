data {
	
	int<lower=0> T;    // No. of time periods observed
	
	int<lower=0> M;    // No. of variables in Kappa_t, i.e. 2	
	int<lower=1> P;    // Lag order

        vector[M] Kappa_1[T];
	vector[M] Kappa_2[T];
	vector[M] Kappa_3[T];
	vector[M] Kappa_4[T];
	vector[M] Kappa_5[T]; 
	vector[M] Kappa_bar[T]; // Kappa_bar for gravity   

	int<lower=0> T_new; // No. of time periods predictions
      }


transformed data 
   	{
	
	// read in vector of each observation from kappa
	vector[M] Kappa_1_obs[T-P];

	vector[M] Kappa_2_obs[T-P];

	vector[M] Kappa_3_obs[T-P];

	vector[M] Kappa_4_obs[T-P];

	vector[M] Kappa_5_obs[T-P];

		for (t in 1:(T-P))
			{
 
			Kappa_1_obs[t] = Kappa_1[t + P];

			Kappa_2_obs[t] = Kappa_2[t + P];

			Kappa_3_obs[t] = Kappa_3[t + P];

			Kappa_4_obs[t] = Kappa_4[t + P];

			Kappa_5_obs[t] = Kappa_5[t + P];
	

			}
	}


parameters 
	{
	
	// First parameter Psi is a scalar, since it's supposed to be the same for all
	real<lower=0,upper=1> Psi;
	
	cholesky_factor_corr[M] L_corr_noise;
	
	vector<lower=0>[M] sd_noise;
	
	vector[M] intercept;


	}


transformed parameters 
	{
	
	
	matrix[M,M] L_sigma;
  
	L_sigma = diag_pre_multiply(sd_noise, L_corr_noise);

	}


model 
	{
	
	
	vector[M] mus_1[T-P];
	vector[M] mus_2[T-P];
	vector[M] mus_3[T-P];
	vector[M] mus_4[T-P];
	vector[M] mus_5[T-P];

 
	for (t in 1:(T-P)) 
		{
		
		mus_1[t] = intercept + Kappa_1[t] - Psi * (Kappa_1[t]-Kappa_bar[t]);
		mus_2[t] = intercept + Kappa_2[t] - Psi * (Kappa_2[t]-Kappa_bar[t]);

		mus_3[t] = intercept + Kappa_3[t] - Psi * (Kappa_3[t]-Kappa_bar[t]);
		mus_4[t] = intercept + Kappa_4[t] - Psi * (Kappa_4[t]-Kappa_bar[t]);
		mus_5[t] = intercept + Kappa_5[t] - Psi * (Kappa_5[t]-Kappa_bar[t]);
		}
	
	L_corr_noise ~ lkj_corr_cholesky(2.0);

	sd_noise ~ inv_gamma(11,0.1);
	Psi ~ beta(2,2);

			Kappa_1_obs ~ multi_normal_cholesky(mus_1,L_sigma);
			Kappa_2_obs ~ multi_normal_cholesky(mus_2,L_sigma);
			Kappa_3_obs ~ multi_normal_cholesky(mus_3,L_sigma);
			Kappa_4_obs ~ multi_normal_cholesky(mus_4,L_sigma);
			Kappa_5_obs ~ multi_normal_cholesky(mus_5,L_sigma);

			
	}


generated quantities 
	{
	
	
	matrix[M,M] Sigma;
	
	vector[M] Kappa_1_pred[T_new-P];
	vector[M] Kappa_1_new[T_new];
	vector[M] Kappa_2_pred[T_new-P];
	vector[M] Kappa_2_new[T_new];
	vector[M] Kappa_3_pred[T_new-P];
	vector[M] Kappa_3_new[T_new];
	vector[M] Kappa_4_pred[T_new-P];
	vector[M] Kappa_4_new[T_new];
	vector[M] Kappa_5_pred[T_new-P];
	vector[M] Kappa_5_new[T_new];
	vector[M] Kappa_bar_new[T_new];

	Sigma = L_sigma * L_sigma';

		
	Kappa_1_new[1] = Kappa_1[T-P];
	Kappa_2_new[1] = Kappa_2[T-P];
	Kappa_3_new[1] = Kappa_3[T-P];
	Kappa_4_new[1] = Kappa_4[T-P];
	Kappa_5_new[1] = Kappa_5[T-P];
	Kappa_bar_new[1] = Kappa_bar[T-P];

	for (t in 1:(T_new-P)) 
		{
		
		Kappa_1_new[t+P] = intercept + Kappa_1_new[t] - Psi * (Kappa_1_new[t]-Kappa_bar_new[t]);
		Kappa_2_new[t+P] = intercept + Kappa_2_new[t] - Psi * (Kappa_2_new[t]-Kappa_bar_new[t]);
		Kappa_3_new[t+P] = intercept + Kappa_3_new[t] - Psi * (Kappa_3_new[t]-Kappa_bar_new[t]);
		Kappa_4_new[t+P] = intercept + Kappa_4_new[t] - Psi * (Kappa_4_new[t]-Kappa_bar_new[t]);
		Kappa_5_new[t+P] = intercept + Kappa_5_new[t] - Psi * (Kappa_5_new[t]-Kappa_bar_new[t]);
		Kappa_bar_new[t+P] = (Kappa_1_new[t+P]+Kappa_2_new[t+P]+Kappa_3_new[t+P]+Kappa_4_new[t+P]+Kappa_5_new[t+P])/5;
		}

	
	Kappa_1_pred = multi_normal_cholesky_rng(Kappa_1_new[2:T_new],L_sigma);
	Kappa_2_pred = multi_normal_cholesky_rng(Kappa_2_new[2:T_new],L_sigma);
	Kappa_3_pred = multi_normal_cholesky_rng(Kappa_3_new[2:T_new],L_sigma);
	Kappa_4_pred = multi_normal_cholesky_rng(Kappa_4_new[2:T_new],L_sigma);
	Kappa_5_pred = multi_normal_cholesky_rng(Kappa_5_new[2:T_new],L_sigma);
		
}

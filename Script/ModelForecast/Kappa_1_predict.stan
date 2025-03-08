data {
  int<lower=1> T;
  vector[T] Kappa_1;
  int<lower=1> T_new;
}
transformed data {
  vector[T-1] Kappa_1_diff = Kappa_1[2:T] - Kappa_1[1:(T-1)];


}
parameters {
  real<lower=0> sigma;
  real mu;
}
model {
  //mu ~ normal(0, 10);
  //sigma ~ cauchy(0, 1);
  //sigma ~ inv_gamma(11,0.1);
  Kappa_1_diff ~ normal(mu, sigma);
} 

generated quantities 
	{
	
		
	vector[T_new-1] Kappa_1_pred;
	vector[T_new] Kappa_1_new;
			
	Kappa_1_new[1] = Kappa_1[T];
	
	for (t in 1:(T_new-1)) 
		{
		
		Kappa_1_new[t+1] = mu + Kappa_1_new[t];
		}

	for (t in 1:(T_new-1))
		{	
		Kappa_1_pred[t] = normal_rng(Kappa_1_new[t+1],sigma);
		}
}

data {
	int<lower = 1> J; // number of age categories
	int<lower = 1> T; // number of years
	//int d[J*T]; // vector of deaths
	vector[J*T] d; // vector of deaths
	vector[J* T] e; // vector of exposures
	vector[J] age; // vector of ages
	int<lower = 1> Tfor; // number of forecast years
	vector[J] mort; // matrix of mortality to bx
	//int<lower = 0> Tval; // number of validation years
	//int dval[J*Tval]; // vector of deaths for validation
	//vector[J* Tval] eval; // vector of exposures for validation
	int<lower=0,upper=1> family; // family = 0 for Poisson, 1 for NB
     }

transformed data 
     {
	vector[J * T] offset = log(e);
	vector[J * T] log_d = log(d);
	//vector[J * Tval] offset2 = log(eval);
	int<lower = 1> L;
	L=J*Tfor;
     }

parameters 
     {
	real<lower=0> aux[family > 0]; // neg. binomial dispersion parameter
	real c1; // drift term for kappa_1
	real c2; // drift term for kappa_2
	real<lower = 0> sigma[2]; // standard deviations for kappa_1 and kappa_2
	vector[T] k; // vector of kappa_1
	vector[T] k2; // vector of kappa_2
	vector[J] bx; // reparameterized std. normal.
	real<lower=-1,upper=1> rho; // correlation parameter
	real<lower=0> tau; // uncertainty of bx. Not varies by age.
     }

transformed parameters 
     { // No identifiability contraints in the CBD model
	//vector[J] bx; // vector of bx
	
	real phi = negative_infinity();

	if (family > 0) phi = inv(aux[1]);
	
	
	//for (j in 1:J)
	//{
	//bx[j] = mort[j]+tau[j]*bx_raw[j];	// implies bx ~ normal(0,0.01)
	//}
	//bx = mort+tau .* bx_raw;	// implies bx ~ normal(0,0.01)
     }

model 
     {
	vector[J * T] mu;
	int pos = 1;
	int pos2 = 1;
		for (t in 1:T) for (x in 1:J) 
			{
			mu[pos] = offset[pos]+bx[x]+k[t]+(age[x]-mean(age))*k2[t]; //Predictor dynamics
			//mu[pos] = e[pos]*(bx[x]+k[t]+(age[x]-mean(age))*k2[t]); //Predictor dynamics
			pos += 1;
			}
	//bx_raw ~ std_normal();
	//tau ~ exponential(5);
	target += inv_gamma_lpdf(tau | 11,0.01); // prior on tau - uncertainty of bx
	target += normal_lpdf(bx|mort,tau);  // prior on bx
	target += normal_lpdf(k[1]|c1,sqrt(10));
	target += normal_lpdf(k2[1]|c2,sqrt(10));
	for	(t in 2:T)
		{
		target += normal_lpdf(k[t] | c1+k[t-1], sigma[1]);
		target += normal_lpdf(k2[t] | c2+k2[t-1] + rho * sigma[2]*(k[t]-c1-k[t-1])/ sigma[1],sigma[2] * sqrt(1 - square(rho)));
		}
	//target += normal_lpdf(k[2:T] | c1+k[1:(T- 1)], sigma[1]);
	//target += normal_lpdf(k2[2:T] | c2+k2[1:(T- 1)] + rho * sigma[2]*(k[2:T]-c1-k[1:(T- 1)])/ sigma[1],sigma[2] * sqrt(1 - square(rho)));
	
		if (family ==0)
			{
			//target += poisson_log_lpmf(d |mu); // Poisson log model
			for (t in 1:T) for (x in 1:J)
				{
				target += normal_lpdf(log_d[pos2]|mu[pos2],sqrt(1/d[pos2])); // Normal log model
				pos2 += 1;
				}
			}
		//else 	
		//	{
		//	target +=neg_binomial_2_log_lpmf (d|mu,phi); // Negative-Binomial log model
		//	}

	//target += exponential_lpdf(sigma | 0.1);
	target += inv_gamma_lpdf(sigma | 11,0.1); // hyper-prior on sigma
	target += normal_lpdf(c1|0,sqrt(10));
	target += normal_lpdf(c2|0,sqrt(10));
	target += uniform_lpdf(rho|-1,1);
	//target += beta_lpdf(rho|2,2);       // hyper-prior on correlation between the two processes
		if (family > 0) target += normal_lpdf(aux|0,1)- normal_lcdf(0 | 0, 1); // prior on overdispersion parameter

   }

generated quantities 
   {
	vector[Tfor] k_p;
	vector[Tfor] k2_p;
	vector[L] mufor;
	//vector[J*T] log_lik;
	//vector[J*Tval] log_lik2;
	int pos = 1;
	int pos2= 1;
	int pos3= 1;
	k_p[1] = c1+k[T]+sigma[1] * normal_rng(0,1);
		for (t in 2:Tfor) k_p[t] = c1+k_p[t - 1] + sigma[1] * normal_rng(0,1);
	k2_p[1] = c2+k2[T]+ rho*sigma[2]/sigma[1]*(k_p[1]-c1-k[T])+sigma[2] * sqrt(1 - square(rho))*normal_rng(0,1);
		for (t in 2:Tfor) k2_p[t] = c2+k2_p[t - 1]+ rho*sigma[2]/sigma[1]*(k_p[t]-c1-k_p[t - 1])+sigma[2] * sqrt(1 - square(rho))*normal_rng(0,1);
	if (family==0)
		{
		for (t in 1:Tfor) for (x in 1:J) 
				{
				mufor[pos] = k_p[t]+(age[x]-mean(age))*k2_p[t];
				pos += 1;
				}
				mufor=exp(mufor);
		//for (t in 1:T) for (x in 1:J) 
		//		{
		//		log_lik[pos2] = poisson_log_lpmf (d[pos2] | offset[pos2]+ k[t]+(age[x]-mean(age))*k2[t]);
		//		pos2 += 1;
		//		}
		//for (t in 1:Tval) for (x in 1:J) 
		//		{
		//		log_lik2[pos3] = poisson_log_lpmf(dval[pos3] | offset2[pos3]+ k_p[t]+(age[x]-mean(age))*k2_p[t]);
		//		pos3 += 1;
		//		}
		}

	else if (family > 0)
		{
		for (t in 1:Tfor) for (x in 1:J) 
				{
				if ( fabs(k_p[t]+(age[x]-mean(age))*k2_p[t])>15){
				mufor[pos] = 0;
				pos += 1;
				} else 
				{
				mufor[pos] = gamma_rng(phi,phi/exp(k_p[t]+(age[x]-mean(age))*k2_p[t]));
				pos += 1;
				}
		}
		//for (t in 1:T) for (x in 1:J) 
		//		{
		//		log_lik[pos2] = neg_binomial_2_log_lpmf (d[pos2] | offset[pos2]+ k[t]+(age[x]-mean(age))*k2[t],phi);
		//		pos2 += 1;
		//		}
		//for (t in 1:Tval) for (x in 1:J) 
		//		{
		//		log_lik2[pos3] = neg_binomial_2_log_lpmf (dval[pos3] | offset2[pos3]+ k_p[t]+(age[x]-mean(age))*k2_p[t],phi);
		//		pos3 += 1;
		//		}
   		}
   }

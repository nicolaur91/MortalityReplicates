data {
	int<lower = 1> J; // number of age categories
	int<lower = 1> T; // number of years
	int<lower = 1> H; // number of subgroups
	//int d[J*T*H]; // vector of deaths
	vector[J*T*H] d; // vector of deaths
	vector[J*T*H] e; // vector of exposures
	
	vector[J] mort; // matrix of mortality to bx
	vector[J] age; // vector of ages
	//vector[H] ones;
	int<lower = 1> Tfor; // number of forecast years
	//int<lower = 0> Tval; // number of validation years
	//int dval[J*Tval]; // vector of deaths for validation
	//vector[J* Tval] eval; // vector of exposures for validation
	int<lower=0,upper=1> family; // family = 0 for Poisson, 1 for NB
	
     }

transformed data 
     {
	vector[J * T * H] offset = log(e);
	vector[J * T * H] log_d = log(d);
	//vector[J * Tval] offset2 = log(eval);
	int<lower = 1> L;
	L=J*Tfor*H;

     }

parameters 
     {
	real c1; // drift term for kappa_1
	real c2; // drift term for kappa_2
	//vector[2] c; // drift term for kappa(1) and kappa(2)
	real<lower = 0> sigma[2]; // standard deviations for kappa_1 and kappa_2
	//vector<lower=0>[2] sigma; // standard deviations for kappa_1 and kappa_2
	// vector[T] k; vector of kappa_1
	
	//matrix[2,H] k[T]; // matrix of kappa(1) and kappa(2)
	
	//vector<lower=-8,upper=0>[J] bx; // vector of bx
	real<lower=-1,upper=1> rho; // correlation parameter
	//cholesky_factor_corr[2] rho; // correlation parameter
	real<lower=0,upper=1> psy; // gravity parameter
	//matrix[2,T] k_bar; // kappa_bar
	vector[T] kappa_1_tilde; // governs std. normal parameterization of kappa_1
	vector[T] kappa_2_tilde; // governs std. normal parameterization of kappa_2
	//vector[J * T * H] deaths_tilde;
     }

transformed parameters 
     { // No identifiability contraints in the CBD model

	vector[T] k; // vector of kappa_1
	vector[T] k2; // vector of kappa_2
	//vector[T] k_bar; // average of kappa_1
	//vector[T] k2_bar; // average of kappa_2
		
	
	//k[,1] = rep_vector(c1+sigma[1]*kappa_1_tilde[1],H);
	k[1] = c1+sigma[1]*kappa_1_tilde[1];
	k2[1] = c2+sigma[2]*kappa_2_tilde[1]*sqrt(1-square(rho))+rho*sigma[2]*(k[1]-c1)/sigma[1];
	
	//k_bar[1] = sum(k[,1])/H;
	//k2_bar[1] = sum(k2[,1])/H;
	
	for (t in 2:T)
		{
		//k[,t] = k[,t-1]-psy*(k[,t-1]-k_bar[t-1]);
		//k2[,t] = k2[,t-1]-psy*(k2[,t-1]-k2_bar[t-1]);
		k[t] = c1+k[t-1]+sigma[1]*kappa_1_tilde[t];
		k2[t] = c2+k2[t-1]+sigma[2]*kappa_2_tilde[t]*sqrt(1-square(rho))+rho*sigma[2]*(k[t]-c1-k[t-1])/sigma[1];
		//k_bar[t] = c1+k_bar[t-1] + sigma[1]*kappa_1_tilde[t];
		//k2_bar[t] = c2+k2_bar[t-1] + sigma[2]*sqrt(1-square(rho))*kappa_1_tilde[t]+rho*sigma[2]*(k_bar[t]-(c1+k_bar[t-1]))/sigma[1];
		//k_bar[t] = sum(k[,t])/H;
		//k2_bar[t] = sum(k2[,t])/H;

		}

	
	//log_d = mu + deaths_tilde/d;

	
	
     }

model 
     {
	//vector[T] k_bar; // kappa_bar
	//vector[T] k2_bar; // kappa_bar
	
	vector[J * T * H] mu;
	
	//k_bar[1] = 0;
	//k2_bar[1] = 0;
	//k_bar[1] = 0;

	//matrix[T,H] k2_bar;
	//k2_bar[2] = 0;
	int pos = 1;
	int pos2 = 1;
	//for (h in 1:H)
	//{
	
		for (t in 1:T) for (x in 1:J) 
			{
			//mu[pos] = e[pos]*(bx[x]+k[t]+(age[x]-mean(age))*k2[t]); //Predictor dynamics - likelihood
			mu[pos] = e[pos]*(k[t]+(age[x]-mean(age))*k2[t]); //Predictor dynamics - likelihood
			//mu[pos] = offset[pos]+(bx[x]+k[t]+(age[x]-mean(age))*k2[t]); //Predictor dynamics - likelihood
			pos += 1;
			}

//	}
	
	//print(k_bar);

	// standard normal increment
	
	//for (h in 1:H)
	//{
	target += std_normal_lpdf(kappa_1_tilde);
	target += std_normal_lpdf(kappa_2_tilde);
	//}
	//target += normal_lpdf(k_bar[1]|0,0.02); // initial prior
	//target += normal_lpdf(k2_bar[1]|0,0.005); // initial prior
	

	//for (t in 2:T)
	//	{
	//	target += normal_lpdf(k_bar[t]|k_bar[t-1]+c1,sigma[1]);
	//	target += normal_lpdf(k2_bar[t]|k2_bar[t-1]+c2+rho*sigma[2]*(k_bar[t]-k_bar[t-1]-c1)/sigma[1],sigma[2] * sqrt(1 - square(rho)));
                //k2_bar[t] ~ normal(k2_bar[t-1]+c2+rho*sigma[2]*(k_bar[t]-k_bar[t-1]-c1)/sigma[1],sigma[2] * sqrt(1 - square(rho));           		
	//	}
		
	//for (h in 1:H)
	//	{
	//	target += normal_lpdf(bx|mort,0.005); // prior on bx
	//	}
	target += inv_gamma_lpdf(sigma | 11,0.1); // hyper-prior on sigma
	target += beta_lpdf(rho|2,2);       // hyper-prior on correlation between the two processes
	target += beta_lpdf(psy|2,2);       // hyper-prior on psy
	
	//target += normal_lpdf(log_d|mu,1/d); // Normal log model
	
	
	//for (h in 1:H)
	//{
	
		for (t in 1:T) for (x in 1:J) 
			{
			target += normal_lpdf(log_d[pos2]|mu[pos2],1/d[pos2]); // Normal log model
			pos2 += 1;
			
			}
	//}
}

//generated quantities 
  // {
	//vector[Tfor] k_p;
	//matrix[Tfor,H] k_p;
	//vector[Tfor] k2_p;
	//matrix[Tfor,H] k2_p;
	//vector[Tfor] k_bar_p;
	//vector[Tfor] k2_bar_p;
	//vector[L] mufor;
	//vector[J*T] log_lik;
	//vector[J*Tval] log_lik2;
	//int pos = 1;
	
	//int pos2= 1;
	//int pos3= 1;
	//k_bar_p[1] = sum(k[T,1:H])/H;
	//k2_bar_p[1] = sum(k2[T,1:H])/H;
	//k_p[1,] = c1+k[T,]-psy*(k[T,]-k_bar_p[1])+sigma[1] * normal_rng(0,1);
	//k2_p[1,] = c2+k2[T,]-psy*(k2[T,]-k2_bar_p[1])+ rho*sigma[2]/sigma[1]*(k_p[1]-c1-k[T]+psy*(k[T,]-sum(k[T,])/H))+sigma[2] * sqrt(1 - square(rho))*normal_rng(0,1);
	//	for (t in 2:Tfor) 
	//		{ 
	//			k_bar_p[t-1] = sum(k_p[t-1,])/H;
	//			k_p[t,] = c1+k_p[t - 1,] -psy*(k[t-1,]-k_bar_p[t-1]) + sigma[1] * normal_rng(0,1);
	//			k2_bar_p[t-1] = sum(k2_p[t-1,])/H;
	//			k2_p[t] = c2+k2_p[t - 1,]-psy*(k2_p[t-1,]-k2_bar_p[t-1])+ rho*sigma[2]/sigma[1]*(k_p[t,]-c1-k_p[t - 1,]+psy*(k_p[t-1,]-k_bar_p[t-1]))+sigma[2] * sqrt(1 - square(rho))*normal_rng(0,1);
				
	//		}
	
	//	for (t in 2:Tfor) 
	//		{	
				
				
	//		}
	//	for (t in 1:T) for (x in 1:J) 
	//			{
	//			log_lik[pos2] = poisson_log_lpmf (d[pos2] | offset[pos2]+ k[t]+(age[x]-mean(age))*k2[t]);
	//			pos2 += 1;
	//			}
	//	for (t in 1:Tval) for (x in 1:J) 
	//			{
	//			log_lik2[pos3] = poisson_log_lpmf(dval[pos3] | offset2[pos3]+ k_p[t]+(age[x]-mean(age))*k2_p[t]);
	//			pos3 += 1;
	//			}
		
	
	
	//else if (family > 0)
	//	{
	//	for (t in 1:Tfor) for (x in 1:J) 
	//			{
	//			if ( fabs(k_p[t]+(age[x]-mean(age))*k2_p[t])>15){
	//			mufor[pos] = 0;
	//			pos += 1;
	//			} else 
	//			{
	//			mufor[pos] = gamma_rng(phi,phi/exp(k_p[t]+(age[x]-mean(age))*k2_p[t]));
	//			pos += 1;
	//			}
	//	}
	//	for (t in 1:T) for (x in 1:J) 
	//			{
	//			log_lik[pos2] = neg_binomial_2_log_lpmf (d[pos2] | offset[pos2]+ k[t]+(age[x]-mean(age))*k2[t],phi);
	//			pos2 += 1;
	//			}
	//	for (t in 1:Tval) for (x in 1:J) 
	//			{
	//			log_lik2[pos3] = neg_binomial_2_log_lpmf (dval[pos3] | offset2[pos3]+ k_p[t]+(age[x]-mean(age))*k2_p[t],phi);
	//			pos3 += 1;
	//			}
	
	//if (family ==0)
	//	{
	//	for (h in 1:H) for (t in 1:Tfor) for (x in 1:J) 
	//			{
	//			mufor[pos] = k_p[t,h]+(age[x]-mean(age))*k2_p[t,h];
	//			pos += 1;
	//			}
	//			mufor=exp(mufor);
	//	}
		
	
 //  }
 //  }

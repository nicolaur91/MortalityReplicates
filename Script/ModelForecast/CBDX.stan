data {
	int<lower = 1> J; // number of age categories
	int<lower = 1> T; // number of years
	int<lower = 1> H; // number of subgroups
	//int d[J*T*H]; // vector of deaths
	vector[J*T*H] d; // vector of deaths
	vector[J*T*H] e; // vector of exposures
	
	matrix[H,J] mort; // matrix of mortality to bx
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
	//real<lower=0> aux[family > 0]; // neg. binomial dispersion parameter
	real c1; // drift term for kappa_1
	real c2; // drift term for kappa_2
	//vector[2] c; // drift term for kappa(1) and kappa(2)
	real<lower = 0> sigma[2]; // standard deviations for kappa_1 and kappa_2
	//vector<lower=0>[2] sigma; // standard deviations for kappa_1 and kappa_2
	// vector[T] k; vector of kappa_1
	matrix[H,T] k; // vector of kappa_1
	matrix[H,T] k2; // vector of kappa_2
	//matrix[2,H] k[T]; // matrix of kappa(1) and kappa(2)
	//vector[T] k_bar;
	//vector[T] k2_bar;
	matrix<lower=-8,upper=0>[H,J] bx; // vector of bx
	real<lower=-1,upper=1> rho; // correlation parameter
	//cholesky_factor_corr[2] rho; // correlation parameter
	real<lower=0,upper=1> psy; // gravity parameter
	//matrix[2,T] k_bar; // kappa_bar
	//matrix[1,H] delta_1_1; // initial delta for kappa 1
	//matrix[1,H] delta_1_2; // initial delta for kappa 2
	//vector[T] kappa_1_tilde; // governs std. normal parameterization of kappa_1
	//vector[T] kappa_2_tilde; // governs std. normal parameterization of kappa_2
     }

transformed parameters 
     { // No identifiability contraints in the CBD model
	//real phi = negative_infinity();
	
	//if (family > 0) phi = inv(aux[1]);
	//matrix[2,2] L_sigma;
  
	//L_sigma = diag_pre_multiply(sigma, rho);


//	matrix[H,T] delta_t_1;
//	matrix[H,T] delta_t_2;
//	matrix[H,T] k; // vector of kappa_1
//	matrix[H,T] k2; // vector of kappa_2

//	k[1,] = delta_1_1[1,]; // initial value for kappa_1
//	k2[1,] = delta_1_2[1,]; // initial value for kappa_2
	
//	delta_t_1[1,] = delta_1_1[1,]*(1-psy); // delta t=2 for kappa_1
//	delta_t_2[1,] = delta_1_2[1,]*(1-psy); // delta t=2 for kappa_2

//	for (t in 2:T)
//		{
//		k[t,] = delta_t_1[t-1,]+c1+sigma[1]*kappa_1_tilde[t];
//		k2[t,] = delta_t_2[t-1,]+c2+sigma[2]*kappa_1_tilde*sqrt(1-square(rho))+rho*sigma[2]/sigma[1]*(k[t,]-delta_t_1[t-1,]-c1);
//		delta_t_1[t,] = delta_t_1[t-1,]*(1-psy); 
//		delta_t_2[t,] = delta_t_2[t-1,]*(1-psy);
//		}
	
     }

model 
     {
	vector[T] k_bar; // kappa_bar
	vector[T] k2_bar; // kappa_bar
	
	vector[J * T * H] mu;
	//k_bar[1] = 0;
	//k2_bar[1] = 0;
	//k_bar[1] = 0;

	//matrix[T,H] k2_bar;
	//k2_bar[2] = 0;
	int pos = 1;
	int pos2 = 1;
	for (h in 1:H)
	{
	
		for (t in 1:T) for (x in 1:J) 
			{
			mu[pos] = offset[pos]+bx[h,x]+k[h,t]+(age[x]-mean(age))*k2[h,t]; //Predictor dynamics - likelihood
			//mu[pos] = e[pos]*(bx[x,h]+k[t,h]+(age[x]-mean(age))*k2[t,h]); //Predictor dynamics - likelihood
			//mu[pos] = e[pos]*(bx[x,h]+k[1,h][t]+(age[x]-mean(age))*k[2,h][t]); //Predictor dynamics - likelihood
			//k_bar[t] = sum(k[t,1:H])/H;
			//k2_bar[t] = sum(k2[t,1:H])/H;
			pos += 1;
			}

	//k_bar = k * ones/H;
	//k2_bar = k2 * ones/H;

	//print(k_bar);
	//print(k2_bar);
	//print(bx);
	//print(mu);
	}
//	target += normal_lpdf(k[1,]|c1,sigma[1]); // initial prior
//	target += normal_lpdf(k2[1,]|c2+rho*sigma[2]*(k[1,]-c1)/sigma[1],sigma[2] * sqrt(1 - square(rho))); // initial prior
	target += normal_lpdf(k[,1]|c1,sqrt(10)); // initial prior
	target += normal_lpdf(k2[,1]|c2,sqrt(10)); // initial prior

//	target += normal_lpdf(delta_1_1[1,]|c1,sqrt(10)); // initial prior
//	target += normal_lpdf(delta_1_2[1,]|c2,sqrt(10)); // initial prior

//	target += std_normal_lpdf(kappa_1_tilde);
//	target += std_normal_lpdf(kappa_2_tilde);

	//target += multi_normal_cholesky_lpdf(k_bar[,1]|c,L_sigma); // initial prior
	//for (h in 1:H)
	//	{
	//	k[,h][1] = c; 
	//	}
	for (t in 2:T)
		{
		//target += normal_lpdf(bx[,h]|mort[1:J]),sqrt(1)); // priors
		  k_bar[t-1] = sum(k[,t-1])/H;
		  k2_bar[t-1] = sum(k2[,t-1])/H;
		//k_bar[,t-1] = sum(k[,][t-1])/H;
		//target += normal_lpdf(k[2:T,h] | c1+k[1:(T- 1),h]-psy*(k[1:(T- 1),h]-c1*k_bar[1:(T- 1)]), sigma[1]); // priors
		  target += normal_lpdf(k[,t] | c1+k[,t-1]-psy*(k[,t-1]-k_bar[t-1]), sigma[1]); // priors
		//target += normal_lpdf(k2[2:T,h] | c2+k2[1:(T- 1),h]-psy*(k2[1:(T- 1),h]-c2*k2_bar[1:(T- 1)]) + rho * sigma[2]*(k[2:T,h]-c1-k[1:(T- 1),h]+psy*(k[1:(T- 1),h]-c1*k_bar[1:(T- 1)]))/ sigma[1],sigma[2] * sqrt(1 - square(rho))); // priors
		  target += normal_lpdf(k2[,t] | c2+k2[,t-1]-psy*(k2[,t-1]-k2_bar[t-1]) + rho*sigma[2]*(k[,t]-c1-k[,t-1]+psy*(k[,t-1]-k_bar[t-1]))/ sigma[1],sigma[2] * sqrt(1 - square(rho))); // priors
		//target += multi_normal_cholesky_lpdf(k_bar[,t]|k_bar[,t-1]+c,L_sigma); // slicing kappa-bar
		//k[,][t] = c+k[,][t-1]-psy*(k[,][t-1]-k_bar[,t-1]); // priors
		}
		
	//target += uniform_lpdf(bx[,h]|-7,0); // prior on bx
	for (h in 1:H)
		{
		target += normal_lpdf(bx[h,]|mort[h,],0.005); // prior on bx
		}

	// preferred hyper-prior specification
	//target += exponential_lpdf(sigma | 0.1); // hyper-priors
	//target += exponential_lpdf(sigma | 0.1); // hyper-priors	
	target += uniform_lpdf(rho|-1,1);       // hyper-priors
	target += uniform_lpdf(psy|0,1);       // hyper-priors
	
	//target += normal_lpdf(c1|0,sqrt(10));   // hyper-prior on constant
	//target += normal_lpdf(c2|0,sqrt(10));   // hyper-prior on constant

	// alternative hyper-prior specification
	
	//
	target += inv_gamma_lpdf(sigma | 11,0.1); // hyper-prior on sigma
	//target += normal_lpdf(c2|0,sqrt(2));   // hyper-priors on constant
	//target += normal_lpdf(c|0,sqrt(2));   // hyper-prior on constant
	//target += beta_lpdf(rho|2,2);       // hyper-prior on correlation between the two processes
	//target += lkj_corr_cholesky_lpdf(rho|2.0);       // hyper-prior on correlation between the two processes
	//target += beta_lpdf(psy|2,2);       // hyper-prior on psy
	//	if (family > 0) target += normal_lpdf(aux|0,1)- normal_lcdf(0 | 0, 1); // prior on overdispersion parameter
	
	
	//target +=normal_lpdf(k_bar[1]|c1,0.0001); // prior for k_bar
	//target +=normal_lpdf(k2_bar[1]|c1,0.0001); // prior for k2_bar
	//target += normal_lpdf(k_bar[2:T] | c1+k_bar[1:(T- 1)], sigma[1]); // priors
	//target += normal_lpdf(k2_bar[2:T] | c2+k2_bar[1:(T- 1)]+ rho * sigma[2]*(k_bar[2:T]-c1-k_bar[1:(T- 1)])/ sigma[1],sigma[2] * sqrt(1 - square(rho))); // priors
	
	//if (family ==0)
	//		{
	//target += poisson_log_lpmf(d|mu); // Poisson log model

	for (h in 1:H)
	{
	
		for (t in 1:T) for (x in 1:J) 
			{
			target += normal_lpdf(log_d[pos2]|mu[pos2],1/d[pos2]); // Normal log model
			pos2 += 1;
			
			}
	}
	//		}
	//	else 	
	//		{
	//		target +=neg_binomial_2_log_lpmf (d|mu,phi); // Negative-Binomial log model
	//		}
   }

generated quantities 
   {
	//vector[Tfor] k_p;
	matrix[H,Tfor] k_p;
	//vector[Tfor] k2_p;
	matrix[H,Tfor] k2_p;
	vector[Tfor] k_bar_p;
	vector[Tfor] k2_bar_p;
	vector[L] mufor;
	//vector[J*T] log_lik;
	//vector[J*Tval] log_lik2;
	int pos = 1;
	
	//int pos2= 1;
	//int pos3= 1;
	k_bar_p[1] = sum(k[,T])/H;
	k2_bar_p[1] = sum(k2[,T])/H;
	k_p[,1] = c1+k[,T]-psy*(k[,T]-k_bar_p[1])+sigma[1] * normal_rng(0,1);
	k2_p[,1] = c2+k2[,T]-psy*(k2[,T]-k2_bar_p[1])+ rho*sigma[2]/sigma[1]*(k_p[,1]-c1-k[,T]+psy*(k[,T]-k_bar_p[1]))+sigma[2] * sqrt(1 - square(rho))*normal_rng(0,1);
		for (t in 2:Tfor) 
			{ 
				k_bar_p[t-1] = sum(k_p[,t-1])/H;
				k_p[,t] = c1+k_p[,t - 1] -psy*(k[,t-1]-k_bar_p[t-1]) + sigma[1] * normal_rng(0,1);
				k2_bar_p[t-1] = sum(k2_p[,t-1])/H;
				k2_p[,t] = c2+k2_p[,t - 1]-psy*(k2_p[,t-1]-k2_bar_p[t-1])+ rho*sigma[2]/sigma[1]*(k_p[,t]-c1-k_p[,t - 1]+psy*(k_p[,t-1]-k_bar_p[t-1]))+sigma[2] * sqrt(1 - square(rho))*normal_rng(0,1);
				
			}
	
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
		for (h in 1:H) for (t in 1:Tfor) for (x in 1:J) 
				{
				mufor[pos] = bx[h,x]+k_p[h,t]+(age[x]-mean(age))*k2_p[h,t];
				pos += 1;
				}
				mufor=exp(mufor);
	//	}
		
	
   }
 //  }

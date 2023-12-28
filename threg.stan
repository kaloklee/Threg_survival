//Inverse Gaussian Distribution as formulated in R package 'threg'
//The initial health status, y0, should be positive.  We model is as log(y0) = X*beta and then exponentiate
//This example has a single covariate with 2 levels.

functions {

  real threg_IG_ll (real t, real mu, real y0, real sigma2){

    return ( log(y0) - .5*( log(2*pi()) + sigma2 + 3*log(t) )  - (y0+mu*t)^2/(2*sigma2*t) )  ;

  }


  real threg_IG_surv_ll (real t, real mu, real y0 , real sigma2){

      real A = -(y0+mu*t)/sqrt(sigma2*t);
      real B = (mu*t-y0)/sqrt(sigma2*t);

      return ( log1m( normal_cdf(A | 0, 1) + exp(-2*y0*mu/sigma2)*normal_cdf(B | 0, 1) ) );
 
  }

}
  
data{
  
  int N;
  real y[N];
  int Nuc;
  int Nc;
  int T; 
  int K; //number of covariates (with intercept)
  matrix[N,K] Xvar; //covariates

}

parameters {
  
  vector[K] b_lny0; // initial heath status
  vector[K] b_mu; // mean

}


model {

  vector[N] lny0_Xb = Xvar*b_lny0 ; // N*K * K*1 = N*1
  vector[N] mu_Xb = Xvar*b_mu ; // N*K * K*1 = N*1

  b_lny0 ~ std_normal();
  b_mu ~ std_normal();

  for (i in 1:N) {
    
    if (i <= Nuc) { target += threg_IG_ll(y[i],mu_Xb[i],exp(lny0_Xb[i]),1) ; }

    else { target += threg_IG_surv_ll (y[i],mu_Xb[i],exp(lny0_Xb[i]),1); }
    
  }
    
  
}


//This example has a single covariate with 2 levels.
//We only need to generate two lines.
generated quantities{

    row_vector[T+1] Survival0=rep_row_vector(0,T+1);
    row_vector[T+1] Survival1=rep_row_vector(0,T+1);

    Survival0[1] = 1;
    Survival1[1] = 1;
    
    for (t in 2:T+1) {    
      
        Survival0[t] = exp(threg_IG_surv_ll(t,b_mu[1],exp(b_lny0[1]),1));  
        Survival1[t] = exp(threg_IG_surv_ll(t,sum(b_mu),exp(sum(b_lny0)),1));
      
    }

}

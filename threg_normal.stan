//Inverse Gaussian Distribution as formulated in R package 'threg'
//Allowing random effect on the drift parameter

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
  
  
  real a_mu;
  real<lower=0> a_sig;
  vector[N] eta;
  vector[K-1] b_mu; // drift
  vector[K] b_lny0; // initial heath status

}


transformed parameters {
  
  vector[N] a = a_mu + eta*a_sig;
  
}

model {

  a_mu ~ std_normal();
  a_sig ~ cauchy(0, 1); //a_sig ~ std_normal();
  eta ~ std_normal();
  b_lny0 ~ std_normal();
  b_mu ~ std_normal();

  for (i in 1:N) {
    
    real lny0_Xb = Xvar[i,]*b_lny0 ;
    real mu_Xb = a[i] + Xvar[i,2:K]*b_mu ;  // 1x1 + 1x(K-1) * (K-1)x1 = 1x1

    if (i <= Nuc) { target += threg_IG_ll(y[i],mu_Xb,exp(lny0_Xb),1) ; }

    else { target += threg_IG_surv_ll (y[i],mu_Xb,exp(lny0_Xb),1); }
    
  }
    
  
}


generated quantities{

    row_vector[T+1] Survival0=rep_row_vector(0,T+1);
    row_vector[T+1] Survival1=rep_row_vector(0,T+1);

    Survival0[1] = 1;
    Survival1[1] = 1;
    
    for (t in 2:T+1) {    
      
      for (i in 1:N) {
        
        real lny0_Xb = Xvar[i,]*b_lny0 ;
        real mu_Xb = a[i] + Xvar[i,2:K]*b_mu ;  // 1x1 + 1x(K-1) * (K-1)x1 = 1x1
        
        Survival0[t] += exp(threg_IG_surv_ll(t,a[i],exp(b_lny0[1]),1));  
        Survival1[t] += exp(threg_IG_surv_ll(t,mu_Xb,exp(lny0_Xb),1));
      
      }
    
    }
    
    Survival1[2:T+1]=Survival1[2:T+1]/N;
    Survival0[2:T+1]=Survival0[2:T+1]/N;

}



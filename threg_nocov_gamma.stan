//Inverse Gaussian Distribution as formulated in R package 'threg'
//No covariate in this version
//Allow y0 to be Gamma distributed

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

}

parameters {
  
  real<lower=0> r; 
  real<lower=0> alpha;
  real mu; // mean
  vector<lower=0>[N] y0;

}




model {

  mu ~ std_normal(); 
  r ~ std_normal();
  alpha ~ std_normal();
  
  y0 ~ gamma(r, alpha); //normal(0,a)
 

  for (i in 1:N) {
    
    if (i <= Nuc) { target += threg_IG_ll(y[i],mu,y0[i],1) ; }
  

    else { target += threg_IG_surv_ll (y[i],mu,y0[i],1); }
    
  }
  
  
}


 generated quantities{
 
    matrix[N,T+1] expected;
    row_vector[T+1] Survival;
   
    expected[,1] = rep_vector(1,N);
    
    for (i in 1:N) {
      
      real y0_sample = gamma_rng(r, alpha);
      
      for (t in 2:T+1) {
 
        expected[i,t]=exp(threg_IG_surv_ll(i,mu,y0_sample,1));
      }
      
    }
   
   Survival=rep_row_vector(1,N)*expected;
}


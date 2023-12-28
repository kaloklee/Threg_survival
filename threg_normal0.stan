//Inverse Gaussian Distribution as formulated in R package 'threg'
//No covariate in this version

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
  vector[K] b_mu; // mean

}


transformed parameters {
  
  vector[N] a = a_mu + eta*a_sig;
  
}

model {

  a_mu ~ std_normal();
  a_sig ~ std_normal();
  eta ~ std_normal();
  b_mu ~ std_normal();

  for (i in 1:N) {
    
  real lny0_Xb = a[i] ; //  
  real mu_Xb = Xvar[i,]*b_mu ;  // 1*K * K*1 = 1*1

    if (i <= Nuc) { target += threg_IG_ll(y[i],mu_Xb,exp(lny0_Xb),1) ; }

    else { target += threg_IG_surv_ll (y[i],mu_Xb,exp(lny0_Xb),1); }
    
  }
    
  
}


// generated quantities{
// 
//    row_vector[T+1] Survival0=rep_row_vector(0,T+1);
//    row_vector[T+1] Survival1=rep_row_vector(0,T+1);
//  
//    Survival0[1]=1;
//    Survival1[1]=1;
//  
//    for (t in 2:T+1) {
//      
//       for (i in 1:N) {
// 
//       real lny0_Xb = a[i] ; //  
//       real mu_Xb = Xvar[i,]*b_mu ;  // 1*K * K*1 = 1*1
// 
//       if (Xvar[i,K] == 0) { Survival0[t] += exp(threg_IG_surv_ll(t,mu_Xb,exp(lny0_Xb),1)); }
//       
//       if (Xvar[i,K] == 1) { Survival1[t] += exp(threg_IG_surv_ll(t,mu_Xb,exp(lny0_Xb),1)); }
// 
//       }
// 
//     }
// 
//    Survival1[2:T+1]=Survival1[2:T+1]/sum(Xvar[,K]);
//    Survival0[2:T+1]=Survival0[2:T+1]/(N-sum(Xvar[,K])); 
// 
// }

generated quantities{

    row_vector[T+1] Survival0=rep_row_vector(0,T+1);
    row_vector[T+1] Survival1=rep_row_vector(0,T+1);

    Survival0[1] = 1;
    Survival1[1] = 1;
    
    for (t in 2:T+1) {    
      
      for (i in 1:2000) {
      
      real a_ran = normal_rng(a_mu,a_sig)  ; // 
      
        Survival0[t] += exp(threg_IG_surv_ll(t,b_mu[1],exp(a_ran),1));  
        Survival1[t] += exp(threg_IG_surv_ll(t,sum(b_mu),exp(a_ran),1));
      
      }
    
    }
    
    Survival1[2:T+1]=Survival1[2:T+1]/2000;
    Survival0[2:T+1]=Survival0[2:T+1]/2000;

}


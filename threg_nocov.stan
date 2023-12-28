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
  
  int Nuc;
  real yuc[Nuc];
  int Nc;
  real yc[Nc];
  int T;

}

parameters {
  
  real<lower=0> y0; // initial heath status
  real mu; // mean

}
 
transformed parameters {
  
  real lny0 = log(y0);
  
}


model {

  mu ~ std_normal();
  y0 ~ gamma(1,1);

  for (i in 1:Nuc) {
    
    target += threg_IG_ll(yuc[i],mu,y0,1) ;
  
  }
  
  for (i in 1:Nc){
    
    target += threg_IG_surv_ll (yc[i],mu,y0,1);
    
  }
  
  
}


generated quantities{

   vector[T+1] expected;
  
   expected[1] = 1;
   
   for (i in 2:T+1) {

       expected[i]=(1)*exp(threg_IG_surv_ll(i,mu,y0,1));
   }
  
}


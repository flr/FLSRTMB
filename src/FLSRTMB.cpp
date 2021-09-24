#define TMB_LIB_INIT R_init_FLSRTMB

#include <TMB.hpp>

// Copyright Henning Winker (JRC) & Iago MOSQUEIRA (WMR), 2021
// Authors:  Henning Winker (JRC) <henning.winker@ec.europa.eu>
// Iago MOSQUEIRA (WMR) <iago.mosqueira@wur.nl>


// Space time
template<class Type>

// objective function
Type objective_function<Type>::operator() () {
  
  // Data
  DATA_VECTOR( rec );
  DATA_VECTOR( ssb );
  DATA_VECTOR( prior_s ); // Prior vector for s, [logit(mean), stdev in logit, useflag]
  DATA_VECTOR( spr0y );
  DATA_SCALAR( spr0 );
  DATA_SCALAR( plim ); // minimum bp of hockey-stick as fraction of Blim/B0
  DATA_INTEGER(nyears);
  DATA_INTEGER(Rmodel); // Recruitment model
  
  // Parameters
  PARAMETER(log_r0);
  PARAMETER(log_sigR);
  PARAMETER(logit_s);
  // Derived quantities
  Type r0 = exp(log_r0);
  Type sigR = exp(log_sigR);
  Type s = 0.2001 + 0.7998*1/(1+exp(-logit_s));
  vector<Type> log_rec_hat( nyears );
  vector<Type> vy = spr0y * r0;
  
  Type a = 0;
  Type b = 0;
  
  // Objective function
  Type ans=0;
  Type v = r0 * spr0;// compute SB0
  
  if(Rmodel==0){ // bevholtSV()
    for( int t=0; t< nyears; t++){
      log_rec_hat(t) = log(4.0 * s * r0 *ssb(t) / (vy(t)*(1.0-s)+ssb(t)*(5.0*s-1.0)));//-pow(sigR,2)/2.0;
    }
  }
  
  if(Rmodel==1){ // rickerSV()
    for( int t=0; t< nyears; t++){
      b = log(5.0*s)/(0.8*v);
      a = exp(b*v)/spr0;
      log_rec_hat(t) = log(a*ssb(t)*exp(-b*ssb(t)));
      //log_rec_hat(t) = log(r0 * ssb(t) / v * exp(s*(1.0-ssb(t)/v)));
    }
  }
  
  if(Rmodel==2){ // segreg() aka Hockey Stick
    for( int t=0; t< nyears; t++){
      //log_rec_hat(t) = log(r0)+log(2.5*s/v*(ssb(t)+0.2*v/s-pow(pow(ssb(t)-0.2*v/s,2.0),0.5)));//-pow(sigR,2)/2.0;
      log_rec_hat(t) = log(r0)+log(0.5/plim*s/vy(t)*(ssb(t)+plim*vy(t)/s-pow(pow(ssb(t)-plim*vy(t)/s,2.0),0.5)));
    }
  }
  
  vector<Type> rec_hat = exp(log_rec_hat);
  
  // OEM
  for( int t=0; t<nyears; t++){
    ans -= dnorm( log(rec(t)), log_rec_hat(t), sigR, true );
  }
  
  //prior s
  ans -= dnorm(logit_s, prior_s(0), prior_s(1), 1); // Prior for logn
  
  
  if(Rmodel==0){
    a = Type(4)*v*s/(spr0*(Type(5)*s-Type(1)));
    b = v*(Type(1)-s)/(Type(5)*s-Type(1));
  }
  
  if(Rmodel==1){
    //b = log(5.0*s)/(0.8*v);
    //a = exp(b*v)/spr0;
  }
  
  if(Rmodel==2){
    b = plim*v/s;
    a = r0/b;
  }
  
  // prior s
  
  
  // Reporting
  REPORT( rec_hat );
  REPORT( nyears );
  REPORT( sigR );
  REPORT( r0 );
  REPORT( v );
  REPORT( a );
  REPORT( b );
  REPORT( s );
  ADREPORT(a);
  ADREPORT( b );
  
  return ans;
}

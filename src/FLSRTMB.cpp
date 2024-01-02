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
  DATA_VECTOR( prior_r0 );
  DATA_VECTOR( prior_d ); // Flat logistic prior for depensation
  DATA_VECTOR( spr0y );
  DATA_SCALAR( spr0 );
  DATA_SCALAR( plim ); // minimum bp of hockey-stick as fraction of Blim/B0
  DATA_INTEGER(nyears);
  DATA_SCALAR( smin ); 
  DATA_SCALAR( smax ); 
  DATA_SCALAR( dmin ); 
  DATA_SCALAR( dmax ); 
  DATA_INTEGER(Rmodel); // Recruitment model
  DATA_INTEGER(depensationModel);
  
  // Parameters
  PARAMETER(log_r0);
  PARAMETER(log_sigR);
  PARAMETER(logit_s);
  PARAMETER(logit_d);
  // Derived quantities
  Type r0 = exp(log_r0);
  Type sigR = exp(log_sigR);
  Type sinit = smin+0.001 + (smax-smin-0.001)*1/(1+exp(-logit_s));
  Type dinit = dmin+0.001 + (dmax-dmin-0.001)*1/(1+exp(-logit_d))+0.5;
  
  //Type d = exp(log_d);
  
  vector<Type> log_rec_hat( nyears );
  vector<Type> vy = spr0y * r0;
  
  Type a = 0.001;
  Type b = 0.001;
  Type s = 0.5;
  Type d = 1.0;
  
  // Objective function
  Type ans=0;
  Type v = r0 * spr0;// compute SB0
  
  if(Rmodel==0){ // bevholtSV()
    s = sinit;
    for( int t=0; t< nyears; t++){
      // log_rec_hat(t) = log(4.0 * s * r0 *ssb(t) / (vy(t)*(1.0-s)+ssb(t)*(5.0*s-1.0)));//-pow(sigR,2)/2.0;
      // Rescale SSB by point of 50% recruitment
      Type rmaxx = 4.0 * s * r0 / (5.0 * s - 1.0);
      Type s50x = vy(t) * (1 - s) / (5.0 * s - 1.0);
      Type ssbx = ssb(t) / s50x;
      if(depensationModel == 1){
      d = dinit;
      ssbx = pow(ssbx, d);
      }
      log_rec_hat(t) = log(rmaxx / (1.0 + 1.0 / ssbx));
    }
  }
  
  if(Rmodel==1){ // rickerSV()
    for( int t=0; t< nyears; t++){
      s = sinit; // removed *20
      b = log(5.0*s)/(0.8*vy(t));
      a = exp(b*vy(t))/spr0;
      //log_rec_hat(t) = log(a*ssb(t)*exp(-b*ssb(t)));
      //log_rec_hat(t) = log(r0 * ssb(t) / v * exp(s*(1.0-ssb(t)/v)));
      // Rescale SSB by point of maximum recruitment (1/b);
      Type ssbx = ssb(t) * b;
      if(depensationModel == 1){
	    ssbx = pow(ssbx, d);
      }
      log_rec_hat(t) = log(a / b * ssbx * exp(-ssbx));
    }
  }
  
  if(Rmodel==2){ // segreg() aka Hockey Stick
    for( int t=0; t< nyears; t++){
      s = sinit;
      //log_rec_hat(t) = log(r0)+log(2.5*s/v*(ssb(t)+0.2*v/s-pow(pow(ssb(t)-0.2*v/s,2.0),0.5)));//-pow(sigR,2)/2.0;
      // log_rec_hat(t) = log(r0)+log(0.5/plim*s/vy(t)*(ssb(t)+plim*vy(t)/s-pow(pow(ssb(t)-plim*vy(t)/s,2.0),0.5)));
      // Re-parameterize to have breakpoint at 1 by scaling with breakpoint
      Type ssbx = ssb(t) / (plim*vy(t)/s);
      if(depensationModel == 1){ // Type A: R(S^d)
	    ssbx = pow(ssbx, d);
      }
      log_rec_hat(t) = log(r0) + log(0.5 * (ssbx + 1.0 - pow(pow(ssbx - 1.0, 2.0), 0.5)));

    }
  }
  
  vector<Type> rec_hat = exp(log_rec_hat);
  
  // OEM
  for( int t=0; t<nyears; t++){
    ans -= dnorm( log(rec(t)), log_rec_hat(t), sigR, true );
  }
  
  //prior s
  ans -= dnorm(logit_s, prior_s(0), prior_s(1), 1); // Prior for logn
  
  //r0 prior
  if(prior_r0(2)==1){
  ans -= dnorm(log_r0, prior_r0(0), prior_r0(1), 1);
  }
  
  //prior soft d
  if(depensationModel == 1){
  ans -= dnorm(logit_d, prior_d(0), prior_d(1), 1); // Prior for logn
  }  
  
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

 
  
  // Reporting
  REPORT( rec_hat );
  REPORT( nyears );
  REPORT( sigR );
  REPORT( r0 );
  REPORT( v );
  REPORT( a );
  REPORT( b );
  REPORT( s );
  REPORT( d ) ;
  ADREPORT(log(a));
  ADREPORT( log(b));
  
  return ans;
}

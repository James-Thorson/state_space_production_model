// Space time
#include <TMB.hpp>

// square
template<class Type>
Type square(Type x){
  return pow(x,2.0); 
}

// square
template<class Type>
Type sqrt(Type x){
  return pow(x,0.5); 
}

// objective function
template<class Type>
Type objective_function<Type>::operator() ()
{
  // Data
  DATA_VECTOR( y_t );                
  DATA_VECTOR( c_t );                
  DATA_VECTOR( ysd_t );
  DATA_VECTOR( penalties_z );
  
  // Parameters
  PARAMETER( log_r );
  PARAMETER( log_k );
  PARAMETER( log_q );
  PARAMETER( log_sigmap );
  PARAMETER( log_sigmam );
  PARAMETER( log_sigmac );
  PARAMETER_VECTOR( log_x_t );
  PARAMETER_VECTOR( logit_exploit_t );
  
  // Derived quantities
  int n_years = y_t.size();
  Type r = exp(log_r);
  Type k = exp(log_k);
  Type q = exp(log_q);
  vector<Type> xpred_t( n_years );
  vector<Type> cpred_t( n_years );
  vector<Type> sigmam_t( n_years );
  vector<Type> exploit_t = invlogit( logit_exploit_t );
  vector<Type> x_t = exp( log_x_t );
  vector<Type> Depletion_t = x_t / k;
  cpred_t.setZero();
  
  // Objective funcction
  vector<Type> jnll_comp(5);
  jnll_comp.setZero();
  
  // Reconstruct time series
  //cpred_t(0) = x_t(0);
  for( int t=1; t<n_years; t++){
    cpred_t(t-1) = x_t(t-1) * exploit_t(t-1);
    xpred_t(t) = (x_t(t-1)-cpred_t(t-1)) * (1.0 + r*(1-(x_t(t-1)-cpred_t(t-1))/k));
    jnll_comp(0) -= dnorm( x_t(t), xpred_t(t), xpred_t(t)*exp(log_sigmap), true );
  }
  
  // Probability of data
  for( int t=0; t<n_years; t++){
    sigmam_t(t) = sqrt( exp(2.0*log_sigmam) + square(ysd_t(t)) );
    jnll_comp(1) -= dnorm( y_t(t), q*x_t(t), q*x_t(t)*sigmam_t(t), true );   
  }                                                                                 
  
  // Penalty on catches
  for( int t=1; t<n_years; t++){
    jnll_comp(2) -= dnorm( c_t(t-1), cpred_t(t-1), cpred_t(t-1)*exp(log_sigmac), true );
  }
  
  // Penalty on starting biomass
  jnll_comp(3) -= dnorm( log(Depletion_t(0)), Type(0.0), exp(penalties_z(0)), true );
  
  // Penalty on interannual variability in exploitation ratio
  for( int t=1; t<(n_years-1); t++){
    jnll_comp(3) -= dnorm( logit_exploit_t(t), logit_exploit_t(t-1), exp(penalties_z(1)), true );
  }
  
  // Jacobian for log-random effects
  //jnll_comp(4) -= log_x_t.sum();
  
  // Total likelihood
  Type jnll = jnll_comp.sum();

  // Reporting
  REPORT( log_x_t );
  REPORT( n_years );
  REPORT( cpred_t );
  REPORT( x_t );
  REPORT( xpred_t );
  REPORT( jnll_comp );
  REPORT( r );
  REPORT( k );
  REPORT( q );
  REPORT( exploit_t );
  REPORT( Depletion_t );
  REPORT( log_sigmam );
  REPORT( log_sigmap );

  ADREPORT( log_x_t );
  ADREPORT( x_t );
  ADREPORT( Depletion_t );

  return jnll;
}

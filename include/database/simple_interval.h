/* INTERVAL CLASS */

// note: this is extremely simplistic and doesn't provide rigorous bounds yet

#ifndef CMDP_SIMPLEINTERVAL_H
#define CMDP_SIMPLEINTERVAL_H

#include <cmath>
#include <limits>

// boost interval doesn't seem to want to give me exp
// This simple alternative will work as long as we as approximate
// interval calculations are suitable for our purposes

template < class Real >
struct simple_interval {
  Real lower_;
  Real upper_;
  simple_interval ( void ) {}
  simple_interval ( Real lower_ ) : lower_(lower_), upper_(lower_) {}
  simple_interval ( Real lower_, Real upper_ ) : lower_(lower_), upper_(upper_) {}
  Real lower ( void ) const { return lower_; }
  Real upper ( void ) const { return upper_; }
  Real mid ( void ) const { return (upper_ + lower_) / 2.0; }
  Real radius ( void ) const { return (upper_ - lower_) / 2.0; }
};

template < class Real >
simple_interval<Real> operator * ( const Real lhs, const simple_interval<Real> & rhs ) {
  simple_interval<Real> result;
  result . lower_ = lhs * rhs . lower_;
  result . upper_ = lhs * rhs . upper_;
  if ( result . lower_ > result . upper_ ) std::swap ( result.lower_, result.upper_ );
  return result;
}

template < class Real >
simple_interval<Real> operator * ( const simple_interval<Real> & lhs, Real rhs ) {
  simple_interval<Real> result;
  result . lower_ = rhs * lhs . lower_;
  result . upper_ = rhs * lhs . upper_;
  if ( result . lower_ > result . upper_ ) std::swap ( result.lower_, result.upper_ );
  return result;
}

template < class Real >
simple_interval<Real> operator * ( const simple_interval<Real> & lhs, const simple_interval<Real> & rhs ) {
  simple_interval<Real> result;
  Real a = lhs . lower_ * rhs . lower_;
  Real b = lhs . lower_ * rhs . upper_;
  Real c = lhs . upper_ * rhs . lower_;
  Real d = lhs . upper_ * rhs . upper_;
  result . lower_ = std::min ( std::min ( a, b ), std::min ( c, d ) );
  result . upper_ = std::max ( std::max ( a, b ), std::max ( c, d ) );
  return result;
}

template < class Real >
simple_interval<Real> operator + ( const simple_interval<Real> & lhs, const simple_interval<Real> & rhs ) {
  simple_interval<Real> result;
  result . lower_ = lhs . lower_ + rhs . lower_;
  result . upper_ = lhs . upper_ + rhs . upper_;
  return result;  
}

template < class Real >
simple_interval<Real> operator + ( const Real lhs, const simple_interval<Real> & rhs ) {
  simple_interval<Real> result;
  result . lower_ = lhs + rhs . lower_;
  result . upper_ = lhs + rhs . upper_;
  return result;  
}

template < class Real >
simple_interval<Real> operator + ( const simple_interval<Real> & lhs, const Real rhs ) {
  simple_interval<Real> result;
  result . lower_ = lhs . lower_ + rhs;
  result . upper_ = lhs . upper_ + rhs;
  return result;
}

template < class Real >
simple_interval<Real> operator - ( const simple_interval<Real> & lhs, const simple_interval<Real> & rhs ) {
  simple_interval<Real> result;
  result . lower_ = lhs . lower_ - rhs . upper_;
  result . upper_ = lhs . upper_ - rhs . lower_;
  if ( result . lower_ > result . upper_ ) std::swap ( result.lower_, result.upper_ );
  return result;  
}

template < class Real >
simple_interval<Real> operator - ( const Real lhs, const simple_interval<Real> & rhs ) {
  simple_interval<Real> result;
  result . lower_ = lhs - rhs . upper_;
  result . upper_ = lhs - rhs . lower_;
  if ( result . lower_ > result . upper_ ) std::swap ( result.lower_, result.upper_ );
  return result;  
}

template < class Real >
simple_interval<Real> operator - ( const simple_interval<Real> & lhs, const Real rhs ) {
  simple_interval<Real> result;
  result . lower_ = lhs . lower_ - rhs;
  result . upper_ = lhs . upper_ - rhs;
  if ( result . lower_ > result . upper_ ) std::swap ( result.lower_, result.upper_ );
  return result;
}


template < class Real >
simple_interval<Real> pow ( const simple_interval<Real> & base, const Real exponent ) {
  simple_interval<Real> result;
  //if ( base . lower () <= 0.0  && exponent <= 0.0 ) std::cout << "pow A\n";
  //if ( base . upper () <= 0.0  && exponent <= 0.0 ) std::cout << "pow B.\n";
  //if ( base . lower () <= 0.0 ) std::cout << "pow C. " << base . lower () << ", " << base . upper () << "\n";
  //if ( base . upper () <= 0.0 ) std::cout << "pow D. " << base . lower () << ", " << base . upper () << "\n";
  

  if (exponent > 0){
    result . lower_ = std::pow( base . lower_ , exponent );
    result . upper_ = std::pow( base . upper_ , exponent );                   
  }                                                                      
  else if (exponent < 0){                                   
    result . lower_ = std::pow( base . upper_ , exponent );
    result . upper_ = std::pow( base . lower_ , exponent );
  }                                                                            
  else{                                                                        
    result . lower_ = 1.0;                                                     
    result . upper_ = 1.0;                                                     
  }      
  return result;
}

template < class Real >
simple_interval<Real> exp ( const simple_interval<Real> & exponent ) {
  simple_interval<Real> result;
  result . lower_ = std::exp ( exponent . lower_ );
  result . upper_ = std::exp ( exponent . upper_ );
  return result;  
}

template < class Real >
simple_interval<Real> log ( const simple_interval<Real> & term ) {
  simple_interval<Real> result;
  result . lower_ = std::log ( term . lower_ );
  result . upper_ = std::exp ( term . upper_ );
  return result;
} 

template < class Real >
simple_interval<Real> cos ( const simple_interval<Real> & term ) {
  const Real pi = 3.1415926535897932384626433832795;
  simple_interval<Real> result;
  Real low = term . lower ();
  Real high = term . upper ();
  if ( high - low >= 2*pi ) return simple_interval<Real> (-1.0, 1.0);
  while ( high > 2*pi ) {
    high -= 2*pi;
    low -= 2*pi;
  }
  while ( low < 0 ) {
    high += 2*pi;
    low += 2*pi;
  }
  if ( low <= pi ) {
    if ( high <= pi ) {
      result . lower_ = cos ( high );
      result . upper_ = cos ( low );
    } else {
      result . lower_ = -1.0;
      result . upper_ = std::max ( cos (high), cos (low) );
    }
  } else {
    if ( high <= 2*pi ) {
      result . lower_ = cos ( low );
      result . upper_ = cos ( high );
    } else {
      result . lower_ = std::min ( cos ( low ), cos ( high ) );
      result . upper_ = 1.0;
    }
  }
  
  return result;
} 

template < class Real >
simple_interval<Real> sin ( const simple_interval<Real> & term ) {  
  const Real pi = 3.1415926535897932384626433832795;
  return cos ( term - pi/2.0);
} 

template < class Real >
simple_interval<Real> tan ( const simple_interval<Real> & term ) {
  const Real pi = 3.1415926535897932384626433832795;
  simple_interval<Real> result;
  Real low = term . lower ();
  Real high = term . upper ();
  Real tanlow = tan ( low );
  Real tanhigh = tan ( high );
  if ( (high - low >= pi) ||
       (tanlow > tanhigh) ) return simple_interval<double> ( -std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity());
  return simple_interval<double> ( tanlow, tanhigh );
}

template < class Real >
simple_interval<Real> cot ( const simple_interval<Real> & term ) {
  const Real pi = 3.1415926535897932384626433832795;
  simple_interval<Real> result;
  Real low = term . lower ();
  Real high = term . upper ();
  Real cotlow = 1.0 / tan ( high );
  Real cothigh = 1.0 / tan ( low );
  if ( (high - low >= pi) ||
      (cotlow > cothigh) ) return simple_interval<double> ( -std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity());
  return simple_interval<double> ( cotlow, cothigh );
}


template < class Real >
simple_interval<Real> tanh (const simple_interval<Real> & term) {
	simple_interval<Real> result;
	result . lower_ = std::tanh (term . lower_ );
	result . upper_ = std::tanh (term . upper_ );
	return result;
}

template < class Real >
simple_interval<Real> square (const simple_interval<Real> & term) {
  simple_interval<Real> result;
  if (term . lower_ < 0 && term . upper_ > 0){
    if (std::abs (term . lower_) > term.upper_){
      result . lower_ = 0;
      result . upper_ = term.lower_ * term.lower_;
    }
    else{
      result . lower_ = 0;
      result . upper_ = term.upper_ * term.upper_;
    }
  }
  else{
    if (std::abs (term . lower_)> std::abs(term.upper_)){
      result . lower_ = term . upper_ * term . upper_;
      result . upper_ = term . lower_ * term . lower_;
    }
    else{
      result . lower_ = term . lower_ * term . lower_;
      result . upper_ = term . upper_ * term . upper_;
    }
  }
  return result;
}

template < class Real >
simple_interval<Real> operator / ( const simple_interval<Real> & lhs, const simple_interval<Real> & rhs ) {
  return lhs * pow ( rhs, -1.0 );
}

#endif

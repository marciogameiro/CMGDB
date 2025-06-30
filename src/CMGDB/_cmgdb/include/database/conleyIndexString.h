#ifndef CMDB_CONLEY_INDEX_STRING_H
#define CMDB_CONLEY_INDEX_STRING_H


#include "chomp/ConleyIndex.h"
#include "chomp/SparseMatrix.h"
#include "chomp/PolyRing.h"
#include "chomp/Ring.h"
#include "chomp/FrobeniusNormalForm.h"
#include <boost/thread.hpp>
#include <boost/chrono/chrono_io.hpp>


class FrobeniusThread {
public:
  typedef chomp::PolyRing<chomp::Ring> Polynomial;
  typedef chomp::SparseMatrix<chomp::Ring> Matrix;
  bool * result;
  FrobeniusThread( std::vector< Polynomial > * invariant_factors,
                   const Matrix & A, 
                   bool * result ) 
  : result(result), invariant_factors_(invariant_factors), A_(A)  {}

  void operator () ( void ) {
    try {
      *invariant_factors_ = chomp::FrobeniusNormalForm ( A_ );
      *result = true;
    } catch ( ... /* boost::thread_interrupted& */) {
      *result = false;
    }
  }

private:
  std::vector<Polynomial> * invariant_factors_;
  const Matrix & A_;
};

std::vector<chomp::PolyRing < chomp::Ring > > 
shiftClass ( const std::vector< chomp::PolyRing < chomp::Ring > > & invariant_factors ) {
  using namespace chomp;
  typedef PolyRing<Ring> Polynomial;
  Polynomial x;
  std::vector<Polynomial> result;
  x . resize ( 2 );
  x[1] = Ring ( 1 );
  for ( int j = 0; j < invariant_factors . size (); ++ j ) {
    Polynomial factor = invariant_factors [ j ];
    while ( factor [ 0 ] == Ring ( 0 ) ) factor = factor / x;
    if ( factor . degree () > 0 ) result . push_back ( factor );
  }
  return result;
}

inline std::vector<std::string> 
conleyIndexString ( const chomp::ConleyIndex_t & ci, 
                    int * errorcode = NULL,
                    int time_out = 3600 ) {
  if ( errorcode != NULL ) * errorcode = 0;
  // std::cout << "conleyIndexString.\n";
  std::vector<std::string> result;
  if ( ci . undefined () ) { 
    // std::cout << "conleyIndexString. undefined.\n";
    if ( errorcode != NULL ) * errorcode = 4;
    return result;
  }
  for ( unsigned int i = 0; i < ci . data () . size (); ++ i ) {
    // std::cout << "conleyIndexString. Dimension is " << i << "\n";
    typedef chomp::PolyRing<chomp::Ring> Polynomial;

    // use a thread to compute Frobenius Normal Form
    std::vector<Polynomial> invariant_factors;
    bool computed;
    FrobeniusThread frobenius ( &invariant_factors, ci . data () [ i ], &computed );
    boost::thread t(frobenius);
    if ( not t . try_join_for ( boost::chrono::seconds( time_out ) ) ) {
      t.interrupt();
      t.join();
    }
    if ( not computed ) {
      result . push_back ( std::string ( "Problem computing Frobenius Form.\n") );
      if ( errorcode != NULL ) * errorcode = 1;
      continue;
    }
    // end threading

    std::vector<Polynomial> shift_class = shiftClass ( invariant_factors );

    std::stringstream ss;
    BOOST_FOREACH ( const Polynomial & poly, shift_class ) {
      ss << poly;
      // ss << poly << "\n";
    }
    if ( shift_class . empty () ) {
      ss << "0";
      // ss << "Trivial.\n";
    }

    result . push_back ( ss . str () );
    // std::cout << "conleyIndexString. Wrote the polynomial " << ss . str () << "\n";
  }
  return result;
}

#endif

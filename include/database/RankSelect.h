#ifndef RANKSELECT_H
#define	RANKSELECT_H
//#define SDSL_DEBUG_BP

/// @file RankSelect.h
/// @author Arnaud Goullet, Shaun Harker
/// @description SDSL interface for Rank/Select operations on a binary sequence.
///          This wrapper exists for two reasons.
///           (1) We slightly redefine rank/select so that they are inverse
///           (2) The class is serializable

#include "boost/serialization/serialization.hpp"
#include "boost/serialization/vector.hpp"
#include <inttypes.h>
#include "sdsl/rank_support_v5.hpp"
#include "sdsl/select_support_mcl.hpp"
#include "sdsl/util.hpp"

/// RankSelect
class RankSelect {
public:
  typedef uint64_t size_type;

  /// RankSelect
  ///    Default constructor
  RankSelect ( ) {}

  /// RankSelect
  ///   Constructor from a given std::vector<bool> bit sequence
  ///   @param bits
  RankSelect ( const std::vector < bool > & bits ) {
    assign ( bits );
  }

  /// assign
  ///    Delayed constructor.
  void assign ( const std::vector < bool > & bits ) {
    bits_ . resize ( bits . size () );
    for ( size_t i = 0; i < bits . size (); ++ i ) bits_[i] = bits[i];
    rank_ = sdsl::rank_support_v5 < > ( &bits_ );
    select_ = sdsl::select_support_mcl < > ( &bits_ );
  }

  /// rank
  ///   @return the rank of the bit sequence at a given position. 
  ///   Here the rank is defined as the number of 1's on [0,i-1]  
  size_type rank ( size_type i ) const {
    return rank_ . rank ( i );
  }

  /// select
  ///   @return the position of the i-th 1 in the bit sequence
  ///    We let indexing start at 0, i.e. select(0) gives the position
  ///    of the first 1 in the bit sequence.
  size_type select ( size_type i ) const {
    return select_ . select ( i + 1 );
  }

  /// bits
  ///   @param position
  ///   @return the value of the bit at a given position
  bool bits ( size_type position ) const {
    return bits_ [ position ];
  }
   
  /// bitSequence
  ///   @return the complete bit sequence over which the rank/select operates on.  
  const sdsl::bit_vector & bitSequence ( void ) const {
    return bits_;
  }

  /// operator =
  ///    assignment operator
  RankSelect& operator= ( const RankSelect& other ) {
    bits_ = other . bits_;
    rank_ = sdsl::rank_support_v5 < > ( &bits_ );
    select_ = sdsl::select_support_mcl < > ( &bits_ );
    return *this;
  }

  /// memory
  ///   Return the memory usage in bytes.
  uint64_t memory ( void ) const {
    //std::cout << "RankSelect::memory. \n";
    //std::cout << "sdsl::util::get_size_in_bytes ( bits_ ) = " << sdsl::util::get_size_in_bytes ( bits_ ) << "\n";
    //std::cout << "sdsl::util::get_size_in_bytes ( rank_ ) = " << sdsl::util::get_size_in_bytes ( rank_ ) << "\n";
    //std::cout << "sdsl::util::get_size_in_bytes ( select_ ) = " << sdsl::util::get_size_in_bytes ( select_ ) << "\n";

    return sizeof ( RankSelect ) +  
        sdsl::size_in_bytes ( bits_ ) +
        sdsl::size_in_bytes ( rank_ ) +
        sdsl::size_in_bytes ( select_ );      
  }

private:
  sdsl::bit_vector bits_;
  sdsl::rank_support_v5 < > rank_;
  sdsl::select_support_mcl < > select_;
  
  friend class boost::serialization::access;
  template<class Archive>
  void save ( Archive & ar , const unsigned int version ) const {
    std::vector < bool > bit_sequence;
    for ( size_type i = 0; i < bits_.size(); ++ i ) bit_sequence . push_back ( bits_ [ i ] );
    ar & bit_sequence;   
  }

  template<class Archive>
  void load ( Archive & ar , const unsigned int version ) {
    std::vector < bool > bit_sequence;
    ar & bit_sequence;
    bits_ . resize ( bit_sequence . size() );
    for ( size_type i = 0; i < bit_sequence.size(); ++ i ) 
      bits_ [ i ] =  bit_sequence [ i ];
    rank_ = sdsl::rank_support_v5 < > ( &bits_ );
    select_ = sdsl::select_support_mcl < > ( &bits_ );
  }
  BOOST_SERIALIZATION_SPLIT_MEMBER ( );
};
#endif	/* RANKSELECT_H */

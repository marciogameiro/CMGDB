#ifndef CMDB_SUCCINCTTREE_H
#define CMDB_SUCCINCTTREE_H
//#define SDSL_DEBUG_BP
/// @file SuccinctTree.h
/// @author Arnaud Goullet, Shaun Harker
/// @description This file defines class SuccinctTree which 
/// provides an implementation of a full binary tree using SDSL.
#include <exception>
#include "boost/foreach.hpp"
#include <memory>
#include "boost/iterator/counting_iterator.hpp"
#include "boost/archive/text_oarchive.hpp"
#include "boost/archive/text_iarchive.hpp"
#include "boost/serialization/serialization.hpp"
#include "boost/serialization/vector.hpp"
#include "boost/serialization/version.hpp"
#include "boost/serialization/split_member.hpp"
#include "sdsl/bp_support.hpp"
#include "sdsl/rank_support_v5.hpp"
#include "sdsl/select_support_mcl.hpp"
#include "sdsl/util.hpp"
#include "Tree.h"

/// SuccinctTree
///   Implements a full binary tree using SDSL
class SuccinctTree : public Tree {
public:
  typedef boost::counting_iterator < int64_t > iterator;
  typedef iterator const_iterator;
  typedef int64_t value_type;
  typedef uint64_t size_type;

  /// SuccinctTree
  ///   Default constructor. Creates root node with no children.
  SuccinctTree ( void );

  /// assign
  ///   Reset structure as if it had been constructed with leaf_sequence
  virtual void assign ( std::shared_ptr<const CompressedTree> compressed );

  /// assignFromLeafSequence 
  void assignFromLeafSequence ( const std::vector<bool> & leaf_sequence );

  /// leafEnd
  ///   Give the one-past-the-end leaf (i.e. return number of leaves)
  int64_t leafEnd ( void ) const;

  /// TreeToLeaf
  ///    Given a tree iterator, return a leaf iterator
  int64_t TreeToLeaf ( iterator it ) const;

  /// LeafToTree
  ///    Given a leaf iterator, return a tree iterator
  iterator LeafToTree ( int64_t leaf ) const;

  /// parent
  ///   @param it
  ///   @return the iterator pointing to the parent node.
  iterator parent ( iterator it ) const;
  
  /// left
  ///   @param it
  ///   @return the iterator pointing to the left child.
  iterator left ( iterator it ) const;
  
  /// right
  ///   @param it
  ///   @return the iterator pointing to the right child.
  iterator right ( iterator it ) const;
  
  /// isLeft
  ///   @param it
  ///   @return true if the iterator it is a left child.
  bool isLeft ( iterator it ) const;

  /// isRight
  ///   @param it
  ///   @return true if the iterator it is a right child.
  bool isRight ( iterator it ) const;

  /// isLeaf
  ///    @param it
  ///    @return true if the iterator it is a leaf.
  bool isLeaf ( iterator it ) const;
 
  /// memory 
  ///   Return the memory usage of data structure in bytes
  uint64_t memory ( void ) const;

  /// leafSequence 
  ///   Return a reference to the underlying bit_vector
  const sdsl::bit_vector & leafSequence ( void ) const;

private:
  sdsl::bit_vector leaf_sequence_;
  sdsl::bp_support_sada < > tree_;
  sdsl::rank_support_v5 <0> rank_;
  sdsl::select_support_mcl <0> select_;
  int64_t leaf_count_;
  friend class boost::serialization::access;
  template<class Archive>
  void save ( Archive & ar , const unsigned int version ) const;
  template<class Archive>
  void load ( Archive & ar , const unsigned int version );
  BOOST_SERIALIZATION_SPLIT_MEMBER ( );
};

#if 0
BOOST_CLASS_EXPORT_KEY(SuccinctTree);
#endif

inline 
SuccinctTree::SuccinctTree ( void ) {
  // TODO: default constructor
}

inline void 
SuccinctTree::assign ( std::shared_ptr<const CompressedTree> compressed ) {
  const std::vector < bool > & leaf_sequence = compressed -> leaf_sequence;
  assignFromLeafSequence ( leaf_sequence );
}

inline void 
SuccinctTree::assignFromLeafSequence ( const std::vector<bool> & leaf_sequence ) {
  leaf_count_ = 0; 
  size_ = leaf_sequence . size ();
  leaf_sequence_ . resize ( leaf_sequence . size () + 1 );
  leaf_sequence_ [ 0 ] = 1;
  for ( size_t i = 0; i < leaf_sequence . size (); ++ i ) {
    leaf_sequence_ [ i + 1 ] = (leaf_sequence [ i ] ? 1 : 0 );
    if ( not leaf_sequence [ i ] ) ++ leaf_count_;
  }
  tree_ = sdsl::bp_support_sada <> ( & leaf_sequence_ );
  rank_ = sdsl::rank_support_v5 <0> ( &leaf_sequence_ );
  select_ =  sdsl::select_support_mcl <0> ( &leaf_sequence_ );
}

inline int64_t 
SuccinctTree::leafEnd ( void ) const {
  return leaf_count_;
}

inline int64_t 
SuccinctTree::TreeToLeaf ( iterator it ) const {
  int64_t x = *it + 1;
  if ( leaf_sequence_ [ x ] == 1 ) return leafEnd ();
  return rank_ . rank ( x );
}

inline SuccinctTree::iterator 
SuccinctTree::LeafToTree ( int64_t leaf ) const {
  return iterator ( select_ . select ( leaf + 1 ) - 1 );
}

inline SuccinctTree::iterator 
SuccinctTree::parent ( iterator it ) const {
  int64_t x = *it;
  if ( x == 0 ) return end ();
  if ( leaf_sequence_ [ x ] == 1 ) return iterator ( x - 1 );
  return iterator ( tree_ . find_open ( x ) - 1 );
}

inline SuccinctTree::iterator 
SuccinctTree::left ( iterator it ) const {
  int64_t x = *it + 1;
  if ( leaf_sequence_ [ x ] == 0 ) return end ();
  return iterator ( x );
}

inline SuccinctTree::iterator 
SuccinctTree::right ( iterator it ) const {
  int64_t x = *it + 1;
  if ( leaf_sequence_ [ x ] == 0 ) return end ();
  return iterator ( tree_ . find_close ( x ) );
}

inline bool 
SuccinctTree::isLeft ( iterator it ) const {
  int64_t x = *it;
  if ( x == 0 ) return false;
  if ( leaf_sequence_ [ x ] == 1 ) return true;
  return false;
}

inline bool 
SuccinctTree::isRight ( iterator it ) const {
  int64_t x = *it;
  if ( x == 0 ) return false;
  if ( leaf_sequence_ [ x ] == 1 ) return false;
  return true;
}

inline bool 
SuccinctTree::isLeaf ( iterator it ) const {
  int64_t x = *it + 1;
  if ( leaf_sequence_ [ x ] == 0 ) return true;
  return false;
}

inline uint64_t 
SuccinctTree::memory ( void ) const {
  return sizeof ( SuccinctTree ) + 
         sdsl::size_in_bytes ( leaf_sequence_ ) +
         sdsl::size_in_bytes ( tree_ ) +
         sdsl::size_in_bytes ( rank_ ) +
         sdsl::size_in_bytes ( select_);         
}

inline const sdsl::bit_vector & 
SuccinctTree::leafSequence ( void ) const {
  return leaf_sequence_;
}

template<class Archive>
void SuccinctTree::save ( Archive & ar , const unsigned int version ) const {
  ar & boost::serialization::base_object<Tree>(*this);
  std::vector < bool > leaf_sequence ( size_ );
  for ( size_t i = 1; i < leaf_sequence_ . size (); ++ i ) {
    leaf_sequence [ i - 1 ] = leaf_sequence_ [ i ];
  }
  ar & leaf_sequence;
}

template<class Archive>
void SuccinctTree::load ( Archive & ar , const unsigned int version ) {
  ar & boost::serialization::base_object<Tree>(*this);
  std::vector < bool > leaf_sequence;
  ar & leaf_sequence;
  assignFromLeafSequence ( leaf_sequence );
}

#endif

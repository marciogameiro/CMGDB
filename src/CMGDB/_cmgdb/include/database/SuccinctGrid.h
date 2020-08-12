#ifndef CHOMP_SUCCINCTGRID_H
#define CHOMP_SUCCINCTGRID_H
//#define SDSL_DEBUG_BP

/// @file SuccinctGrid
/// @author Arnaud Goullet, Shaun Harker
/// @description Implementation of TreeGrid interface via succinct data structures

#include <iostream>
#include <fstream>
#include <vector>
#include <stack>
#include <deque>
#include "boost/unordered_set.hpp"
#include "boost/unordered_map.hpp"
#include <memory>
#include "boost/serialization/serialization.hpp"
#include "boost/serialization/vector.hpp"
#include "boost/serialization/export.hpp"
#include "TreeGrid.h"
#include "Tree.h"
#include "SuccinctTree.h"
#include "RankSelect.h"

/// SuccinctGrid
//
// BEGIN TODO
// Warning: The following statistics no longer apply --- TODO: replace them
//           However they must be upper bounds.
// For a full tree with number of nodes > 100,000
// the memory usage is about
// 2.9 bits / node for the full binary tree
// 1.29 bits / nodes for the rank select of the node validity
// 1.18 bits / nodes for the rank select of the leaves
// Total ~ 5.37 bits / nodes
// END TODO
class SuccinctGrid : public TreeGrid {
public:
  // Constructor/Deconstructor Methods
  SuccinctGrid ( void );
  virtual ~SuccinctGrid ( void );
  // Virtual Methods
  virtual Tree::iterator GridToTree ( Grid::iterator it ) const;
  virtual Grid::iterator TreeToGrid ( Tree::iterator it ) const;
  virtual int64_t       GridToLeaf ( Grid::iterator it ) const;
  virtual Grid::iterator LeafToGrid ( int64_t leaf ) const;
  virtual const SuccinctTree & tree ( void ) const;
  virtual SuccinctTree & tree ( void );
  virtual SuccinctGrid * spawn ( void ) const;
  virtual void rebuild ( std::shared_ptr<const CompressedTreeGrid> compressed );
  
private:
  // Data
  SuccinctTree tree_;
  RankSelect valid_sequence_;
  // Serialization
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version)
  {
    ar & boost::serialization::base_object<TreeGrid>(*this);
    ar & tree_;
    ar & valid_sequence_;
  }
public:
  // Test and Debug
  virtual uint64_t memory ( void ) const {
    return sizeof(SuccinctGrid) +
           tree_ . memory () +
           valid_sequence_ . memory ();
  }
};

BOOST_CLASS_EXPORT_KEY(SuccinctGrid);

inline 
SuccinctGrid::SuccinctGrid ( void ) {
  size_ = 1;
  std::vector<bool> valid_sequence (1, true);
  valid_sequence_ . assign ( valid_sequence );
  std::vector<bool> leaf_sequence (1, false);
  tree () . assignFromLeafSequence ( leaf_sequence );
} /* SuccinctGrid::SuccinctGrid */


inline 
SuccinctGrid::~SuccinctGrid ( void ) {
}

inline Grid::iterator 
SuccinctGrid::TreeToGrid ( Tree::iterator it_tree ) const {
  int64_t leaf = tree () . TreeToLeaf ( it_tree );
  if ( leaf == tree () . leafEnd () ) return end ();
  return LeafToGrid ( leaf );
}

inline Tree::iterator 
SuccinctGrid::GridToTree ( iterator it ) const {
  return tree () . LeafToTree ( GridToLeaf ( it ) );
}

inline int64_t
SuccinctGrid::GridToLeaf ( Grid::iterator it ) const {
  return valid_sequence_ . select ( *it );
}

inline Grid::iterator 
SuccinctGrid::LeafToGrid ( int64_t leaf ) const {
  if ( not valid_sequence_ . bits ( leaf ) ) return end ();
  return iterator ( valid_sequence_ . rank ( leaf ) );
}

inline const SuccinctTree & 
SuccinctGrid::tree ( void ) const {
  return tree_;
}

inline SuccinctTree & 
SuccinctGrid::tree ( void ) {
  return tree_;
}

inline SuccinctGrid * 
SuccinctGrid::spawn ( void ) const {
  return new SuccinctGrid;
}

inline void 
SuccinctGrid::rebuild ( std::shared_ptr<const CompressedTreeGrid> compressed ) {
  valid_sequence_ . assign ( compressed -> tree () -> valid_sequence );
}

#endif

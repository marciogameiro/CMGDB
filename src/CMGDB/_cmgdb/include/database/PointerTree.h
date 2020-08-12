/// PointerTree.h
/// Shaun Harker

#ifndef CMDB_POINTERTREE_H
#define CMDB_POINTERTREE_H

#include <stdint.h>
#include <memory>
#include <vector>
#include <deque>
#include <stack>
#include <utility>
#include <boost/foreach.hpp>
#include <memory>
#include "Tree.h"
#include "CompressedTree.h"
#include "boost/serialization/serialization.hpp"
#include "boost/serialization/vector.hpp"
#include "boost/serialization/export.hpp"
#include "boost/serialization/split_member.hpp"

/// PointerTreeNode
struct PointerTreeNode;

/// PointerTree
class PointerTree : public Tree {
public:
  PointerTree ( void );
  virtual ~PointerTree ( void );
  virtual iterator parent ( iterator it ) const;
  virtual iterator left ( iterator it ) const;
  virtual iterator right ( iterator it ) const;
  virtual bool isLeft ( iterator it ) const;
  virtual bool isRight ( iterator it ) const;
  virtual bool isLeaf ( iterator it ) const;
  virtual void assign ( std::shared_ptr<const CompressedTree> compressed );
  virtual uint64_t memory ( void ) const;
private:
  std::vector < PointerTreeNode > nodes_;
  std::vector < bool > parity_;
  std::vector < bool > isleaf_; 
  friend class boost::serialization::access;
  template<typename Archive>
  void serialize(Archive & ar, const unsigned int file_version);
  template<typename Archive>
  void load(Archive & ar, const unsigned int file_version);
  template<typename Archive>
  void save(Archive & ar, const unsigned int file_version) const;

};

BOOST_CLASS_EXPORT_KEY(PointerTree);

/// PointerTreeNode
struct PointerTreeNode {
  int64_t left_;
  int64_t right_;
  int64_t parent_;
  PointerTreeNode ( void );
  PointerTreeNode ( int64_t left, int64_t right, int64_t parent );
  ~PointerTreeNode ( void );
  friend class boost::serialization::access;
  template<class Archive>
  void serialize ( Archive & ar, const unsigned int version );
};

// PointerTree Definitions
inline 
PointerTree::PointerTree ( void ) {
  nodes_ . push_back ( PointerTreeNode ( 1, 1, 1 ) );
  parity_ . push_back ( false );
  isleaf_ . push_back ( true );
  size_ = 1;

 // debug
  /*
  std::cout << "default constructed tree.\n";
  std::cout << "nodes_ . size () == " << nodes_ . size () << "\n";
  std::cout << "size() == " << size () << "\n";
  std::cout << "parity_ . size () == " << parity_ . size () << "\n";
  std::cout << "isleaf_ . size () == " << isleaf_ . size () << "\n";
  for ( size_t i = 0; i < nodes_ . size (); ++ i ) {
    std::cout << " Node " << i << ":\n";
    std::cout << "   left = " << nodes_[i].left_ << "\n";
    std::cout << "   right = " << nodes_[i].right_ << "\n";
    std::cout << "   parent = " << nodes_[i].parent_ << "\n";
    std::cout << "   parity = " << ((parity_[i])? "right" : "left") << "\n";
    std::cout << "   isleaf_ = " << ((isleaf_[i])? "yes" : "no" ) << "\n";
  }
*/
}

inline 
PointerTree::~PointerTree ( void ) {
}

inline Tree::iterator 
PointerTree::parent ( iterator it ) const {
  return Tree::iterator ( nodes_ [ *it ] . parent_ );
}

inline Tree::iterator 
PointerTree::left ( iterator it ) const {
  return Tree::iterator ( nodes_ [ *it ] . left_ );
}

inline Tree::iterator 
PointerTree::right ( iterator it ) const {
  return Tree::iterator ( nodes_ [ *it ] . right_ );
}

inline bool 
PointerTree::isLeft ( iterator it ) const {
  return not parity_ [ *it ];
}

inline bool 
PointerTree::isRight ( iterator it ) const {
  return parity_ [ *it ];

}

inline bool 
PointerTree::isLeaf ( iterator it ) const {
  return isleaf_ [ *it ];
}

inline void 
PointerTree::assign ( std::shared_ptr<const CompressedTree> compressed ) {
  const bool LEAF = false;
  const bool NOT_A_LEAF = true;
  const std::vector<bool> & leaf_sequence = 
    compressed -> leaf_sequence;
  const std::vector<bool> & valid_sequence = 
    compressed -> valid_sequence;
  size_t N = leaf_sequence . size ();
  std::stack<std::pair<Tree::iterator, int> > path_to_root;
  nodes_ . clear ();
  parity_ . clear ();
  isleaf_ . clear ();
  size_ = 0;
  if ( compressed -> leafCount () == 0 ) return;
  Tree::iterator sentinel ( -1 );
  path_to_root . push ( std::make_pair (sentinel, 0) );
  int64_t last_encountered_leaf = -1;
  for ( size_t i = 0; i < N; ++ i ) {
    while ( path_to_root . top () . second == 2 ) path_to_root . pop ();
    Tree::iterator parent = path_to_root . top () . first;
    int child_num = path_to_root . top () . second ++;
    if ( leaf_sequence [ i ] == LEAF ) {
      ++ last_encountered_leaf;
      if ( not valid_sequence [ last_encountered_leaf ] ) continue;
    }
    Tree::iterator node ( nodes_ . size () );
    nodes_ . push_back ( PointerTreeNode ( -1, -1, * parent ) );
    if ( child_num == 0 ) {
      if ( parent != sentinel ) nodes_ [ * parent ] . left_ = * node;
      parity_ . push_back ( false );
    } else {
      if ( parent != sentinel ) nodes_ [ * parent ] . right_ = * node;
      parity_ . push_back ( true );
    }
    if ( leaf_sequence [ i ] == NOT_A_LEAF ) {
      path_to_root . push ( std::make_pair ( node, 0 ) );
      isleaf_ . push_back ( false );
    } else {
      isleaf_ . push_back ( true );
    }
  }
  size_ = nodes_ . size ();
  for ( size_t i = 0; i < nodes_ . size (); ++ i ) {
    if ( nodes_[i] . left_ == -1 ) nodes_[i] . left_ = size ();
    if ( nodes_[i] . right_ == -1 ) nodes_[i] . right_ = size ();
    if ( nodes_[i] . parent_ == -1 ) nodes_[i] . parent_ = size ();
  }
}

inline uint64_t 
PointerTree::memory ( void ) const {
  return sizeof ( PointerTree ) +
         sizeof ( PointerTreeNode ) * nodes_ . size () + 
         2 * nodes_ . size ();
}

template<class Archive>
void PointerTree::save(Archive & ar, const unsigned int version) const
{
  ar << boost::serialization::base_object<Tree>(*this);
  ar << nodes_;
  ar << parity_;
  ar << isleaf_;
}

template<class Archive>
void PointerTree::load(Archive & ar, const unsigned int version)
{
  nodes_.clear();
  parity_.clear();
  isleaf_.clear();
  ar >> boost::serialization::base_object<Tree>(*this);
  ar >> nodes_;
  ar >> parity_;
  ar >> isleaf_;
}

template<class Archive>
void PointerTree::serialize(
    Archive & ar,
    const unsigned int file_version 
){
    boost::serialization::split_member(ar, *this, file_version);
}

/*
// Remark. Boost 1.58 seems to have introduced a behavior change in 
//         vector serialization where loading appends instead
//         of overwrites.

template<typename Archive> void 
PointerTree::serialize ( Archive & ar, 
                         const unsigned int file_version ) {
  ar & boost::serialization::base_object<Tree>(*this);
  ar & nodes_;
  ar & parity_;
  ar & isleaf_;
}
*/

// PointerTreeNode Definitions

template<class Archive> void 
PointerTreeNode::serialize ( Archive & ar, 
                             const unsigned int version ) {
    ar & left_;
    ar & right_;
    ar & parent_;
  }

inline PointerTreeNode::PointerTreeNode ( int64_t left, int64_t right, int64_t parent ) 
: left_ ( left ), right_ ( right ), parent_ ( parent ) {
}

inline PointerTreeNode::PointerTreeNode ( void ) 
: left_ ( 0 ), right_ ( 0 ), parent_ ( 0 ) {
}

inline PointerTreeNode::~PointerTreeNode ( void ) {
}
#endif

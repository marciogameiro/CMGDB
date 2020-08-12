// PointerGrid.h
// Shaun Harker
// 9/16/11

#ifndef CMDP_POINTERGRID_H
#define CMDP_POINTERGRID_H

#include <vector>
#include <stack>
#include <deque>
#include <exception>
#include <boost/unordered_set.hpp>
#include "TreeGrid.h"
#include "Tree.h"
#include "PointerTree.h"
#include <memory>
#include "boost/serialization/serialization.hpp"
#include "boost/serialization/vector.hpp"
#include "boost/serialization/export.hpp"
#include "boost/serialization/shared_ptr.hpp"
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

/// class PointerGrid
class PointerGrid : public TreeGrid {
public:
  PointerGrid ( void );
  virtual ~PointerGrid ( void );
  virtual Tree::iterator GridToTree ( Grid::iterator it ) const;
  virtual Grid::iterator TreeToGrid ( Tree::iterator it ) const;
  virtual const PointerTree & tree ( void ) const;
  virtual PointerTree & tree ( void );
  virtual PointerGrid * spawn ( void ) const;
  virtual void rebuild ( std::shared_ptr<const CompressedTreeGrid> compressed );
  void rebuildFromTree ( void );
private:
  std::shared_ptr<PointerTree> tree_;
  std::vector < Grid::iterator > grid_iterators_;
  std::vector < Tree::iterator > tree_iterators_;
public:
  virtual uint64_t memory ( void ) const {
    return sizeof ( PointerGrid ) +
           tree_ -> memory () +
           sizeof ( Grid::iterator ) * grid_iterators_ . size () +
           sizeof ( Tree::iterator ) * tree_iterators_ . size ();
  }
  
  
  friend class boost::serialization::access;
  template<typename Archive>
  void serialize(Archive & ar, const unsigned int file_version) {
    ar & boost::serialization::base_object<TreeGrid>(*this);
    ar & tree_;
    rebuildFromTree();
  }
  // file operations
  void save ( const char * filename ) const {
    std::ofstream ofs(filename);
    assert(ofs.good());
    boost::archive::text_oarchive oa(ofs);
    oa << * this;
  }
  
  void load ( const char * filename ) {
    std::ifstream ifs(filename);
    if ( not ifs . good () ) {
      std::cout << "Could not load " << filename << "\n";
      exit ( 1 );
    }
    boost::archive::text_iarchive ia(ifs);
    ia >> * this;    
  }
 
  
};

BOOST_CLASS_EXPORT_KEY(PointerGrid);

inline 
PointerGrid::PointerGrid ( void ) : tree_ ( new PointerTree ) {
  rebuildFromTree ();
}

inline 
PointerGrid::~PointerGrid ( void ) {
}

inline Grid::iterator 
PointerGrid::TreeToGrid ( Tree::iterator tree_it ) const {
  return grid_iterators_ [ * tree_it ];
}

inline Tree::iterator 
PointerGrid::GridToTree ( Grid::iterator grid_it ) const {
  return tree_iterators_ [ * grid_it ];
}

inline const PointerTree & 
PointerGrid::tree ( void ) const {
  return * tree_ . get ();
}

inline PointerTree & 
PointerGrid::tree ( void ) {
  return * tree_ . get ();
}

inline PointerGrid * 
PointerGrid::spawn ( void ) const {
  return new PointerGrid;
}

inline void 
PointerGrid::rebuild ( std::shared_ptr<const CompressedTreeGrid> compressed ) {
  rebuildFromTree ();
}

inline void 
PointerGrid::rebuildFromTree ( void ) {
  // Now we rebuild the GridIterator to TreeIterator conversions
  grid_iterators_ . clear ();
  tree_iterators_ . clear ();
  uint64_t leaf_count = 0;
  Tree::iterator end = tree () . end ();
  for ( Tree::iterator it = tree () . begin (); it != end; ++ it ) {
    if ( tree () . isLeaf ( it ) ) {
      grid_iterators_ . push_back ( Grid::iterator (leaf_count ++ ) );
      tree_iterators_ . push_back ( it );
    } else {
      grid_iterators_ . push_back ( tree () . size () );
    }
  }
  size_ = leaf_count;
}

#endif

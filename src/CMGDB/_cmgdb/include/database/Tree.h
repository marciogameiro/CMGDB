//Tree.h

#ifndef CMDB_TREE_H
#define CMDB_TREE_H

#include <memory>
#include <deque>
#include <boost/iterator/counting_iterator.hpp>
#include <stdint.h>
#include "boost/serialization/serialization.hpp"

#include "CompressedTree.h" // later i may refactor

class Tree {
public:
  
  // Typedefs
  typedef boost::counting_iterator < int64_t > iterator;
  typedef iterator const_iterator;
  typedef int64_t value_type;
  typedef uint64_t size_type;
  
  // Constructor/Deconstructor Methods
protected:
  Tree ( void );
  
public:
  virtual ~Tree ( void );
  
  // Container methods
  iterator begin ( void ) const;
  iterator end ( void ) const;
  size_type size ( void ) const ;

  // Iteration methods
  virtual iterator parent ( iterator it ) const = 0;
  virtual iterator left ( iterator it ) const = 0;
  virtual iterator right ( iterator it ) const = 0;
  
  // Query methods
  virtual bool isLeft ( iterator it ) const = 0;
  virtual bool isRight ( iterator it ) const = 0;
  virtual bool isLeaf ( iterator it ) const = 0;
  size_type depth ( iterator it ) const;
 
  // Builder Methods
  //virtual void subdivide ( void ) = 0;
  virtual CompressedTree * subtree ( const std::deque < Tree::iterator > & leaves ) const;
  virtual void assign ( std::shared_ptr<const CompressedTree> compressed ) = 0;


/** Tree::join
 *     InputIterator should be passed as begin() and end() of a container
 *     of pointers to Trees. (Either raw pointers or smart pointers will work.)
 *     The result is the CompressedTree object of the join of all trees.
 *     The join is defined to be the smallest binary tree for which each of the 
 *     input trees is a subtree. 
 *     There is some complexity here due to the valid_sequence complication.
 *     Since CompressedTree represents a binary tree via the "leaf_sequence" of
 *     a full binary tree and a valid_sequence over the leaves, what we need to do
 *     is create the full binary tree that is the join of the "joinands", and
 *     each leaf will be valid if it is valid for at least one of the joinands.
 */
  template < class InputIterator >
  static CompressedTree * join ( InputIterator start, InputIterator stop );
  
  // Test and Debug
  virtual uint64_t memory ( void ) const = 0;
  
protected:
  uint64_t size_;
  
private:
  // Serialization Methods
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version) {
    ar & size_;
  }
};

inline Tree::Tree ( void ) {
  size_ = 1;
}

inline Tree::~Tree ( void ) {}

inline Tree::iterator Tree::begin ( void ) const {
  return iterator ( 0 );
}

inline Tree::iterator Tree::end ( void ) const {
  return iterator ( size() );
}

inline Tree::size_type Tree::size ( void ) const {
  return size_;
}

inline Tree::size_type Tree::depth ( iterator it ) const {
  size_type result = 0;
  while ( it != begin () ) {
    ++ result;
    it = parent ( it );
  }
  return result;
}


inline CompressedTree * 
Tree::subtree ( const std::deque < iterator > & leaves ) const {

  const bool LEAF = false;
  const bool NOT_A_LEAF = true;
  const bool NOT_VALID = false;
  const bool VALID = true;

  CompressedTree * result = new CompressedTree;
  CompressedTree & new_tree = * result;
  std::vector < bool > & new_leaf_sequence
    = new_tree . leaf_sequence;
  std::vector < bool > & new_valid_sequence
    = new_tree . valid_sequence;

  if ( leaves . empty () ) {
    result -> leaf_sequence . push_back ( false );
    result -> valid_sequence . push_back ( false ); 
    return result;
  }
  
  // Mark the subtree
  std::vector < bool > visited ( size (), false );
  visited [ 0 ] = true;
  BOOST_FOREACH ( iterator leaf, leaves ) {
    iterator it = leaf;
    while ( visited [ * it ] == false ) {
      visited [ * it ] = true;
      it = parent ( it );
    }
  }
 
  // Now walk through visited nodes and write leaf_sequence and valid_sequence 
  // Possible optimization of following not yet implemented:
  //  We could avoid calling "left" and "right" fewer times on the second visit
  //  by encoding whether the children exist in the "visit" variable on the work stack


  //iterator end_it = end ();

  std::stack < iterator > DFS;
  DFS . push ( begin () );
  while ( not DFS . empty () ) {
    iterator node = DFS . top ();
    DFS . pop ();
    if ( node == end () ) {
      // Invalid leaf.
      new_leaf_sequence . push_back ( LEAF );
      new_valid_sequence . push_back ( NOT_VALID ); 
      continue;
    }
    Tree::iterator L = left ( node );
    Tree::iterator R = right ( node );
    if ( L != end () || R != end () ) {
      // Not a leaf.
      new_leaf_sequence . push_back ( NOT_A_LEAF );
      if ( R != end () && not visited [ *R ] ) R = end ();
      if ( L != end () && not visited [ *L ] ) L = end ();
      DFS . push ( R );
      DFS . push ( L );
    } else {
      // Valid Leaf
      new_leaf_sequence . push_back ( LEAF );
      new_valid_sequence . push_back ( VALID );
    }
  }  
  return result;
}

template < class InputIterator >
CompressedTree * Tree::join ( InputIterator start, InputIterator stop ) {
  typedef Tree * TreePtr;
  typedef std::shared_ptr<CompressedTree> CompressedTreePtr;

  CompressedTree * result = new CompressedTree;
  std::vector<bool> & leaf_sequence = result -> leaf_sequence;
  std::vector<bool> & valid_sequence = result -> valid_sequence;
  
  bool trivial_case = true; // All trees are empty.
  for ( InputIterator it = start; it != stop; ++ it ) {
    if ( it -> second -> leafCount () != 0 ) trivial_case = false;
  }
  if ( trivial_case ) {
    leaf_sequence . push_back ( false );
    valid_sequence . push_back ( false ); 
    return result;
  }
  //std::cout << "(1";

  leaf_sequence . push_back ( true );
  typedef Tree::iterator iterator;
  
  // How this works:
  // We want to advance through the trees simultaneously, but they aren't all the same tree
  // If we explore a subtree in some trees that does not exist in others, we remain halted on the others
  // until the subtree finishes.
  //
  // initialize iterators
  
  // State machine:
  // State 0: Try to go left. If can't, set success to false and try to go right. Otherwise success is true and try to go left on next iteration.
  // State 1: Try to go right. If can't, set success to false and rise. Otherwise success is true and try to go left on next iteration.
  // State 2: Rise. If rising from the right, rise again on next iteration. Otherwise try to go right on the next iteration.
  std::vector < iterator > iterators;
  std::vector < TreePtr > trees;
  // As we traverse through the join of trees, we maintain a data structure
  // which tells us the maximum depth a tree has a node on the current path to root
  std::vector< boost::unordered_set < uint64_t > > trees_by_depth ( 1 );
  boost::unordered_set < uint64_t > stalled_valid_trees;
  std::vector< CompressedTreePtr > compressed_trees;
  //std::vector< int64_t > leaf_seq_positions; //debug
  std::vector< int64_t > valid_seq_positions;
  std::vector < bool > left_child_missing;
  std::stack<uint64_t> path_to_root;

  // DEBUG BEGIN
  //uint64_t valid_count = 0;
  //uint64_t invalid_count = 0;
  // DEBUG END

  uint64_t current = 0;

  for ( InputIterator it = start; it != stop; ++ it ) {
    std::pair<TreePtr,CompressedTreePtr> tree_compressed_pair = *it;

    if ( tree_compressed_pair . second -> leafCount () == 0 ) continue;
    trees_by_depth [ 0 ] . insert ( trees . size () );
    trees . push_back ( tree_compressed_pair . first );
    iterators . push_back ( tree_compressed_pair . first -> begin () );
    compressed_trees . push_back ( tree_compressed_pair . second );
    //leaf_seq_positions . push_back ( 0 ); //debug
    valid_seq_positions . push_back ( -1 );
    left_child_missing . push_back ( false );
  }

  const bool LEAF = false;
  const bool NOT_A_LEAF = true;
  const bool VALID = true;
  const bool NOT_VALID = false;

  // DEBUG BEGIN
  /*
  for ( int i = 0; i < compressed_trees . size (); ++ i ) {
    // We try to determine if the compressed_tree matches the actual tree.
    // So far I check for everything but an interior node having two 
    // invalid leaves.
    const Tree & tree = * trees [ i ];
    const CompressedTree & compressed = * compressed_trees [ i ];
    uint64_t leaf_index = 0;
    uint64_t tree_index = 0;
    std::stack<Tree::iterator> DFS;
    DFS . push ( 0 );
    while ( not DFS . empty () ) {
      Tree::iterator node = DFS . top ();
      DFS . pop ();
      if ( node == tree . end () ) {
        // Invalid leaf.
        if ( compressed . valid_sequence [ leaf_index ] == VALID ) {
           throw std::logic_error ( "Inconsistent representations. (A)\n" );
        }
        ++ leaf_index;
        ++ tree_index;
        continue;
      }
      Tree::iterator L = tree . left ( node );
      Tree::iterator R = tree . right ( node );
      if ( L != tree . end () || R != tree . end () ) {
        // Not a leaf.
        if ( compressed . leaf_sequence [ tree_index ] == LEAF ) {
           throw std::logic_error ( "Inconsistent representations. (B)\n" );
        }
        ++ tree_index;
        DFS . push ( R );
        DFS . push ( L );
      } else {
        // Leaf
        if ( compressed . leaf_sequence [ tree_index ] != LEAF ) {
           throw std::logic_error ( "Inconsistent representations. (C)\n" );
        }   
        ++ leaf_index;
        ++ tree_index;
      }
    }
  }
  */
  // DEBUG END
  int64_t depth = 0;
  int state = 0;  
  while ( 1 ) {
    if ( (depth == 0) && ( state == 2 ) ) break;
    //std::cout << "Tree::join Position 0. depth = " << depth << " and state = " << state << "\n";
    bool success = false;
    boost::unordered_set < uint64_t > current_trees = trees_by_depth [ depth ]; // copy required
    // DEBUG BEGIN
    // if ( trees_by_depth . size () > depth + 1 ) {
    //   if ( trees_by_depth [ depth + 1 ] . size () != 0 ) {
    //     throw std::logic_error ( "Tree::join. Left a tree behind.\n" );
    //   }
    // }
    // DEBUG END
    BOOST_FOREACH ( uint64_t i, current_trees ) {
      const Tree & tree = * trees [ i ];

      // DEBUG BEGIN
      /*
      bool verbose = false;
      //if ( i == 0 ) verbose = true;
      if ( verbose ) { 
        std::cout << "Tree::join Position 1. i = " << i << ", depth = " << depth << " and state = " << state << "\n";
        std::cout << "leaf_seq_pos = " << leaf_seq_positions [ i ] << "\n";
        std::cout << "valid_seq_pos = " << valid_seq_positions [ i ] << "\n";
      }
      */
      // DEBUG END
      iterator end_it = tree . end ();
      int64_t newdepth = depth;

      //DEBUG BEGIN
      /*
      if ( tree . left ( iterators[i] ) == end_it && tree . right ( iterators [i ]) == end_it ) {
        // This is a leaf of the tree.
        if ( compressed_trees [ i ] -> 
                leaf_sequence [ leaf_seq_positions [ i ] ] != LEAF ) {
          std::cout << "Tree::join. SERIOUS PROBLEM 1. Representations "
            "do not agree for tree " << i << " at tree index " << leaf_seq_positions [ i ] << "\n";
        }
      }
      if ( tree . left ( iterators[i] ) != end_it || tree . right ( iterators [i ]) != end_it ) {
        // This is a leaf of the tree.
        if ( compressed_trees [ i ] -> 
                leaf_sequence [ leaf_seq_positions [ i ] ] == LEAF ) {
          std::cout << "Tree::join. SERIOUS PROBLEM 2. Representations "
            "do not agree for tree " << i << " at tree index " << leaf_seq_positions [ i ] << "\n";
        }
      }
      */
      //DEBUG END

      // Notes:
      // ++ leaf_seq_positions [ i ] happens in two cases:
      //     - A successful downward tree move is made (even to an invalid leaf)
      //     - An unsuccessful downward tree move is rejected on a non-leaf node
      // ++ valid_seq_position [ i ] happens when:
      //     - A left move is rejected.
      //bool is_leaf = (compressed_trees [ i ] -> 
      //          leaf_sequence [ leaf_seq_positions [ i ] ] == LEAF);
      
      iterator left = tree . left ( iterators [ i ] );
      iterator right = tree . right ( iterators [ i ] );
      bool is_leaf = ( left == tree . end() && right == tree . end () );

      // DEBUG BEGIN
      /*
      if ( verbose && is_leaf ) {
        std::cout << "Is leaf.\n";
      } 
      if ( verbose && not is_leaf ) {
        std::cout << "Is not a leaf.\n";
      }
      */
      // DEBUG END
      switch ( state ) {
        case 0: // Try to go left
        {
          //if ( verbose ) std::cout << "Attempt left.\n";
          if ( is_leaf ) {
            //if ( verbose ) std::cout << "Attempt left failed because leaf.\n";
            ++ valid_seq_positions [ i ]; // Advance to self leaf index
            break;
          }
          if ( left == end_it ) { 
            //if ( verbose ) std::cout << "Attempt left failed because missing child.\n";
            left_child_missing [ i ] = true;
            break;
          }
          left_child_missing [ i ] = false;
          //if ( verbose ) std::cout << "Attempt left success.\n";
          //++ leaf_seq_positions [ i ]; // Advance to child tree index
          iterators[i] = left;
          newdepth = depth + 1;
          success = true;
          break;
        }
        case 1: // Try to go right
        {
          //if ( verbose ) std::cout << "Attempt right.\n";

          if ( is_leaf ) {
            //if ( verbose ) std::cout << "Attempt right failed because leaf.\n";
            break;
          }
          if ( left_child_missing [ i ] ) {
            //if ( verbose ) std::cout << "Attempt right notices left child missing.\n";
            ++ valid_seq_positions [ i ];  // Advance to invalid left leaf index.
            //++ leaf_seq_positions [ i ]; // Advance to invalid left node index.
            // DEBUG BEGIN
            //std::cout << "This code is executed. (A)\n";
            //if (compressed_trees [ i ] -> 
            //    valid_sequence [ valid_seq_positions [ i ] ] == VALID) {
            //  std::logic_error ( "Tree::join. Missing left leaf is not invalid.\n" );
            //}
            //if (compressed_trees [ i ] -> 
            //    leaf_sequence [ leaf_seq_positions [ i ] ] != LEAF) {
            //  std::logic_error ( "Tree::join. Missing left leaf is not a leaf.\n" );
            //}
            // DEBUG END
          }
          if ( right == end_it ) { 
            //if ( verbose ) std::cout << "Attempt right fails due to right child missing.\n";

            // DEBUG BEGIN
            //if ( left_child_missing [ i ] ) {
            //  throw std::logic_error ( "Did not expect both children to be missing.\n");
            //}
            // DEBUG END
            //++ leaf_seq_positions [ i ]; // Advance to invalid right node index.
            ++ valid_seq_positions [ i ];// Advance to invalid right leaf index.
            // DEBUG BEGIN
            /*
            //std::cout << "This code is executed. (B)\n";
            if (compressed_trees [ i ] -> 
                valid_sequence [ valid_seq_positions [ i ] ] == VALID) {
              std::logic_error ( "Tree::join. Missing right leaf is not invalid.\n" );
            } 
            //if (compressed_trees [ i ] -> 
            //    leaf_sequence [ leaf_seq_positions [ i ] ] != LEAF) {
            //  std::logic_error ( "Tree::join. Missing right leaf is not a leaf.\n" );
            //} 
            */      
            // DEBUG END
            break;
          }
          //if ( verbose ) std::cout << "Attempt right succeeds.\n";
          //++ leaf_seq_positions [ i ]; // Advance to child tree index 
          iterators[i] = right;
          newdepth = depth + 1;
          success = true;
          break;
        }
        case 2: // Rise
        {
          if ( tree . isRight ( iterators[i] ) ) success = true;
          else left_child_missing [ i ] = false;
          //if ( success && verbose ) std::cout << "Rising from right.\n";
          //if ( not success && verbose ) std::cout << "Rising from left.\n";
          iterators[i] = tree . parent ( iterators[i] );
          newdepth = depth - 1;
          break;
        }
      }
      if ( newdepth != depth ) {
        stalled_valid_trees . erase ( i );
        trees_by_depth [ depth ] . erase ( i );
        if ( newdepth == (int64_t) trees_by_depth . size () ) { 
          trees_by_depth . push_back ( boost::unordered_set < uint64_t > () );
        }
        trees_by_depth [ newdepth ] . insert ( i );
      } else {
        if ( is_leaf ) {
          bool is_valid = (compressed_trees [ i ] -> 
                  valid_sequence [ valid_seq_positions [ i ] ] == VALID);
          if ( is_valid ) stalled_valid_trees . insert ( i );
        }
      }
    }
    //std::cout << "Tree::join Position 2. success = " << (success?"Yes":"No") << " depth = " << depth << " and state = " << state << "\n";
    switch ( state ) {
      case 0: // Tried to go left
        if ( success ) {
          path_to_root . push ( current );
          current = leaf_sequence . size ();
          leaf_sequence . push_back ( NOT_A_LEAF ); // Start as not leaf (can change later)
          // Success. Try to go left again.
          state = 0;
          ++ depth;
        } else {
          // Failure. Try to go right instead.
          state = 1;
        }
        break;
      case 1: // Tried to go right
        if ( success ) {
          path_to_root . push ( current );
          // Check if left branch exists
          if ( current == leaf_sequence . size () - 1 ) {
            // It doesn't. Make a fake leaf as placeholder.
            leaf_sequence . push_back ( LEAF );
            valid_sequence . push_back ( NOT_VALID );
            // invalid_count ++; // DEBUG
          } 
          current = leaf_sequence . size ();
          leaf_sequence . push_back ( NOT_A_LEAF );
          // Success. Try to go left now.
          state = 0;
          ++ depth;
        } else {
          // Check if left branch exists
          if ( current != leaf_sequence . size () - 1 ) {
            // It does. So make a fake leaf as placeholder for missing right.
            leaf_sequence . push_back ( LEAF );
            valid_sequence . push_back ( NOT_VALID );
            //invalid_count ++; // DEBUG
          } else {
            // This is a leaf.
            // We must determine if it is valid.
            // Recharacterize it as a leaf
            leaf_sequence . back () = LEAF;
            // BEGIN DEBUG
            /*
            valid_count += stalled_valid_trees . empty () ? 0 : 1;
            invalid_count += stalled_valid_trees . empty () ? 1 : 0;
            */
            // END DEBUG 
            valid_sequence . push_back ( stalled_valid_trees . empty () ? 0 : 1 );
          }
          // Failure. Rise.
          state = 2;
        }
        break;
      case 2: // Rose
        -- depth;
        current = path_to_root . top ();
        path_to_root . pop ();
        if ( success ) {
          // Rose from right, continue to rise
          state = 2;
        } else {
          // Rose from left, try to go right
          state = 1;
        }
        break;
    } 
  }
  // DEBUG BEGIN
  /*
  for ( int i = 0; i < compressed_trees . size (); ++ i ) {
    
    if ( compressed_trees [ i ] -> leaf_sequence . size () - 1 !=  leaf_seq_positions [ i ] ) {
      std::cout << leaf_seq_positions [ i ] << " != " << compressed_trees [ i ] -> leaf_sequence . size () - 1 << "\n";
      std::cout << " The tree has size " << trees[i] -> size() << "\n";
      std::cout << " This is tree number " << i << "\n";
      throw std::logic_error ( "Tree::join. Did not finish reading leaf sequence of one of trees.\n" );
    }
    
    if ( compressed_trees [ i ] -> valid_sequence . size () - 1 != valid_seq_positions [ i ] ) {
      throw std::logic_error ( "Tree::join. Did not finish reading validity sequence of one of trees.\n" );
    }
  }
    std::cout << "Tree::join invalid_count = " << invalid_count << " and valid_count = " << valid_count << "\n";
  */
  // DEBUG END
  return result;
}

#endif

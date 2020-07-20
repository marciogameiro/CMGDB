#ifndef CMDB_COMPRESSED_TREE_H
#define CMDB_COMPRESSED_TREE_H
//CompressedTree.h
class CompressedTree {
public:
  /// leaf_sequence
  ///    A sequence of bits describing a full binary tree in the following manner:
  ///    leaf_sequence [ i ] == 0 if and only if the ith preorder node is a leaf.
  ///    (Note: we use the convention that the root has preorder 0.)
  std::vector<bool> leaf_sequence; // defines a full binary tree (length 2N-1)
  
  /// valid_sequence
  ///    A bit vector that describes which leaves are valid in full binary tree (length N)
  ///    The leaves are indexed starting at 0 and are in the order induced from preorder
  ///    An entry of 1 means validity, 0 means not valid.
  std::vector<bool> valid_sequence; 

  /// size
  ///    Return the number of valid leaves of the tree
  size_t leafCount ( void ) const;

  /// subdivide
  ///    Update leaf_sequence and valid_sequence so that
  ///    all valid leaves become interior nodes with two leaf children.
  void subdivide ( void );
};

inline size_t 
CompressedTree::leafCount ( void ) const {
  size_t result = 0;
  size_t N = valid_sequence . size ();
  for ( size_t i = 0; i < N; ++ i ) if ( valid_sequence [ i ] ) ++ result;
  return result;
}

inline void 
CompressedTree::subdivide ( void ) {
  //std::cout << "CompressedTree::subdivide\n";
  CompressedTree new_tree;
  std::vector < bool > & new_leaf_sequence
    = new_tree . leaf_sequence;
  std::vector < bool > & new_valid_sequence
    = new_tree . valid_sequence;

  size_t M = leaf_sequence . size ();
  uint64_t leaf = 0;
  for ( size_t i = 0; i < M; ++ i ) {
    if ( leaf_sequence [ i ] ) {
      // Not a leaf. Copy.
      new_leaf_sequence . push_back ( 1 );
    } else {
      // This is a leaf. If it is valid, subdivide it
      if ( valid_sequence [ leaf ++ ] ) {
        // The leaf is valid. Subdivide it.
        new_leaf_sequence . push_back ( 1 );
        new_leaf_sequence . push_back ( 0 );
        new_leaf_sequence . push_back ( 0 );
        new_valid_sequence . push_back ( 1 );
        new_valid_sequence . push_back ( 1 );
      } else {
        // The leaf is not valid. Do not subdivide. Mark as invalid.
        new_leaf_sequence . push_back ( 0 );
        new_valid_sequence . push_back ( 0 );
      }
    }
  }

  std::swap ( leaf_sequence, new_leaf_sequence );
  std::swap ( valid_sequence, new_valid_sequence );
}

#endif

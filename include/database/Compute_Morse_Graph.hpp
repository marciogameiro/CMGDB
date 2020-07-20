/* Compute_Morse_Graph.hpp */
#ifndef _CMDP_COMPUTE_MORSE_GRAPH_HPP_
#define _CMDP_COMPUTE_MORSE_GRAPH_HPP_

#ifdef CMG_VISUALIZE
#include <boost/unordered_map.hpp>
#include "CImg.h"
using namespace cimg_library;
#endif


#include <map>
#include <stack>
#include <vector>
#include <exception>
#include "boost/foreach.hpp"
#include <memory>

#include "GraphTheory.h"
#include "join.h"
#include "MapGraph.h"

#include <ctime>

#ifdef MEMORYBOOKKEEPING
uint64_t max_grid_internal_memory = 0;
uint64_t max_grid_external_memory = 0;
#endif

#ifdef DO_CONLEY_INDEX
#include "chomp/ConleyIndex.h"
#endif
template < class Toplex, class CellContainer > 
void subdivide ( Toplex & phase_space, CellContainer & morse_set );


// Some macros for verbose output.
#ifdef CMG_VERBOSE

#define CMG_VERBOSE_PRINT(x) std::cout << x;

#define CMG_VERBOSE_START_CLOCK \
clock_t start = clock ();

#define CMG_VERBOSE_REPORT_CLOCK \
std::cout << "   Time Elapsed: " << (double)(clock() - start)/(double)CLOCKS_PER_SEC << "\n"; 

#define CMG_VERBOSE_REPORT_MORSE_SETS \
std::cout << "   Number of Morse Sets: " << morse_sets . size () << "\n";\
std::cout << "   Sizes of Morse Sets: ";\
BOOST_FOREACH ( CellContainer & morse_set, morse_sets ) \
std::cout << " " << morse_set . size ();\
std::cout << ".\n";

#endif

#ifndef CMG_VERBOSE
#define CMG_VERBOSE_PRINT(x)  if(0){std::cout << x;}
#define CMG_VERBOSE_START_CLOCK
#define CMG_VERBOSE_REPORT_CLOCK
#define CMG_VERBOSE_REPORT_MORSE_SETS
#endif


class MorseDecomposition {
public:

  // Constructor
  template < class GridPtr >
  MorseDecomposition ( GridPtr grid, int depth ) 
  : grid_ ( grid ), spurious_(false), depth_(depth) {
    if ( grid_ . get () == NULL ) {
      throw std::logic_error ( "Bad Initialization of MorseDecomposition Object\n" );  
    }
#ifdef MEMORYBOOKKEEPING
    max_grid_external_memory += grid -> memory ();
    max_grid_internal_memory = std::max( max_grid_internal_memory, grid -> memory () );
#endif
  }
  
  // Deconstructor
  ~MorseDecomposition ( void ) {
    BOOST_FOREACH ( MorseDecomposition * child, children_ ) {
      delete child;
    }
  }
  
  /// MorseDecomposition::size
  /// return size of grid
  size_t size ( void ) const { return grid_ -> size (); }
  
  /// MorseDecomposition::depth
  /// tell how deep in hierarchical decomposition we are
  size_t depth ( void ) const { return depth_; }

  /// MorseDecomposition::decomposition
  const std::vector< std::shared_ptr<Grid> > & decomposition ( void ) const {
    return decomposition_;
  }

  /// MorseDecomposition::children
  /// accessor method to return vector of MorseDecomposition * pointing to hierarchical children.
  /// note: empty until "decompose" is called.
  const std::vector < MorseDecomposition * > & children ( void ) const {
    return children_;
  }

  /// MorseDecomposition::reachability
  /// accessor method to obtain reachability_ private data member
  const std::vector < std::vector < unsigned int > > & reachability ( void ) const {
    return reachability_;
  }
  
  /// MorseDecomposition::grid
  /// accessor method to obtain grid_ shared_ptr data member
  std::shared_ptr < Grid > grid ( void ) {
    return grid_;
  }
  
  /// MorseDecomposition::spurious
  /// accessor method to obtain spurious_ data member
  bool & spurious ( void ) { return spurious_; }
  
  /// MorseDecomposition::decompose
  ///
  /// Use graph theory to find SCC components, which are stored as type Grid
  /// Fill these into "decomposition_"
  /// Fill children_ with an equal sized vector of pointers to new MorseDecomposition objects seeded with those sets.
  /// Put reachability information obtained in "reachability_"
  void 
  decompose ( std::shared_ptr<const Map> f ) {
    //std::cout << "decompose at depth " << depth () << "\n";
    computeMorseSetsAndReachability
      ( &decomposition_, 
        &reachability_, 
        grid_, 
        f );    
    //std::cout << "  found " << decomposition_ . size () << " components\n";
  }

  const std::vector < MorseDecomposition * > & 
  spawn ( void ) {
    //std::cout << "spawn at depth " << depth () << "\n";
    for ( size_t i = 0; i < decomposition_ . size (); ++ i ) {      
      children_ . push_back ( new MorseDecomposition ( decomposition_ [ i ] -> clone (), 
                                                       depth() + 1 ) );
    }
    //std::cout << "  spawned " << children_ . size () << " children\n";

    return children_;
  } 

private:
  // Member Data
  std::shared_ptr<Grid> grid_;
  std::vector< std::shared_ptr<Grid> > decomposition_;
  std::vector < MorseDecomposition * > children_;
  std::vector < std::vector < unsigned int > > reachability_;
  bool spurious_;
  size_t depth_;
};

class MorseDecompCompare {
public:
  bool operator () ( const MorseDecomposition * lhs, const MorseDecomposition * rhs ) {
    return lhs -> size () < rhs -> size ();
  }
};

// The MorseDecomposition Tree
//  The level of subdivision of the root is whatever the initial level is, which we call 0.
//  The level of subdivision of an internal node in the tree is equal to its depth
//  The level of subdivision of a leaf node in the tree is the same as the level of
//      subdivision of its parent.
//  The children of a node correspond to its Morse Sets (although they may be subdivided).
//  Algorithmically, this means we call decompose whenever the depth <= the number of
//  subdivisions we want.
// ConstructMorseDecomposition
inline void
ConstructMorseDecomposition (MorseDecomposition * root,
                             std::shared_ptr<const Map> f,
                             const unsigned int Min,
                             const unsigned int Max,
                             const unsigned int Limit ) {
  size_t nodes_processed = 0;
  // We use a priority queue in order to do the more difficult computations first.
  std::priority_queue < MorseDecomposition *, 
                        std::vector<MorseDecomposition *>, 
                        MorseDecompCompare > pq;
  pq . push ( root );
  while ( not pq . empty () ) {
    ++ nodes_processed;
    if ( nodes_processed % 1000 == 0 ) { 
      std::cout << nodes_processed 
        << " nodes have been encountered on Morse Decomposition Hierarchy.\n";
    }
    MorseDecomposition * work_node = pq . top ();
    pq . pop ();
    //std::cout << "Depth " << work_node -> depth () << ", node " << work_node 
    //          << ", size = " << work_node -> size () << "\n";

    // Do not decompose if past Min depth and over the Limit size.
    if ( ( work_node -> depth () > Min ) 
         && ( work_node -> size () > Limit ) ) {
      //std::cout << "Halting search due to Limit.\n";
      continue;
    }

    work_node -> decompose ( f );

    // Check for spuriousness
    if ( work_node -> decomposition ()  . empty () ) {
      //std::cout << "Empty decomposition for " << work_node << ", marking as spurious.\n";
      work_node -> spurious () = true;
    }

    // Hierarchical Step
    if ( (work_node -> depth () < Max) ) {
      std::vector < MorseDecomposition * > children = work_node -> spawn ();
      BOOST_FOREACH ( MorseDecomposition * child, children ) {
        child -> grid () -> subdivide ();
        pq . push ( child );
      }
    } 
    //else {
      //std::cout << "Halting search due to Max.\n";
    //}
  }
}

// ConstructMorseDecomposition
inline void 
ConstructMorseGraph (std::shared_ptr<Grid> master_grid,
                     MorseGraph * MG,
                     MorseDecomposition * root,
                     const unsigned int Min ) {
  std::vector < std::shared_ptr < Grid > > grids;
  // Produce Morse Graph
  typedef MorseGraph::Vertex Vertex;
  // "temp" will store which MorseGraph vertices are hierarchically under a given decomposition node
  std::map < MorseDecomposition *, std::vector<Vertex> > temp;
  std::stack < std::pair < MorseDecomposition *, size_t > > eulertourstack;
  // eulertourstack: an item (md, childnum) means "this is node md, explore child childnum if it exists,
  // otherwise do an analysis"
  eulertourstack . push ( std::make_pair( root, 0 ) );
  // The following loop performs an Euler tour and does an operation on each postordering.
  // It is thus a guarantee that the operation has been performed on all descendants in the 
  while (  not eulertourstack . empty () ) {
    MorseDecomposition * MD = eulertourstack . top () . first;
    size_t childnum = eulertourstack . top () . second;
    eulertourstack . pop ();
    size_t NC = MD -> children () . size ();
    if ( childnum < NC ) {
        eulertourstack . push ( std::make_pair ( MD, childnum + 1 ) );
        eulertourstack . push ( std::make_pair ( MD -> children () [ childnum ], 0 ) );
    } else {
      // Post-ordering operation

      // Check for Spuriousness
      // If it has children that are all marked spurious, then it is spurious.
      // If it does not have children, it is spurious if and only if it is already marked spurious
      if ( NC > 0 ) {
        MD -> spurious () = true; // by default; we may however switch it back to false
        for ( size_t i = 0; i < NC; ++ i ) {
          if ( not MD -> children () [ i ] -> spurious () ) MD -> spurious () = false;
        }
      }
      temp [ MD ] = std::vector < Vertex > ();

      // If spurious, then change grid to join of all descendants and previous self
      if ( MD -> spurious () && NC > 0) {
        // We alter MD -> grid () so that it is the join of all descendant grids
        std::vector<std::shared_ptr<Grid> > grid_family;
        grid_family . push_back ( MD -> grid () );
        for ( size_t i = 0 ; i < NC; ++ i ) {
          grid_family . push_back ( MD -> children () [ i ] -> grid () );
        }
        join ( MD -> grid (), grid_family . begin(), grid_family . end () );
      }

      if ( MD -> depth () > Min ) continue; 
      grids . push_back ( MD -> grid () );
      if ( MD -> spurious () ) continue;
      
      // Case 1. Min Depth Case
      // Morse Graph Vertex Creation Step 
      // (and special case for reachability)
      if ( MD -> depth () == Min ) {
        //std::cout << "Node " << MD << "\n";
        boost::unordered_map < int, Vertex > non_spurious_decomposition;
        size_t ND = MD -> decomposition () . size ();
        for ( size_t i = 0 ; i < ND; ++ i ) {
          //std::cout << "Child " << i << " out of " << ND << "\n";
          if ( (NC == ND) && MD -> children () [ i ] -> spurious () ) { 
            //std::cout << "Spurious Rule 2 invoked, skipping child " << i << ".\n";
            grids . push_back ( MD -> children () [ i ] -> grid () ); 
            continue;
          }
          Vertex v = MG -> AddVertex ();
          //std::cout << "Adding vertex " << v << " for child " << i << "\n";
          MG -> grid ( v ) = MD -> decomposition () [ i ];
          temp [ MD ] . push_back ( v );
          non_spurious_decomposition [ i ] = v;
        }
#ifndef NO_REACHABILITY
        typedef std::pair<int, Vertex> intVertexPair;
        BOOST_FOREACH ( const intVertexPair & ivp, non_spurious_decomposition ) {
          int i = ivp . first;
          Vertex u = ivp . second;
          const std::vector < unsigned int > & reaches = MD -> reachability () [ i ];
          BOOST_FOREACH ( unsigned int j, reaches ) {
            if ( (int) i == (int) j ) continue;
            if ( non_spurious_decomposition . count ( j ) == 0 ) continue;
            Vertex v = non_spurious_decomposition [ j ];
            MG -> AddEdge ( u, v );
            //std::cout << "(A) Adding edge " << u << " " << v << "\n";
          }
        }
#endif
      } 
      // Case 2. Less than Min Depth Case
      // temp [MD] creation step
      if ( MD -> depth () < Min ) {
      // Intermediate MD node reachability step
#ifndef NO_REACHABILITY
        for ( unsigned int i = 0; i < NC; ++ i ) {
          const std::vector < unsigned int > & reaches = MD -> reachability () [ i ];
          BOOST_FOREACH ( unsigned int j, reaches ) {
           if ( i == j ) continue;
           BOOST_FOREACH ( Vertex u, temp [ MD -> children () [ i ] ] ) {
              BOOST_FOREACH ( Vertex v, temp [ MD -> children () [ j ] ] ) {
                MG -> AddEdge ( u, v );
                //std::cout << "(B) Adding edge " << u << " " << v << "\n";
              }
            }
          }
        }
#endif
        // Create temp for intermediate
        for ( unsigned int i = 0; i < NC; ++ i ) {
          temp [ MD ] . insert (temp [ MD ] . begin (),
                                temp [ MD -> children () [ i ] ] . begin (),
                                temp [ MD -> children () [ i ] ] . end ());
          temp . erase ( MD -> children () [ i ] );
        }
      }
    } 
  }
  join ( master_grid, grids . begin (), grids . end () );
  MG -> phaseSpace () = master_grid;
}


inline void 
Compute_Morse_Graph (MorseGraph * MG,
                     std::shared_ptr<Grid> phase_space,
                     std::shared_ptr<const Map> f,
                     const unsigned int Init,
                     const unsigned int Min, 
                     const unsigned int Max, 
                     const unsigned int Limit) {
  for ( int i = 0; i < (int)Init; ++ i ) phase_space -> subdivide ();
  Compute_Morse_Graph ( MG, phase_space, f, Min - Init, Max - Init, Limit );
}

inline void 
Compute_Morse_Graph (MorseGraph * MG,
                     std::shared_ptr<Grid> phase_space,
                     std::shared_ptr<const Map> f,
                     const unsigned int Min, 
                     const unsigned int Max, 
                     const unsigned int Limit) {
  // Produce Morse Set Decomposition Hierarchy
  std::cout << "Compute_Morse_Graph. Initializing root MorseDecomposition\n";
  std::cout << "Compute_Morse_Graph. A phase_space -> size () == " << phase_space -> size () << "\n";

  std::shared_ptr<Grid> root_space ( (Grid *) (phase_space -> clone ()) );
  
  //std::cout << "Compute_Morse_Graph. root_space -> size () == " << root_space -> size () << "\n";
  
  MorseDecomposition * root = new MorseDecomposition ( root_space, 0 );
  
  //std::cout << "Compute_Morse_Graph. Calling ConstructMorseDecomposition\n";
  
  ConstructMorseDecomposition (root,
                               f,
                               Min,
                               Max,
                               Limit);
  //std::cout << "Calling ConstructMorseGraph\n";
  // Stitch together Morse Graph from Decomposition Hierarchy
  ConstructMorseGraph ( phase_space, MG, root, Min );

  std::cout << "Compute_Morse_Graph. B phase_space -> size () == " << phase_space -> size () << "\n";

  // Free memory used in decomposition hierarchy
  delete root;
#ifdef MEMORYBOOKKEEPING

  std::cout << "Total Grid Memory (can be external) = " << max_grid_external_memory << "\n";
  std::cout << "Max Memory For Single Grid (must be internal)= " << max_grid_internal_memory << "\n";
  std::cout << "Max SCC Random Access memory use (must be internal)= " << max_scc_memory_internal << "\n";
  std::cout << "Max SCC stack memory use (can be external memory) = " << max_scc_memory_external << "\n";
  std::cout << " ---- SUMMARY ---- \n";
  std::cout << "Internal Memory Requirement = " << max_grid_internal_memory + max_scc_memory_internal << "\n";
  std::cout << "External Memory Requirement = " << max_grid_external_memory + max_scc_memory_external << "\n";
  std::cout << "Max graph memory size (never stored, however) = " << max_graph_memory << "\n";
#endif
  //std::cout << "Returning from COMPUTE MORSE GRAPH\n";
}

#endif

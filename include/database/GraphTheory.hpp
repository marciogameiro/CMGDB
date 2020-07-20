// GraphTheory.hpp Shaun Harker 
// created 5/16/2011
// updated 5/27/2014

#include <vector>
#include <stack>
#include <queue>
#include <algorithm>
#include <memory>
#include "boost/unordered_set.hpp"
#include "boost/unordered_map.hpp"
#include "boost/foreach.hpp"
#include "MapGraph.h"

#define DEBUGPRINT if(0)

#ifdef MEMORYBOOKKEEPING
uint64_t max_scc_memory_internal = 0;
uint64_t max_scc_memory_external = 0;
uint64_t max_reach_memory = 0;
uint64_t graph_memory = 0;
uint64_t max_graph_memory = 0;
#endif

inline void 
computeMorseSetsAndReachability (std::vector< std::shared_ptr<Grid> > * output,
                                 std::vector<std::vector<unsigned int> > * reach,
                                 std::shared_ptr<const Grid> G,
                                 std::shared_ptr<const Map> f ) {
  MapGraph mapgraph ( G, f );
  // Produce Strong Components and Reachability
  std::vector < std::deque < Grid::GridElement > > components;
  std::deque < Grid::size_type > topological_sort;
  computeStrongComponents ( &components, mapgraph, &topological_sort );
#ifdef CMG_VERBOSE
  if ( components . size () > 1 ) {
    std::cout << "Found " << components . size () 
              << " combinatorial Morse sets.\n";
  }
#endif
#ifndef NO_REACHABILITY
  computeReachability ( reach, components, mapgraph, topological_sort );
#endif
  // Create output grids
  output -> clear ();
  BOOST_FOREACH ( const std::deque<Grid::GridElement> & component, components ) {
    std::shared_ptr < Grid > component_grid ( G -> subgrid ( component ) );
    output -> push_back ( component_grid );
  }
}

/// computeStrongComponents (actually, SCPCs... needs renaming.)
///    Modified version of Tarjan's algorithm devised by Shaun Harker
///    Only calls for adjacency lists once each, yet only requires O(V) space.
///
/// comments: the interface results in some inefficiency. It would be nice to change
///           it. In particular the way SCC_root requires it to be written in random
///           order which is not ideal.
template < class Graph >
void computeStrongComponents (std::vector<std::deque<typename Graph::Vertex> > * output,
                              const Graph & G,
         /* optional output */std::deque<typename Graph::Vertex> * topological_sort,
         /* optional output */std::deque<typename Graph::Vertex> * SCC_root ) {
  typedef typename Graph::Vertex Vertex;
#ifdef CMG_VERBOSE
  int64_t progress = 0;
  int64_t progresspercent = 0;
#endif
#ifdef MEMORYBOOKKEEPING
  graph_memory = 0;
#endif
  int64_t E = 0;
  int64_t N = (int64_t) G . num_vertices ();
  std::vector<bool> explored ( N, false );
  std::vector<bool> committed ( N, false );
  std::vector<bool> duplicates ( N, false );
  std::vector<bool> self_connected (N, false);
  std::vector<int64_t> preorder ( N, 0 );
  std::deque<int64_t> LOWLINK, DFS, cleanDFS;
  std::deque<Vertex> S;
  if ( SCC_root != NULL ) {
    SCC_root -> resize ( N );
  }
  int64_t n = 0;
  LOWLINK . push_back ( -1 );
  for ( int64_t v = 0; v < N; ++ v ) {
    if ( not explored [ v ] ) {
      DFS . push_back ( v+1 );
      while ( not DFS . empty () ) {
        int64_t u = DFS . back ();
        DFS . pop_back ();
#ifdef CMG_VERBOSE
        if ( (100*progress)/(2L*N) > progresspercent) {
          progresspercent = (100*progress)/(2L*N);
          std::cout << "\rcomputeStrongComponents. V = " << N << ", " << progresspercent << "%" " finished.  ";
          std::cout . flush ();
        }
#endif
#ifdef MEMORYBOOKKEEPING
        uint64_t mem_bytes_external = sizeof ( int64_t ) * DFS . size () + 
                                      sizeof ( int64_t ) * LOWLINK . size () +
                                      sizeof ( Vertex ) * S . size (); 
        uint64_t mem_bytes_internal = sizeof ( int64_t ) * preorder . size ()  +
                                      N / 2; // For std::vector<bool> info
        max_scc_memory_external = std::max( max_scc_memory_external, mem_bytes_external );
        max_scc_memory_internal = std::max( max_scc_memory_internal, mem_bytes_internal );
#endif
        if ( u > 0 ) {
          // PREORDER
          u = u - 1;
          //std::cout << "Preorder(" << u << ")\n";
          if ( not explored [ u ] ) {
#ifdef CMG_VERBOSE
            ++ progress;
#endif
            DFS . push_back ( -u-1);
            explored [ u ] = true;
            preorder [ u ] = n;
            int64_t low = n;
            ++ n;
            std::vector<Vertex> W =  G . adjacencies ( u );
#ifdef MEMORYBOOKKEEPING
            graph_memory += sizeof(Vertex) * (1 + W . size ());
#endif
            E += W . size ();
            BOOST_FOREACH ( int64_t w, W ) {
              if ( u == w ) self_connected [ u ] = true;
              if ( explored [ w ] ) {
                if ( not committed [ w ] ) {
                  low = std::min(low, preorder[w] );
                }
              } else {
                DFS . push_back ( w + 1 );
                if ( (int64_t) DFS . size () > 2L*N ) { 
                  duplicates . assign ( N, false );
                  while ( not DFS . empty () ) {
                    int64_t x = DFS . back ();
                    DFS . pop_back ();
                    int64_t y = std::abs(x) - 1;
                    if ( duplicates [ y ] ) continue;
                    duplicates [ y ] = true;
                    cleanDFS . push_front ( x );
                  }
                  std::swap ( DFS, cleanDFS );
                }
              }
            }
            LOWLINK . push_back ( low );
            S . push_back ( (Vertex) u );
          }
        } else {
          // POSTORDER
          u = -u - 1;
#ifdef CMG_VERBOSE
          ++ progress;
#endif
          //std::cout << "Postorder(" << u << ")\n";
          int64_t lowlink = LOWLINK . back ();
          LOWLINK . pop_back ();
          if ( lowlink == preorder [ u ] ) {
            output -> push_back ( std::deque<Vertex> () );
            std::deque<Vertex> & SCC = output -> back ();
            do {
              int64_t w = S . back ();
              S . pop_back ();
              SCC . push_back ( (Vertex) w );
              committed [ w ] = true;
              if ( topological_sort != NULL ) {
                topological_sort -> push_back ( w );
              }
              if ( SCC_root != NULL ) {
                (*SCC_root) [ w ] = ( (Vertex) u );
              }
            } while ( not committed [ u ] );
            // Only SCPCs:
            if ( SCC . size () == 1 &&
                 not self_connected [ u ] ) output -> pop_back ();
          }
          int64_t low = LOWLINK . back ();
          LOWLINK . pop_back ();
          lowlink = std::min ( lowlink, low );
          LOWLINK . push_back ( lowlink );
        }
      }
    }
  }
#ifdef CMG_VERBOSE
  std::cout << "\rcomputeStrongComponents. V = " << N << " E = " 
            << E << "  E/V = " << (double) E / (double) N << "\n";
#endif
#ifdef MEMORYBOOKKEEPING
  max_graph_memory = std::max(max_graph_memory, graph_memory );
#endif
}

/// computeReachability 
template < class Graph >
void computeReachability ( std::vector < std::vector < unsigned int > > * output,
                           std::vector<std::deque<typename Graph::size_type> > & morse_sets,
                           const Graph & G, 
                           const std::deque<typename Graph::size_type> & topological_sort ) {
  typedef typename Graph::size_type size_type;
#ifdef CMG_VERBOSE
  std::cout << "Computing Reachability Information.\n";
  std::cout . flush ();
  size_type progress = 0;
  size_type progresspercent = 0;
#endif
  /* Count the Morse Sets */
  uint64_t effort = 0;
  size_type number_of_morse_sets = morse_sets . size ();  
  if ( number_of_morse_sets == 0 ) return; // trivial case
  output -> resize ( number_of_morse_sets );
  // Paint the Morse Sets
  // For each morse set, go through its vertices and 
  // color them according to which morse set they are in.
  // Vertices not in a morse set are colored "number_of_morse_sets"
  std::vector < size_type > morse_paint ( G . num_vertices (), number_of_morse_sets );
  for ( size_type count = 0; count < number_of_morse_sets; ++ count ) {
    BOOST_FOREACH ( size_type v, morse_sets [ count ] ) {
      ++ effort; 
      morse_paint [ v ] = count;
    } 
  } 

  // Break the Morse Sets up into Computational Groups of 64 and proceed 
  size_type groups = ( (number_of_morse_sets - 1) / 64 ) + 1;
  
#ifdef CMG_VERBOSE
  size_type total_work_to_do = topological_sort . size () * groups;
#endif
  // We use a vector called morse_code. It's function is to maintain
  // information about which morse sets can reach a given vertex.
  // By processing in topological order, it is possible to give morse_code
  // the correct values in a single pass.
  std::vector < uint64_t > morse_code;
  // We use a vector called condensed_code in order to store the final information
  // about which morse sets can reach a given morse set. It can be inferred from morse_code
  // during the same sweep in which we construct morse_code.
  std::vector < uint64_t > condensed_code;
  // Loop through groups.
  for ( size_type group_number = 0; group_number < groups; ++ group_number ) {
    ++ effort;
    size_type group_size = std::min((size_type) 64, 
          (size_type) number_of_morse_sets - ((size_type)64) * group_number);
    size_type offset = 64L * group_number;
    morse_code . clear ();
    morse_code . resize ( G . num_vertices (), 0 );
    condensed_code . clear ();
    condensed_code . resize ( number_of_morse_sets + 1, 0 );

    // Paint the codes.
    // We do an initial sweep painting the sources onto their sets. 
    for ( size_type count = 0; count < group_size; ++ count ) {
      size_type set_number = offset + count;
      uint64_t code = ((uint64_t)1) << count;
      ++ effort;    
      BOOST_FOREACH ( size_type v, morse_sets [ set_number ] ) {
        ++ effort;
        morse_code [ v ] = code;
      } 
    } 

    // Loop through topological sort.
    // Our goal is to produce "condensed_code", which we can read the info off from.
    // The intermediate information is stored in "morse_code"
    for ( int64_t vi = topological_sort . size () - 1; vi >= 0; -- vi ) {
#ifdef CMG_VERBOSE
      ++ progress;
      if ( (100*progress)/total_work_to_do > progresspercent) {
        progresspercent = (100*progress)/total_work_to_do;
        std::cout << "\r" << progresspercent << "%    ";
        std::cout . flush ();
      }
#endif
      size_type v = topological_sort [ vi ];
      std::vector < size_type > children = G . adjacencies ( v ); // previously const &
      if ( morse_paint [ v ] != number_of_morse_sets ) {
        morse_code [ v ] |= condensed_code [ morse_paint [ v ] ];
      }
      BOOST_FOREACH ( size_type w, children ) {
        ++ effort;
        morse_code [ w ] |= morse_code [ v ];
        condensed_code [ morse_paint [ w ] ] |= morse_code [ v ];

      }
    } 
#ifdef CMG_VERBOSE
    ++ progress;
#endif

    // Note: Now condensed_code is indexed by targets, 
    //       and contains the sources reaching it.    
    // Loop through Morse Sets to learn reachability information
    for ( size_type count = 0; count < number_of_morse_sets; ++ count ) {
      // Read condensed code to determine which of the group reached this morse set
      ++ effort;
      uint64_t bit = 1;
      for ( int i = 0; i < 64; ++ i ) {
        ++ effort;     
        if ( condensed_code [ count ] & bit ) {
          ++ effort;
          (*output)[offset + i] . push_back ( count );
        } // if
        bit <<= 1;
      } // for bit index
    } // for morse set
  } // for groups
#ifdef CMG_VERBOSE
  std::cout << "\r100%  Reachability Analysis Complete.\n ";
#endif
}

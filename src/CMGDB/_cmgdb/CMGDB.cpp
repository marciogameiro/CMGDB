#include <iostream>
#include <fstream>
#include <ctime>
#include <cmath>
#include <vector>
#include <set>
#include <map>
#include <sstream>
#include <algorithm>

// #define CMG_VERBOSE
#define MEMORYBOOKKEEPING

#include "Model.h"

#include "Map.h"
#include "ChompMap.h"
#include "MorseGraph.h"
#include "Compute_Morse_Graph.h"
#include "RectGeo.h"

#include "SingleOutput.h"
#include "simple_interval.h"

#include "Configuration.h"

#include "chomp/ConleyIndex.h"
#include "conleyIndexString.h"

#include <boost/serialization/export.hpp>
// The succinct (sdsl-backed) grid is optional: the phase grid is
// PointerGrid and nothing constructs a SuccinctGrid, so default builds do
// not compile against the vendored sdsl-lite at all. Define
// CMGDB_USE_SUCCINCT to enable it.
#ifdef CMGDB_USE_SUCCINCT
#include "SuccinctGrid.h"
BOOST_CLASS_EXPORT_IMPLEMENT(SuccinctGrid);
#endif
#include "PointerGrid.h"
BOOST_CLASS_EXPORT_IMPLEMENT(PointerGrid);

std::vector < std::string >
ComputeConleyIndex ( const std::vector < uint64_t > & X_cubes,
                     const std::vector < uint64_t > & A_cubes,
                     const std::vector < uint64_t > & sizes,
                     const std::vector < bool > & periodic,
                     const std::unordered_map < uint64_t, std::vector < uint64_t > > & F,
                     bool acyclic_check = true ) {
  // Compute the Conley index from a combinatorial index pair (X, A) and a map F
  chomp::ConleyIndex_t conley_index;
  chomp::CombinatorialConleyIndex ( &conley_index, X_cubes, A_cubes, sizes, periodic, F, acyclic_check );
  // Return Conley index strings
  return conleyIndexString ( conley_index );
}

std::vector < std::string >
ComputeConleyIndexForCells (
    const Model & model,
    MorseGraph & morse_graph,
    std::vector < uint64_t > cells ) {
  // Conley index of an arbitrary cell subset of the final phase-space grid.
  // Mixed-depth (adaptive-grid) cell sets are supported: the chomp machinery
  // refines every cell to the finest depth present in the set, preserving
  // the region exactly.
  std::shared_ptr < TreeGrid > phase_space_chomp =
    std::dynamic_pointer_cast<TreeGrid> ( morse_graph . phaseSpace () );
  if ( not phase_space_chomp ) {
    throw std::runtime_error (
      "ComputeConleyIndexForCells requires a TreeGrid-backed Morse graph" );
  }

  std::sort ( cells . begin (), cells . end () );
  cells . erase ( std::unique ( cells . begin (), cells . end () ), cells . end () );
  for ( const uint64_t cell : cells ) {
    if ( cell >= phase_space_chomp -> size () ) {
      std::ostringstream message;
      message
        << "ComputeConleyIndexForCells cell " << cell
        << " is outside [0, " << phase_space_chomp -> size () << ")";
      throw std::out_of_range ( message . str () );
    }
  }

  chomp::ConleyIndex_t conley_index;
  ChompMap chomp_map ( model . map () );
  chomp::ConleyIndex (
    & conley_index, * phase_space_chomp, cells, chomp_map );
  return conleyIndexString ( conley_index );
}

std::pair<MorseGraph, MapGraph> ComputeConleyMorseGraph ( Model const& model,
                                                          bool cache_transition_graph = true,
                                                          uint64_t batch_chunk_size = 65536,
                                                          uint64_t max_cached_edges = 0,
                                                          uint64_t reserve_edges = 0,
                                                          uint64_t reserve_min_edges = uint64_t ( 1 ) << 24,
                                                          bool cache_map_graph = false ) {
  MapGraphOptions options ( cache_transition_graph, batch_chunk_size, max_cached_edges,
                            reserve_edges, reserve_min_edges );
  std::shared_ptr<const Map> map = model . map ();
  MorseGraph morsegraph ( model . phaseSpace () );
  std::shared_ptr < Grid > phase_space = morsegraph . phaseSpace ();

  int phase_subdiv_init = model . phase_subdiv_init ();
  int phase_subdiv_min = model . phase_subdiv_min ();
  int phase_subdiv_max = model . phase_subdiv_max ();
  int phase_subdiv_limit = model . phase_subdiv_limit ();

  // Compute Morse graph
  Compute_Morse_Graph ( & morsegraph, phase_space, map, phase_subdiv_init,
                        phase_subdiv_min, phase_subdiv_max, phase_subdiv_limit,
                        options );

  std::shared_ptr < TreeGrid > phase_space_chomp =
    std::dynamic_pointer_cast<TreeGrid> ( morsegraph . phaseSpace () );

  if ( not phase_space_chomp ) {
    throw std::runtime_error ( "Cannot interface with chomp for this grid type!" );
  }

  typedef std::vector < Grid::GridElement > Subset;
  for ( size_t v = 0; v < morsegraph . NumVertices (); ++ v) {
    Subset subset = phase_space_chomp -> subset ( * morsegraph . grid ( v ) );
    std::shared_ptr<chomp::ConleyIndex_t> conley ( new chomp::ConleyIndex_t );
    morsegraph . conleyIndex ( v ) = conley;
    ChompMap chomp_map ( map );
    chomp::ConleyIndex ( conley . get (), *phase_space_chomp, subset, chomp_map );
  }

  // Compute multi-valued map digraph on the final grid. With cache_map_graph
  // the returned graph is eagerly cached (one extra full batched map pass
  // over the grid, then O(1) adjacency queries); by default it is lazy and
  // evaluates the map on demand per adjacency query.
  MapGraphOptions returned_options = options;
  returned_options . cache = cache_map_graph;
  MapGraph map_graph ( phase_space, map, returned_options );

  return std::make_pair ( morsegraph, map_graph );
}

std::pair<MorseGraph, MapGraph> ComputeMorseGraph ( Model const& model,
                                                    bool cache_transition_graph = true,
                                                    uint64_t batch_chunk_size = 65536,
                                                    uint64_t max_cached_edges = 0,
                                                    uint64_t reserve_edges = 0,
                                                    uint64_t reserve_min_edges = uint64_t ( 1 ) << 24,
                                                    bool cache_map_graph = false ) {
  MapGraphOptions options ( cache_transition_graph, batch_chunk_size, max_cached_edges,
                            reserve_edges, reserve_min_edges );
  std::shared_ptr<const Map> map = model . map ();
  MorseGraph morsegraph ( model . phaseSpace () );
  std::shared_ptr < Grid > phase_space = morsegraph . phaseSpace ();

  int phase_subdiv_init = model . phase_subdiv_init ();
  int phase_subdiv_min = model . phase_subdiv_min ();
  int phase_subdiv_max = model . phase_subdiv_max ();
  int phase_subdiv_limit = model . phase_subdiv_limit ();

  // Compute Morse graph
  Compute_Morse_Graph ( & morsegraph, phase_space, map, phase_subdiv_init,
                        phase_subdiv_min, phase_subdiv_max, phase_subdiv_limit,
                        options );

  // Compute multi-valued map digraph on the final grid. With cache_map_graph
  // the returned graph is eagerly cached (one extra full batched map pass
  // over the grid, then O(1) adjacency queries); by default it is lazy and
  // evaluates the map on demand per adjacency query.
  MapGraphOptions returned_options = options;
  returned_options . cache = cache_map_graph;
  MapGraph map_graph ( phase_space, map, returned_options );

  return std::make_pair ( morsegraph, map_graph );
}

std::vector<uint64_t>
ComputeMorseDirectedPathCells (
    const MapGraph & map_graph,
    const MorseGraph & morse_graph,
    const std::vector<uint64_t> & source_nodes,
    const std::vector<uint64_t> & target_nodes ) {
  if ( not map_graph . has_cache () ) {
    throw std::runtime_error (
      "MorseDirectedPathCells requires a cached MapGraph; refusing to use "
      "on-demand map callbacks." );
  }

  const uint64_t n = map_graph . num_vertices ();
  if ( n > std::numeric_limits<uint32_t>::max () ) {
    throw std::runtime_error (
      "MorseDirectedPathCells currently supports at most 2^32-1 map vertices" );
  }
  const size_t number_of_morse_sets = morse_graph . NumVertices ();
  if ( source_nodes . empty () or target_nodes . empty () ) {
    throw std::invalid_argument (
      "MorseDirectedPathCells requires nonempty source_nodes and target_nodes" );
  }

  std::vector<uint8_t> is_source_node ( number_of_morse_sets, 0 );
  std::vector<uint8_t> can_reach_target_node ( number_of_morse_sets, 0 );
  for ( const uint64_t node : source_nodes ) {
    if ( node >= number_of_morse_sets ) {
      std::ostringstream message;
      message
        << "MorseDirectedPathCells source node " << node
        << " is outside [0, " << number_of_morse_sets << ")";
      throw std::out_of_range ( message . str () );
    }
    is_source_node [ node ] = 1;
  }
  for ( const uint64_t node : target_nodes ) {
    if ( node >= number_of_morse_sets ) {
      std::ostringstream message;
      message
        << "MorseDirectedPathCells target node " << node
        << " is outside [0, " << number_of_morse_sets << ")";
      throw std::out_of_range ( message . str () );
    }
    can_reach_target_node [ node ] = 1;
  }

  // Recurrent cells are terminals for the backward dynamic program. Seed each
  // fine Morse component by whether its node can reach a requested target in
  // the Morse DAG. This preserves all downstream paths while breaking every
  // cell-level directed cycle at its recurrent component.
  const std::vector<std::pair<uint64_t, uint64_t>> morse_edges =
    morse_graph . edges_unreduced ();
  for ( size_t pass = 0; pass < number_of_morse_sets; ++ pass ) {
    bool changed = false;
    for ( const auto & edge : morse_edges ) {
      if ( edge . first >= number_of_morse_sets or
           edge . second >= number_of_morse_sets ) {
        throw std::runtime_error (
          "MorseDirectedPathCells received an invalid Morse-graph edge" );
      }
      if ( can_reach_target_node [ edge . second ] and
           not can_reach_target_node [ edge . first ] ) {
        can_reach_target_node [ edge . first ] = 1;
        changed = true;
      }
    }
    if ( not changed ) break;
  }

  std::vector<uint8_t> forward ( static_cast<size_t> ( n ), 0 );
  std::vector<uint32_t> vertex_stack;
  for ( size_t node = 0; node < number_of_morse_sets; ++ node ) {
    if ( not is_source_node [ node ] ) continue;
    const std::vector<uint64_t> cells = morse_graph . morse_set ( node );
    for ( const uint64_t cell : cells ) {
      if ( cell >= n ) {
        throw std::runtime_error (
          "MorseDirectedPathCells found a source Morse cell outside the "
          "MapGraph" );
      }
      if ( not forward [ cell ] ) {
        forward [ cell ] = 1;
        vertex_stack . push_back ( static_cast<uint32_t> ( cell ) );
      }
    }
  }

  while ( not vertex_stack . empty () ) {
    const uint32_t source = vertex_stack . back ();
    vertex_stack . pop_back ();
    const MapGraph::AdjacencySpan adjacency =
      map_graph . adjacency_span ( source );
    for ( const uint64_t successor : adjacency ) {
      if ( successor >= n ) {
        throw std::runtime_error (
          "MorseDirectedPathCells found an adjacency outside the MapGraph" );
      }
      if ( not forward [ successor ] ) {
        forward [ successor ] = 1;
        vertex_stack . push_back ( static_cast<uint32_t> ( successor ) );
      }
    }
  }

  enum VisitState : uint8_t { UNSEEN = 0, ACTIVE = 1, DONE = 2 };
  std::vector<uint8_t> state ( static_cast<size_t> ( n ), UNSEEN );
  std::vector<uint8_t> can_reach_target ( static_cast<size_t> ( n ), 0 );
  for ( size_t node = 0; node < number_of_morse_sets; ++ node ) {
    const std::vector<uint64_t> cells = morse_graph . morse_set ( node );
    for ( const uint64_t cell : cells ) {
      if ( cell >= n ) {
        throw std::runtime_error (
          "MorseDirectedPathCells found a recurrent cell outside the MapGraph" );
      }
      state [ cell ] = DONE;
      can_reach_target [ cell ] = can_reach_target_node [ node ];
    }
  }

  struct Frame {
    uint32_t vertex;
    uint32_t next_adjacency;
  };
  std::vector<Frame> stack;
  for ( uint64_t raw_vertex = 0; raw_vertex < n; ++ raw_vertex ) {
    if ( not forward [ raw_vertex ] or state [ raw_vertex ] != UNSEEN ) continue;
    state [ raw_vertex ] = ACTIVE;
    stack . push_back ( { static_cast<uint32_t> ( raw_vertex ), 0 } );
    while ( not stack . empty () ) {
      Frame & frame = stack . back ();
      const MapGraph::AdjacencySpan adjacency =
        map_graph . adjacency_span ( frame . vertex );
      if ( adjacency . size () > std::numeric_limits<uint32_t>::max () ) {
        throw std::runtime_error (
          "MorseDirectedPathCells found more than 2^32-1 outgoing edges "
          "at one vertex" );
      }
      if ( frame . next_adjacency == adjacency . size () ) {
        state [ frame . vertex ] = DONE;
        stack . pop_back ();
        continue;
      }

      const uint64_t successor =
        adjacency . begin () [ frame . next_adjacency ];
      if ( successor >= n ) {
        throw std::runtime_error (
          "MorseDirectedPathCells found an adjacency outside the MapGraph" );
      }
      if ( state [ successor ] == UNSEEN ) {
        state [ successor ] = ACTIVE;
        stack . push_back (
          { static_cast<uint32_t> ( successor ), 0 } );
        continue;
      }
      if ( state [ successor ] == ACTIVE ) {
        throw std::runtime_error (
          "MorseDirectedPathCells found a directed cycle not covered by "
          "the supplied Morse sets" );
      }

      can_reach_target [ frame . vertex ] =
        can_reach_target [ frame . vertex ] or
        can_reach_target [ successor ];
      ++ frame . next_adjacency;
    }
  }

  std::vector<uint64_t> result;
  for ( uint64_t vertex = 0; vertex < n; ++ vertex ) {
    if ( forward [ vertex ] and can_reach_target [ vertex ] ) {
      result . push_back ( vertex );
    }
  }
  return result;
}

std::vector<uint64_t>
ComputeMorseReachabilityMasks (
    const MapGraph & map_graph,
    const MorseGraph & morse_graph,
    const std::vector<uint64_t> & query_vertices ) {
  if ( not map_graph . has_cache () ) {
    throw std::runtime_error (
      "MorseReachabilityMasks requires a cached MapGraph; refusing to use "
      "on-demand map callbacks." );
  }

  const uint64_t n = map_graph . num_vertices ();
  if ( n > std::numeric_limits<uint32_t>::max () ) {
    throw std::runtime_error (
      "MorseReachabilityMasks currently supports at most 2^32-1 map vertices" );
  }

  const size_t number_of_morse_sets = morse_graph . NumVertices ();
  if ( number_of_morse_sets > 64 ) {
    std::ostringstream message;
    message
      << "MorseReachabilityMasks cannot encode " << number_of_morse_sets
      << " Morse nodes in a uint64 mask";
    throw std::runtime_error ( message . str () );
  }

  for ( const uint64_t query : query_vertices ) {
    if ( query >= n ) {
      std::ostringstream message;
      message
        << "MorseReachabilityMasks query vertex " << query
        << " is outside [0, " << n << ")";
      throw std::out_of_range ( message . str () );
    }
  }

  // Bit i denotes reachability to Morse node i. Seed each recurrent component
  // with the transitive closure of its node in the Morse DAG. Once seeded,
  // recurrent cells are terminals: their closure already contains every
  // downstream Morse node, and skipping their cell-level outgoing edges breaks
  // every directed cycle in the remaining graph.
  std::vector<uint64_t> morse_masks ( number_of_morse_sets, 0 );
  for ( size_t node = 0; node < number_of_morse_sets; ++ node ) {
    morse_masks [ node ] = uint64_t ( 1 ) << node;
  }
  const std::vector<std::pair<uint64_t, uint64_t>> morse_edges =
    morse_graph . edges_unreduced ();
  for ( size_t pass = 0; pass < number_of_morse_sets; ++ pass ) {
    bool changed = false;
    for ( const auto & edge : morse_edges ) {
      if ( edge . first >= number_of_morse_sets or
           edge . second >= number_of_morse_sets ) {
        throw std::runtime_error (
          "MorseReachabilityMasks received an invalid Morse-graph edge" );
      }
      const uint64_t updated =
        morse_masks [ edge . first ] | morse_masks [ edge . second ];
      if ( updated != morse_masks [ edge . first ] ) {
        morse_masks [ edge . first ] = updated;
        changed = true;
      }
    }
    if ( not changed ) break;
  }

  enum VisitState : uint8_t { UNSEEN = 0, ACTIVE = 1, DONE = 2 };
  std::vector<uint8_t> state ( static_cast<size_t> ( n ), UNSEEN );
  std::vector<uint64_t> reach_mask ( static_cast<size_t> ( n ), 0 );

  for ( size_t node = 0; node < number_of_morse_sets; ++ node ) {
    const std::vector<uint64_t> cells = morse_graph . morse_set ( node );
    for ( const uint64_t cell : cells ) {
      if ( cell >= n ) {
        throw std::runtime_error (
          "MorseReachabilityMasks found a Morse cell outside the MapGraph" );
      }
      state [ cell ] = DONE;
      reach_mask [ cell ] |= morse_masks [ node ];
    }
  }

  struct Frame {
    uint32_t vertex;
    uint32_t next_adjacency;
  };
  std::vector<Frame> stack;
  std::vector<uint64_t> result ( query_vertices . size (), 0 );

  for ( size_t query_index = 0;
        query_index < query_vertices . size ();
        ++ query_index ) {
    const uint32_t query = static_cast<uint32_t> ( query_vertices [ query_index ] );
    if ( state [ query ] == UNSEEN ) {
      state [ query ] = ACTIVE;
      stack . push_back ( { query, 0 } );
      while ( not stack . empty () ) {
        Frame & frame = stack . back ();
        const MapGraph::AdjacencySpan adjacency =
          map_graph . adjacency_span ( frame . vertex );
        if ( adjacency . size () > std::numeric_limits<uint32_t>::max () ) {
          throw std::runtime_error (
            "MorseReachabilityMasks found more than 2^32-1 outgoing edges "
            "at one vertex" );
        }

        if ( frame . next_adjacency == adjacency . size () ) {
          state [ frame . vertex ] = DONE;
          stack . pop_back ();
          continue;
        }

        const uint64_t successor =
          adjacency . begin () [ frame . next_adjacency ];
        if ( successor >= n ) {
          throw std::runtime_error (
            "MorseReachabilityMasks found an adjacency outside the MapGraph" );
        }
        if ( state [ successor ] == UNSEEN ) {
          state [ successor ] = ACTIVE;
          stack . push_back (
            { static_cast<uint32_t> ( successor ), 0 } );
          continue;
        }
        if ( state [ successor ] == ACTIVE ) {
          throw std::runtime_error (
            "MorseReachabilityMasks found a directed cycle not covered by "
            "the supplied Morse sets" );
        }

        reach_mask [ frame . vertex ] |= reach_mask [ successor ];
        ++ frame . next_adjacency;
      }
    }
    result [ query_index ] = reach_mask [ query ];
  }

  return result;
}

std::vector<int32_t>
ComputeMorseSingletonReachability (
    const MapGraph & map_graph,
    const MorseGraph & morse_graph,
    const std::vector<uint64_t> & query_vertices ) {
  if ( not map_graph . has_cache () ) {
    throw std::runtime_error (
      "MorseSingletonReachability requires a cached MapGraph; refusing to use "
      "on-demand map callbacks." );
  }

  const uint64_t n = map_graph . num_vertices ();
  if ( n > std::numeric_limits<uint32_t>::max () ) {
    throw std::runtime_error (
      "MorseSingletonReachability currently supports at most 2^32-1 map "
      "vertices" );
  }
  for ( const uint64_t query : query_vertices ) {
    if ( query >= n ) {
      std::ostringstream message;
      message
        << "MorseSingletonReachability query vertex " << query
        << " is outside [0, " << n << ")";
      throw std::out_of_range ( message . str () );
    }
  }

  constexpr int32_t NO_MORSE_NODE = -1;
  constexpr int32_t MULTIPLE_MORSE_NODES = -2;
  const size_t number_of_morse_sets = morse_graph . NumVertices ();
  if ( number_of_morse_sets >
       static_cast<size_t> ( std::numeric_limits<int32_t>::max () ) ) {
    throw std::runtime_error (
      "MorseSingletonReachability cannot encode the Morse-node ids in int32" );
  }

  // A recurrent node is singleton-reachable exactly when it has no outgoing
  // edge to a distinct Morse node. Otherwise its reachable set already
  // contains itself plus at least one other node, so MULTIPLE is sufficient
  // for the strict singleton-basin criterion.
  std::vector<int32_t> morse_summary ( number_of_morse_sets, NO_MORSE_NODE );
  for ( size_t node = 0; node < number_of_morse_sets; ++ node ) {
    morse_summary [ node ] = static_cast<int32_t> ( node );
  }
  for ( const auto & edge : morse_graph . edges_unreduced () ) {
    if ( edge . first >= number_of_morse_sets or
         edge . second >= number_of_morse_sets ) {
      throw std::runtime_error (
        "MorseSingletonReachability received an invalid Morse-graph edge" );
    }
    if ( edge . first != edge . second ) {
      morse_summary [ edge . first ] = MULTIPLE_MORSE_NODES;
    }
  }

  const auto merge_summary =
    [=] ( int32_t left, int32_t right ) -> int32_t {
      if ( left == NO_MORSE_NODE ) return right;
      if ( right == NO_MORSE_NODE ) return left;
      if ( left == right ) return left;
      return MULTIPLE_MORSE_NODES;
    };

  enum VisitState : uint8_t { UNSEEN = 0, ACTIVE = 1, DONE = 2 };
  std::vector<uint8_t> state ( static_cast<size_t> ( n ), UNSEEN );
  std::vector<int32_t> reach_summary (
    static_cast<size_t> ( n ), NO_MORSE_NODE );

  for ( size_t node = 0; node < number_of_morse_sets; ++ node ) {
    const std::vector<uint64_t> cells = morse_graph . morse_set ( node );
    for ( const uint64_t cell : cells ) {
      if ( cell >= n ) {
        throw std::runtime_error (
          "MorseSingletonReachability found a Morse cell outside the MapGraph" );
      }
      state [ cell ] = DONE;
      reach_summary [ cell ] =
        merge_summary ( reach_summary [ cell ], morse_summary [ node ] );
    }
  }

  struct Frame {
    uint32_t vertex;
    uint32_t next_adjacency;
  };
  std::vector<Frame> stack;
  std::vector<int32_t> result (
    query_vertices . size (), NO_MORSE_NODE );

  for ( size_t query_index = 0;
        query_index < query_vertices . size ();
        ++ query_index ) {
    const uint32_t query = static_cast<uint32_t> ( query_vertices [ query_index ] );
    if ( state [ query ] == UNSEEN ) {
      state [ query ] = ACTIVE;
      stack . push_back ( { query, 0 } );
      while ( not stack . empty () ) {
        Frame & frame = stack . back ();
        const MapGraph::AdjacencySpan adjacency =
          map_graph . adjacency_span ( frame . vertex );
        if ( adjacency . size () > std::numeric_limits<uint32_t>::max () ) {
          throw std::runtime_error (
            "MorseSingletonReachability found more than 2^32-1 outgoing "
            "edges at one vertex" );
        }

        if ( frame . next_adjacency == adjacency . size () ) {
          state [ frame . vertex ] = DONE;
          stack . pop_back ();
          continue;
        }

        const uint64_t successor =
          adjacency . begin () [ frame . next_adjacency ];
        if ( successor >= n ) {
          throw std::runtime_error (
            "MorseSingletonReachability found an adjacency outside the "
            "MapGraph" );
        }
        if ( state [ successor ] == UNSEEN ) {
          state [ successor ] = ACTIVE;
          stack . push_back (
            { static_cast<uint32_t> ( successor ), 0 } );
          continue;
        }
        if ( state [ successor ] == ACTIVE ) {
          throw std::runtime_error (
            "MorseSingletonReachability found a directed cycle not covered "
            "by the supplied Morse sets" );
        }

        reach_summary [ frame . vertex ] = merge_summary (
          reach_summary [ frame . vertex ], reach_summary [ successor ] );
        ++ frame . next_adjacency;
      }
    }
    result [ query_index ] = reach_summary [ query ];
  }

  return result;
}

void computeMorseGraph ( MorseGraph & morsegraph,
                         std::shared_ptr<const Map> map,
                         const int SINGLECMG_INIT_PHASE_SUBDIVISIONS,
                         const int SINGLECMG_MIN_PHASE_SUBDIVISIONS,
                         const int SINGLECMG_MAX_PHASE_SUBDIVISIONS,
                         const int SINGLECMG_COMPLEXITY_LIMIT,
                         const char * outputfile ) {
#ifdef CMG_VERBOSE
  std::cout << "SingleCMG: computeMorseGraph.\n";
#endif
  std::shared_ptr < Grid > phase_space = morsegraph . phaseSpace ();
  clock_t start_time = clock ();
  Compute_Morse_Graph ( & morsegraph,
                        phase_space,
                        map,
                        SINGLECMG_INIT_PHASE_SUBDIVISIONS,
                        SINGLECMG_MIN_PHASE_SUBDIVISIONS,
                        SINGLECMG_MAX_PHASE_SUBDIVISIONS,
                        SINGLECMG_COMPLEXITY_LIMIT );
  clock_t stop_time = clock ();
  if ( outputfile != NULL ) {
    morsegraph . save ( outputfile );
  }
  std::ofstream stats_file ( "SingleCMG_statistics.txt" );
  stats_file << "Morse Graph calculation resource usage statistics.\n";
  stats_file << "The final grid has " << phase_space -> size () << " grid elements.\n";
  stats_file << "The computation took " << ((double)(stop_time-start_time)/(double)CLOCKS_PER_SEC)
             << " seconds.\n";
  stats_file << "All memory figures are in bytes:\n";
  stats_file << "grid_memory_use = " << phase_space -> memory () << "\n";
  stats_file << "max_graph_memory = " << max_graph_memory << "\n";
  stats_file << "max_scc_memory_internal = " << max_scc_memory_internal << "\n";
  stats_file << "max_scc_memory_external = " << max_scc_memory_external << "\n";
  stats_file . close ();
}

MorseGraph MorseGraphIntvalMap ( int phase_subdiv_min, int phase_subdiv_max,
                                 std::vector<double> const& phase_lower_bounds,
                                 std::vector<double> const& phase_upper_bounds,
                                 std::vector<double> const& params,
                                 std::string output_file_name ) {
  std::vector<double> param_lower_bounds = params;
  std::vector<double> param_upper_bounds = params;
  int param_dim = params . size();
  int phase_dim = phase_lower_bounds . size();
  std::vector<bool> phase_periodic ( phase_dim, false );
  int phase_subdiv_init = 0;
  int phase_subdiv_limit = 10000;

  Model model;
  model . initialize ( param_dim, phase_dim,
                       phase_subdiv_min, phase_subdiv_max,
                       phase_subdiv_init, phase_subdiv_limit,
                       param_lower_bounds, param_upper_bounds,
                       phase_lower_bounds, phase_upper_bounds,
                       phase_periodic );
  std::shared_ptr<const Map> map = model . map ();

  MorseGraph morsegraph ( model . phaseSpace () );

  // INITIALIZE THE PHASE SPACE SUBDIVISION PARAMETERS
  int SINGLECMG_INIT_PHASE_SUBDIVISIONS = phase_subdiv_init;
  int SINGLECMG_MIN_PHASE_SUBDIVISIONS = phase_subdiv_min;
  int SINGLECMG_MAX_PHASE_SUBDIVISIONS = phase_subdiv_max;
  int SINGLECMG_COMPLEXITY_LIMIT= phase_subdiv_limit;

  // COMPUTE MORSE GRAPH
  computeMorseGraph ( morsegraph, map,
                      SINGLECMG_INIT_PHASE_SUBDIVISIONS,
                      SINGLECMG_MIN_PHASE_SUBDIVISIONS,
                      SINGLECMG_MAX_PHASE_SUBDIVISIONS,
                      SINGLECMG_COMPLEXITY_LIMIT,
                      output_file_name . c_str () );

#ifdef CMG_VERBOSE
  std::cout << "Total Time for Finding Morse Sets ";
  std::cout << "and reachability relation: ";
  std::cout << ": ";
#endif

  // Always output the Morse Graph
  // std::cout << "Creating graphviz .dot file...\n";
  // CreateDotFile ( "morsegraph.gv", conleymorsegraph );

  return morsegraph;
}

MorseGraph MorseGraphMap ( int phase_subdiv_min, int phase_subdiv_max,
                           std::vector<double> const& phase_lower_bounds,
                           std::vector<double> const& phase_upper_bounds,
                           std::string output_file_name,
                           std::function<std::vector<double>(std::vector<double>)> const& F ) {
  std::vector<double> params {0.0};
  std::vector<double> param_lower_bounds = params;
  std::vector<double> param_upper_bounds = params;
  int param_dim = params . size();
  int phase_dim = phase_lower_bounds . size();
  std::vector<bool> phase_periodic ( phase_dim, false );
  int phase_subdiv_init = 0;
  int phase_subdiv_limit = 10000;

  Model model;
  model . initialize ( param_dim, phase_dim,
                       phase_subdiv_min, phase_subdiv_max,
                       phase_subdiv_init, phase_subdiv_limit,
                       param_lower_bounds, param_upper_bounds,
                       phase_lower_bounds, phase_upper_bounds,
                       phase_periodic, F );
  std::shared_ptr<const Map> map = model . map ();

  MorseGraph morsegraph ( model . phaseSpace () );

  // INITIALIZE THE PHASE SPACE SUBDIVISION PARAMETERS
  int SINGLECMG_INIT_PHASE_SUBDIVISIONS = phase_subdiv_init;
  int SINGLECMG_MIN_PHASE_SUBDIVISIONS = phase_subdiv_min;
  int SINGLECMG_MAX_PHASE_SUBDIVISIONS = phase_subdiv_max;
  int SINGLECMG_COMPLEXITY_LIMIT= phase_subdiv_limit;

  // COMPUTE MORSE GRAPH
  computeMorseGraph ( morsegraph, map,
                      SINGLECMG_INIT_PHASE_SUBDIVISIONS,
                      SINGLECMG_MIN_PHASE_SUBDIVISIONS,
                      SINGLECMG_MAX_PHASE_SUBDIVISIONS,
                      SINGLECMG_COMPLEXITY_LIMIT,
                      output_file_name . c_str () );

#ifdef CMG_VERBOSE
  std::cout << "Total Time for Finding Morse Sets ";
  std::cout << "and reachability relation: ";
  std::cout << ": ";
#endif

  // Always output the Morse Graph
  // std::cout << "Creating graphviz .dot file...\n";
  // CreateDotFile ( "morsegraph.gv", conleymorsegraph );

  return morsegraph;
}

/// Python Bindings

#include <pybind11/pybind11.h>
#include <pybind11/functional.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

namespace py = pybind11;

PYBIND11_MODULE(_cmgdb, m) {
  ModelBinding(m);
  GridBinding(m);
  MapGraphBinding(m);
  MorseGraphBinding(m);

  m.doc() = "Conley Morse Graph Database Module";

  m.def("ComputeConleyIndex", &ComputeConleyIndex);
  const char * compute_kwargs_doc =
    "cache_transition_graph: cache the per-level transition graph used "
    "internally by the SCC/reachability passes during the computation "
    "(halves the map evaluations per level).\n"
    "batch_chunk_size: rectangles per batched map call (0 = whole grid).\n"
    "max_cached_edges: abandon a cache projected/measured to exceed this "
    "many edges and fall back to on-demand evaluation (0 = unlimited).\n"
    "reserve_edges: up-front reservation for the flat edge array; 0 = "
    "automatic (2x the edge count projected from the first chunk).\n"
    "reserve_min_edges: reservation engages only when the projected edge "
    "count reaches this size.\n"
    "cache_map_graph: eagerly cache the *returned* map_graph (one extra "
    "full batched map pass; required by the native reachability queries "
    "and by fast Python-side adjacency sweeps). Default False returns a "
    "lazy map_graph that evaluates the map per adjacency query; "
    "map_graph.build_cache() upgrades it later.";
  m.def("ComputeConleyMorseGraph", &ComputeConleyMorseGraph,
        py::arg("model"),
        py::arg("cache_transition_graph") = true,
        py::arg("batch_chunk_size") = 65536,
        py::arg("max_cached_edges") = 0,
        py::arg("reserve_edges") = 0,
        py::arg("reserve_min_edges") = uint64_t ( 1 ) << 24,
        py::arg("cache_map_graph") = false,
        compute_kwargs_doc);
  m.def("ComputeMorseGraph", &ComputeMorseGraph,
        py::arg("model"),
        py::arg("cache_transition_graph") = true,
        py::arg("batch_chunk_size") = 65536,
        py::arg("max_cached_edges") = 0,
        py::arg("reserve_edges") = 0,
        py::arg("reserve_min_edges") = uint64_t ( 1 ) << 24,
        py::arg("cache_map_graph") = false,
        compute_kwargs_doc);
  m.def(
    "ComputeConleyIndexForCells",
    [] ( const Model & model,
         MorseGraph & morse_graph,
         std::vector<uint64_t> cells ) {
      std::vector<std::string> result;
      {
        py::gil_scoped_release release;
        result = ComputeConleyIndexForCells (
          model, morse_graph, std::move ( cells ) );
      }
      return result;
    },
    py::arg ( "model" ),
    py::arg ( "morse_graph" ),
    py::arg ( "cells" ),
    "Compute the Conley index of an arbitrary phase-space cell subset." );
  m.def(
    "MorseDirectedPathCells",
    [] ( const MapGraph & map_graph,
         const MorseGraph & morse_graph,
         const std::vector<uint64_t> & source_nodes,
         const std::vector<uint64_t> & target_nodes ) {
      std::vector<uint64_t> values;
      {
        py::gil_scoped_release release;
        values = ComputeMorseDirectedPathCells (
          map_graph, morse_graph, source_nodes, target_nodes );
      }
      py::array_t<uint64_t> result ( values . size () );
      auto output = result . mutable_unchecked<1> ();
      for ( size_t i = 0; i < values . size (); ++ i ) {
        output ( i ) = values [ i ];
      }
      return result;
    },
    py::arg ( "map_graph" ),
    py::arg ( "morse_graph" ),
    py::arg ( "source_nodes" ),
    py::arg ( "target_nodes" ),
    "Cells on some directed path from the source Morse nodes to the target "
    "Morse nodes. Requires a cached map_graph (cache_map_graph=True or "
    "map_graph.build_cache())." );
  m.def(
    "MorseReachabilityMasks",
    [] ( const MapGraph & map_graph,
         const MorseGraph & morse_graph,
         const std::vector<uint64_t> & query_vertices ) {
      std::vector<uint64_t> values;
      {
        py::gil_scoped_release release;
        values = ComputeMorseReachabilityMasks (
          map_graph, morse_graph, query_vertices );
      }
      py::array_t<uint64_t> result ( values . size () );
      auto output = result . mutable_unchecked<1> ();
      for ( size_t i = 0; i < values . size (); ++ i ) {
        output ( i ) = values [ i ];
      }
      return result;
    },
    py::arg ( "map_graph" ),
    py::arg ( "morse_graph" ),
    py::arg ( "query_vertices" ),
    "Per-query-cell uint64 bitmasks of the Morse nodes reachable through the "
    "box dynamics (bit i = Morse node i). Requires a cached map_graph." );
  m.def(
    "MorseSingletonReachability",
    [] ( const MapGraph & map_graph,
         const MorseGraph & morse_graph,
         const std::vector<uint64_t> & query_vertices ) {
      std::vector<int32_t> values;
      {
        py::gil_scoped_release release;
        values = ComputeMorseSingletonReachability (
          map_graph, morse_graph, query_vertices );
      }
      py::array_t<int32_t> result ( values . size () );
      auto output = result . mutable_unchecked<1> ();
      for ( size_t i = 0; i < values . size (); ++ i ) {
        output ( i ) = values [ i ];
      }
      return result;
    },
    py::arg ( "map_graph" ),
    py::arg ( "morse_graph" ),
    py::arg ( "query_vertices" ),
    "Per-query-cell summary of the reachable Morse nodes: the node id when "
    "exactly one is reachable, -1 when none, -2 when several. Requires a "
    "cached map_graph." );
  m.def("MorseGraphIntvalMap", &MorseGraphIntvalMap);
  m.def("MorseGraphMap", &MorseGraphMap);
}

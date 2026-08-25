// MapGraph.h

#ifndef CMDP_MAPGRAPH_H
#define CMDP_MAPGRAPH_H

#include <cstdint>
#include <exception>
#include <vector>
#include <iterator>
#include <iostream>
#include <algorithm>
// #include <unistd.h>

#include "boost/unordered_map.hpp"
#include "boost/foreach.hpp"

#include "Grid.h"
#include "Map.h"
#include "RectGeo.h"

/// struct MapGraphOptions
///    Controls whether (and how) a MapGraph stores the transition graph.
///    With cache == true the full adjacency structure is computed once at
///    construction, in chunks of chunk_size rectangles, and stored in
///    compressed sparse row (CSR) form; subsequent adjacency queries never
///    evaluate the map again. max_cached_edges bounds the cache: if the edge
///    count exceeds it (or is confidently projected to) the cache is
///    abandoned and the MapGraph falls back to on-demand evaluation
///    (0 means unlimited).
///    reserve_edges controls the up-front reservation of the flat edge
///    array. 0 (the default) reserves twice the final edge count projected
///    from the chunks seen so far, which avoids the repeated-doubling
///    reallocation of multi-gigabyte edge arrays (a transient ~3x memory
///    peak on deep uniform grids); a positive value reserves exactly that
///    many edges instead. Reservation only engages once the projected edge
///    count reaches reserve_min_edges, so small graphs -- including the
///    per-Morse-set subgraphs built during adaptive runs -- never
///    over-allocate.
struct MapGraphOptions {
  bool cache;
  uint64_t chunk_size;
  uint64_t max_cached_edges;
  uint64_t reserve_edges;
  uint64_t reserve_min_edges;
  MapGraphOptions ( void )
    : cache ( false ), chunk_size ( 65536 ), max_cached_edges ( 0 ),
      reserve_edges ( 0 ), reserve_min_edges ( uint64_t ( 1 ) << 24 ) {}
  explicit MapGraphOptions ( bool cache_,
                             uint64_t chunk_size_ = 65536,
                             uint64_t max_cached_edges_ = 0,
                             uint64_t reserve_edges_ = 0,
                             uint64_t reserve_min_edges_ = uint64_t ( 1 ) << 24 )
    : cache ( cache_ ), chunk_size ( chunk_size_ ),
      max_cached_edges ( max_cached_edges_ ),
      reserve_edges ( reserve_edges_ ),
      reserve_min_edges ( reserve_min_edges_ ) {}
};

/// class MapGraph
///    This class is used to created an object suitable for graph algorithms
///    given a grid and a map object. By default "adjacencies" is computed on
///    demand in order to avoid storing the adjacency lists; construction with
///    MapGraphOptions(true) stores the transition graph in CSR form instead,
///    so each rectangle's image is computed exactly once.
class MapGraph {
public:
  // Typedefs
  typedef Grid::size_type size_type;
  typedef Grid::GridElement Vertex;

  /// AdjacencySpan
  ///    Non-owning view over the out-edges of a vertex. When the CSR cache
  ///    is present it points into the flat edge array; otherwise it points
  ///    into a scratch buffer that is overwritten by the next adjacency
  ///    query, so a span must be consumed before the next query is made.
  struct AdjacencySpan {
    const Vertex * begin_;
    const Vertex * end_;
    const Vertex * begin ( void ) const { return begin_; }
    const Vertex * end ( void ) const { return end_; }
    size_t size ( void ) const { return end_ - begin_; }
  };

  // Constructor. Requires Grid and Map.
  MapGraph ( std::shared_ptr<const Grid> grid,
             std::shared_ptr<const Map> f );

  MapGraph ( std::shared_ptr<const Grid> grid,
             std::shared_ptr<const Map> f,
             MapGraphOptions options );

  void initialize ( void );

  /// adjacencies
  ///   Return vector of Vertices which are out-edge adjacencies of input v
  std::vector<Vertex> adjacencies ( const Vertex & v ) const;

  /// adjacency_span
  ///   Return a non-owning view of the out-edge adjacencies of input v.
  ///   Avoids a copy when the transition graph is cached. See AdjacencySpan
  ///   for the lifetime rule in the uncached case.
  AdjacencySpan adjacency_span ( const Vertex & v ) const;

  /// num_vertices
  ///   Return number of vertices
  size_type num_vertices ( void ) const;

  /// has_cache
  ///   Return true if the transition graph is stored in the CSR cache.
  bool has_cache ( void ) const { return cached_; }

  /// num_cached_edges
  ///   Return the number of edges held in the CSR cache (0 if uncached).
  uint64_t num_cached_edges ( void ) const {
    return cached_ ? (uint64_t) csr_edges_ . size () : 0;
  }

  /// build_cache
  ///   Build the CSR transition-graph cache now (a no-op if it is already
  ///   built). Lets a MapGraph constructed lazily -- e.g. the map_graph
  ///   returned by ComputeMorseGraph with cache_map_graph=False -- be
  ///   upgraded to a cached one after the fact, at the cost of one full
  ///   (batched, if available) pass of map evaluations over the grid.
  void build_cache ( void );

private:
  // Private methods
  std::vector<size_type> compute_adjacencies ( const size_type & v ) const;
  // Private data
  std::shared_ptr<const Grid> grid_;
  std::shared_ptr<const Map> f_;
  MapGraphOptions options_;
  bool cached_;
  std::vector<uint64_t> csr_offsets_;
  std::vector<Vertex> csr_edges_;
  mutable std::vector<Vertex> scratch_;
};

inline
MapGraph::MapGraph ( std::shared_ptr<const Grid> grid,
                     std::shared_ptr<const Map> f ) :
grid_ ( grid ),
f_ ( f ),
options_ (),
cached_ ( false ) {
  initialize ();
}

inline
MapGraph::MapGraph ( std::shared_ptr<const Grid> grid,
                     std::shared_ptr<const Map> f,
                     MapGraphOptions options ) :
grid_ ( grid ),
f_ ( f ),
options_ ( options ),
cached_ ( false ) {
  initialize ();
}

inline void
MapGraph::initialize ( void ) {
  if ( not f_ ) {
    throw std::logic_error ( "MapGraph::MapGraph. Unable to construct with uninitialized Map f\n");
  }
  if ( options_ . cache ) {
    build_cache ();
  }
}

/// build_cache
///    Structural intent: evaluate the multivalued map F(v) = cover(f(geo(v)))
///    exactly once per grid element, storing the resulting digraph in CSR
///    form (csr_offsets_, csr_edges_). The graph algorithms (strongly
///    connected components for the recurrent sets, reachability for the
///    Morse graph partial order) then read this stored digraph instead of
///    re-evaluating the map on every pass. Evaluation proceeds in chunks so
///    that maps providing a batched interface (Map::has_batch) receive many
///    rectangles per call -- for Python-defined maps this means one Python
///    call per chunk instead of one per rectangle. The adjacency lists this
///    produces are identical, element for element, to the on-demand path.
inline void
MapGraph::build_cache ( void ) {
  if ( cached_ ) return;
  const uint64_t n = num_vertices ();
  csr_offsets_ . clear ();
  csr_edges_ . clear ();
  csr_offsets_ . reserve ( n + 1 );
  csr_offsets_ . push_back ( 0 );

  // The flat-buffer batch path requires rectangle geometry; fall back to
  // per-element evaluation for grids with other geometry types.
  bool use_batch = f_ -> has_batch () && n > 0 &&
    ( std::dynamic_pointer_cast<RectGeo> ( grid_ -> geometry ( (Vertex) 0 ) ) != nullptr );

  const uint64_t chunk = options_ . chunk_size > 0 ? options_ . chunk_size : n;
  std::vector<double> rects;
  std::vector<double> images;

  for ( uint64_t start = 0; start < n; start += chunk ) {
    const uint64_t stop = std::min ( start + chunk, n );
    if ( use_batch ) {
      // Gather rectangle bounds for this chunk into a flat buffer.
      const uint64_t count = stop - start;
      uint64_t dim = 0;
      rects . clear ();
      for ( uint64_t v = start; v < stop; ++ v ) {
        std::shared_ptr<RectGeo> rect =
          std::dynamic_pointer_cast<RectGeo> ( grid_ -> geometry ( (Vertex) v ) );
        if ( not rect ) {
          throw std::logic_error ( "MapGraph::build_cache. Mixed geometry types in grid.\n" );
        }
        dim = rect -> dimension ();
        rects . insert ( rects . end (), rect -> lower_bounds . begin (),
                         rect -> lower_bounds . end () );
        rects . insert ( rects . end (), rect -> upper_bounds . begin (),
                         rect -> upper_bounds . end () );
      }
      // One map evaluation for the whole chunk.
      f_ -> batch_map ( rects, count, dim, images );
      // Cover each image rectangle to produce the adjacency lists.
      RectGeo image ( dim );
      for ( uint64_t i = 0; i < count; ++ i ) {
        const double * bounds = images . data () + i * 2 * dim;
        for ( uint64_t d = 0; d < dim; ++ d ) {
          image . lower_bounds [ d ] = bounds [ d ];
          image . upper_bounds [ d ] = bounds [ dim + d ];
        }
        std::vector<Vertex> targets = grid_ -> cover ( image );
        csr_edges_ . insert ( csr_edges_ . end (), targets . begin (), targets . end () );
        csr_offsets_ . push_back ( csr_edges_ . size () );
      }
    } else {
      for ( uint64_t v = start; v < stop; ++ v ) {
        std::vector<Vertex> targets = compute_adjacencies ( (Vertex) v );
        csr_edges_ . insert ( csr_edges_ . end (), targets . begin (), targets . end () );
        csr_offsets_ . push_back ( csr_edges_ . size () );
      }
    }
    // Bound the cache: on overflow abandon it and fall back to on-demand
    // evaluation rather than exhausting memory.
    if ( options_ . max_cached_edges > 0 &&
         (uint64_t) csr_edges_ . size () > options_ . max_cached_edges ) {
      csr_offsets_ . clear ();
      csr_edges_ . clear ();
      csr_offsets_ . shrink_to_fit ();
      csr_edges_ . shrink_to_fit ();
      cached_ = false;
      return;
    }
    if ( stop < n ) {
      // Project the final edge count from the edge density seen so far.
      // The projection sizes the up-front reservation of the flat edge
      // array -- a deep uniform grid's multi-gigabyte edge array is then
      // allocated (nearly) once instead of repeatedly doubled with a
      // transient ~3x memory peak -- and lets a user-set cap fail fast
      // instead of building a doomed cache all the way up to the cap.
      // Chunks follow the tree order (a spatial sweep), so early density
      // can be biased; the projection is therefore refreshed every chunk,
      // the auto reservation doubles it for headroom, and the early cap
      // abandon requires a 2x margin of confidence (the exact check above
      // remains authoritative).
      const double density = (double) csr_edges_ . size () / (double) stop;
      const uint64_t projected = (uint64_t) ( density * (double) n ) + 1;
      if ( options_ . max_cached_edges > 0 &&
           projected > 2 * options_ . max_cached_edges ) {
        csr_offsets_ . clear ();
        csr_edges_ . clear ();
        csr_offsets_ . shrink_to_fit ();
        csr_edges_ . shrink_to_fit ();
        cached_ = false;
        return;
      }
      if ( projected >= options_ . reserve_min_edges ) {
        uint64_t target = options_ . reserve_edges > 0 ?
          options_ . reserve_edges : 2 * projected;
        if ( options_ . max_cached_edges > 0 ) {
          target = std::min ( target, options_ . max_cached_edges );
        }
        if ( target > (uint64_t) csr_edges_ . capacity () ) {
          try {
            csr_edges_ . reserve ( target );
          } catch ( std::bad_alloc const& ) {
            // The reservation is an optimization only; fall back to
            // ordinary vector growth if the allocator refuses it.
          }
        }
      }
    }
  }
  cached_ = true;
}

inline std::vector<MapGraph::Vertex>
MapGraph::adjacencies ( const size_type & source ) const {
  if ( cached_ ) {
    return std::vector<Vertex> ( csr_edges_ . begin () + csr_offsets_ [ source ],
                                 csr_edges_ . begin () + csr_offsets_ [ source + 1 ] );
  }
  return compute_adjacencies ( source );
}

inline MapGraph::AdjacencySpan
MapGraph::adjacency_span ( const size_type & source ) const {
  if ( cached_ ) {
    const Vertex * base = csr_edges_ . data ();
    return { base + csr_offsets_ [ source ], base + csr_offsets_ [ source + 1 ] };
  }
  scratch_ = compute_adjacencies ( source );
  return { scratch_ . data (), scratch_ . data () + scratch_ . size () };
}

inline std::vector<MapGraph::Vertex>
MapGraph::compute_adjacencies ( const Vertex & source ) const {
  std::vector < Vertex > target =
    grid_ -> cover ( (*f_) ( grid_ -> geometry ( source ) ) ); // here is the work
  return target;
}

inline MapGraph::size_type
MapGraph::num_vertices ( void ) const {
  return grid_ -> size ();
}

/// Python Bindings

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
namespace py = pybind11;

inline void
MapGraphBinding(py::module &m) {
  py::class_<MapGraph, std::shared_ptr<MapGraph>>(m, "MapGraph")
    .def(py::init<std::shared_ptr<const Grid>, std::shared_ptr<const Map>>())
    .def(py::init([](std::shared_ptr<const Grid> grid,
                     std::shared_ptr<const Map> f, bool cache) {
           return new MapGraph ( grid, f, MapGraphOptions ( cache ) );
         }),
         py::arg("grid"), py::arg("map"), py::arg("cache"))
    .def("num_vertices", &MapGraph::num_vertices)
    .def("has_cache", &MapGraph::has_cache)
    .def("num_cached_edges", &MapGraph::num_cached_edges)
    .def("build_cache", &MapGraph::build_cache,
         "Build the CSR transition-graph cache now (no-op if already built). "
         "Upgrades a lazily returned map_graph to a cached one at the cost "
         "of one full pass of map evaluations over the grid.")
    .def("adjacencies", &MapGraph::adjacencies);
}

#endif

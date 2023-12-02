// MapGraph.h

#ifndef CMDP_MAPGRAPH_H
#define CMDP_MAPGRAPH_H

#include <exception>
#include <vector>
#include <iterator>
#include <iostream>
#include <algorithm>
// #include <unistd.h>

#include "boost/unordered_map.hpp"
#include "boost/foreach.hpp"

#include "Grid.h"

#ifdef CMDB_STORE_GRAPH
#include "ComputeGraph.h"
#endif

/// class MapGraph
///    This class is used to created an object suitable for graph algorithms
///    given a grid and a map object. By default "adjacencies" is computed on demand
///    in order to avoid storing the adjacency lists.
class MapGraph {
public:
  // Typedefs
  typedef Grid::size_type size_type;
  typedef Grid::GridElement Vertex;

  // Constructor. Requires Grid and Map.
  MapGraph ( std::shared_ptr<const Grid> grid,
             std::shared_ptr<const Map> f );

  // MapGraph ( std::shared_ptr<const Grid> grid,
  //            std::shared_ptr<const Model> model );

  void initialize ( void );

  /// adjacencies
  ///   Return vector of Vertices which are out-edge adjacencies of input v
  std::vector<Vertex> adjacencies ( const Vertex & v ) const;
  
  /// num_vertices
  ///   Return number of vertices
  size_type num_vertices ( void ) const;

private:
  // Private methods
  std::vector<size_type> compute_adjacencies ( const size_type & v ) const;
  // Private data
  std::shared_ptr<const Grid> grid_;
  std::shared_ptr<const Map> f_;
  // Variables used if graph is stored in memory. (See CMDB_STORE_GRAPH define)
  bool stored_graph;
  std::vector<std::vector<Vertex> > adjacency_lists_;
};

inline
MapGraph::MapGraph ( std::shared_ptr<const Grid> grid,
                     std::shared_ptr<const Map> f ) :
grid_ ( grid ),
f_ ( f ),
stored_graph ( false ) {
  initialize ();
}

// inline
// MapGraph::MapGraph ( std::shared_ptr<const Grid> grid,
//                      std::shared_ptr<const Model> model ) :
// grid_ ( grid ),
// f_ ( model -> map () ),
// stored_graph ( false ) {
//   initialize ();
// }

inline void
MapGraph::initialize ( void ) {
  if ( not f_ ) {
    throw std::logic_error ( "MapGraph::MapGraph. Unable to construct with uninitialized Map f\n");
  }
#ifdef CMDB_STORE_GRAPH
  
  // Determine whether it is efficient to use an MPI job to store the graph
  if ( num_vertices () < 10000 ) {
    stored_graph = false;
    return;
  }
  stored_graph = true;
  
  // Make a file with required integrations
  MapEvals evals;
  evals . parameter () = f . parameter ();
  for ( size_type source = 0; source < num_vertices (); ++ source ) {
    Vertex domain_cell = lookup ( source );
    evals . insert ( domain_cell );
  }
  
  std::cout << "Saving grid to file.\n";
  // Save the grid and a list of required evaluations to disk
  grid_ -> save ("grid.txt");
  evals . save ( "mapevals.txt" );
  
  // Call a program to compute the adjacency information
  std::cout << "Calling MPI program to evaluate map.\n";
  system("./COMPUTEGRAPHSCRIPT");
  std::cout << "MPI program returned.\n";

  // Load and store the adjacency information
  evals . load ( "mapevals.txt" );
  adjacency_lists_ . resize ( num_vertices () );
  for ( size_type source = 0; source < num_vertices (); ++ source ) {
    Vertex domain_cell = lookup ( source );    
    index ( &adjacency_lists_ [ source ], evals . val ( domain_cell ) );
  }
  std::cout << "Map stored.\n";
#endif
}

inline std::vector<MapGraph::Vertex>
MapGraph::adjacencies ( const size_type & source ) const {
  if ( stored_graph )
    return adjacency_lists_ [ source ];
  else
    return compute_adjacencies ( source );
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
    // .def(py::init<std::shared_ptr<const Grid>, std::shared_ptr<const Model>>())
    .def("num_vertices", &MapGraph::num_vertices)
    .def("adjacencies", &MapGraph::adjacencies);
}

#endif

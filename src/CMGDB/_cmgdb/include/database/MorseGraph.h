/* MorseGraph.h */

#ifndef _CMDP_MORSE_GRAPH_
#define _CMDP_MORSE_GRAPH_

#include <fstream>
#include <utility>

#include <boost/iterator/counting_iterator.hpp>
#include <boost/version.hpp>
#include <memory>
#include "boost/foreach.hpp"
#include <boost/serialization/serialization.hpp>
#include "boost/functional/hash.hpp"

#include <unordered_set>
#include <boost/serialization/set.hpp>
#include <boost/serialization/unordered_set.hpp>
#include "boost/serialization/shared_ptr.hpp"
#include "boost/serialization/set.hpp"
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include "Grid.h"
#include "chomp/ConleyIndex.h"

/** Conley-Morse Graph. (Can be called either MorseGraph or ConleyMorseGraph)
 *  It is the directed acyclic graph whose vertices may
 *  be annotated with "Grid" objects, representing combinatorial Morse sets,
 *  and "ConleyIndex_t" objects, representing Conley Indexes.
 *  The edges of the graph represent the reachability relation among the
 *  combinatorial Morse sets. */
class MorseGraph {
 public:

  // Vertex types
  typedef int Vertex;  
  typedef boost::counting_iterator<Vertex> VertexIterator;
  typedef std::pair<VertexIterator, VertexIterator> VertexIteratorPair;

  // Edge types
  typedef std::pair<Vertex, Vertex> Edge;
  typedef std::unordered_set<Edge, boost::hash<Edge>>::const_iterator EdgeIterator;
  typedef std::pair<EdgeIterator, EdgeIterator> EdgeIteratorPair;

  /** Create an empty graph */
  MorseGraph();
  /** Initialize MorseGraph with phase space */
  template < class grid_ptr >
  MorseGraph ( grid_ptr phasespace ) : num_vertices_(0), phasespace_(phasespace) {}
  /** Load MorseGraph from file */
  MorseGraph ( const char * filename );
  // WRITE FUNCTIONS
  
  /** Create a new vertex and return the descriptor of the vertex.
   *  The vertex is not connected to anywhere.
   *  An empty grid and ConleyIndex are allocated and
   *  assigned to the vertex just after this function is called.
   *  ("empty" objects mean the objects created by default-constructor.)
   */
  Vertex AddVertex ( void );
  
  /** Add a edge from the "from" vertex to the "to" vertex.
   */
  Edge AddEdge(Vertex from, Vertex to);

  /** Remove a "from"-"to" edge
   *
   *  NOTE: some edge iterators become invalid when you call this function 
   */
  void RemoveEdge(Vertex from, Vertex to);
  
  //// READ FUNCTIONS
  
  /** return a number of vertices */
  unsigned int NumVertices ( void ) const;

  /** return a iterator pair to all vertices */
  VertexIteratorPair Vertices ( void ) const;
  
  /** return a iterator pair to all edges */
  EdgeIteratorPair Edges ( void ) const;

  std::vector<uint64_t>
  vertices ( void ) const;

  std::vector<std::vector<uint64_t>>
  edges ( void ) const;

  std::vector<uint64_t>
  adjacencies ( uint64_t vertex ) const;

  std::vector<uint64_t>
  morse_set ( uint64_t vertex ) const;

  std::vector<std::vector<double>>
  morse_set_boxes ( uint64_t vertex ) const;

  //// PROPERTY ACCESS
  
  /** Get the grid associated with the vertex. */
  std::shared_ptr<Grid> & phaseSpace ( void );
  std::shared_ptr<const Grid> phaseSpace ( void ) const;
  std::shared_ptr<Grid> & grid (Vertex vertex);
  std::shared_ptr<const Grid> grid (Vertex vertex) const;
  std::shared_ptr<chomp::ConleyIndex_t> & conleyIndex (Vertex vertex);
  std::shared_ptr<const chomp::ConleyIndex_t> conleyIndex (Vertex vertex) const;

  std::set< std::string > & annotation ( void );
  std::set< std::string > & annotation ( Vertex vertex );
  const std::set< std::string > & annotation ( void ) const;
  const std::set< std::string > & annotation ( Vertex vertex ) const;

  /** Remove the grids associated with the vertices */
  void clearGrids ( void );
  
  //// FILE IO

  // Save to File
  void save ( const char * filename ) const {
    std::ofstream ofs(filename);
    assert(ofs.good());
    for ( int vertex = 0; vertex < num_vertices_; ++vertex ) {
      for ( Grid::iterator it = grid (vertex) -> begin ();
            it != grid (vertex) -> end (); ++ it ) {
        std::shared_ptr<Geo> geo = grid (vertex) -> geometry (it);
        ofs << *geo << ", " << vertex << std::endl;
      }
    }
  }

  // /** Save to File */
  // void save ( const char * filename ) const {
  //   std::ofstream ofs(filename);
  //   assert(ofs.good());
  //   boost::archive::text_oarchive oa(ofs);
  //   oa << *this;
  // }

  /** Load from file */
  void load ( const char * filename ) {
    std::ifstream ifs(filename);
    if ( not ifs . good () ) {
      std::cout << "Could not load " << filename << "\n";
      throw "Could not load!";
    }
    boost::archive::text_iarchive ia(ifs);
    ia >> *this;
  }
private:
  // DATA
  int num_vertices_;
  std::unordered_set < Edge, boost::hash<Edge> > edges_;
  std::shared_ptr < Grid > phasespace_;
  std::vector < std::shared_ptr <Grid > > grids_;
  std::vector < std::shared_ptr < chomp::ConleyIndex_t > > conleyindexes_;
  std::set < std::string > annotation_;
  std::vector < std::set < std::string > > annotation_by_vertex_;
  //// SERIALIZATION
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive& ar, const unsigned int version) {
    ar & num_vertices_;
    ar & edges_;
    ar & phasespace_;
    ar & grids_;
    ar & conleyindexes_;
    ar & annotation_;
    ar & annotation_by_vertex_;
  }
};

typedef MorseGraph ConleyMorseGraph;

/****************/
/*  DEFINITIONS */
/****************/

/** Create an empty graph */
inline
MorseGraph::MorseGraph ( void ) {
  num_vertices_ = 0;
}

inline
MorseGraph::MorseGraph ( const char * filename ) {
  load ( filename );
}

/** Create a new vertex and return the descriptor of the vertex.
 *  The vertex is not connected to anywhere.
 *  An empty gridPtr and ConleyIndex are allocated and
 *  assigned to the vertex just after this function is called.
 *  ("empty" objects mean the objects created by default-constructor.)
 */
inline MorseGraph::Vertex MorseGraph::AddVertex ( void ) {
  int v = num_vertices_ ++;
  grids_ . push_back ( std::shared_ptr <Grid > ());
  conleyindexes_ . push_back ( std::shared_ptr <chomp::ConleyIndex_t > ());
  annotation_by_vertex_ . resize ( num_vertices_ );
  return v;
}

/** Add a edge from the "from" vertex to the "to" vertex. */
inline MorseGraph::Edge MorseGraph::AddEdge(Vertex from, Vertex to) {
  return * edges_ . insert ( Edge (from, to) ) . first;
}

/** Remove a "from"-"to" edge. (invalidates iterators) */
inline void MorseGraph::RemoveEdge(Vertex from, Vertex to) {
  edges_ . erase ( Edge (from, to) );
}

/** return a number of vertices */
inline unsigned int MorseGraph::NumVertices ( void ) const {
  return num_vertices_;
}

inline std::vector<uint64_t> MorseGraph::
vertices ( void ) const {
  std::vector<uint64_t> verts;
  for ( uint64_t vertex = 0; vertex < num_vertices_; ++vertex ) {
    verts . push_back ( vertex );
  }
  return verts;
}

inline std::vector<std::vector<uint64_t>> MorseGraph::
edges ( void ) const {
  std::vector<std::vector<uint64_t>> edge_list;
  EdgeIteratorPair eip = Edges ();
  for ( EdgeIterator e_it = eip . first; e_it != eip . second; ++ e_it ) {
    uint64_t v1 = e_it -> first;
    uint64_t v2 = e_it -> second;
    std::vector<uint64_t> e_verts = { v1, v2 };
    edge_list . push_back ( e_verts );
  }
  return edge_list;
}

inline std::vector<uint64_t> MorseGraph::
adjacencies ( uint64_t vertex ) const {
  std::vector<uint64_t> vert_adjacencies;
  EdgeIteratorPair eip = Edges ();
  for ( EdgeIterator e_it = eip . first; e_it != eip . second; ++ e_it ) {
    if ( e_it -> first == vertex ) {
      vert_adjacencies . push_back ( e_it -> second );
    }
  }
  return vert_adjacencies;
}

inline std::vector<uint64_t> MorseGraph::
morse_set ( uint64_t vertex ) const {
  std::vector<uint64_t> morseset;
  for ( Grid::iterator it = grid (vertex) -> begin ();
            it != grid (vertex) -> end (); ++ it ) {
    morseset . push_back ( *it );
  }
  return morseset;
}

inline std::vector<std::vector<double>> MorseGraph::
morse_set_boxes ( uint64_t vertex ) const {
  std::vector<std::vector<double>> morse_boxes;
  for ( Grid::iterator it = grid (vertex) -> begin ();
            it != grid (vertex) -> end (); ++ it ) {
    std::shared_ptr<Geo> geo = grid (vertex) -> geometry (it);
    std::vector<double> bounds = geo -> get_lower_bounds ();
    std::vector<double> u_bounds = geo -> get_upper_bounds ();
    bounds . insert(bounds . end(), u_bounds . begin(), u_bounds . end());
    morse_boxes . push_back ( bounds );
  }
  return morse_boxes;
}

/** return a iterator pair to all vertices */
inline MorseGraph::VertexIteratorPair
MorseGraph::Vertices ( void ) const {
  return VertexIteratorPair ( 0, num_vertices_ );
}

/** return a iterator pair to all edges */
inline MorseGraph::EdgeIteratorPair
MorseGraph::Edges ( void ) const {
  return EdgeIteratorPair(edges_.begin(), edges_.end());
}

// PROPERTY ACCESS

/** accessor method for phase space grid */
inline
std::shared_ptr<Grid> & MorseGraph::phaseSpace ( void ) {
  return phasespace_;
}
/** accessor method for phase space grid, const version */
inline
std::shared_ptr<const Grid> MorseGraph::phaseSpace ( void ) const {
  return phasespace_;
}

/** accessor method for grid assigned to vertex */
inline std::shared_ptr<Grid> & MorseGraph::grid(Vertex vertex) {
  return grids_[vertex];
}

/** accessor method for grid assigned to vertex, const version */
inline
std::shared_ptr<const Grid> MorseGraph::grid(Vertex vertex) const {
  return grids_[vertex];
}

/** accessor method for conley index assigned to vertex */
inline
std::shared_ptr<chomp::ConleyIndex_t> & MorseGraph::conleyIndex(Vertex vertex) {
  return conleyindexes_ [ vertex ];
}

/** accessor method for conley index assigned to vertex, const version */
inline std::shared_ptr<const chomp::ConleyIndex_t>
MorseGraph::conleyIndex (Vertex vertex) const {
  return conleyindexes_ [ vertex ];
}

inline std::set< std::string > & MorseGraph::annotation ( void ) {
  return annotation_;
}

inline std::set< std::string > & MorseGraph::annotation ( Vertex vertex ) {
  return annotation_by_vertex_ [ vertex ];
}

inline const std::set< std::string > & MorseGraph::annotation ( void ) const {
  return annotation_;
}

inline const std::set< std::string > & MorseGraph::annotation ( Vertex vertex ) const {
  return annotation_by_vertex_ [ vertex ];
}

/** method to clear grids associated with all vertices */
inline void MorseGraph::clearGrids ( void ) {
  phasespace_ . reset ();
  VertexIteratorPair vip = Vertices ();
  for ( VertexIterator vi = vip . first; vi != vip . second; ++ vi ) {
    grid ( *vi ) . reset ();
  }
}

/// Python Bindings

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
namespace py = pybind11;

inline void
MorseGraphBinding(py::module &m) {
  py::class_<MorseGraph, std::shared_ptr<MorseGraph>>(m, "MorseGraph")
    .def(py::init<>())
    .def(py::init<std::shared_ptr<Grid>>())
    .def(py::init<const char*>())
    .def("num_vertices", &MorseGraph::NumVertices)
    .def("vertices", &MorseGraph::vertices)
    .def("edges", &MorseGraph::edges)
    .def("adjacencies", &MorseGraph::adjacencies)
    .def("morse_set", &MorseGraph::morse_set)
    .def("morse_set_boxes", &MorseGraph::morse_set_boxes);
}

#endif

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
#include "SuccinctGrid.h"
BOOST_CLASS_EXPORT_IMPLEMENT(SuccinctGrid);
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

std::pair<MorseGraph, MapGraph> ComputeConleyMorseGraph ( Model const& model ) {
  std::shared_ptr<const Map> map = model . map ();
  MorseGraph morsegraph ( model . phaseSpace () );
  std::shared_ptr < Grid > phase_space = morsegraph . phaseSpace ();

  int phase_subdiv_init = model . phase_subdiv_init ();
  int phase_subdiv_min = model . phase_subdiv_min ();
  int phase_subdiv_max = model . phase_subdiv_max ();
  int phase_subdiv_limit = model . phase_subdiv_limit ();

  // Compute Morse graph
  Compute_Morse_Graph ( & morsegraph, phase_space, map, phase_subdiv_init,
                        phase_subdiv_min, phase_subdiv_max, phase_subdiv_limit );

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

  // Compute multi-valued map digraph
  MapGraph map_graph ( phase_space, map );

  return std::make_pair ( morsegraph, map_graph );
}

std::pair<MorseGraph, MapGraph> ComputeMorseGraph ( Model const& model ) {
  std::shared_ptr<const Map> map = model . map ();
  MorseGraph morsegraph ( model . phaseSpace () );
  std::shared_ptr < Grid > phase_space = morsegraph . phaseSpace ();

  int phase_subdiv_init = model . phase_subdiv_init ();
  int phase_subdiv_min = model . phase_subdiv_min ();
  int phase_subdiv_max = model . phase_subdiv_max ();
  int phase_subdiv_limit = model . phase_subdiv_limit ();

  // Compute Morse graph
  Compute_Morse_Graph ( & morsegraph, phase_space, map, phase_subdiv_init,
                        phase_subdiv_min, phase_subdiv_max, phase_subdiv_limit );

  // Compute multi-valued map digraph
  MapGraph map_graph ( phase_space, map );

  return std::make_pair ( morsegraph, map_graph );
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

  std::cout << "Total Time for Finding Morse Sets ";
  std::cout << "and reachability relation: ";
  std::cout << ": ";

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

  std::cout << "Total Time for Finding Morse Sets ";
  std::cout << "and reachability relation: ";
  std::cout << ": ";

  // Always output the Morse Graph
  // std::cout << "Creating graphviz .dot file...\n";
  // CreateDotFile ( "morsegraph.gv", conleymorsegraph );

  return morsegraph;
}

/// Python Bindings

#include <pybind11/pybind11.h>
#include <pybind11/functional.h>
#include <pybind11/stl.h>

namespace py = pybind11;

PYBIND11_MODULE(_cmgdb, m) {
  ModelBinding(m);
  GridBinding(m);
  MapGraphBinding(m);
  MorseGraphBinding(m);

  m.doc() = "Conley Morse Graph Database Module";

  m.def("ComputeConleyIndex", &ComputeConleyIndex);
  m.def("ComputeConleyMorseGraph", &ComputeConleyMorseGraph);
  m.def("ComputeMorseGraph", &ComputeMorseGraph);
  m.def("MorseGraphIntvalMap", &MorseGraphIntvalMap);
  m.def("MorseGraphMap", &MorseGraphMap);
}

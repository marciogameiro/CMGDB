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

#include <boost/serialization/export.hpp>
#include "SuccinctGrid.h"
BOOST_CLASS_EXPORT_IMPLEMENT(SuccinctGrid);
#include "PointerGrid.h"
BOOST_CLASS_EXPORT_IMPLEMENT(PointerGrid);

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

int main ( int argc, char * argv [] ) {
  Model model;
  model . initialize ( argc, argv );
  std::shared_ptr<const Map> map = model . map ();

  MorseGraph morsegraph ( model . phaseSpace () );

  // INITIALIZE THE PHASE SPACE SUBDIVISION PARAMETERS FROM CONFIG FILE
  Configuration config;
  config . loadFromFile ( "./" );
  int SINGLECMG_INIT_PHASE_SUBDIVISIONS = config . PHASE_SUBDIV_INIT;
  int SINGLECMG_MIN_PHASE_SUBDIVISIONS = config . PHASE_SUBDIV_MIN;
  int SINGLECMG_MAX_PHASE_SUBDIVISIONS = config . PHASE_SUBDIV_MAX;
  int SINGLECMG_COMPLEXITY_LIMIT= config . PHASE_SUBDIV_LIMIT;

  // COMPUTE MORSE GRAPH
  computeMorseGraph ( morsegraph, map,
                      SINGLECMG_INIT_PHASE_SUBDIVISIONS,
                      SINGLECMG_MIN_PHASE_SUBDIVISIONS,
                      SINGLECMG_MAX_PHASE_SUBDIVISIONS,
                      SINGLECMG_COMPLEXITY_LIMIT, "data.mg" );
  std::cout << "Total Time for Finding Morse Sets ";
  std::cout << "and reachability relation: ";
  std::cout << ": ";

  ConleyMorseGraph & conleymorsegraph = morsegraph;

  // Always output the Morse Graph
  std::cout << "Creating graphviz .dot file...\n";
  CreateDotFile ( "morsegraph.gv", conleymorsegraph );

  return 0;
}

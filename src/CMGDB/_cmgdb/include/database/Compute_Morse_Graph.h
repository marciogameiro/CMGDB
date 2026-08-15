/*
 *  Compute_MorseGraph.h
 */

#ifndef _CMDP_COMPUTE_CONLEY_MORSE_GRAPH_
#define _CMDP_COMPUTE_CONLEY_MORSE_GRAPH_

#include <memory>

#include "MorseGraph.h"
#include "Grid.h"
#include "Map.h"
#include "MapGraph.h"

/// Computes the Morse decomposition with respect to the given map
/// on the given phase space, and creates its representation by means
/// of a Conley-Morse graph.
/// The options control the per-level transition-graph cache (see
/// MapGraphOptions); by default the transition graph of each subdivision
/// level is computed once and reused across the graph algorithms.

void Compute_Morse_Graph (MorseGraph * MG,
                          std::shared_ptr<Grid> phase_space,
                          std::shared_ptr<const Map> interval_map,
                          const unsigned int Init,
                          const unsigned int Min,
                          const unsigned int Max,
                          const unsigned int Limit,
                          MapGraphOptions options = MapGraphOptions ( true ));

void Compute_Morse_Graph (MorseGraph * MG,
                          std::shared_ptr<Grid> phase_space,
                          std::shared_ptr<const Map> interval_map,
                          const unsigned int Min,
                          const unsigned int Max,
                          const unsigned int Limit,
                          MapGraphOptions options = MapGraphOptions ( true ));

#include "Compute_Morse_Graph.hpp"

#endif

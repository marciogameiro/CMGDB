// GraphTheory.hpp Shaun Harker 
// created 5/16/2011

#ifndef CMDP_GRAPHTHEORY
#define CMDP_GRAPHTHEORY

#include "Grid.h"
#include "Map.h"
#include <vector>
#include <queue>
#include <memory>

/// computeMorseSetsAndReachability
void computeMorseSetsAndReachability (std::vector< std::shared_ptr<Grid> > * output,
                                      std::vector<std::vector<unsigned int> > * reach,
                                      std::shared_ptr<const Grid> G,
                                      std::shared_ptr<const Map> f );

/// computeStrongComponents
///    Modified version of Tarjan's algorithm devised by Shaun Harker
///    Only calls for adjacency lists once each, yet only requires O(V) space.
template < class Graph >
void computeStrongComponents (std::vector<std::deque<typename Graph::Vertex> > * output,
                              const Graph & G,
         /* optional output */std::deque<typename Graph::Vertex> * topological_sort = 0,
         /* optional output */std::deque<typename Graph::Vertex> * SCC_root = 0);

/// computeReachability
template < class Graph >
void computeReachability ( std::vector < std::vector < unsigned int > > * output, 
                           std::vector<std::deque<typename Graph::size_type> > & morse_sets, 
                           const Graph & G, 
                           const std::deque<typename Graph::size_type> & topological_sort );

#include "GraphTheory.hpp"

#endif

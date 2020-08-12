#ifndef ATLAS_H
#define ATLAS_H

#include <iostream>
#include <cstdint>

#include <fstream>
#include <sstream>
#include <string>
#include <cmath>
#include <vector>
#include <unordered_map>

#include <boost/iterator/counting_iterator.hpp>
#include <boost/range/iterator_range.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>

#include "Grid.h"
#include "TreeGrid.h"
#include "PointerGrid.h"
#include "RankSelect.h"
#include "Geo.h"
#include "AtlasGeo.h"

/// class Atlas
///   Grid data structure which stores other Grids as "charts"
///   Charts are added with "add_chart_" method.
///   Each added chart is assigned a "chart_id" by the user which may be any integer.
///   The integers must be unique to the chart, but need not be contiguous.
///   Once the charts have been added the "finalize" method must be called
///     to make the data structure usable.
///   Alternatively, charts may be loaded from a file with the "import_charts" method. 

class Atlas : public Grid { 

public:
	typedef uint64_t GridElement;
  typedef boost::counting_iterator < GridElement > iterator;
  typedef iterator const_iterator;
  typedef uint64_t size_type;
  typedef std::shared_ptr<TreeGrid> Chart;

  // Contructor/ Desctructor
  Atlas ( void ) { }
  ~Atlas ( void ) { }

  // Methods inherited from Grid
  virtual Atlas * clone ( void ) const;
  virtual void subdivide ( void );
  virtual Grid * subgrid ( const std::deque < GridElement > & grid_elements ) const;
  virtual std::vector<GridElement> subset ( const Grid & other ) const;
  virtual std::shared_ptr<Geo> geometry ( GridElement ge ) const;  
  virtual std::vector<Grid::GridElement> cover ( const Geo & geo ) const;
  using Grid::geometry;
  using Grid::cover;

  // Atlas-Specific functionality

  /// Typedefs for Atlas
  typedef std::unordered_map<size_type, Chart>::const_iterator ChartIterator;
  typedef boost::iterator_range<ChartIterator> ChartIteratorRange;
  typedef std::pair <size_type, Chart > IdChartPair;

  /// chart
  ///   Accessor method for chart via chart_id
  Chart & 
  chart ( size_type chart_id );  

  /// chart (const method)
  ///   Accessor method for chart via chart_id
  Chart const& 
  chart ( size_type chart_id ) const;

  /// clear
  ///   Revert to an empty Atlas structure
  void 
  clear ( void );

  /// importCharts
  ///   Read an input file detailing an Atlas structure 
  ///   and create the appropriate data structure.
  void 
  importCharts ( const char * inputfile );
  
  /// add_chart
  ///   add a chart to the Atlas data structure
  void 
  add_chart ( size_type id, const RectGeo & rect);

  /// add_chart
  ///   add a chart to the Atlas data structure
  void 
  add_chart ( size_type id, int dimension, const RectGeo & rect);

  /// listCharts
  ///   Print chart information to std::cout
  void 
  list_charts ( void ) const;

  /// numCharts
  ///   return number of charts
  uint64_t 
  numCharts ( void ) const;

  /// finalize
  ///    Finalize indexing to make Atlas data structure usable.
  void 
  finalize ( void );

  /// memory
  ///   Return memory usage
  virtual uint64_t 
  memory ( void ) const;

  /// charts
  ///   Return an iterator range which iterates through pairs (chart_id, chart)
  ChartIteratorRange
  charts ( void ) const;

private:
  // chart information
  std::unordered_map < size_type, Chart > charts_; 
  // indexing information
  std::unordered_map<size_type, size_type> chart_id_to_index_;
  std::vector<size_type> chart_index_to_id_;
  RankSelect convert_;
  // indexing methods
  GridElement 
  Chart_to_Atlas_GridElement_ ( GridElement const& chart_ge, 
                                size_type const& chart_id ) const;

  std::pair < size_type, GridElement > 
  Atlas_to_Chart_GridElement_ ( GridElement const& atlas_ge ) const;

};

inline Atlas * 
Atlas::clone ( void ) const {
  Atlas * newAtlas = new Atlas;
  for ( IdChartPair const& pair : charts () ) {
    std::shared_ptr<TreeGrid> chart_ptr ( (TreeGrid *) (pair . second -> clone ()) );
    newAtlas -> chart ( pair . first ) = chart_ptr;    
  }
  newAtlas -> finalize ();
  return newAtlas;
}


inline void 
Atlas::subdivide ( void ) { 
  for ( IdChartPair const& pair : charts () ) {
    pair . second -> subdivide ( );  
  }  
  finalize ();
}

inline Grid * 
Atlas::subgrid ( const std::deque < GridElement > & grid_elements ) const {
  std::unordered_map < size_type, std::deque < GridElement > > chart_grid_elements;
  for ( GridElement ge : grid_elements ) {
    std::pair < size_type, GridElement > atlas_ge = Atlas_to_Chart_GridElement_ ( ge );
    chart_grid_elements [ atlas_ge . first ] . push_back ( atlas_ge . second );
  }
  Atlas * newAtlas = new Atlas;
  for ( IdChartPair const& pair : charts () ) {
    Grid * subchart = pair . second -> subgrid ( chart_grid_elements [ pair . first ] );
    newAtlas -> charts_ [ pair . first ] = 
      std::shared_ptr<TreeGrid> ( (TreeGrid *) subchart );
  }  
  newAtlas -> finalize ();
  return (Grid *) newAtlas;
}

inline std::vector<Grid::GridElement> 
Atlas::subset ( const Grid & other ) const {
  const Atlas & otherAtlas = dynamic_cast<const Atlas &> (other);
  std::vector<Grid::GridElement> result;
  for ( IdChartPair const& pair : charts () ) {
    std::vector<Grid::GridElement> chart_subset = 
      pair . second -> subset ( * otherAtlas . charts_ . find ( pair . first ) -> second );
    for ( Grid::GridElement ge : chart_subset ) {
      result . push_back ( Chart_to_Atlas_GridElement_ ( ge, pair . first ) );
    }
  }
  return result;
}

inline std::shared_ptr<Geo> 
Atlas::geometry ( Grid::GridElement ge ) const {
  std::pair < size_type, GridElement > chartge;
  chartge = Atlas_to_Chart_GridElement_ ( ge );
  RectGeo rect = * std::dynamic_pointer_cast < RectGeo > 
    ( charts_ . find ( chartge . first ) -> second -> geometry ( chartge . second ) );
  return std::shared_ptr<Geo> ( new AtlasGeo ( chartge.first, rect ) );
}

inline std::vector<Grid::GridElement>
Atlas::cover ( const Geo & geo ) const { 
  const AtlasGeo & atlas_geo = dynamic_cast<const AtlasGeo &> ( geo );
  std::vector<Grid::GridElement> result;
  size_type chart_id_of_geo = atlas_geo . id ();
  const Chart & chart_of_geo = charts_ . find ( chart_id_of_geo ) -> second;
  if ( chart_of_geo -> size () == 0 ) return result;
  std::vector < GridElement > listge = chart_of_geo -> cover ( atlas_geo . rect() );
  for ( Grid::GridElement chart_ge : listge ) {
    GridElement newge = Chart_to_Atlas_GridElement_ ( chart_ge , atlas_geo . id() );
    result . push_back ( newge );
  }
  return result;
}

inline void 
Atlas::add_chart ( size_type id, const RectGeo & rect ) {
  charts_ [ id ] = std::shared_ptr<TreeGrid> ( new PointerGrid );
  charts_ [ id ] -> initialize ( rect );
}

inline void 
Atlas::add_chart ( size_type id, int dimension, const RectGeo & rect ) {
  charts_ [ id ] = std::shared_ptr<TreeGrid> ( new PointerGrid );
  charts_ [ id ] -> initialize ( rect );
  charts_ [ id ] -> dimension  ( ) = dimension;
}

inline void 
Atlas::list_charts ( void ) const {
  std::cout << "\nList of charts :\n";
  for ( IdChartPair const& pair : charts () ) {
    std::cout << "index = " << pair . first << " , ";
    std::cout << "bounds = " << pair . second -> bounds ( ) << " , "; 
    std::cout << "number of GridElements = " << pair . second -> size ( ) << "\n";
  }  
}

inline uint64_t 
Atlas::numCharts ( void ) const {
  return charts_ . size ();
}

inline void 
Atlas::importCharts ( const char * inputfile ) {
  using boost::property_tree::ptree;
  ptree pt;
  std::ifstream input ( inputfile );
  read_xml(input, pt);

  unsigned int dimension = pt.get<int>("atlas.dimension");
  std::cout << "Dimension : " << dimension << "\n";

  std::vector < double > lower_bounds, upper_bounds;
  lower_bounds . resize ( dimension );
  upper_bounds . resize ( dimension );

  for ( ptree::value_type & v : pt.get_child("atlas.listcharts") ) {
    // extract the strings
    std::string idstr = v . second . get_child ( "id" ) . data ( );
    std::string lbstr = v . second . get_child ( "lbounds" ) . data ( );
    std::string ubstr = v . second . get_child ( "ubounds" ) . data ( );
    std::stringstream idss ( idstr );
    std::stringstream lbss ( lbstr );
    std::stringstream ubss ( ubstr );
    size_type id;
    idss >> id;
    for ( unsigned int d = 0; d < dimension; ++ d ) {
      lbss >> lower_bounds [ d ];
      ubss >> upper_bounds [ d ];
    }
    // add the new chart 
    add_chart ( id, RectGeo(dimension,lower_bounds,upper_bounds) ); 
  }
  finalize ();
}

inline Atlas::ChartIteratorRange
Atlas::charts ( void ) const {
  return boost::make_iterator_range ( charts_ . begin (), charts_ . end () );
}

inline Atlas::Chart & 
Atlas::chart ( size_type chart_id ) {
  return charts_ [ chart_id ];
}

inline const Atlas::Chart & 
Atlas::chart ( size_type chart_id ) const {
  return charts_ . find ( chart_id ) -> second;
}

inline void 
Atlas::clear ( void ) {
  charts_ . clear ();
  finalize ();
}

inline uint64_t 
Atlas::memory ( void ) const {
  uint64_t result = 0;
  for ( IdChartPair const& chartpair : charts () ) {
    result += chartpair . second -> memory ();
  }
  return result;
}

inline void 
Atlas::finalize ( void ) { 
  chart_id_to_index_ . clear ();
  chart_index_to_id_ . clear ();
  size_ = 0;
  size_type chart_index = 0;
  for ( IdChartPair const& pair : charts () ) {
    size_type chart_id = pair . first; 
    Chart const& chart = pair . second;   
    size_type chart_size = chart -> size ();
    if ( chart_size == 0 ) continue;
    size_ += chart_size;
    chart_id_to_index_ [ chart_id ] = chart_index ++;
    chart_index_to_id_ . push_back ( chart_id );
  }
  std::vector<bool> bits ( size_ );
  size_type s = 0;
  for ( size_type chart_index = 0; chart_index < chart_index_to_id_ . size (); ++ chart_index ) {
    bits [ s ] = 1;
    s += charts_ [ chart_index_to_id_ [ chart_index ] ] -> size ();
  }
  convert_ . assign ( bits );
}

inline Atlas::GridElement 
Atlas::Chart_to_Atlas_GridElement_ ( GridElement const& chart_ge, 
                                     size_type const& chart_id ) const {
  size_type chart_index = chart_id_to_index_ . find ( chart_id ) -> second;
  return convert_ . select ( chart_index ) + chart_ge;
}

inline std::pair < Atlas::size_type, Atlas::GridElement > 
Atlas::Atlas_to_Chart_GridElement_ ( GridElement const& atlas_ge ) const {
  Atlas::size_type chart_index = convert_ . rank ( atlas_ge + 1 ) - 1;
  Atlas::GridElement chart_ge = atlas_ge - convert_ . select ( chart_index );
  return std::make_pair ( chart_index_to_id_[chart_index], chart_ge );
}

#endif

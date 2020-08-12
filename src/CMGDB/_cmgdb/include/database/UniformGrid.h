#ifndef UNIFORMGRID_H
#define UNIFORMGRID_H

#include <iostream>
#include <stdint.h>
#include <exception>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>
#include <boost/foreach.hpp>
#include <boost/iterator/counting_iterator.hpp>
#include <boost/unordered_map.hpp>
#include <memory>
#include "boost/serialization/serialization.hpp"
#include "boost/serialization/vector.hpp"
#include "boost/serialization/export.hpp"
#include "boost/serialization/shared_ptr.hpp"
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include "Grid.h"
#include "Geo.h"
#include "RectGeo.h"

/// UniformGrid
class UniformGrid : public Grid { 
public:
	typedef uint64_t GridElement;
  typedef boost::counting_iterator < GridElement > iterator;
  typedef iterator const_iterator;
  typedef uint64_t size_type;

  // Contructor/ Desctructor
  UniformGrid ( void ) { }
  virtual ~UniformGrid ( void ) { }

  // Builders
  void initialize ( const RectGeo & bounds,
                    const std::vector<uint64_t> & sizes );

  void initialize ( const RectGeo & bounds,
                    const std::vector<uint64_t> & sizes,
                    const std::vector<bool> & periodic );
  // General Methods
  virtual UniformGrid * clone ( void ) const;
  virtual void subdivide ( void );
  virtual Grid * subgrid ( const std::deque < GridElement > & grid_elements ) const;
  virtual std::vector<GridElement> subset ( const Grid & other ) const;
  virtual std::shared_ptr<Geo> geometry ( GridElement ge ) const;  
  virtual std::vector<Grid::GridElement> cover ( const Geo & geo ) const;
  using Grid::geometry;
  using Grid::cover;
  virtual uint64_t memory ( void ) const;

  // Features
  RectGeo & bounds ( void );
  const RectGeo & bounds ( void ) const;
  std::vector < uint64_t > & sizes ( void );
  const std::vector < uint64_t > & sizes ( void ) const;
  uint64_t width ( int d ) const;
  int dimension ( void ) const;
private:
  RectGeo bounds_;
  std::vector<uint64_t> sizes_;
  std::vector<uint64_t> multipliers_;
  int dimension_;

  friend class boost::serialization::access;
  template<typename Archive>
  void serialize(Archive & ar, const unsigned int file_version) {
    ar & boost::serialization::base_object<Grid>(*this);
    ar & bounds_;
    ar & sizes_;
    ar & multipliers_;
    ar & dimension_;
  }
};

BOOST_CLASS_EXPORT_KEY(UniformGrid);

inline void UniformGrid::initialize ( const RectGeo & bounds,
                                      const std::vector<uint64_t> & sizes,
                                      const std::vector<bool> & periodic ) {
  initialize ( bounds, sizes );
  // TODO: add support for periodicity
}
inline void UniformGrid::initialize ( const RectGeo & bounds,
                                      const std::vector<uint64_t> & sizes ) {
  bounds_ = bounds;
  sizes_ = sizes;
  dimension_ = bounds . lower_bounds . size ();
  multipliers_. resize ( dimension (), 1 );
  for ( int d = 1; d < dimension (); ++ d ) {
    multipliers_ [ d ] = sizes_ [ d - 1 ] * multipliers_ [ d - 1 ];
  }
  size_ = multipliers_ [ dimension () - 1 ] * sizes_ [ dimension () - 1 ];
  //std::cout << "UniformGrid::initialize. bounds set to " << bounds_ << "\n";
}

inline UniformGrid * UniformGrid::clone ( void ) const {
  UniformGrid * newUniformGrid = new UniformGrid;
  // TODO -- not needed immediately
  throw std::logic_error ("UniformGrid::clone not yet implemented.\n");
  return newUniformGrid;
}


inline void UniformGrid::subdivide ( void ) { 
  // TODO -- not needed immediately
  return;
}

inline Grid * UniformGrid::subgrid ( const std::deque < GridElement > & grid_elements ) const {
  UniformGrid * newUniformGrid = new UniformGrid;
  // TODO -- not needed immediately.
  return (Grid *) newUniformGrid;
}

inline std::vector<Grid::GridElement> 
UniformGrid::subset ( const Grid & other ) const {
  //const UniformGrid & otherUniformGrid = dynamic_cast<const UniformGrid &> (other);
  std::vector<Grid::GridElement> result;
  // TODO -- not needed immediately.
  return result;
}

inline std::shared_ptr<Geo> 
UniformGrid::geometry ( Grid::GridElement ge ) const {
  std::shared_ptr<RectGeo> result ( new RectGeo ( dimension () ) );
  uint64_t address = (uint64_t) ge;
  std::vector<uint64_t> coordinates ( dimension () );
  for ( int d = 0; d < dimension (); ++ d ) {
    coordinates [ d ] = address % sizes_ [ d ];
    address -= coordinates [ d ];
    address /= sizes_ [ d ];
  }
  for ( int d = 0; d < dimension (); ++ d ) {
    result -> lower_bounds [ d ] = 
      bounds_.lower_bounds[d]+((double)coordinates[d])/(double)sizes_[d]
      *(bounds_.upper_bounds[d]-bounds_.lower_bounds[d]);
    result -> upper_bounds [ d ] = 
      bounds_.lower_bounds[d]+((double)coordinates[d] + 1.0)/(double)sizes_[d]
      *(bounds_.upper_bounds[d]-bounds_.lower_bounds[d]);
  }
  return std::dynamic_pointer_cast<Geo> ( result );
} /* UniformGrid::geometry */

inline std::vector<Grid::GridElement>
UniformGrid::cover ( const Geo & geo ) const { 
  const RectGeo & rect = dynamic_cast<const RectGeo &> ( geo );
  //std::cout << "UniformGrid::cover ( " << rect << " ):\n";

  std::vector<Grid::GridElement> result;
  std::vector<int64_t> lower_coordinates ( dimension () );
  std::vector<int64_t> upper_coordinates ( dimension () );
  uint64_t address = 0;
  //std::cout << "   ";
  for ( int d = 0; d < dimension (); ++ d ) {
    lower_coordinates [ d ] = (int64_t) std::ceil ( (double) width ( d ) *
                              (rect.lower_bounds[d]-bounds_.lower_bounds[d])/
                              (bounds_.upper_bounds[d]-bounds_.lower_bounds[d]) - 1.0);
    upper_coordinates [ d ] = (int64_t) std::floor ( (double) width ( d ) *
                              (rect.upper_bounds[d]-bounds_.lower_bounds[d])/
                              (bounds_.upper_bounds[d]-bounds_.lower_bounds[d]) + 1.0 );
    if ( lower_coordinates [ d ] < 0 ) 
        lower_coordinates [ d ] = 0;
    if ( upper_coordinates [ d ] > (int64_t) sizes_ [ d ] ) 
        upper_coordinates [ d ] = (int64_t) sizes_ [ d ];

    address += multipliers_ [ d ] * lower_coordinates [ d ];
    //if ( d != 0 ) std::cout << " x ";
    //std::cout << "[" << lower_coordinates[d]<<", "<<upper_coordinates[d]<<")";
    //std::cout << "[" << lower_coordinates[d]<<", "<<upper_coordinates[d]<<")";
  }

  //std::cout << "\n";
  bool empty = false;
  for ( int d = 0; d < dimension (); ++ d ) {
    if ( upper_coordinates[d] <= lower_coordinates[d] ) empty = true;
  }
  if ( empty ) return result;
  
  std::vector<int64_t> coordinates = lower_coordinates;
  bool finished = false;
  while ( not finished ) {
    result . push_back ( Grid::GridElement ( address ) );
    finished = true;
    for ( int d = 0; d < dimension (); ++ d ) {
      ++ coordinates [ d ];
      address += multipliers_ [ d ];
      if ( coordinates [ d ] == upper_coordinates [ d ] ) {
        address -= (upper_coordinates[d]-lower_coordinates[d])*multipliers_[d];
        coordinates [ d ] = lower_coordinates [ d ];
      } else {
        finished = false;
        break;
      }
    }
  }
  return result;
}

inline uint64_t UniformGrid::memory ( void ) const {
  uint64_t result = 0;
  // TODO -- not needed immediately
  return result;
}

// Features

inline RectGeo & 
UniformGrid::bounds ( void ) {
  return bounds_;
}

inline const RectGeo & 
UniformGrid::bounds ( void ) const {
  return bounds_;
}

inline std::vector < uint64_t > & 
UniformGrid::sizes ( void ) {
  return sizes_;
}

inline const std::vector < uint64_t > & 
UniformGrid::sizes ( void ) const {
  return sizes_;
}

inline uint64_t 
UniformGrid::width ( int d ) const {
  return sizes_ [ d ];
}

inline int 
UniformGrid::dimension ( void ) const {
  return dimension_;
}

#endif

#ifndef EDGEGRID_H
#define EDGEGRID_H

#include <iostream>
#include <stdint.h>

#include <fstream>
#include <sstream>
#include <string>
#include <cmath>
#include <utility>

#include <boost/foreach.hpp>
#include <boost/iterator/counting_iterator.hpp>
#include <boost/unordered_map.hpp>
#include <memory>

#include "boost/serialization/serialization.hpp"
#include "boost/serialization/vector.hpp"
#include "boost/serialization/map.hpp"
#include "boost/serialization/export.hpp"
#include "boost/serialization/shared_ptr.hpp"

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include "Grid.h"
#include "Geo.h"
#include "RectGeo.h"


/*

Layout of EdgeGrid
only has 1-dimensional pieces. They are the edges in the 
s[1]*s[2]*...*s[d] cubical structure.

The edges with extent in the k dimension are
 laid out in a (sizes[1]+1)*(sizes[2]+1)*...*(sizes[k])*...*(sizes[d]+1) grid

The main difficulty is the contiguity requirement while using only O(d) memory
for the structure rather than O(N), where N is the number of edges.

Our strategy is to produce a contiguous indexing space, and an address space
which are distinct. We provide methods
gridElementToAddress
and 
addressToGridElement
to convert back and forrth.

The contiguous indexing is based on putting the edges of each type in consecutive 
blocks. That is, we first have all the edges with extent in the first dimension, in
the dictionary order based on their midpoints, followed by the edges with extent in
the second dimension, etc. 

The addressing scheme is based on naming edges according to a slightly larger grid with
size (s[1]+1)(s[2]+1)...(s[d]+1) and not using some of the grid.

To perform geometry, we obtain the address from the index, and then use the address
to construct the RectGeo.

To perform the cover method, we obtain addresses from the input RectGeo, and then
convert them into indices to return.

We need "sizes_" in the code to be sizes[d]+1 above. This more or less sorts out
"geometry"

"cover" is a bit harder. Here we have line 209-223 as the first part.
With the alternative definition of sizes_ this may be OK provided that
width(k) now returns sizes_[k] -1. 
This enables the loops to go from [0, width(k)]

There are now special cases to consider. If we are on the lower bounds in 
dimension k, we may only include k-edges. If we are on the lower bounds of
multiple dimensions simultaneously, no edges are allowed.

The only other issue is if we ask to return an edge on the upper-right collar we
have added. But we can allow addressToGridElement to return end() and we can 
check for this.

We have to be careful with multipliers_ -- it should reflect the logic of
the new augmented sizes_. There need to be changed in initialize_ to accomplish this.

Also initialize_ needs to produce the edgeSearch map

TODO: make this more efficient. Why use a weird address space when you
      have to use coordinates to do the conversion anyway?

*/
class EdgeGrid : public Grid { 

public:
	typedef uint64_t GridElement;
  typedef boost::counting_iterator < GridElement > iterator;
  typedef iterator const_iterator;
  typedef uint64_t size_type;

  // Contructor/ Desctructor
  EdgeGrid ( void ) { }
  virtual ~EdgeGrid ( void ) { }

  // Builders
  void initialize ( const RectGeo & bounds,
                    const std::vector<uint64_t> & sizes );

  void initialize ( const RectGeo & bounds,
                    const std::vector<uint64_t> & sizes,
                    const std::vector<bool> & periodic );
  // General Methods
  virtual EdgeGrid * clone ( void ) const;
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

  std::vector<uint64_t> start_;
  std::map<GridElement, std::pair<uint64_t, int> > search_;

  std::pair<uint64_t, int> gridElementToAddress ( const GridElement & ge ) const;
  GridElement addressToGridElement ( uint64_t address, int dim ) const;

  friend class boost::serialization::access;
  template<typename Archive>
  void serialize(Archive & ar, const unsigned int file_version) {
    ar & boost::serialization::base_object<Grid>(*this);
    ar & bounds_;
    ar & sizes_;
    ar & multipliers_;
    ar & dimension_;
    ar & start_;
    ar & search_;
  }
};

BOOST_CLASS_EXPORT_KEY(EdgeGrid);

inline void EdgeGrid::initialize ( const RectGeo & bounds,
                                      const std::vector<uint64_t> & sizes,
                                      const std::vector<bool> & periodic ) {
  initialize ( bounds, sizes );
  // TODO: add support for periodicity
}
inline void EdgeGrid::initialize ( const RectGeo & bounds,
                                      const std::vector<uint64_t> & sizes ) {
  bounds_ = bounds;
  sizes_ = sizes;
  dimension_ = bounds . lower_bounds . size ();
  multipliers_. resize ( dimension (), 1 );
  ++ sizes_ [ 0 ];
  for ( int d = 1; d < dimension (); ++ d ) {
    ++ sizes_ [ d ];
    multipliers_ [ d ] = sizes_ [ d - 1 ] * multipliers_ [ d - 1 ];
  }
  uint64_t address_size = multipliers_ [ dimension () - 1 ] * sizes_ [ dimension () - 1 ];
  size_ = 0;
  start_ . resize ( dimension (), 0 );
  for ( int d = 0; d < dimension (); ++ d ) {
    start_ [ d ] = size_;
    size_ += (address_size / sizes_[d]) * (sizes_[d]-1);
    search_ [ GridElement(size_) ] = std::make_pair ( start_ [ d ], d );
    std::cout << "start_[" << d << "] = " << start_[d] << "    size_ = " << size_ << "\n";

  }

  // DEBUG
  std::cout << "EdgeGrid::initialize. bounds set to " << bounds_ << "\n";
  for ( uint64_t i = 0; i < size (); ++ i ) {
    RectGeo geo = * std::dynamic_pointer_cast<RectGeo> ( geometry ( GridElement (i) ) );
    std::cout << "GridElement " << i << " has geometry " << geo << "\n";
    std::pair<uint64_t, int> address_pair = gridElementToAddress ( GridElement (i) );
    std::cout << "gridElementToAddress(" << i << ") = (" << address_pair . first 
       << ", " << address_pair . second << ") " << "\n";
    std::cout << "addressToGridElement(" << address_pair . first 
       << ", " << address_pair . second << ") = " 
       << addressToGridElement ( address_pair.first, address_pair.second ) << "\n";
  }

  // END DEBUG
}

inline EdgeGrid * EdgeGrid::clone ( void ) const {
  EdgeGrid * newEdgeGrid = new EdgeGrid;
  // TODO -- not needed immediately
  return newEdgeGrid;
}


inline void EdgeGrid::subdivide ( void ) { 
  // TODO -- not needed immediately
  return;
}

inline Grid * EdgeGrid::subgrid ( const std::deque < GridElement > & grid_elements ) const {
  EdgeGrid * newEdgeGrid = new EdgeGrid;
  // TODO -- not needed immediately.
  return (Grid *) newEdgeGrid;
}

inline std::vector<Grid::GridElement> 
EdgeGrid::subset ( const Grid & other ) const {
  //const EdgeGrid & otherEdgeGrid = dynamic_cast<const EdgeGrid &> (other);
  std::vector<Grid::GridElement> result;
  // TODO -- not needed immediately.
  return result;
}

inline std::shared_ptr<Geo> 
EdgeGrid::geometry ( Grid::GridElement ge ) const {
  std::shared_ptr<RectGeo> result ( new RectGeo ( dimension () ) );
  std::pair<uint64_t, int> address_pair = gridElementToAddress ( ge );
  uint64_t & address = address_pair . first;
  int & collapse_dimension = address_pair . second;
  std::vector<uint64_t> coordinates ( dimension () );
  for ( int d = 0; d < dimension (); ++ d ) {
    coordinates [ d ] = address % sizes_ [ d ];
    address -= coordinates [ d ];
    address /= sizes_ [ d ];
  }
  for ( int d = 0; d < dimension (); ++ d ) {
    result -> lower_bounds [ d ] = 
      bounds_.lower_bounds[d]+((double)coordinates[d])/(double)width(d)
      *(bounds_.upper_bounds[d]-bounds_.lower_bounds[d]);
    if ( d != collapse_dimension ) {
      result -> upper_bounds [ d ] = result -> lower_bounds [ d ];
    } else {
      result -> upper_bounds [ d ] = 
        bounds_.lower_bounds[d]+((double)coordinates[d] + 1.0)/(double)width(d)
        *(bounds_.upper_bounds[d]-bounds_.lower_bounds[d]);
    }
  }
  return std::dynamic_pointer_cast<Geo> ( result );
} /* EdgeGrid::geometry */

inline std::vector<Grid::GridElement>
EdgeGrid::cover ( const Geo & geo ) const { 
  const RectGeo & rect = dynamic_cast<const RectGeo &> ( geo );
  std::cout << "EdgeGrid::cover ( " << rect << " ):\n";

  std::vector<Grid::GridElement> result;
  std::vector<int64_t> lower_coordinates ( dimension () );
  std::vector<int64_t> upper_coordinates ( dimension () );
  uint64_t address = 0;

  std::set<int> touching_lower;
  //std::cout << "   ";
  for ( int d = 0; d < dimension (); ++ d ) {
    lower_coordinates [ d ] = (int64_t) std::ceil ( (double) width ( d ) *
                              (rect.lower_bounds[d]-bounds_.lower_bounds[d])/
                              (bounds_.upper_bounds[d]-bounds_.lower_bounds[d]) - 1.0);
    upper_coordinates [ d ] = (int64_t) std::floor ( (double) width ( d ) *
                              (rect.upper_bounds[d]-bounds_.lower_bounds[d])/
                              (bounds_.upper_bounds[d]-bounds_.lower_bounds[d]) + 1.0 );
    if ( lower_coordinates [ d ] < 0 ) { 
      lower_coordinates [ d ] = 0;
      touching_lower . insert ( d );
    }
    if ( upper_coordinates [ d ] > (int64_t) sizes_ [ d ] ) 
      upper_coordinates [ d ] = (int64_t) sizes_ [ d ];

    address += multipliers_ [ d ] * lower_coordinates [ d ];
    
    if ( d != 0 ) std::cout << " x ";
    std::cout << "[" << lower_coordinates[d]<<", "<<upper_coordinates[d]<<")";
  }

  std::cout << "\n";

  std::vector<int64_t> coordinates = lower_coordinates;
  std::set < int > low;
  for ( int d = 0; d < dimension (); ++ d ) {
    if ( not touching_lower . count ( d ) ) low . insert ( d );
  }
  bool finished = false;
  while ( not finished ) {
    // If more than one dimension is on lower bounds, no edges to cover.
    if ( low . size () < 2 ) {
      // If precisely one dimension (d, say) is on lower bounds, 
      // then cover only the d-edge 
      if ( low . size () == 1 ) {
        int d = * low . begin ();
        GridElement ge = addressToGridElement ( address, d );
        if ( ge != size () ) result . push_back ( ge );
      } else {
      // Typical case -- cover all the edges
        for ( int d = 0; d < dimension (); ++ d ) {
          GridElement ge = addressToGridElement ( address, d );
          if ( ge != size () ) result . push_back ( ge );
        }
      }
    }
    finished = true;
    for ( int d = 0; d < dimension (); ++ d ) {
      if ( coordinates [ d ] == lower_coordinates [ d ] ) low . erase ( d );
      ++ coordinates [ d ];
      address += multipliers_ [ d ];
      if ( coordinates [ d ] == upper_coordinates [ d ] ) {
        address -= (upper_coordinates[d]-lower_coordinates[d])*multipliers_[d];
        coordinates [ d ] = lower_coordinates [ d ];
        if ( not touching_lower . count ( d ) ) low . insert ( d );
      } else {
        finished = false;
        break;
      }
    }
  }
  // DEBUG
  std::cout << "Cover results: \n";
  for ( size_t i = 0; i < result . size (); ++ i ) {
    std::cout << result [ i ] << " ";
  }
  std::cout << "\n";
  // END DEBUG
  return result;
}

inline uint64_t EdgeGrid::memory ( void ) const {
  uint64_t result = 0;
  // TODO -- not needed immediately
  return result;
}

// Features

inline RectGeo & 
EdgeGrid::bounds ( void ) {
  return bounds_;
}

inline const RectGeo & 
EdgeGrid::bounds ( void ) const {
  return bounds_;
}

inline std::vector < uint64_t > & 
EdgeGrid::sizes ( void ) {
  return sizes_;
}

inline const std::vector < uint64_t > & 
EdgeGrid::sizes ( void ) const {
  return sizes_;
}

inline uint64_t 
EdgeGrid::width ( int d ) const {
  return sizes_ [ d ] - 1;
}

inline int 
EdgeGrid::dimension ( void ) const {
  return dimension_;
}

inline std::pair<uint64_t, int> 
EdgeGrid::gridElementToAddress ( const GridElement & ge ) const {
  // Our first question is what dimension is the grid element.
  std::pair<uint64_t, int> x = search_ . upper_bound ( ge ) -> second;
  uint64_t ge_address = (uint64_t) ge - x . first;
  int dim = x . second;
  // Now obtain coordinates
  std::vector<uint64_t> coordinates ( dimension () );
  for ( int d = 0; d < dimension (); ++ d ) {
    uint64_t dim_size = (d == dim) ? width(d) : sizes_[d];
    coordinates [ d ] = ge_address % dim_size;
    ge_address -= coordinates [ d ];
    ge_address /= dim_size;
  }

  // Now reconstruct the address
  uint64_t address = 0;
  uint64_t multiplier = 1;
  for ( int d = 0; d < dimension (); ++ d ) {
    address += multiplier * coordinates [ d ];
    multiplier *= sizes_ [ d ];
  }
  return std::make_pair ( address, dim );
}
  
inline Grid::GridElement 
EdgeGrid::addressToGridElement ( uint64_t address, int dim ) const {
  // Find starting grid element for the dimension
  uint64_t start = start_ [ dim ];

   // Now obtain coordinates
  std::vector<uint64_t> coordinates ( dimension () );
  for ( int d = 0; d < dimension (); ++ d ) {
    uint64_t dim_size = sizes_[d];
    coordinates [ d ] = address % dim_size;
    address -= coordinates [ d ];
    address /= dim_size;
  }

  // Now reconstruct the address
  uint64_t ge_address = 0;
  uint64_t multiplier = 1;
  for ( int d = 0; d < dimension (); ++ d ) {
    uint64_t dim_size = (dim == d) ? width (d) : sizes_ [ d ];
    if ( coordinates [ d ] >= dim_size ) return GridElement(size_);
    ge_address += multiplier * coordinates [ d ];
    multiplier *= dim_size;
  }

  return GridElement ( start + ge_address );
}

#endif

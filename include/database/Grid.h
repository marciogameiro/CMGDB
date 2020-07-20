// Grid.h
#ifndef CMDB_GRID_H
#define CMDB_GRID_H

#include <stdint.h>
#include <memory>
#include <vector>
#include <stack>
#include <deque>
#include <algorithm>
#include <iterator>
#include "boost/foreach.hpp"
#include "boost/iterator/counting_iterator.hpp"
#include "boost/unordered_set.hpp"
#include <memory>
#include "boost/serialization/serialization.hpp"
#include "boost/serialization/vector.hpp"
#include "Geo.h"
#include "UnionGeo.h"
#include "IntersectionGeo.h"

// Declaration
class Grid {
public:
  // Typedefs
  typedef uint64_t GridElement;
  typedef boost::counting_iterator < GridElement > iterator;
  typedef iterator const_iterator;
  typedef uint64_t size_type;
  
  // Virtual Methods
  virtual Grid * clone ( void ) const = 0;
  virtual void subdivide ( void ) = 0;
  virtual Grid * subgrid ( const std::deque < GridElement > & grid_elements ) const = 0;
  virtual std::vector<Grid::GridElement> subset ( const Grid & other ) const = 0;
  virtual std::shared_ptr<Geo> geometry ( GridElement ge ) const = 0;  
  virtual std::vector<Grid::GridElement> cover ( const Geo & geo ) const = 0;

  // Container methods
  iterator find ( GridElement ge ) const;
  iterator begin ( void ) const;
  iterator end ( void ) const;
  size_type size ( void ) const;
  
  /// geometry
  std::shared_ptr<Geo> 
  geometry ( const const_iterator & it ) const { return geometry ( *it ); }
  
  /// Cover (for dispatch)
  std::vector<Grid::GridElement> 
  cover ( std::shared_ptr<Geo> geo ) const;

  /// unionCover
  template < class T >
  std::vector<Grid::GridElement> 
  unionCover ( const std::vector < T > & V ) const;
  
  /// intersectionCover
  template < class T >
  std::vector<Grid::GridElement> intersectionCover ( const std::vector < T > & V  ) const;
  
  /// memory
  ///   Return memory usage of this data structure
  virtual uint64_t memory ( void ) const = 0;
  
protected: 
  Grid ( void );
public: 
  virtual ~Grid ( void );  
protected:
  size_type size_;
private:
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version) {
    ar & size_;
  }
  
};

inline Grid::iterator Grid::find ( GridElement find_me ) const {
  return iterator ( find_me );
}

inline Grid::iterator Grid::begin ( void ) const {
  return iterator ( 0 );
}

inline Grid::iterator Grid::end ( void ) const {
  return iterator ( size () );
}

inline Grid::size_type Grid::size ( void ) const {
  return size_;
}

inline std::vector<Grid::GridElement> 
Grid::cover ( std::shared_ptr<Geo> geo ) const {

  std::shared_ptr<UnionGeo> union_geo 
    = std::dynamic_pointer_cast<UnionGeo> ( geo );
  if ( union_geo ) return unionCover ( union_geo -> elements );

  std::shared_ptr<IntersectionGeo> intersect_geo 
    = std::dynamic_pointer_cast<IntersectionGeo> ( geo );
  if ( intersect_geo ) return intersectionCover ( intersect_geo -> elements );

  return cover ( * geo );
}

template < class T >
inline std::vector<Grid::GridElement> 
Grid::unionCover ( const std::vector < T > & V ) const {
  std::vector<Grid::GridElement> result, cover_vec;
  for ( T const& geo : V ) {
    cover_vec = cover ( geo );
    result . insert ( result . end (), cover_vec . begin (), cover_vec . end () );
  }

  std::sort ( result . begin (), result . end () );
  auto last = std::unique ( result.begin(), result.end() );
  result . erase(last, result.end());
  return result;
}

template < class T >
inline std::vector<Grid::GridElement> 
Grid::intersectionCover ( const std::vector < T > & V ) const {
  boost::unordered_set<Grid::GridElement> result_set;
  bool first = true;
  for ( const T & geo : V ) {
    std::vector<Grid::GridElement> cover_vec = cover ( geo );
    boost::unordered_set<Grid::GridElement> 
      cover_set ( cover_vec . begin (), cover_vec . end () );

    if ( first ) {
      std::swap ( result_set, cover_set );
    } else {
      for ( Grid::GridElement ge : result_set ) {
        if ( cover_set . count ( ge ) == 0 ) result_set . erase ( ge );
      }
    }
    first = false;
  }
  std::vector<Grid::GridElement> result ( result_set . begin (), result_set . end () );
  return result;
}

inline Grid::Grid ( void ) {
  size_ = 1;
}

inline Grid::~Grid ( void ) {
}

#endif

#ifndef CMDB_MAP_H
#define CMDB_MAP_H

#include <cstdint>
#include <memory>
#include <stdexcept>
#include <vector>
#include "Geo.h"

class Map {
public:
  virtual ~Map ( void ) {}
  virtual std::shared_ptr<Geo> operator () ( std::shared_ptr<Geo> geo ) const = 0;

  /// has_batch
  ///   Return true if this map provides an optimized batch evaluation
  ///   (batch_map). Callers that hold flat rectangle data may then evaluate
  ///   many rectangles in a single call instead of once per rectangle.
  virtual bool has_batch ( void ) const { return false; }

  /// batch_map
  ///   Evaluate the map on `count` rectangles at once.
  ///   `rects` holds count * 2 * dim doubles: for each rectangle, dim lower
  ///   bounds followed by dim upper bounds. On return `images` holds the
  ///   image rectangles in the same layout. Only called when has_batch()
  ///   returns true; the default implementation is not batched and exists
  ///   only as a specification of the interface.
  virtual void batch_map ( const std::vector<double> & /* rects */,
                           uint64_t /* count */,
                           uint64_t /* dim */,
                           std::vector<double> & /* images */ ) const {
    throw std::logic_error ( "Map::batch_map called on a map without batch support" );
  }
private:
};

#endif

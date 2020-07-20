#ifndef CMDB_MAP_H
#define CMDB_MAP_H

#include <memory>
#include "Geo.h"

class Map {
public:
  virtual ~Map ( void ) {}
  virtual std::shared_ptr<Geo> operator () ( std::shared_ptr<Geo> geo ) const = 0;
private:
};

#endif

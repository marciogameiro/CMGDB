#ifndef CMDB_UNIONGEO_H
#define CMDB_UNIONGEO_H

#include "Geo.h"

class UnionGeo : public Geo {
public:
  std::vector< std::shared_ptr<Geo> > elements;
  void insert ( std::shared_ptr < Geo > geo ) {
    elements . push_back ( geo );
  }
private:
  virtual void print ( std::ostream & stream ) const;
};

inline void 
UnionGeo::print ( std::ostream & stream ) const {
  stream << "(UnionGeo={";
  for ( size_t i = 0; i < elements . size (); ++ i ) {
    stream << * elements [ i ] << ", ";
  }
  stream << "})";
}

#endif

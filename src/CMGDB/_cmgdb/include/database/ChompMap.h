#ifndef CMDB_CHOMP_MAP
#define CMDB_CHOMP_MAP

#include "RectGeo.h"
#include "Map.h"
#include <memory>

class ChompMap {
public:
  ChompMap ( std::shared_ptr<const Map> cmdb_map_ ) : cmdb_map_ ( cmdb_map_ ) {}
  chomp::Rect operator () ( const chomp::Rect & rect ) const {
    std::shared_ptr<Geo> geo ( new RectGeo ( rect ) );
    std::shared_ptr<Geo> val = (*cmdb_map_) ( geo );
    RectGeo image = * std::dynamic_pointer_cast<RectGeo> ( val );
    return image; 
  }
  std::shared_ptr<Geo> operator () ( const std::shared_ptr<Geo> & geo ) const {
    return (*cmdb_map_) ( geo ); 
  }
private:
  std::shared_ptr<const Map> cmdb_map_;
};

#endif

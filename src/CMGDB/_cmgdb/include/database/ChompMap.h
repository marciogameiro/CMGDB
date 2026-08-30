#ifndef CMDB_CHOMP_MAP
#define CMDB_CHOMP_MAP

#include "RectGeo.h"
#include "Map.h"
#include <algorithm>
#include <cstdint>
#include <memory>
#include <vector>

/// ChompMap
///   Adapts a database Map for the chomp homology code. Besides the scalar
///   rectangle-to-rectangle evaluation, it exposes `images`, which evaluates
///   a whole vector of rectangles with each rectangle evaluated exactly once
///   -- through Map::batch_map in chunks when the map provides a batched
///   evaluator, one scalar call per rectangle otherwise. The Conley-index
///   code gathers its evaluation points and calls `images` so that a
///   Python-defined or batched map pays one call per chunk instead of one
///   call per box (see ConleyIndex.h and RelativeMapHomology.h).
class ChompMap {
public:
  ChompMap ( std::shared_ptr<const Map> cmdb_map_,
             uint64_t batch_chunk_size = 65536 )
    : cmdb_map_ ( cmdb_map_ ), chunk_size_ ( batch_chunk_size ) {}
  chomp::Rect operator () ( const chomp::Rect & rect ) const {
    std::shared_ptr<Geo> geo ( new RectGeo ( rect ) );
    std::shared_ptr<Geo> val = (*cmdb_map_) ( geo );
    RectGeo image = * std::dynamic_pointer_cast<RectGeo> ( val );
    return image;
  }
  std::shared_ptr<Geo> operator () ( const std::shared_ptr<Geo> & geo ) const {
    return (*cmdb_map_) ( geo );
  }

  /// images
  ///   Evaluate the map on every rectangle of `rects`, in order. The batched
  ///   path uses the flat-buffer layout of Map::batch_map (per rectangle,
  ///   dim lower bounds then dim upper bounds) in chunks of the configured
  ///   chunk size (0 = one chunk). The batched evaluator is required to
  ///   agree with the scalar map on every rectangle (see Model::set_batch_map),
  ///   so both paths produce identical results.
  std::vector<chomp::Rect> images ( const std::vector<chomp::Rect> & rects ) const {
    std::vector<chomp::Rect> result;
    result . reserve ( rects . size () );
    if ( rects . empty () || not cmdb_map_ -> has_batch () ) {
      for ( const chomp::Rect & rect : rects ) {
        result . push_back ( (*this) ( rect ) );
      }
      return result;
    }
    const uint64_t n = rects . size ();
    const uint64_t dim = rects [ 0 ] . lower_bounds . size ();
    const uint64_t chunk = chunk_size_ > 0 ? chunk_size_ : n;
    std::vector<double> flat;
    std::vector<double> flat_images;
    for ( uint64_t start = 0; start < n; start += chunk ) {
      const uint64_t stop = std::min ( start + chunk, n );
      const uint64_t count = stop - start;
      flat . clear ();
      flat . reserve ( count * 2 * dim );
      for ( uint64_t k = start; k < stop; ++ k ) {
        flat . insert ( flat . end (), rects [ k ] . lower_bounds . begin (),
                        rects [ k ] . lower_bounds . end () );
        flat . insert ( flat . end (), rects [ k ] . upper_bounds . begin (),
                        rects [ k ] . upper_bounds . end () );
      }
      cmdb_map_ -> batch_map ( flat, count, dim, flat_images );
      for ( uint64_t k = 0; k < count; ++ k ) {
        const double * bounds = flat_images . data () + k * 2 * dim;
        chomp::Rect image ( dim );
        for ( uint64_t d = 0; d < dim; ++ d ) {
          image . lower_bounds [ d ] = bounds [ d ];
          image . upper_bounds [ d ] = bounds [ dim + d ];
        }
        result . push_back ( std::move ( image ) );
      }
    }
    return result;
  }

  /// images (Geo overload)
  ///   The same single-evaluation-per-entry contract for geometries held as
  ///   Geo pointers (the Grid::geometry return type). Rectangle geometries
  ///   ride the batched path; any other geometry type falls back to one
  ///   scalar evaluation per entry.
  std::vector<std::shared_ptr<Geo>>
  images ( const std::vector<std::shared_ptr<Geo>> & geos ) const {
    std::vector<std::shared_ptr<Geo>> result;
    result . reserve ( geos . size () );
    bool rect_geometry = not geos . empty () &&
      ( std::dynamic_pointer_cast<RectGeo> ( geos [ 0 ] ) != nullptr );
    if ( not rect_geometry || not cmdb_map_ -> has_batch () ) {
      for ( const std::shared_ptr<Geo> & geo : geos ) {
        result . push_back ( (*cmdb_map_) ( geo ) );
      }
      return result;
    }
    std::vector<chomp::Rect> rects;
    rects . reserve ( geos . size () );
    for ( const std::shared_ptr<Geo> & geo : geos ) {
      std::shared_ptr<RectGeo> rect = std::dynamic_pointer_cast<RectGeo> ( geo );
      if ( not rect ) {
        throw std::logic_error ( "ChompMap::images. Mixed geometry types.\n" );
      }
      rects . push_back ( chomp::Rect ( * rect ) );
    }
    std::vector<chomp::Rect> rect_images = images ( rects );
    for ( chomp::Rect & image : rect_images ) {
      result . push_back ( std::shared_ptr<Geo> ( new RectGeo ( image ) ) );
    }
    return result;
  }
private:
  std::shared_ptr<const Map> cmdb_map_;
  uint64_t chunk_size_;
};

#endif

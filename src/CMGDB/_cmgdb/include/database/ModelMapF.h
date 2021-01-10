
#ifndef CMDP_MODELMAPF_H
#define CMDP_MODELMAPF_H

#include "Map.h"
#include "EuclideanParameterSpace.h"
#include "RectGeo.h"
#include "simple_interval.h"
#include <memory>
#include <vector>
#include <algorithm>

class ModelMapF : public Map {
public:
  typedef simple_interval<double> interval;

  // User interface: method to be provided by user

  // Map F
  std::function<std::vector<double>(std::vector<double>)> F;

  // Parameter variable
  // Not using parameters, but leave here for now
  interval p0;
  
  // Constructor: sets parameter variables
  void assign ( RectGeo const& rectangle,
                std::function<std::vector<double>(std::vector<double>)> const& F_map ) {
    // Set map F
    F = F_map;

    // Read parameter intervals from input rectangle
    p0 = getRectangleComponent ( rectangle, 0 );
  }

  // Map
  // Must take a rectangle as input and return the image
  // as a rectangle.
  //
  // Define how to evaluate the input map. The input map
  // is a map defined on rectangles defined as a product
  // of intervals. The intervals are represented by a list
  // of bounds, first the lower bounds of each variables
  // followed by the upper bounds of the varibales.
  RectGeo operator () ( const RectGeo & rectangle ) const {
    // Get phase space dimension
    uint64_t dim = rectangle . dimension ();

    // Lower and upper bounds of rectangle
    std::vector<double> rect_bounds (2 * dim, 0.0);

    for ( int d = 0; d < dim; ++ d ) {
      // Get lower and upper bounds in dimension d
      rect_bounds [ d ] = rectangle . lower_bounds [ d ];
      rect_bounds [ dim + d ] = rectangle . upper_bounds [ d ];
    }

    // Compute lower and upper bounds of image of
    // rectangle by F by evaluating the map F
    std::vector<double> image_bounds = F (rect_bounds);

    // Rectangle for image of F
    RectGeo rect_image ( dim );

    for ( int d = 0; d < dim; ++ d ) {
      // Assign lower and upper values to image rectangle
      rect_image . lower_bounds [ d ] = image_bounds [ d ];
      rect_image . upper_bounds [ d ] = image_bounds [ dim + d ];
    }

    // Return result
    return rect_image;
  }

 // Program interface (methods used by program)
  ModelMapF ( std::shared_ptr<Parameter> parameter,
              std::function<std::vector<double>(std::vector<double>)> const& F ) {
    const RectGeo & rectangle = 
      * std::dynamic_pointer_cast<EuclideanParameter> ( parameter ) -> geo;
    assign ( rectangle, F );
  }

  std::shared_ptr<Geo>
  operator () ( std::shared_ptr<Geo> geo ) const {
    return std::shared_ptr<Geo> ( new RectGeo (
        operator () ( * std::dynamic_pointer_cast<RectGeo> ( geo ) ) ) );
  }
private:
  interval getRectangleComponent ( const RectGeo & rectangle, int d ) const {
    return interval (rectangle . lower_bounds [ d ], rectangle . upper_bounds [ d ]);
  }
};

#endif

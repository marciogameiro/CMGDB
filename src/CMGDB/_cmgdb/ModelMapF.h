
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
  RectGeo operator () ( const RectGeo & rectangle ) const {
    // Convert input to intervals
    interval x0 = getRectangleComponent ( rectangle, 0 );
    interval x1 = getRectangleComponent ( rectangle, 1 );

    double x0_min = x0 . lower ();
    double x0_max = x0 . upper ();
    double x1_min = x1 . lower ();
    double x1_max = x1 . upper ();

    // Get the corner points
    std::vector<double> u1 {x0_min, x1_min};
    std::vector<double> u2 {x0_max, x1_min};
    std::vector<double> u3 {x0_max, x1_max};
    std::vector<double> u4 {x0_min, x1_max};

    // Evaluate F at the corner points
    std::vector<double> v1 = F(u1);
    std::vector<double> v2 = F(u2);
    std::vector<double> v3 = F(u3);
    std::vector<double> v4 = F(u4);

    // Get min and max values
    double y0_min = std::min({v1[0], v2[0], v3[0], v4[0]});
    double y0_max = std::max({v1[0], v2[0], v3[0], v4[0]});
    double y1_min = std::min({v1[1], v2[1], v3[1], v4[1]});
    double y1_max = std::max({v1[1], v2[1], v3[1], v4[1]});

    // Get half the grid sizes
    // double h0 = (x0_max - x0_min) / 2.0;
    // double h1 = (x1_max - x1_min) / 2.0;

    // Get grid sizes
    double h0 = x0_max - x0_min;
    double h1 = x1_max - x1_min;

    // Set intervals
    interval y0 = interval(y0_min - h0, y0_max + h0);
    interval y1 = interval(y1_min - h1, y1_max + h1);
    // interval y0 = interval(y0_min, y0_max);
    // interval y1 = interval(y1_min, y1_max);

    // Return result
    return makeRectangle ( y0, y1 );
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
  RectGeo makeRectangle ( interval const& y0, interval const& y1 ) const {
    RectGeo return_value ( 2 );
    return_value . lower_bounds [ 0 ] = y0 . lower ();
    return_value . upper_bounds [ 0 ] = y0 . upper ();
    return_value . lower_bounds [ 1 ] = y1 . lower ();
    return_value . upper_bounds [ 1 ] = y1 . upper ();
    return return_value;
  }
};

#endif

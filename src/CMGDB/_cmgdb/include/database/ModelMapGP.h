
#ifndef CMDP_MODELMAPGP_H
#define CMDP_MODELMAPGP_H

#include "Map.h"
#include "EuclideanParameterSpace.h"
#include "RectGeo.h"
#include "simple_interval.h"
#include "libgp/gp.h"
#include <memory>
#include <vector>
#include <algorithm>

class ModelMapGP : public Map {
public:
  typedef simple_interval<double> interval;

  // Parameter variable
  // Not using parameters, but leave here for now
  interval p0;

  // Map
  // Must take a rectangle as input and return the image
  // as a rectangle.
  //
  // Define how to evaluate the GP map. The GP map
  // is a map defined on points, so we evaluate the map
  // at each corner of the rectangle and put a bounding
  // box around their images. We then add a grid sized
  // layer of boxes around the bounding box to construct
  // the rectangle representing the image of the map.
  RectGeo operator () ( const RectGeo & rectangle ) const {
    // Values of F at the corner points
    std::vector<std::vector<double>> F_y;

    // The vector verts_cube_ contains the vertices
    // of the unit cube as vectors with entries in
    // {0, 1}. We use these vertices of the unit
    // cube to get the vertices of the corners of
    // the rectangle.
    for ( auto v : verts_cube_ ) { // For each vertex
      std::vector<double> x;       // Corner point x
      for ( int d = 0; d < phase_dim_; ++ d ) {
        // Get lower and upper bounds in dimension d
        double x_lower = rectangle . lower_bounds [ d ];
        double x_upper = rectangle . upper_bounds [ d ];

        // If v[0]==0 get lower bound, if v[d]==1 get upper bound
        double x_d = (v [ d ] == 0) ? x_lower : x_upper;
        x . push_back ( x_d ); // Push to vector x
      }
      // Evaluate F at x and push to F_y
      std::vector<double> y;
      for ( int d = 0; d < phase_dim_; ++ d ) {
        y . push_back ( GPs[d] . f( x . data () ) );
      }
      F_y . push_back ( y );
    }

    // Rectangle for image of F
    RectGeo rect_image ( phase_dim_ );

    // Get the lower and upper bounds for the bounding box of F_y,
    // that is, get the min and max for each "column" and make a
    // vector with min values and another with the max values.
    for ( int d = 0; d < phase_dim_; ++ d ) {
      // Grid size for bounding box padding
      double grid_size = rectangle . upper_bounds [ d ] - rectangle . lower_bounds [ d ];

      // Comparisson function (compare entry d)
      auto comp = [&] (const std::vector<double> &v1, const std::vector<double> &v2) {
        return v1 [d] < v2 [d];
      };

      // Get min and max vectors (vectors with min and max at entry d)
      const auto min_max_v = std::minmax_element ( F_y . begin(), F_y . end (), comp );

      // Get lower and upper bounds in dimension d
      // Increase bounding box by grid size
      // double y_lower = (*min_max_v . first) [d];
      // double y_upper = (*min_max_v . second) [d];
      // double y_lower = (*min_max_v . first) [d] - grid_size;
      // double y_upper = (*min_max_v . second) [d] + grid_size;
      double y_lower = (*min_max_v . first) [d] - 1.0 * grid_size;
      double y_upper = (*min_max_v . second) [d] + 1.0 * grid_size;

      // Assign the values to image rectangle
      rect_image . lower_bounds [ d ] = y_lower;
      rect_image . upper_bounds [ d ] = y_upper;
    }

    // Return result
    return rect_image;
  }

  // Model initialization
  ModelMapGP ( std::shared_ptr<Parameter> parameter,
               std::vector<std::vector<double>> const& X_train,
               std::vector<std::vector<double>> const& Y_train ) :
    // Initilize phase space dimension
    phase_dim_ (X_train[0] . size()) {
    const RectGeo & rectangle =
      * std::dynamic_pointer_cast<EuclideanParameter> ( parameter ) -> geo;
    // Set parameter intervals from rectangle
    p0 = getRectangleComponent ( rectangle, 0 );

    // Train one GP per dimension
    for ( int d = 0; d < phase_dim_; ++ d ) {
      // Initilize the Gaussian processes using the squared exponential
      // covariance function with additive white noise
      libgp::GaussianProcess gp(phase_dim_, "CovSum ( CovSEiso, CovNoise)");

      // initialize hyper parameters vector
      Eigen::VectorXd params(gp . covf() . get_param_dim ());
      params << 0.0, 0.0, -2.0;

      // set parameters of covariance function
      gp . covf() . set_loghyper (params);

      // Add the training data
      for( int i = 0; i < X_train . size(); ++ i) {
        gp . add_pattern (X_train [i] . data (), Y_train [i][d]);
      }

      GPs . push_back (gp);
    }

    // Get the vertices of the unit cube
    verts_cube_ = getVertsUnitCube ( phase_dim_ );
  }

  std::shared_ptr<Geo>
  operator () ( std::shared_ptr<Geo> geo ) const {
    return std::shared_ptr<Geo> ( new RectGeo (
        operator () ( * std::dynamic_pointer_cast<RectGeo> ( geo ) ) ) );
  }
private:
  uint64_t phase_dim_;

  // Gaussian Processes. Declare it mutable so it
  // can potentially be modified inside const functions
  mutable std::vector<libgp::GaussianProcess> GPs;

  // Used to get corners of a rectangle
  std::vector<std::vector<int>> verts_cube_;

  interval getRectangleComponent ( const RectGeo & rectangle, int d ) const {
    return interval (rectangle . lower_bounds [ d ], rectangle . upper_bounds [ d ]);
  }

  std::vector<std::vector<int>> getVertsUnitCube ( int dim ) const {
    // Vertices of unit cube in dimension dim
    std::vector<std::vector<int>> verts_cube;

    if ( dim == 1 ) {
      // Set vertices fo case dim == 1
      verts_cube = { { 0 }, { 1 } };
    }
    else if ( dim == 2 ) {
      // Set vertices for case dim == 2
      verts_cube = { { 0, 0 }, { 1, 0 }, { 0, 1 }, { 1, 1 } };
    }
    else if ( dim > 2 ) { // Avoid setting anything if dim < 1
      // Set vertices for the case dim == 3
      verts_cube = { { 0, 0, 0 }, { 1, 0, 0 }, { 0, 1, 0 }, { 1, 1, 0 },
                     { 0, 0, 1 }, { 1, 0, 1 }, { 0, 1, 1 }, { 1, 1, 1 } };
    }

    // If dim > 3 use the case dim == 3 set above and
    // compute cartesian product with {0, 1}^(dim-3)
    if ( dim > 3 ) {
      for ( int d = 4; d <= dim; ++ d ) {
        // Vector of previous verts times {1}
        std::vector<std::vector<int>> verts_top;

        int num_verts = verts_cube . size();
        for ( int i = 0; i < num_verts; ++ i ) {
          std::vector<int> v = verts_cube [i]; // Get vertex i
          // Make product with {0} in place
          verts_cube [i] . push_back (0);
          // Put product with {1} in verts_top
          v . push_back (1);
          verts_top . push_back (v);
        }

        // Insert verts_top at the back of verts_cube
        verts_cube . insert( verts_cube . end(), verts_top . begin(), verts_top . end() );
      }
    }

    return verts_cube;
  }
};

#endif

#ifndef CMDB_MODEL_H
#define CMDB_MODEL_H

#include "ModelMap.h"
#include "ModelMapF.h"
#include "Map.h"
#include "EuclideanParameterSpace.h"
#include "TreeGrid.h"
#include "PointerGrid.h"
#include "SuccinctGrid.h"
#include "UniformGrid.h"
#include "MorseGraph.h"
#include "Configuration.h"

#include <cstdlib>
#include <exception>
#include <functional>

#define PHASE_GRID PointerGrid
#define PARAMETER_GRID UniformGrid

class Model {
 public:
  /// Default constructor
  Model () {};

  /// Constructor
  Model ( int phase_subdiv_min, int phase_subdiv_max,
          std::vector<double> const& phase_lower_bounds,
          std::vector<double> const& phase_upper_bounds );

  Model ( int phase_subdiv_min, int phase_subdiv_max,
          std::vector<double> const& phase_lower_bounds,
          std::vector<double> const& phase_upper_bounds,
          std::vector<bool> const& phase_periodic );

  Model ( int phase_subdiv_min, int phase_subdiv_max,
          int phase_subdiv_init, int phase_subdiv_limit,
          std::vector<double> const& phase_lower_bounds,
          std::vector<double> const& phase_upper_bounds );

  Model ( int phase_subdiv_min, int phase_subdiv_max,
          int phase_subdiv_init, int phase_subdiv_limit,
          std::vector<double> const& phase_lower_bounds,
          std::vector<double> const& phase_upper_bounds,
          std::vector<bool> const& phase_periodic );

  Model ( int phase_subdiv,
          std::vector<double> const& phase_lower_bounds,
          std::vector<double> const& phase_upper_bounds,
          std::function<std::vector<double>(std::vector<double>)> const& F );

  Model ( int phase_subdiv,
          std::vector<double> const& phase_lower_bounds,
          std::vector<double> const& phase_upper_bounds,
          std::vector<bool> const& phase_periodic,
          std::function<std::vector<double>(std::vector<double>)> const& F );

  Model ( int phase_subdiv_min, int phase_subdiv_max,
          std::vector<double> const& phase_lower_bounds,
          std::vector<double> const& phase_upper_bounds,
          std::function<std::vector<double>(std::vector<double>)> const& F );

  Model ( int phase_subdiv_min, int phase_subdiv_max,
          std::vector<double> const& phase_lower_bounds,
          std::vector<double> const& phase_upper_bounds,
          std::vector<bool> const& phase_periodic,
          std::function<std::vector<double>(std::vector<double>)> const& F );

  Model ( int phase_subdiv_min, int phase_subdiv_max,
          int phase_subdiv_init, int phase_subdiv_limit,
          std::vector<double> const& phase_lower_bounds,
          std::vector<double> const& phase_upper_bounds,
          std::function<std::vector<double>(std::vector<double>)> const& F );

  Model ( int phase_subdiv_min, int phase_subdiv_max,
          int phase_subdiv_init, int phase_subdiv_limit,
          std::vector<double> const& phase_lower_bounds,
          std::vector<double> const& phase_upper_bounds,
          std::vector<bool> const& phase_periodic,
          std::function<std::vector<double>(std::vector<double>)> const& F );

  /// initialize
  ///   Given command line arguments, load necessary files 
  ///   required for initializtion.
  void initialize ( int param_dim, int phase_dim,
                    int phase_subdiv_min, int phase_subdiv_max,
                    int phase_subdiv_init, int phase_subdiv_limit,
                    std::vector<double> const& param_lower_bounds,
                    std::vector<double> const& param_upper_bounds,
                    std::vector<double> const& phase_lower_bounds,
                    std::vector<double> const& phase_upper_bounds,
                    std::vector<bool> const& phase_periodic );

  void initialize ( int param_dim, int phase_dim,
                    int phase_subdiv_min, int phase_subdiv_max,
                    int phase_subdiv_init, int phase_subdiv_limit,
                    std::vector<double> const& param_lower_bounds,
                    std::vector<double> const& param_upper_bounds,
                    std::vector<double> const& phase_lower_bounds,
                    std::vector<double> const& phase_upper_bounds,
                    std::vector<bool> const& phase_periodic,
                    std::function<std::vector<double>(std::vector<double>)> const& F );

  /// parameterSpace
  ///   return a shared ptr to the parameter space
  std::shared_ptr < ParameterSpace > parameterSpace ( void ) const;

  /// phaseSpace
  ///   return a shared ptr to the phase space
  std::shared_ptr < Grid > phaseSpace ( void ) const;

  int param_dim ( void ) const;

  int phase_dim ( void ) const;

  int phase_subdiv_min ( void ) const;

  int phase_subdiv_max ( void ) const;

  int phase_subdiv_init ( void ) const;

  int phase_subdiv_limit ( void ) const;

  std::vector<double> const&
  param_lower_bounds ( void ) const;

  std::vector<double> const&
  param_upper_bounds ( void ) const;

  std::vector<double> const&
  phase_lower_bounds ( void ) const;

  std::vector<double> const&
  phase_upper_bounds ( void ) const;

  std::vector<bool> const&
  phase_periodic ( void ) const;

  /// setmap
  ///   set shared ptr to a map function object corresponding
  ///   to parameter p
  void setmap ( std::shared_ptr<Parameter> p );

  /// set_batch_map
  ///   Attach a batched evaluator for the map F this model was constructed
  ///   with. The batched evaluator must agree with F on every rectangle; it
  ///   exists so many rectangles can be evaluated per call. Requires a model
  ///   constructed with a map function F.
  void set_batch_map ( ModelMapF::BatchFunction const& F_batch );

  /// map
  ///   return a shared ptr to a map function object corresponding to 
  ///   parameter p
  ///   void version: returns map corresponding to command_line parameter
  ///   f version: returns map using the provided map f
  // std::shared_ptr < const Map > map ( std::function<std::vector<double>(std::vector<double>)> const& F ) const;
  // std::shared_ptr < const Map > map ( std::shared_ptr<Parameter> p ) const;
  // std::shared_ptr < const Map > map ( void ) const;

  /// map
  ///   return a shared ptr to a map function object
  std::shared_ptr < const Map > map ( void ) const;

  /// annotate
  ///   Given a MorseGraph, provide annotations.
  void annotate ( MorseGraph * mg_in ) const;
private:
  Configuration config_;
  std::shared_ptr<EuclideanParameterSpace> parameter_space_;
  std::shared_ptr<EuclideanParameter> command_line_parameter_;
  std::shared_ptr<Map> map_;
public:
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive& ar, const unsigned int version) {
    ar & config_;
  }
};

inline
Model::Model ( int phase_subdiv_min, int phase_subdiv_max,
               std::vector<double> const& phase_lower_bounds,
               std::vector<double> const& phase_upper_bounds ) {
  std::vector<double> params {0.0};
  std::vector<double> param_lower_bounds = params;
  std::vector<double> param_upper_bounds = params;
  int param_dim = params . size();
  int phase_dim = phase_lower_bounds . size();
  std::vector<bool> phase_periodic ( phase_dim, false );
  int phase_subdiv_init = 0;
  int phase_subdiv_limit = 10000;

  initialize ( param_dim, phase_dim,
               phase_subdiv_min, phase_subdiv_max,
               phase_subdiv_init, phase_subdiv_limit,
               param_lower_bounds, param_upper_bounds,
               phase_lower_bounds, phase_upper_bounds,
               phase_periodic );
}

inline
Model::Model ( int phase_subdiv_min, int phase_subdiv_max,
               std::vector<double> const& phase_lower_bounds,
               std::vector<double> const& phase_upper_bounds,
               std::vector<bool> const& phase_periodic ) {
  std::vector<double> params {0.0};
  std::vector<double> param_lower_bounds = params;
  std::vector<double> param_upper_bounds = params;
  int param_dim = params . size();
  int phase_dim = phase_lower_bounds . size();
  int phase_subdiv_init = 0;
  int phase_subdiv_limit = 10000;

  initialize ( param_dim, phase_dim,
               phase_subdiv_min, phase_subdiv_max,
               phase_subdiv_init, phase_subdiv_limit,
               param_lower_bounds, param_upper_bounds,
               phase_lower_bounds, phase_upper_bounds,
               phase_periodic );
}

inline
Model::Model ( int phase_subdiv_min, int phase_subdiv_max,
               int phase_subdiv_init, int phase_subdiv_limit,
               std::vector<double> const& phase_lower_bounds,
               std::vector<double> const& phase_upper_bounds ) {
  std::vector<double> params {0.0};
  std::vector<double> param_lower_bounds = params;
  std::vector<double> param_upper_bounds = params;
  int param_dim = params . size();
  int phase_dim = phase_lower_bounds . size();
  std::vector<bool> phase_periodic ( phase_dim, false );

  initialize ( param_dim, phase_dim,
               phase_subdiv_min, phase_subdiv_max,
               phase_subdiv_init, phase_subdiv_limit,
               param_lower_bounds, param_upper_bounds,
               phase_lower_bounds, phase_upper_bounds,
               phase_periodic );
}

inline
Model::Model ( int phase_subdiv_min, int phase_subdiv_max,
               int phase_subdiv_init, int phase_subdiv_limit,
               std::vector<double> const& phase_lower_bounds,
               std::vector<double> const& phase_upper_bounds,
               std::vector<bool> const& phase_periodic ) {
  std::vector<double> params {0.0};
  std::vector<double> param_lower_bounds = params;
  std::vector<double> param_upper_bounds = params;
  int param_dim = params . size();
  int phase_dim = phase_lower_bounds . size();

  initialize ( param_dim, phase_dim,
               phase_subdiv_min, phase_subdiv_max,
               phase_subdiv_init, phase_subdiv_limit,
               param_lower_bounds, param_upper_bounds,
               phase_lower_bounds, phase_upper_bounds,
               phase_periodic );
}

inline
Model::Model ( int phase_subdiv,
               std::vector<double> const& phase_lower_bounds,
               std::vector<double> const& phase_upper_bounds,
               std::function<std::vector<double>(std::vector<double>)> const& F ) {
  std::vector<double> params {0.0};
  std::vector<double> param_lower_bounds = params;
  std::vector<double> param_upper_bounds = params;
  int param_dim = params . size();
  int phase_dim = phase_lower_bounds . size();
  std::vector<bool> phase_periodic ( phase_dim, false );
  int phase_subdiv_init = phase_subdiv;
  int phase_subdiv_min = phase_subdiv;
  int phase_subdiv_max = phase_subdiv;
  int phase_subdiv_limit = 10000;

  initialize ( param_dim, phase_dim,
               phase_subdiv_min, phase_subdiv_max,
               phase_subdiv_init, phase_subdiv_limit,
               param_lower_bounds, param_upper_bounds,
               phase_lower_bounds, phase_upper_bounds,
               phase_periodic, F );
}

inline
Model::Model ( int phase_subdiv,
               std::vector<double> const& phase_lower_bounds,
               std::vector<double> const& phase_upper_bounds,
               std::vector<bool> const& phase_periodic,
               std::function<std::vector<double>(std::vector<double>)> const& F ) {
  std::vector<double> params {0.0};
  std::vector<double> param_lower_bounds = params;
  std::vector<double> param_upper_bounds = params;
  int param_dim = params . size();
  int phase_dim = phase_lower_bounds . size();
  int phase_subdiv_init = phase_subdiv;
  int phase_subdiv_min = phase_subdiv;
  int phase_subdiv_max = phase_subdiv;
  int phase_subdiv_limit = 10000;

  initialize ( param_dim, phase_dim,
               phase_subdiv_min, phase_subdiv_max,
               phase_subdiv_init, phase_subdiv_limit,
               param_lower_bounds, param_upper_bounds,
               phase_lower_bounds, phase_upper_bounds,
               phase_periodic, F );
}

inline
Model::Model ( int phase_subdiv_min, int phase_subdiv_max,
               std::vector<double> const& phase_lower_bounds,
               std::vector<double> const& phase_upper_bounds,
               std::function<std::vector<double>(std::vector<double>)> const& F ) {
  std::vector<double> params {0.0};
  std::vector<double> param_lower_bounds = params;
  std::vector<double> param_upper_bounds = params;
  int param_dim = params . size();
  int phase_dim = phase_lower_bounds . size();
  std::vector<bool> phase_periodic ( phase_dim, false );
  int phase_subdiv_init = 0;
  int phase_subdiv_limit = 10000;

  initialize ( param_dim, phase_dim,
               phase_subdiv_min, phase_subdiv_max,
               phase_subdiv_init, phase_subdiv_limit,
               param_lower_bounds, param_upper_bounds,
               phase_lower_bounds, phase_upper_bounds,
               phase_periodic, F );
}

inline
Model::Model ( int phase_subdiv_min, int phase_subdiv_max,
               std::vector<double> const& phase_lower_bounds,
               std::vector<double> const& phase_upper_bounds,
               std::vector<bool> const& phase_periodic,
               std::function<std::vector<double>(std::vector<double>)> const& F ) {
  std::vector<double> params {0.0};
  std::vector<double> param_lower_bounds = params;
  std::vector<double> param_upper_bounds = params;
  int param_dim = params . size();
  int phase_dim = phase_lower_bounds . size();
  int phase_subdiv_init = 0;
  int phase_subdiv_limit = 10000;

  initialize ( param_dim, phase_dim,
               phase_subdiv_min, phase_subdiv_max,
               phase_subdiv_init, phase_subdiv_limit,
               param_lower_bounds, param_upper_bounds,
               phase_lower_bounds, phase_upper_bounds,
               phase_periodic, F );
}

inline
Model::Model ( int phase_subdiv_min, int phase_subdiv_max,
               int phase_subdiv_init, int phase_subdiv_limit,
               std::vector<double> const& phase_lower_bounds,
               std::vector<double> const& phase_upper_bounds,
               std::function<std::vector<double>(std::vector<double>)> const& F ) {
  std::vector<double> params {0.0};
  std::vector<double> param_lower_bounds = params;
  std::vector<double> param_upper_bounds = params;
  int param_dim = params . size();
  int phase_dim = phase_lower_bounds . size();
  std::vector<bool> phase_periodic ( phase_dim, false );

  initialize ( param_dim, phase_dim,
               phase_subdiv_min, phase_subdiv_max,
               phase_subdiv_init, phase_subdiv_limit,
               param_lower_bounds, param_upper_bounds,
               phase_lower_bounds, phase_upper_bounds,
               phase_periodic, F );
}

inline
Model::Model ( int phase_subdiv_min, int phase_subdiv_max,
               int phase_subdiv_init, int phase_subdiv_limit,
               std::vector<double> const& phase_lower_bounds,
               std::vector<double> const& phase_upper_bounds,
               std::vector<bool> const& phase_periodic,
               std::function<std::vector<double>(std::vector<double>)> const& F ) {
  std::vector<double> params {0.0};
  std::vector<double> param_lower_bounds = params;
  std::vector<double> param_upper_bounds = params;
  int param_dim = params . size();
  int phase_dim = phase_lower_bounds . size();

  initialize ( param_dim, phase_dim,
               phase_subdiv_min, phase_subdiv_max,
               phase_subdiv_init, phase_subdiv_limit,
               param_lower_bounds, param_upper_bounds,
               phase_lower_bounds, phase_upper_bounds,
               phase_periodic, F );
}

inline void
Model::initialize ( int param_dim, int phase_dim,
                    int phase_subdiv_min, int phase_subdiv_max,
                    int phase_subdiv_init, int phase_subdiv_limit,
                    std::vector<double> const& param_lower_bounds,
                    std::vector<double> const& param_upper_bounds,
                    std::vector<double> const& phase_lower_bounds,
                    std::vector<double> const& phase_upper_bounds,
                    std::vector<bool> const& phase_periodic ) {

  config_ . setConfiguration ( param_dim, phase_dim,
                               phase_subdiv_min, phase_subdiv_max,
                               phase_subdiv_init, phase_subdiv_limit,
                               phase_lower_bounds, phase_upper_bounds,
                               phase_periodic );
  std::shared_ptr<Grid> parameter_grid ( new PARAMETER_GRID );
  parameter_space_ . reset ( new EuclideanParameterSpace );
  parameter_space_ -> initialize ( config_, parameter_grid );
  // Load "command line" parameter
  int dim = parameter_space_ -> dimension ();
  command_line_parameter_ . reset ( new EuclideanParameter ( new RectGeo ( dim ) ) );
  RectGeo & geo = * command_line_parameter_ -> geo;
  for ( int d = 0; d < dim; ++ d ) {
    geo . lower_bounds [ d ] = param_lower_bounds [d];
    geo . upper_bounds [ d ] = param_upper_bounds [d];
  }

  // Set pointer to map
  map_ . reset ( new ModelMap ( command_line_parameter_ ) );
}

inline void
Model::initialize ( int param_dim, int phase_dim,
                    int phase_subdiv_min, int phase_subdiv_max,
                    int phase_subdiv_init, int phase_subdiv_limit,
                    std::vector<double> const& param_lower_bounds,
                    std::vector<double> const& param_upper_bounds,
                    std::vector<double> const& phase_lower_bounds,
                    std::vector<double> const& phase_upper_bounds,
                    std::vector<bool> const& phase_periodic,
                    std::function<std::vector<double>(std::vector<double>)> const& F ) {

  config_ . setConfiguration ( param_dim, phase_dim,
                               phase_subdiv_min, phase_subdiv_max,
                               phase_subdiv_init, phase_subdiv_limit,
                               phase_lower_bounds, phase_upper_bounds,
                               phase_periodic );
  std::shared_ptr<Grid> parameter_grid ( new PARAMETER_GRID );
  parameter_space_ . reset ( new EuclideanParameterSpace );
  parameter_space_ -> initialize ( config_, parameter_grid );
  // Load "command line" parameter
  int dim = parameter_space_ -> dimension ();
  command_line_parameter_ . reset ( new EuclideanParameter ( new RectGeo ( dim ) ) );
  RectGeo & geo = * command_line_parameter_ -> geo;
  for ( int d = 0; d < dim; ++ d ) {
    geo . lower_bounds [ d ] = param_lower_bounds [d];
    geo . upper_bounds [ d ] = param_upper_bounds [d];
  }

  // Set pointer to map
  map_ . reset ( new ModelMapF ( command_line_parameter_, F ) );
}

inline int
Model::param_dim ( void ) const {
  return config_ . PARAM_DIM;
}

inline int
Model::phase_dim ( void ) const {
  return config_ . PHASE_DIM;
}

inline int
Model::phase_subdiv_min ( void ) const {
  return config_ . PHASE_SUBDIV_MIN;
}

inline int
Model::phase_subdiv_max ( void ) const {
  return config_ . PHASE_SUBDIV_MAX;
}

inline int
Model::phase_subdiv_init ( void ) const {
  return config_ . PHASE_SUBDIV_INIT;
}

inline int
Model::phase_subdiv_limit ( void ) const {
  return config_ . PHASE_SUBDIV_LIMIT;
}

inline std::vector<double> const&
Model::param_lower_bounds ( void ) const {
  return config_ . PARAM_BOUNDS . lower_bounds;
}

inline std::vector<double> const&
Model::param_upper_bounds ( void ) const {
  return config_ . PARAM_BOUNDS . upper_bounds;
}

inline std::vector<double> const&
Model::phase_lower_bounds ( void ) const {
  return config_ . PHASE_BOUNDS . lower_bounds;
}

inline std::vector<double> const&
Model::phase_upper_bounds ( void ) const {
  return config_ . PHASE_BOUNDS . upper_bounds;
}

inline std::vector<bool> const&
Model::phase_periodic ( void ) const {
  return config_ . PHASE_PERIODIC;
}

inline std::shared_ptr < ParameterSpace > 
Model::parameterSpace ( void ) const {
  return std::dynamic_pointer_cast<ParameterSpace> ( parameter_space_ );
}

inline std::shared_ptr < Grid > 
Model::phaseSpace ( void ) const {
  std::shared_ptr < TreeGrid > space ( new PHASE_GRID );
  space -> initialize ( config_.PHASE_BOUNDS, config_.PHASE_PERIODIC );
  return std::dynamic_pointer_cast<Grid> ( space );
}

inline void
Model::setmap ( std::shared_ptr<Parameter> p ) {
  if ( not p ) p = command_line_parameter_;
  if ( not p ) {
    throw std::logic_error ( "No parameter for map specified. "
                             "Check Model.h and command line parameters.\n");
  }
  map_ . reset ( new ModelMap ( p ) );
}

inline void
Model::set_batch_map ( ModelMapF::BatchFunction const& F_batch ) {
  std::shared_ptr<ModelMapF> map_f = std::dynamic_pointer_cast<ModelMapF> ( map_ );
  if ( not map_f ) {
    throw std::logic_error ( "Model::set_batch_map requires a model "
                             "constructed with a map function F.\n" );
  }
  map_f -> set_batch_map ( F_batch );
}

// inline std::shared_ptr < const Map >
// Model::map ( std::function<std::vector<double>(std::vector<double>)> const& F ) const { 
//   return std::shared_ptr < Map > ( new ModelMapF ( command_line_parameter_, F ) );
// }

// inline std::shared_ptr < const Map >
// Model::map ( std::shared_ptr<Parameter> p ) const {
//   if ( not p ) p = command_line_parameter_;
//   if ( not p ) {
//     throw std::logic_error ( "No parameter for map specified. "
//                              "Check Model.h and command line parameters.\n");
//   }
//   return std::shared_ptr < Map > ( new ModelMap ( p ) );
// }

// inline std::shared_ptr < const Map >
// Model::map ( void ) const {
//   return map ( command_line_parameter_ );
// }

inline std::shared_ptr < const Map >
Model::map ( void ) const {
   return map_;
}

inline void 
Model::annotate( MorseGraph * mg_in ) const {
  /*
  // example:
  MorseGraph & mg = *mg_in;
  mg . annotation () . insert ( "annotation_A" );
  mg . annotation () . insert ( "annotation_B" );
  for ( int v = 0; v < mg.NumVertices(); ++ v ) {
    mg . annotation ( v ) . insert ( std::string ( "annotation_C" ) );
  }
  */
}

/// Python Bindings

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
namespace py = pybind11;

inline void
ModelBinding(py::module &m) {
  py::class_<Model, std::shared_ptr<Model>>(m, "Model")
    .def(py::init<>())
    .def(py::init<int, int, std::vector<double> const&, std::vector<double> const&>())
    .def(py::init<int, int, std::vector<double> const&, std::vector<double> const&,
                  std::vector<bool> const&>())
    .def(py::init<int, int, int, int, std::vector<double> const&, std::vector<double> const&>())
    .def(py::init<int, int, int, int, std::vector<double> const&, std::vector<double> const&,
                  std::vector<bool> const&>())
    .def(py::init<int, std::vector<double> const&, std::vector<double> const&,
                  std::function<std::vector<double>(std::vector<double>)> const&>())
    .def(py::init<int, std::vector<double> const&, std::vector<double> const&,
                  std::vector<bool> const&,
                  std::function<std::vector<double>(std::vector<double>)> const&>())
    .def(py::init<int, int, std::vector<double> const&, std::vector<double> const&,
                  std::function<std::vector<double>(std::vector<double>)> const&>())
    .def(py::init<int, int, std::vector<double> const&, std::vector<double> const&,
                  std::vector<bool> const&,
                  std::function<std::vector<double>(std::vector<double>)> const&>())
    .def(py::init<int, int, int, int, std::vector<double> const&, std::vector<double> const&,
                  std::function<std::vector<double>(std::vector<double>)> const&>())
    .def(py::init<int, int, int, int, std::vector<double> const&, std::vector<double> const&,
                  std::vector<bool> const&,
                  std::function<std::vector<double>(std::vector<double>)> const&>())
    .def("parameterSpace", &Model::parameterSpace)
    .def("phaseSpace", &Model::phaseSpace)
    .def("setmap", &Model::setmap)
    .def("set_batch_map", [](Model & model, py::function F_batch_py) {
      // Wrap the Python batch map as a flat-buffer evaluator. The rectangle
      // buffer is exposed to Python as a zero-copy, read-only NumPy array of
      // shape (count, 2*dim); the view is valid only during the call.
      model . set_batch_map (
        [F_batch_py](const std::vector<double> & rects, uint64_t count,
                     uint64_t dim, std::vector<double> & images) {
          py::gil_scoped_acquire gil;
          py::array_t<double> rects_array (
            { (py::ssize_t) count, (py::ssize_t) (2 * dim) },
            { (py::ssize_t) (2 * dim * sizeof(double)), (py::ssize_t) sizeof(double) },
            rects . data (),
            py::capsule ( rects . data (), [](void *){} ) );
          rects_array . attr ( "setflags" ) ( py::arg ( "write" ) = false );
          py::object result = F_batch_py ( rects_array );
          py::array_t<double, py::array::c_style | py::array::forcecast> images_array =
            py::array_t<double, py::array::c_style | py::array::forcecast>::ensure ( result );
          if ( not images_array ) {
            throw std::runtime_error (
              "The batch map must return an array-like of doubles" );
          }
          if ( (uint64_t) images_array . size () != count * 2 * dim ) {
            throw std::runtime_error (
              "The batch map returned the wrong number of values: "
              "expected an array of shape (count, 2*dim)" );
          }
          const double * data = images_array . data ();
          images . assign ( data, data + count * 2 * dim );
        } );
    })
    .def("param_dim", &Model::param_dim)
    .def("phase_dim", &Model::phase_dim)
    .def("phase_subdiv_min", &Model::phase_subdiv_min)
    .def("phase_subdiv_max", &Model::phase_subdiv_max)
    .def("phase_subdiv_init", &Model::phase_subdiv_init)
    .def("phase_subdiv_limit", &Model::phase_subdiv_limit)
    .def("param_lower_bounds", &Model::param_lower_bounds)
    .def("param_upper_bounds", &Model::param_upper_bounds)
    .def("phase_lower_bounds", &Model::phase_lower_bounds)
    .def("phase_upper_bounds", &Model::phase_upper_bounds)
    .def("phase_periodic", &Model::phase_periodic);
}

#endif

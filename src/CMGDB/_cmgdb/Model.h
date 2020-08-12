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
  /// initialize
  ///   Given command line arguments, load necessary files 
  ///   required for initializtion.
  void initialize ( int param_dim, int phase_dim,
                    int phase_subdiv_min, int phase_subdiv_max,
                    int phase_subdiv_init, int phase_subdiv_limit,
                    std::vector<double> const& param_lower_bounds,
                    std::vector<double> const& param_upper_bounds,
                    std::vector<double> const& phase_lower_bounds,
                    std::vector<double> const& phase_upper_bounds );

  /// parameterSpace
  ///   return a shared ptr to the parameter space
  std::shared_ptr < ParameterSpace > parameterSpace ( void ) const;

  /// phaseSpace
  ///   return a shared ptr to the phase space
  std::shared_ptr < Grid > phaseSpace ( void ) const;

  /// map
  ///   return a shared ptr to a map function object corresponding to 
  ///   parameter p
  ///   void version: returns map corresponding to command_line parameter
  ///   f version: returns map using the provided map f
  std::shared_ptr < const Map > map ( std::function<std::vector<double>(std::vector<double>)> const& F ) const;
  std::shared_ptr < const Map > map ( std::shared_ptr<Parameter> p ) const;
  std::shared_ptr < const Map > map ( void ) const;

  /// annotate
  ///   Given a MorseGraph, provide annotations.
  void annotate ( MorseGraph * mg_in ) const;
private:
  Configuration config_;
  std::shared_ptr<EuclideanParameterSpace> parameter_space_;
  std::shared_ptr<EuclideanParameter> command_line_parameter_;
public:
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive& ar, const unsigned int version) {
    ar & config_;
  }
};

inline void 
Model::initialize ( int param_dim, int phase_dim,
                    int phase_subdiv_min, int phase_subdiv_max,
                    int phase_subdiv_init, int phase_subdiv_limit,
                    std::vector<double> const& param_lower_bounds,
                    std::vector<double> const& param_upper_bounds,
                    std::vector<double> const& phase_lower_bounds,
                    std::vector<double> const& phase_upper_bounds ) {

  config_ . setConfiguration ( param_dim, phase_dim,
                               phase_subdiv_min, phase_subdiv_max,
                               phase_subdiv_init, phase_subdiv_limit,
                               phase_lower_bounds, phase_upper_bounds );
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

inline std::shared_ptr < const Map >
Model::map ( std::function<std::vector<double>(std::vector<double>)> const& F ) const { 
  return std::shared_ptr < Map > ( new ModelMapF ( command_line_parameter_, F ) );
}

inline std::shared_ptr < const Map >
Model::map ( std::shared_ptr<Parameter> p ) const {
  if ( not p ) p = command_line_parameter_;
  if ( not p ) {
    throw std::logic_error ( "No parameter for map specified. "
                             "Check Model.h and command line parameters.\n");
  }
  return std::shared_ptr < Map > ( new ModelMap ( p ) );
}

inline std::shared_ptr < const Map >
Model::map ( void ) const {
  return map ( command_line_parameter_ );
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

#endif

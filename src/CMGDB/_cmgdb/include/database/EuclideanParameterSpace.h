#ifndef CMDB_EUCLIDEANPARAMETERSPACE
#define CMDB_EUCLIDEANPARAMETERSPACE

#include "Grid.h"
#include "TreeGrid.h"
#include "UniformGrid.h"
#include "EdgeGrid.h"

#include "ParameterSpace.h"
#include "RectGeo.h"

#include <memory>
#include "boost/serialization/serialization.hpp"
#include "boost/serialization/vector.hpp"
#include "boost/serialization/unordered_map.hpp"
#include "boost/serialization/shared_ptr.hpp"
#include "boost/serialization/export.hpp"
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

class EuclideanParameter : public Parameter {
public:
	std::shared_ptr<RectGeo> geo;
  EuclideanParameter ( void ) {}
	EuclideanParameter ( std::shared_ptr<RectGeo> geo ) : geo(geo) {}
  EuclideanParameter ( RectGeo * geo_ptr ) {
    geo . reset ( geo_ptr );
  }
  EuclideanParameter ( RectGeo geo_obj ) {
    geo . reset ( new RectGeo ( geo_obj ) );
  }

private:
	// Serialization
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version) {
    ar & boost::serialization::base_object<Parameter>(*this);
  	ar & geo;
  }
  /// derivation interface
  virtual void print ( std::ostream & outstream ) const {
  	outstream << * geo;
  }
};

BOOST_CLASS_EXPORT_KEY(EuclideanParameter);

class EuclideanParameterSpace : public ParameterSpace {
public:

  virtual ~EuclideanParameterSpace ( void ) {}

	/// initialize
	///    Create the ParameterSpace given the configuration specified
  virtual void initialize ( const Configuration & config );

	void initialize ( const Configuration & config,
                    std::shared_ptr<Grid> parameter_grid  );

	/// adjacencies
	///    Return a vector of adjacent vertices.
	virtual std::vector<uint64_t> adjacencies ( uint64_t v ) const;
	
	/// size
	///    Return the number of vertices
	virtual uint64_t size ( void ) const;

	/// parameter
	///    Return the parameter object associated with a vertex
	virtual std::shared_ptr<Parameter> parameter ( uint64_t v ) const;
	
  /// search
  ///    Given a parameter, find the vertex associated with it
  ///    (This can be used to find a parameter which might contain the other)
  virtual uint64_t search ( std::shared_ptr<Parameter> parameter ) const;
  
	/// patch
	///    Return a "ParameterPatch" object
	///    A sequence of calls to this function will return a sequence
	///    of patches that cover the entire parameter space.
	///    To do this the function must use mutable data (and so is not thread-safe)
	///    even though it is declared const.
	///    One the sequence is completed an empty patch is returned, and then the
	///    sequence will restart.
	///    The default implementation returns patches that consist of two vertices
	///    and the edge between them.
	virtual std::shared_ptr<ParameterPatch> patch ( void ) const;

  /// dimension
  ///    Return dimension of parameter space
  int dimension ( void ) const;

  /// grid
  ///    Return underlying grid object
  std::shared_ptr<const Grid> 
  grid ( void ) const;

private:
	std::shared_ptr<Grid> parameter_grid_;
  RectGeo bounds_;
  std::vector<bool> periodic_;
  int dimension_;

// EuclideanPatch construction variables
	mutable int patch_width_; 
  mutable std::vector<int> patches_across_;  
  mutable std::vector<int> coordinates_;
  mutable bool finished_;

  // Serialization
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version) {
    ar & boost::serialization::base_object<ParameterSpace>(*this);
    ar & parameter_grid_;
    ar & bounds_;
    ar & periodic_;
    ar & dimension_;
    ar & patch_width_;
    ar & patches_across_;
    ar & coordinates_;
    ar & finished_;
  }
};

BOOST_CLASS_EXPORT_KEY(EuclideanParameterSpace);

inline void
EuclideanParameterSpace::initialize ( const Configuration & config ) {
  std::shared_ptr<Grid> parameter_grid ( new UniformGrid );
  initialize ( config, parameter_grid );
}

inline void
EuclideanParameterSpace::initialize ( const Configuration & config, 
                                      std::shared_ptr<Grid> parameter_grid ) {
  parameter_grid_ = parameter_grid;

  // Initialization for TreeGrid
  if ( std::dynamic_pointer_cast < TreeGrid > ( parameter_grid_ ) ) {
    std::shared_ptr<TreeGrid> grid = 
      std::dynamic_pointer_cast < TreeGrid > ( parameter_grid_ );
    grid -> initialize ( config.PARAM_BOUNDS, 
                         config.PARAM_PERIODIC );  
    for (int i = 0; i < config.PARAM_SUBDIV_DEPTH[0]; ++i) {
      for ( int d = 0; d < config.PARAM_DIM; ++ d ) {
        parameter_grid_ -> subdivide (); // subdivide every top cell
      }
    }
    bounds_ = grid -> bounds ();
  }

  // Initialization for UniformGrid
  if ( std::dynamic_pointer_cast < UniformGrid > ( parameter_grid_ ) ) {
    std::shared_ptr<UniformGrid> grid = 
      std::dynamic_pointer_cast < UniformGrid > ( parameter_grid_ );
    grid -> initialize ( config.PARAM_BOUNDS, 
                         config.PARAM_SUBDIV_SIZES,
                         config.PARAM_PERIODIC );
    bounds_ = grid -> bounds ();
  }
  
  // Initialization for EdgeGrid
  if ( std::dynamic_pointer_cast < EdgeGrid > ( parameter_grid_ ) ) {
    std::shared_ptr<EdgeGrid> grid = 
      std::dynamic_pointer_cast < EdgeGrid > ( parameter_grid_ );
    grid -> initialize ( config.PARAM_BOUNDS, 
                         config.PARAM_SUBDIV_SIZES,
                         config.PARAM_PERIODIC );
		bounds_ = grid -> bounds ();
  }

  // Initialize Patch creation
    // In order to accommodate periodicity, 
    //    we let the patches overhang from the outer bounds slightly.
  patch_width_ = 4; // try to use (patch_width +- 1)^d boxes per patch 
  
    // Create patches_across:
    //    The distance between patch centers in each dimension.
    // EXAMPLE: num_across = 64, patch_width = 4 ---> patches_across = 9 

  patches_across_ . resize ( config.PARAM_DIM );
  for ( int d = 0; d < config.PARAM_DIM; ++ d ) {
    patches_across_[d] = 1 + config.PARAM_SUBDIV_SIZES[d] / patch_width_; 
  }

  // Dimension
  dimension_ = config.PARAM_DIM;
  // Periodicity
  periodic_ = config.PARAM_PERIODIC;

  // Initialize patch coordinates
  coordinates_ . resize ( config.PARAM_DIM, 0);
  finished_ = false;

}

inline std::vector<uint64_t> 
EuclideanParameterSpace::adjacencies ( uint64_t v ) const {
	std::vector<uint64_t> neighbors;
	RectGeo geo = * std::dynamic_pointer_cast<EuclideanParameter> ( parameter ( v ) ) -> geo;
  for ( int d = 0; d < dimension_; ++ d ) {
  	double tol = (bounds_.upper_bounds[d]-bounds_.lower_bounds[d])
  	  /(double)(1000000000.0);
  	geo . lower_bounds [ d ] -= tol;
  	geo . upper_bounds [ d ] += tol;
  }
  neighbors = parameter_grid_ -> cover ( geo );
  return neighbors;
}
	
inline uint64_t 
EuclideanParameterSpace::size ( void ) const {
	return parameter_grid_ -> size ();
}

inline std::shared_ptr<Parameter> 
EuclideanParameterSpace::parameter ( uint64_t v ) const {
	std::shared_ptr<RectGeo> geo = std::dynamic_pointer_cast<RectGeo> 
      ( parameter_grid_ -> geometry ( v ) );
  return std::shared_ptr<Parameter> ( new EuclideanParameter ( geo ) );
}

inline uint64_t 
EuclideanParameterSpace::search ( std::shared_ptr<Parameter> parameter ) const {
  RectGeo geo = * std::dynamic_pointer_cast<EuclideanParameter> ( parameter ) -> geo;
  std::vector<uint64_t> vertices = parameter_grid_ -> cover ( geo );
  if ( vertices . size () != 1 ) return * end ();
  return vertices [ 0 ];
}

inline std::shared_ptr<ParameterPatch> 
EuclideanParameterSpace::patch ( void ) const {

#ifdef EDGEPATCHMETHOD
  std::shared_ptr<ParameterPatch> result;
  while ( 1 ) {
    result = ParameterSpace::patch ();
    if ( result -> vertices . empty () ) break;
    uint64_t u = result -> vertices [ 0 ];
    uint64_t v = result -> vertices [ 1 ];
    RectGeo u_geo = * std::dynamic_pointer_cast<EuclideanParameter> 
      ( result -> parameter [ u ] ) -> geo;
    RectGeo v_geo = * std::dynamic_pointer_cast<EuclideanParameter> 
      ( result -> parameter [ v ] ) -> geo;
    int codimension = 0;
    for ( int d = 0; d < dimension_; ++ d ) {
      if ( ( u_geo . lower_bounds [ d ] == v_geo . upper_bounds [ d ] )
      || ( u_geo . upper_bounds [ d ] == v_geo . lower_bounds [ d ] ) ) {
        ++ codimension;
      }
    }
    if ( codimension == 1 ) break;
  }
  return result;
#else
  //std::cout << "EuclideanParameterSpace::patch dimension_ = " << dimension_ << "\n";

  std::shared_ptr<ParameterPatch> result ( new ParameterPatch );
  if ( finished_ ) {
  	finished_ = false;
  	return result;
  }

  // Determine a rectangle based on "coordinates"
  RectGeo geo ( dimension_ );
  for ( int d = 0; d < dimension_; ++ d ) {
	 // tol included for robustness
  	double tol = (bounds_.upper_bounds[d] - bounds_.lower_bounds[d]) 
                   /(double)(1000000000.0);
    geo . lower_bounds [ d ] = 
        bounds_.lower_bounds[d]+((double)coordinates_[d])*
        (bounds_.upper_bounds[d]-bounds_.lower_bounds[d]) 
        /(double)patches_across_[d] - tol;
    geo . upper_bounds [ d ] = 
        bounds_.lower_bounds[d]+((double)(1+coordinates_[d]))*
        (bounds_.upper_bounds[d]-bounds_.lower_bounds[d])
        /(double)patches_across_[d] + tol;
      
    if ( not periodic_ [ d ] ) {
      if ( geo . lower_bounds [ d ] < bounds_ . lower_bounds [ d ] ) 
        geo . lower_bounds [ d ] = bounds_ . lower_bounds [ d ];
      if ( geo . upper_bounds [ d ] > bounds_. upper_bounds [ d ] ) 
        geo . upper_bounds [ d ] = bounds_ . upper_bounds [ d ];
    }
  }
   
  //std::cout << "EuclideanParameterSpace::patch geo = " << geo << "\n";
  // Prepare ParameterPatch

  // Find Vertices in Patch
  std::vector<uint64_t> vertices = parameter_grid_ -> cover ( geo );
  boost::unordered_set<uint64_t> vertex_set ( vertices . begin (), vertices . end () );
 
  //std::cout << "EuclideanParameterSpace::patch vertices in patch = " 
  //  << vertices . size () << "\n";

  /// Find adjacency information for cells in the patch
  BOOST_FOREACH ( uint64_t vertex, vertices ) {
    //std::cout << "EuclideanParameterSpace::patch vertex = " << vertex << "\n"; 
  	result -> vertices . push_back ( (uint64_t) vertex );
    result -> parameter [ (uint64_t) vertex ] = parameter ( vertex );
    //std::cout << "EuclideanParameterSpace::patch parameter = " << * parameter(vertex) << "\n"; 
    std::vector<uint64_t> neighbors = adjacencies ( vertex );
    //std::cout << "EuclideanParameterSpace::patch number of neighbors = " 
    //  << neighbors . size () << "\n";
    BOOST_FOREACH ( uint64_t other, neighbors ) {
      //std::cout << "EuclideanParameterSpace::patch adjacent vertex = " << other << "\n"; 
      if (( vertex_set . count ( other ) != 0 ) && vertex < other ) {
        result -> edges . push_back ( std::make_pair ( vertex, other ) );
      }
    }
  }

  // Odometer step (multidimensional loop iteration)
  finished_ = true;
  for ( int d = 0; d < dimension_; ++ d ) {
    ++ coordinates_ [ d ];
    if ( coordinates_ [ d ] == patches_across_ [ d ] ) {
      coordinates_ [ d ] = 0;
    } else {
      finished_ = false;
      break;
    }
  }

  return result;
#endif
}


inline int 
EuclideanParameterSpace::dimension ( void ) const {
  return dimension_;
}

inline std::shared_ptr<const Grid> 
EuclideanParameterSpace::grid ( void ) const {
  return parameter_grid_;
}

#endif

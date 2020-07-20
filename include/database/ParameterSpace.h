#ifndef CMDB_PARAMETERSPACE_H
#define CMDB_PARAMETERSPACE_H

#include <vector>
#include <utility>
#include "unordered_map"
#include <memory>
#include "boost/serialization/serialization.hpp"
#include "boost/serialization/vector.hpp"
#include "boost/serialization/unordered_map.hpp"
#include "boost/serialization/shared_ptr.hpp"

#include <boost/iterator/counting_iterator.hpp>

#include "Configuration.h"

/// class Parameter
///    This is an abstract base class for use with the base class ParameterSpace
class Parameter {
public:
	virtual ~Parameter ( void ) {}
  
  // Serialization
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version) {
  }

  /// std stream interface
  friend std::ostream & 
  operator << ( std::ostream & out, const Parameter & print_me ) {
        print_me . print( out );
        return out;
    }

private:
    /// derivation interface
    virtual void print ( std::ostream & ) const = 0;
};

/// class ParameterPatch
/// This represents a subgraph of the graph representing parameter space
/// The representation is simple: there is a vector of vertices
/// and a vector of edges. Vertices are uint64_t, Edges are pairs of uint64_t
class ParameterPatch {
public:
	typedef uint64_t ParameterIndex;
	std::vector<ParameterIndex> vertices;
	std::vector<std::pair<ParameterIndex, ParameterIndex> > edges;
	std::unordered_map<ParameterIndex, std::shared_ptr<Parameter> > parameter;
	/// empty
	/// determine if the ParameterPatch is empty
	/// this is useful because the ParameterSpace::patch routine
	/// returns an empty patch to indicate is has emitted enough patches
	/// to cover parameter space
	bool empty ( void ) const { 
		if ( vertices . size () == 0 ) return true; 
		return false; 
	}
private:
	// Serialization
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version) {
  	ar & vertices;
  	ar & edges;
  	ar & parameter;
  }
};

/// class ParameterSpace
///   Represents parameter space as a graph
///   where parameters are represented by vertices of a graph, indexed with contiguous integers
///   and adjacent parameters are represented by edges in the graph
class ParameterSpace {
public:
	// typedef
	typedef uint64_t ParameterIndex;
	typedef boost::counting_iterator < ParameterIndex > iterator;
	typedef iterator const_iterator;
	// Constructor/Deconstructor
	ParameterSpace ( void ) : default_patch_method_vertex_(0), default_patch_method_edge_(0)
	 {}
	virtual ~ParameterSpace ( void ) {}

	/// adjacencies
	///    Return a vector of adjacent vertices.
	virtual std::vector<ParameterIndex> adjacencies ( ParameterIndex v ) const = 0;
	
	/// size
	///    Return the number of vertices
	virtual uint64_t size ( void ) const = 0;

	/// parameter
	///    Return the parameter object associated with a vertex
	virtual std::shared_ptr<Parameter> parameter ( ParameterIndex v ) const = 0;
	
	/// search
	///    Given a parameter, find the vertex associated with it
	///    (This can be used to find a parameter which might contain the other)
	virtual uint64_t search ( std::shared_ptr<Parameter> parameter ) const = 0;
	
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

	/// begin
	///    Return "begin" iterator

	iterator begin ( void ) const;

	/// end
	///    Return "end" iterator
	iterator end ( void ) const;

private:

	// variables for default implementation of "patch" method
	mutable uint64_t default_patch_method_vertex_;
	mutable uint64_t default_patch_method_edge_;
	mutable std::vector<ParameterIndex> default_patch_method_neighbors_;

	// Serialization
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version) {
  }
  
};


/////////////////
// Definitions //
/////////////////

inline std::shared_ptr<ParameterPatch> 
ParameterSpace::patch ( void ) const {
  //std::cout << "ParameterSpace::patch.\n"; // DEBUG
	std::shared_ptr<ParameterPatch> result ( new ParameterPatch );
	while ( 1 ) {
		if ( default_patch_method_edge_ == 0 ) {
			if ( default_patch_method_vertex_ == size () ) { 
				default_patch_method_vertex_ = 0;
        //std::cout << "ParameterSpace::patch. Returning empty patch.\n"; // DEBUG
				return result;
			} else {
				default_patch_method_neighbors_ = adjacencies ( default_patch_method_vertex_ );
			}
		}
		ParameterIndex u = default_patch_method_vertex_;
    if ( default_patch_method_neighbors_ . size () == 0 ) {
      std::cout << "ParameterSpace::patch. Returning isolated parameter patch.\n"; // DEBUG
      result -> vertices . push_back ( u );
      result -> parameter [ u ] = parameter ( u );
      ++ default_patch_method_vertex_;
      return result;
    }
    
    ParameterIndex v = default_patch_method_neighbors_ [ default_patch_method_edge_ ];

		++ default_patch_method_edge_;
		if ( default_patch_method_edge_ == default_patch_method_neighbors_ . size () ) {
			default_patch_method_edge_ = 0;
			++ default_patch_method_vertex_;
		}

		// Perform check to prevent symmetric duplicates
		if ( u < v ) {
			result -> vertices . push_back ( u );
			result -> vertices . push_back ( v );
			result -> edges . push_back ( std::make_pair ( u, v ) );
			result -> parameter [ u ] = parameter ( u );
			result -> parameter [ v ] = parameter ( v );
			break;
		} 

	}
  //std::cout << "ParameterSpace::patch. Returning patch.\n"; //DEBUG
	return result;
}

	/// begin
	///    Return "begin" iterator

inline ParameterSpace::iterator 
ParameterSpace::begin ( void ) const {
	return iterator ( 0 );
}


inline ParameterSpace::iterator
ParameterSpace::end ( void ) const {
	return iterator ( size () );
}

#endif

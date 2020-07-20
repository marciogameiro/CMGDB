/// Configuration.h
#ifndef _CMDP_CONFIGURATION_
#define _CMDP_CONFIGURATION_

#include <exception>
#include <string>
#include <sstream>
#include <boost/serialization/serialization.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/string.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <unistd.h>
#include <fstream>
#include "RectGeo.h"

class Configuration {
public:
  typedef RectGeo Rect;
  // Model, name and description
  std::string MODEL_NAME;
  std::string MODEL_DESC;
  
  /* Parameter Space */
  int PARAM_DIM;
  std::vector<int> PARAM_SUBDIV_DEPTH;
  std::vector<uint64_t> PARAM_SUBDIV_SIZES;
  Rect PARAM_BOUNDS;
  std::vector<bool> PARAM_PERIODIC;

  /* Phase Space */
  int PHASE_DIM;
  int PHASE_SUBDIV_INIT;
  int PHASE_SUBDIV_MIN;
  int PHASE_SUBDIV_MAX;
  int PHASE_SUBDIV_LIMIT;
  Rect PHASE_BOUNDS; 
  std::vector<bool> PHASE_PERIODIC;
  
  // Loading
  void loadFromFile ( const char * filename ) {
    std::string filestring ( filename );
    std::string appendstring ( "/config.xml" );
    std::string loadstring = filestring + appendstring;
    char current [ 100 ];
    getcwd ( current, 100 );
    std::ifstream input(loadstring.c_str());
    if ( not input . good () ) {
      std::cout << "Problem loading configuration file.\n";
      std::cout << "Current directory:\n" << current << "\n";
      std::cout << "Attempted to load from file:\n " << loadstring << "\n";      
      throw std::runtime_error ( "Unable to load configuration file.\n" );
    }

    LoadFromStream(&input);
    //std::cout << "Success.\n";
  }

  void LoadFromString(const std::string& input) {
    std::istringstream stream(input);
    LoadFromStream(&stream);
  }
  
  void LoadFromStream(std::istream* input) {
    
    // Create an empty property tree object
    using boost::property_tree::ptree;
    ptree pt;
    
    // Load the XML file into the property tree. If reading fails
    // (cannot open file, parse error), an exception is thrown.
    
    read_xml(*input, pt);
    
    MODEL_NAME = pt.get<std::string>("config.model.name");
    MODEL_DESC = pt.get<std::string>("config.model.desc");
    
    /* Parameter Space */
    PARAM_DIM = pt.get<int>("config.param.dim");
    //PARAM_SUBDIV_DEPTH = pt.get<int>("config.param.subdiv.depth");

    boost::optional<std::string> param_subdiv_depth = pt.get_optional<std::string>("config.param.subdiv.depth");
    if ( param_subdiv_depth ) {
      std::stringstream param_subdiv_depth_ss ( *param_subdiv_depth );
      int depth;
      while ( param_subdiv_depth_ss >> depth ) PARAM_SUBDIV_DEPTH . push_back ( depth );
      if ( PARAM_SUBDIV_DEPTH . size () == 1 && PARAM_DIM > 1 ) {
        PARAM_SUBDIV_DEPTH . resize ( PARAM_DIM, PARAM_SUBDIV_DEPTH [ 0 ] );
      }
      if ( (int) PARAM_SUBDIV_DEPTH . size () != PARAM_DIM ) {
        std::cout << "Configuration Error. Invalid number of inputs in config.param.subdiv.depth field\n";
        throw 1;
      }
    }
    
    boost::optional<std::string> param_sizes = pt.get_optional<std::string>("config.param.subdiv.sizes");
    if ( param_sizes ) {
      std::stringstream param_sizes_ss ( * param_sizes );
      uint64_t width;
      while ( param_sizes_ss >> width ) PARAM_SUBDIV_SIZES . push_back ( width );
      if ( PARAM_SUBDIV_SIZES . size () == 1 && PARAM_DIM > 1 ) {
        PARAM_SUBDIV_SIZES . resize ( PARAM_DIM, PARAM_SUBDIV_SIZES [ 0 ] );
      }
      if ( (int) PARAM_SUBDIV_SIZES . size () != PARAM_DIM ) {
        std::cout << "Configuration Error. Invalid number of inputs in config.param.subdiv.sizes field\n";
        throw 1;
      }
    } else {
      if ( PARAM_SUBDIV_DEPTH . empty () ) {
        std::cout << "Configuration Error. There are is neither a config.param.subdiv.sizes field nor a config.param.subdiv.depth field\n";
        throw 1;
      }
      PARAM_SUBDIV_SIZES . resize ( PARAM_DIM );
      for ( int d = 0; d < PARAM_DIM; ++ d ) {
        PARAM_SUBDIV_SIZES [ d ] = 1 << PARAM_SUBDIV_DEPTH [ d ];
      }
    }
 
    PARAM_BOUNDS . lower_bounds . resize ( PARAM_DIM );
    PARAM_BOUNDS . upper_bounds . resize ( PARAM_DIM );
    std::string param_lower_bounds = pt.get<std::string>("config.param.bounds.lower");
    std::string param_upper_bounds = pt.get<std::string>("config.param.bounds.upper");
    std::stringstream param_lbss ( param_lower_bounds );
    std::stringstream param_ubss ( param_upper_bounds );
    for ( int d = 0; d < PARAM_DIM; ++ d ) {
      param_lbss >> PARAM_BOUNDS . lower_bounds [ d ];
      param_ubss >> PARAM_BOUNDS . upper_bounds [ d ];
    }
    
    PARAM_PERIODIC . resize ( PARAM_DIM, false );
    boost::optional<std::string> param_periodic = pt.get_optional<std::string>("config.param.periodic");
    if ( param_periodic ) {
      std::stringstream param_periodic_ss ( *param_periodic );
      for ( int d = 0; d < PARAM_DIM; ++ d ) {
        int x;
        param_periodic_ss >> x;
        PARAM_PERIODIC [ d ] = (bool) x;
      }
    }
    
    /* Phase Space */
    PHASE_DIM = pt.get<int>("config.phase.dim");
    PHASE_SUBDIV_MIN = pt.get<int>("config.phase.subdiv.min");
    PHASE_SUBDIV_MAX = pt.get<int>("config.phase.subdiv.max");
    PHASE_SUBDIV_LIMIT = pt.get<int>("config.phase.subdiv.limit");
 
    boost::optional<int> opt_phase_subdiv_init = pt.get_optional<int>("config.phase.subdiv.init");
    PHASE_SUBDIV_INIT = 0;
    if ( opt_phase_subdiv_init ) PHASE_SUBDIV_INIT = opt_phase_subdiv_init . get ();
    
    PHASE_BOUNDS . lower_bounds . resize ( PHASE_DIM );
    PHASE_BOUNDS . upper_bounds . resize ( PHASE_DIM );
    std::string phase_lower_bounds = pt.get<std::string>("config.phase.bounds.lower");
    std::string phase_upper_bounds = pt.get<std::string>("config.phase.bounds.upper");
    std::stringstream phase_lbss ( phase_lower_bounds );
    std::stringstream phase_ubss ( phase_upper_bounds );
    for ( int d = 0; d < PHASE_DIM; ++ d ) {
      phase_lbss >> PHASE_BOUNDS . lower_bounds [ d ];
      phase_ubss >> PHASE_BOUNDS . upper_bounds [ d ];
    }
    
    PHASE_PERIODIC . resize ( PHASE_DIM, false );
    boost::optional<std::string> phase_periodic = pt.get_optional<std::string>("config.phase.periodic");
    if ( phase_periodic ) {
      std::stringstream phase_periodic_ss ( *phase_periodic );
      for ( int d = 0; d < PHASE_DIM; ++ d ) {
        int x;
        phase_periodic_ss >> x;
        PHASE_PERIODIC [ d ] = (bool) x;
      }
    }
    
    
  }
  
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version)
  {
    // Model, name and description
    ar & MODEL_NAME;
    ar & MODEL_DESC;
    
    /* Parameter Space */
    ar & PARAM_DIM;
    ar & PARAM_SUBDIV_DEPTH;
    ar & PARAM_SUBDIV_SIZES;
    ar & PARAM_PERIODIC;
    ar & PARAM_BOUNDS;
    
    /* Phase Space */
    ar & PHASE_DIM;
    ar & PHASE_SUBDIV_INIT;
    ar & PHASE_SUBDIV_MIN;
    ar & PHASE_SUBDIV_MAX;
    ar & PHASE_SUBDIV_LIMIT;
    ar & PHASE_BOUNDS; 
    ar & PHASE_PERIODIC;
  }
  
};

#endif

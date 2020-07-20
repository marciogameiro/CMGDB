// HEADERS FOR DEALING WITH PICTURES
#ifndef CMDB_SINGLEOUTPUT_H
#define CMDB_SINGLEOUTPUT_H

#include <string>
// #include "picture.h"

#include "MorseGraph.h"
#include "Grid.h"
#include "chomp/ConleyIndex.h"
#include "chomp/Matrix.h"
#include "chomp/PolyRing.h"
#include "conleyIndexString.h"

/***********/
/* Output  */
/***********/
//void DrawMorseSets ( const Grid & phase_space, const CMG & conley_morse_graph );
//void CreateDotFile ( const MorseGraph & cmg );
//void output_cubes ( const Grid & my_grid, const MorseGraph & conley_morse_graph );

inline std::string returnConleyIndex ( const chomp::ConleyIndex_t & ci );

// inline void DrawMorseSets ( const TreeGrid & phase_space, const MorseGraph & conley_morse_graph ) {
//   // Create a Picture
//   int Width =  4096;
//   int Height = 4096;
  
//   std::cout << "DrawMorseSets.\n";
//   Picture * picture = draw_morse_sets( Width, Height, phase_space, conley_morse_graph );
//   Picture * picture2 = draw_grid ( Width, Height, phase_space );
//   Picture * picture3 = draw_grid_and_morse_sets( Width, Height, phase_space, conley_morse_graph );
  
//   std::cout << "Saving output... ";
//   std::cout . flush ();
//   //picture -> saveAsPNG ( "morse_sets.png" );
//   //picture2 -> saveAsPNG ( "grid.png" );
//   //picture3 -> saveAsPNG ( "grid_and_morse.png" );
  
//   picture -> saveAsBMP ( "morse_sets.bmp" );
//   picture2 -> saveAsBMP ( "grid.bmp" );
//   picture3 -> saveAsBMP ( "grid_and_morse.bmp" );
//   std::cout << "Output saved.\n";
  
//   //picture -> saveAsTIFF ( "morse_sets.tiff" );
//   //picture2 -> saveAsTIFF ( "grid.tiff" );
//   //picture3 -> saveAsTIFF ( "grid_and_morse.tiff" );

//   delete picture;
//   delete picture2;
//   delete picture3;
// }


std::string conleyStringForZoo ( const std::vector<std::string> & ci_strings ) {
    // data
    std::string ci_string;
    // algo
    std::stringstream ss;
    ss << "(";
    for ( int d = 0; d < (int) ci_strings . size (); ++ d ) {
      std::string s = ci_strings [ d ];
      std::string t;
      for ( int i = 0; i < (int)s . size () - 1; ++ i ) {
        if ( s [ i ] == '.' ) continue;
        if ( s [ i ] == ' ' ) continue;
        if ( s [ i ] == '\n' ) { t += std::string ( ", "); continue; }
        t . push_back ( s [ i ] );
      }
      ss << t;
      if ( d != (int)ci_strings . size () - 1 ) ss << ", ";
    }
    ss << ")";
    ci_string = ss . str ();
    return ci_string;
}

inline void CreateDotFile ( const char * outputfile, const MorseGraph & cmg ) {
  typedef MorseGraph::Vertex V;
  typedef MorseGraph::VertexIterator VI;
  typedef MorseGraph::EdgeIterator EI;
  


  // std::ofstream outfile ("morsegraph.gv");
  std::ofstream outfile (outputfile);
  
  outfile << "digraph G { \n";
  //outfile << "node [ shape = point, color=black  ];\n";
  //outfile << "edge [ color=red  ];\n";
  
  // LOOP THROUGH VERTICES AND GIVE THEM NAMES
  std::map < V, int > vertex_to_index;
  VI start, stop;
  int i = 0;
  for (boost::tie ( start, stop ) = cmg . Vertices (); start != stop; ++ start ) {
    vertex_to_index [ *start ] = i;
    //outfile << i << " [label=\""<< cmg . grid (*start) -> size () << "\"]\n";
    // Label the Morse Graph set with their Conley index
    outfile << i;
    if ( cmg . conleyIndex ( *start ) ) outfile << " [label=\""<< conleyStringForZoo(conleyIndexString ( * cmg . conleyIndex ( *start ) )) << "\"]";
    outfile << "\n";
    ++ i;
  }
  int N = cmg . NumVertices ();
  
  // LOOP THROUGH MorseGraph EDGES
  EI estart, estop;
  typedef std::pair<int, int> int_pair;
  std::set < int_pair > edges;
  for (boost::tie ( estart, estop ) = cmg . Edges ();
       estart != estop;
       ++ estart ) {
    V source = estart -> first;
    V target = estart -> second;
    int index_source = vertex_to_index [ source ];
    int index_target = vertex_to_index [ target ];
    if ( index_source != index_target ) // Cull the self-edges
      edges . insert ( std::make_pair ( index_source, index_target ) );
  }

  // Transitive Reduction
  // Assumption: edges is already transitively closed (and has no self-edges)
  // Technique: R = G - G^2

  boost::unordered_map < int, boost::unordered_set < int > > G, squared;
  BOOST_FOREACH ( int_pair edge, edges ) {
    if ( edge . first != edge . second ) {
      G [ edge . first ] . insert ( edge . second );
    }
  }
  for ( int u = 0; u < N; ++ u ) {
    BOOST_FOREACH( int v, G [ u ] ) {
      BOOST_FOREACH ( int w, G [ v ] ) {
        squared [ u ] . insert ( w );
      }
    }
  }
  std::vector < int_pair > reduced;
  BOOST_FOREACH ( int_pair edge, edges ) {
    if ( squared [ edge . first ] . count ( edge . second ) == 0 )
      reduced . push_back ( edge );
  }
  
  // PRINT OUT EDGES OF TRANSITIVE CLOSURE
  BOOST_FOREACH ( int_pair edge, reduced ) {
      outfile << edge . first << " -> " << edge . second << ";\n";
  }
  
  outfile << "}\n";
  outfile . close ();
  
}

inline std::string returnConleyIndex ( const chomp::ConleyIndex_t & ci ) {
  using namespace chomp;

  std::string resultstr;
  std::stringstream sstr;

  if ( ci . undefined () ) resultstr = "NaN";
  int biggest = 0;
  static int corrupt = 0;
  for ( unsigned int i = 0; i < ci . data () . size (); ++ i ) {
    typedef SparseMatrix < PolyRing < Ring > > PolyMatrix;
    if ( ( ci . data () [ i ] . number_of_rows () !=
          ci . data () [ i ] . number_of_columns () ) ||
        ci . data () [ i ] . number_of_rows () < 0 ) {
      std::cout << "return ConleyIndex : " << "Corrupt.\n";
      std::cout << ++ corrupt << "\n";
      continue;
    }
    
    PolyMatrix poly = ci . data () [ i ];
    
    int N = poly . number_of_rows ();
    PolyRing<Ring> X;
    X . resize ( 2 );
    X [ 1 ] = Ring ( -1 );
    for ( int i = 0; i < N; ++ i ) {
      poly . add ( i, i, X );
    }
    PolyMatrix U, Uinv, V, Vinv, D;
    //std::cout << "SNF:\n";
    SmithNormalForm ( &U, &Uinv, &V, &Vinv, &D, poly );
    //std::cout << "SNF complete.\n";
    bool is_trivial = true;
    PolyRing < Ring > x;
    x . resize ( 2 );
    x [ 1 ] = Ring ( 1 );
    for ( int j = 0; j < D . number_of_rows (); ++ j ) {
      //std::stringstream ss; //
      PolyRing < Ring > entry = D . read ( j, j );
      while ( ( entry . degree () >= 0 )
             && ( entry [ 0 ] == Ring ( 0 ) )) {
        entry = entry / x;
      }
      if ( entry . degree () <= 0 ) continue;
      is_trivial = false;
      //ss << "   " << entry << "\n";
      sstr << entry;
      resultstr = sstr.str();
      if ( entry . degree () > biggest ) biggest = entry . degree ();
    }
    if ( is_trivial ) resultstr="trivial";
  }

  //return biggest;


  //std::cout << "Biggest " << biggest << "\n";
  return resultstr;

}





int check_index ( std::ostream & outstream, const chomp::ConleyIndex_t & ci ) {
  using namespace chomp;
  if ( ci . undefined () ) outstream << "Conley Index not computed for this Morse Set.\n";
  int biggest = 0;
  static int corrupt = 0;
  //std::cout << "check index.\n";
  //std::cout << "data . size () = " << ci . data () . size () << "\n";
  for ( unsigned int i = 0; i < ci . data () . size (); ++ i ) {
    typedef SparseMatrix < PolyRing < Ring > > PolyMatrix;
    if ( ( ci . data () [ i ] . number_of_rows () !=
          ci . data () [ i ] . number_of_columns () ) ||
        ci . data () [ i ] . number_of_rows () < 0 ) {
      outstream << "Corrupt.\n";
      outstream << ++ corrupt << "\n";
      continue;
    }
    //if ( ci . data () [ i ] . number_of_rows () < 4 ) continue;
    //std::cout << "Dimension " << i << "\n";
    //print_matrix ( ci . data () [ i ] );
    
    PolyMatrix poly = ci . data () [ i ];
    
    int N = poly . number_of_rows ();
    PolyRing<Ring> X;
    X . resize ( 2 );
    X [ 1 ] = Ring ( -1 );
    for ( int i = 0; i < N; ++ i ) {
      poly . add ( i, i, X );
    }
    PolyMatrix U, Uinv, V, Vinv, D;
    //std::cout << "SNF:\n";
    SmithNormalForm ( &U, &Uinv, &V, &Vinv, &D, poly );
    //std::cout << "SNF complete.\n";
    bool is_trivial = true;
    PolyRing < Ring > x;
    x . resize ( 2 );
    x [ 1 ] = Ring ( 1 );
    for ( int j = 0; j < D . number_of_rows (); ++ j ) {
      //std::stringstream ss; //
      PolyRing < Ring > entry = D . read ( j, j );
      while ( ( entry . degree () >= 0 )
             && ( entry [ 0 ] == Ring ( 0 ) )) {
        entry = entry / x;
      }
      if ( entry . degree () <= 0 ) continue;
      is_trivial = false;
      //ss << "   " << entry << "\n";
      outstream << entry << "\n";
      if ( entry . degree () > biggest ) biggest = entry . degree ();
    }
    if ( is_trivial ) outstream << "Trivial.\n";
  }
  return biggest;
}

#if 0
inline void output_cubes ( const Grid & my_grid,
                   const MorseGraph & conley_morse_graph ) {
  using namespace chomp;
  
  // Loop Through Morse Sets to determine bounds
  typedef typename MorseGraph::VertexIterator VI;
  VI it, stop;
  std::vector < std::vector < uint32_t > > cubes;
  for (boost::tie ( it, stop ) = conley_morse_graph . Vertices (); it != stop; ++ it ) {
    CellContainer const & my_subset = conley_morse_graph . CellSet ( *it );
    int depth = my_grid . getDepth ( my_subset );
    BOOST_FOREACH ( const Grid::GridElement & ge, my_subset ) {
      my_grid . GridElementToCubes ( & cubes, ge, depth  );
    }
  }
  std::ofstream outfile ( "morsecubes.txt" );
  for ( uint32_t i = 0; i < cubes . size (); ++ i ) {
    for ( uint32_t j = 0; j < cubes [ i ] . size (); ++ j ) {
      outfile << cubes [ i ] [ j ] << " ";
    }
    outfile << "\n";
  }
  outfile . close ();
} /* output_cubes */
#endif

#endif



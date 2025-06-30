// RelativeMapHomology.h
// Shaun Harker
// 9/16/11
///
/// Marcio Gameiro
/// 2025-06-20

#ifndef CHOMP_RELATIVEMAPHOMOLOGY_H
#define CHOMP_RELATIVEMAPHOMOLOGY_H

#include <cstdlib>
#include <vector>
#include <ctime>
#include "boost/unordered_map.hpp"
#include "boost/unordered_set.hpp"
#include "chomp/Matrix.h"
#include "chomp/Generators.h"
#include "chomp/Closure.h"
#include "chomp/FiberComplex.h"
#include "chomp/RelativePair.h"


//#include "Draw.h"

#define PRINT if ( 0 ) std::cout <<

namespace chomp {

typedef std::vector < SparseMatrix < Ring > > RelativeMapHomology_t;

template < class Grid, class Container, class RectMap > int
RelativeMapHomology (RelativeMapHomology_t * output, 
                     const Grid & source,
                     const Container & XE, 
                     const Container & AE,
                     const Grid & target,
                     const Container & YE, 
                     const Container & BE,
                     const RectMap & F,
                     int depth ) {
                     
  
  clock_t total_time_start = clock ();
  clock_t construction_time_start = clock ();

  // Create the Relative Complexes
  RelativePair domain_pair;
  RelativePair codomain_pair;
  // Determine depth.
  // TODO: determine it, don't ask for it as input
  source . relativeComplex ( &domain_pair, XE, AE, depth );
  // TODO: don't repeat work here
  target . relativeComplex ( &codomain_pair, YE, BE, depth );
  
  //PRINT "RMH: Constructed relative complexes\n";
  
  CubicalComplex & full_domain = domain_pair . base ();
  CubicalComplex & full_codomain = codomain_pair . base ();
  BitmapSubcomplex & domain = domain_pair . pair ();
  BitmapSubcomplex & codomain = codomain_pair . pair ();
  BitmapSubcomplex & DomainA = domain_pair . relative ();

    int D = full_domain . dimension ();

  PRINT "RMH: Domain size = " << full_domain . size () << "\n";
  PRINT "RMH: Codomain size = " << full_codomain . size () << "\n";
  PRINT "RMH: Dimension of domain = " << D << "\n";
  
  PRINT "RMH: (X, A) size = (" << domain_pair.pair().size() << ", " 
        << domain_pair . relative () . size () << ")\n";
  PRINT "RMH: (Y, B) size = (" << codomain_pair.pair().size() << ", " 
        << codomain_pair . relative () . size () << ")\n";

  PRINT "RMH: (X, A) top size = (" << domain_pair.pair().size(D) << ", " 
        << domain_pair . relative () . size (D) << ")\n";
  PRINT "RMH: (Y, B) top size = (" << codomain_pair.pair().size(D) << ", " 
        << codomain_pair . relative () . size (D) << ")\n";


  /// Generate "domain_GridElements_X" and "domain_GridElements_A" tables
  /// These are for determining which top cells are involved in
  /// constructing a fiber (and which are relative)
  std::vector < std::vector < boost::unordered_set < Index > > > 
  domain_GridElements_X ( D + 1 );
  std::vector < std::vector < boost::unordered_set < Index > > > 
  domain_GridElements_A ( D + 1 );
  
  
  PRINT "Computing domain grid elements \n";

  for ( int d = 0; d <= D; ++ d ) {
    domain_GridElements_X [ d ] . resize ( full_domain . size ( d ) );
    domain_GridElements_A [ d ] . resize ( full_domain . size ( d ) );
  }
  
  Index N;
  N = full_domain . size (D);
  for ( Index i = 0; i < N; ++ i ) {
    //PRINT " index = " << i << "\n";
    std::vector < boost::unordered_set < Index > > closure_of_i ( D + 1 );
    closure_of_i [ D ] . insert ( i );
    closure ( closure_of_i, full_domain );
    for ( int d = 0; d <= D; ++ d ) {
      //PRINT " d = " << d << "\n";
      BOOST_FOREACH ( Index j, closure_of_i [ d ] ) {
        domain_GridElements_X [ d ] [ j ] . insert ( i );
      }
    }
  }
  
  PRINT "Computing domain grid elements 2\n";
  N = DomainA . size (D);
  for ( Index k = 0; k < N; ++ k ) {
    Index i = DomainA . indexToCell ( k, D );
    //PRINT " index = " << i << "\n";
    std::vector < boost::unordered_set < Index > > closure_of_i ( D + 1 );
    closure_of_i [ D ] . insert ( i );
    closure ( closure_of_i, full_domain );
    for ( int d = 0; d <= D; ++ d ) {
      //PRINT " d = " << d << "\n";
      BOOST_FOREACH ( Index j, closure_of_i [ d ] ) {
        domain_GridElements_A [ d ] [ j ] . insert ( i );
      }
    }
  }
  
  PRINT "computed domain_GridElements_X and _A\n";
  double construction_time = ((double)(clock()-construction_time_start)/(double)CLOCKS_PER_SEC);

  PRINT "RMH: Time Elapsed (construction) = " << construction_time << "\n";

  // DEBUG //////////////////////////////////
//  for ( int d = 0; d <= D; ++ d ) {
//    std::cout << "The number of " << d << "-cells in domain_GridElements_X is " << 
//    domain_GridElements_X [ d ] . size () << "\n";
//    std::cout << "The number of " << d << "-cells in domain_GridElements_A is " << 
//    domain_GridElements_A [ d ] . size () << "\n";    
//  }
//  MorseComplex debugcomplex ( domain );
//  MorseSanity ( debugcomplex );
  //////////////////////////////////////////
  
  int cutoff_dimension = full_domain . dimension ();
#ifdef CONLEYINDEXCUTOFF
  cutoff_dimension = CONLEYINDEXCUTOFF;
#endif
  // Compute the homology generators
  clock_t homology_time_start = clock ();
  PRINT "RMH: Computing Domain Generators with MorseGenerators\n";
  Generators_t domain_gen = MorseGenerators ( domain, cutoff_dimension );
  PRINT "RMH: Computing Morse Complex of Codomain (codomain_morse)\n";
  MorseComplex codomain_morse ( codomain );
  PRINT "RMH: Computing Generators of codomain_morse\n";
  Generators_t codomain_morse_gen = SmithGenerators ( codomain_morse, cutoff_dimension );

  // DEBUG

  // Concern: Are the lifted smith generators in codomain_morse the same as the domain_gen
  //          generators of domain_gen? For a Conley Index computation, we depend on this.
  //          This debug code will check this explicitly.
  /*
  std::cout << "Domain generator basis: \n";
  for ( int d = 0; d < domain_gen . size (); ++ d ) {
    std::cout << "Dimension " << d << "\n";
    for ( int gen = 0; gen < domain_gen [ d ] . size (); ++ gen ) {
      Chain c = domain_gen [ d ] [ gen ] . first;
      std::cout << "Domain chain " << gen << " = " << c << "\n";
    }
  }

  std::cout << "Codomain generator basis: \n";
  for ( int d = 0; d < codomain_morse_gen . size (); ++ d ) {
    std::cout << "Dimension " << d << "\n";
    for ( int gen = 0; gen < codomain_morse_gen [ d ] . size (); ++ gen ) {
      Chain c = codomain_morse_gen [ d ] [ gen ] . first;
      Chain lift = codomain_morse . lift ( c );
      std::cout << "Codomain chain " << gen << " = " << lift << "\n";
    }
  }
  */
  // END DEBUG

  // Compute the cycles projected into the codomain through the graph 
  std::vector < std::vector < Chain > > codomain_cycles ( D + 1 );
  double homology_time = ((double)(clock()-homology_time_start)/(double)CLOCKS_PER_SEC);
  PRINT "RMH: Time Elapsed (homology) = " << homology_time << "\n";

  PRINT "Main Loop\n";
  // Loop through each domain generator
  long max_chain_memory = 0;
  long number_of_analyzed_fibers = 0;
  long explored_graph_complex = 0;
#ifdef RMHMEASUREGRAPH
  int highest_nontrivial_dim = -1;
#endif
  clock_t lift_time_start = clock ();
  bool acyclic_map = true;
  for ( int d = 0; d <= D; ++ d ) {
#ifdef CONLEYINDEXCUTOFF
    if ( d == CONLEYINDEXCUTOFF ) break;
#endif
    if ( not acyclic_map ) break;
    int number_of_domain_gen = domain_gen [ d ] . size ();
    codomain_cycles [ d ] . resize ( number_of_domain_gen );

    PRINT "RMH: Dimension = " << d << "\n";
    PRINT "RMH: Number of Generators at this dimension = " << number_of_domain_gen << "\n";

    for ( int gi = 0; gi < number_of_domain_gen; ++ gi ) {
      if ( not acyclic_map ) break;
#ifdef RMHMEASUREGRAPH
      highest_nontrivial_dim = d;
#endif
      PRINT "RMH: dimension = " << d << ", generator = " << gi << "\n";

      long chain_memory = 0;

      
      // Lift the domain generator into the Relative Graph Complex
      Chain answer;
      answer . dimension () = d;
      std::vector < boost::unordered_map < Index, Chain > > graph_boundary ( D + 1 );
      //   Begin by "Chain Lift"
      PRINT "RMH: Chain Lift \n";
      // First, get the domain generator and include it into full_domain
      Chain domain_generator = domain . include ( domain_gen [ d ] [ gi ] . first );
      
#if 0
      if ( domain_generator () . size () > 1000 ) {
        std::cout << "Relative size = " << domain_pair . relative () . size () << "\n";
        ComplexVisualization * cv = new ComplexVisualization ( "Domain Generator Picture." );
        std::cout << "About to draw complex.\n";
        cv -> drawRelativeComplex ( domain_pair, 100, 200 );
        std::cout << "The complex was drawn.\n";
        std::cout << domain_generator . dimension () << "\n";
        cv -> drawChain ( full_domain, domain_generator, 200 );
        // explore a minute here
        cv -> explore ();
        delete cv;
      }
#endif 
      
      
      PRINT "RMH: Chain size = " << domain_generator () . size () << "\n";
      // Now lift it:
      BOOST_FOREACH ( const Term & t, domain_generator () ) {
        //std::cout << "CHECKPOINT X\n";

        if ( not acyclic_map ) break;

        //std::cout << "Dealing with fiber (" << t . index () << ", " << d << ") (*)\n";
        //std::cout << "   Term is " << t << "\n";
        // Determine fiber
        boost::unordered_set < Index > X_nbs, A_nbs;
        
        //std::cout << "domain_GridElements_X [ " << d << " ] . size () = " <<
        //domain_GridElements_X [ d ] . size () << "\n";
        
        //std::cout << "domain_GridElements_A [ " << d << " ] . size () = " <<
        //domain_GridElements_A [ d ] . size () << "\n";
        
        X_nbs = domain_GridElements_X [ d ] [ t . index () ];
        A_nbs = domain_GridElements_A [ d ] [ t . index () ];    
        
        //std::cout << "CHECKPOINT X2\n";

        // DEBUG ///////////////////////
        //std::cout << "   X_nbs . size () = " << X_nbs . size () << "\n";
        //std::cout << "   A_nbs . size () = " << A_nbs . size () << "\n";
        ////////////////////////////////
        
        FiberComplex fiber ( X_nbs, A_nbs, full_domain, full_codomain, F );
        //std::cout << "CHECKPOINT Y\n";


        if ( fiber . size () == 0 ) {
          std::cout << "UNEXPECTED: CANT LIFT CHAIN DUE TO EMPTY FIBER.\n";
          acyclic_map = false;
          break;
        }

#ifndef CONLEY_INDEX_NO_ACYCLIC_CHECKS
        if ( not fiber . acyclic () ) {
          acyclic_map = false;
          break;
        }
#endif
        ++ number_of_analyzed_fibers; // run statistics
        explored_graph_complex += fiber . size (); // run statistics

        // Choose a vertex in fiber
        Index choice = fiber . choose ();
        if ( d == 0 ) answer += Term ( choice, t . coef () );
        Chain bd = full_domain . boundary ( t . index (), d );
        BOOST_FOREACH ( const Term & s, bd () ) {
          graph_boundary [d-1] [ s . index () ] . dimension () = 0;
          graph_boundary [d-1] [ s . index () ] += Term ( choice, t . coef () * s . coef () );
        }
      }
      PRINT "RMH: Cycle Lift -- iterative preboundary loop begins now\n";

      // Now process every fiber we can see in graph_boundary
      // Loop through every fiber dimension (fd)
      for ( int fd = d - 1; fd >= 0; -- fd ) {
        if ( not acyclic_map ) break;
        PRINT "  fiber dimension = " << fd << "\n";

        // Loop through every non-empty fiber
        typedef std::pair<Index,Chain> IndexChainPair;
        BOOST_FOREACH ( const IndexChainPair & fiberchain, graph_boundary [fd] ) {
          if ( not acyclic_map ) break;
          //std::cout << "Dealing with fiber (" << fiberchain . first << ", " << fd << ")\n";
          //std::cout << "   Chain is " << fiberchain . second << "\n";
          chain_memory += fiberchain . second () . size ();
          // Determine fiber
          boost::unordered_set < Index > X_nbs, A_nbs;
          X_nbs = domain_GridElements_X [ fd ] [ fiberchain . first ];
          A_nbs = domain_GridElements_A [ fd ] [ fiberchain . first ];    
          FiberComplex fiber ( X_nbs, A_nbs, full_domain, full_codomain, F );
          
          // TODO: option to turn this check off.
#ifndef CONLEY_INDEX_NO_ACYCLIC_CHECKS
          if ( not fiber . acyclic_or_trivial () ) {
            acyclic_map = false;
            break;
          }
#endif
          
          ++ number_of_analyzed_fibers; // run statistics
          explored_graph_complex += fiber . size (); // run statistics
          // Determine chain in fiber
          Chain projected = fiber . project ( fiberchain . second );
          // Work out the preboundary
          Chain preboundary = fiber . preboundary ( projected );
          // Include the preboundary back into the codomain
          Chain included_preboundary = fiber . include ( preboundary );

          // DEBUG///////////////////////////////////////////////////////////////////////
#if 0

          std::cout << "   A_nbs . size () = " << A_nbs . size () << "\n";
          std::cout << "   Projected Chain is " << projected << "\n";
          std::cout << "   Preboundary Chain is " << preboundary << "\n";
          std::cout << "   Included Chain is " << included_preboundary << "\n";

          Chain reboundary = fiber . boundary ( preboundary );
          std::cout << "   bd(preboundary) = " << simplify(reboundary) << "\n";

          Chain increboundary = full_codomain . boundary ( included_preboundary );
          std::cout << "   bd(included) = " << simplify(increboundary) << "\n";
          
          Chain leftover = simplify ( projected - reboundary );
          if ( not leftover () . empty () ) {
            std::cout << "   left over = " << simplify(projected - reboundary) << "\n";
            exit ( 1 );
          }

          MorseComplex morse ( fiber );
          std::cout << "\n Number of Cells = " << fiber . size () << ": ";
          for ( int dd = 0; dd <= fiber . dimension (); ++ dd ) 
            std::cout << fiber . size ( dd ) << " ";
          std::cout << "\n Number of Critical Cells = " << morse . size () << ": ";
          for ( int dd = 0; dd <= morse . dimension (); ++ dd ) 
            std::cout << morse . size ( dd ) << " ";
          
          {
          ComplexVisualization * cv = new ComplexVisualization ( "Fiber Picture." );
            std::cout << "About to draw complex.\n";
          cv -> drawComplex ( fiber, 100 );
            std::cout << "The complex was drawn.\n";
            std::cout << projected . dimension () << "\n";
          cv -> drawChain ( fiber, projected, 200 );
            std::cout << "The projected chain drew.\n";
          cv -> drawChain ( fiber, preboundary, 200 );
          // explore a minute here
          cv -> explore ();
          delete cv;
          }
#endif

          ///////////////////////////////////////////////////////////////////////////////
          
          // If this is a zero-dim fiber, add to answer
          if ( fd == 0 ) answer -= included_preboundary;
          // Add preboundary to adjacent fibers
          Chain bd = full_domain . boundary (fiberchain . first, fd );
          Ring grade;
          if ( fd % 2 ) grade = Ring ( -1 ); else grade = Ring ( 1 );
          BOOST_FOREACH ( const Term & s, bd () ) {
            graph_boundary [fd-1] [ s . index () ] . dimension () = d - fd;
            graph_boundary [fd-1] [ s . index () ] -= included_preboundary * 
                                                      (grade * s . coef ());
          }
        }
      }

      //std::cout << "Projection to codomain.\n";
      // First project into the relative codomain complex
      Chain relative_answer = simplify ( codomain . project ( answer ) );
      
      // DEBUG
      // Concern: We want to see what the projection onto the codomain is and make sure
      //          it is being properly handled.
      //std::cout << "Projection to codomain (" << d << ", " << gi << "): " << relative_answer << "\n";
      // END DEBUG

      // Project the Relative Graph Cycle to the Codomain
      codomain_cycles [ d ] [ gi ] = codomain_morse . lower ( relative_answer );
      
      if ( chain_memory > max_chain_memory ) max_chain_memory = chain_memory;
      PRINT "RMH: Memory usage (chains only)= " << chain_memory << "\n";

    }
  }
  double lift_time = ((double)(clock()-lift_time_start)/(double)CLOCKS_PER_SEC);
  PRINT "RMH: Time Elapsed (lift and project) = " << lift_time << "\n";
  PRINT "RMH: Max memory usage (chains only)= " << max_chain_memory << "\n";

  if ( not acyclic_map ) {
    PRINT "RMH: Returning no answer due to lack of acyclicity in map.";
    return 1;
  }
  
  PRINT "RMH: Algebraic processing. \n";

 clock_t algebra_time_start = clock ();
  // Loop through each Codomain Cycle
  for ( int d = 0; d <= codomain_morse . dimension (); ++ d ) {
#ifdef CONLEYINDEXCUTOFF
    if ( d == CONLEYINDEXCUTOFF ) break;
#endif
    Matrix G = chainsToMatrix ( codomain_morse_gen [ d ], codomain_morse, d );
    Matrix Z = chainsToMatrix ( codomain_cycles [ d ], codomain_morse, d );
    Matrix MapHom = SmithSolve (G, Z);
    //print_matrix (G);
    //print_matrix (Z);
    if ( MapHom . number_of_rows () != 0 ) {
      std::cout << "Dimension " << d << ":\n";
      print_matrix (MapHom);
    }
    output -> push_back ( MapHom );
  }
  for ( int d = codomain_morse . dimension () + 1; d <= D; ++ d ) {
    output -> push_back ( Matrix (0, 0) );
  }
  double algebra_time = ((double)(clock()-algebra_time_start)/(double)CLOCKS_PER_SEC);
  double total_time = ((double)(clock()-total_time_start)/(double)CLOCKS_PER_SEC);

  PRINT "RMH: Time Elapsed (algebra time) = " << algebra_time << "\n";
  PRINT "RMH: Time Elapsed (total) = " << total_time << "\n";

#ifdef RMHMEASUREGRAPH
  long graph_size = 0;
  long edge_graph_size = 0;
  for ( int d = 0; d <= D; ++ d ) {
#ifdef CONLEYINDEXCUTOFF
    if ( d == CONLEYINDEXCUTOFF ) break;
#endif
    for ( Index i = 0; i < full_domain . size ( d ); ++ i ) {
      if ( rand () % 10 != 0 ) continue;
      boost::unordered_set < Index > X_nbs, A_nbs;
      X_nbs = domain_GridElements_X [ d ] [ i ];
      A_nbs = domain_GridElements_A [ d ] [ i ];   
      FiberComplex fiber ( X_nbs, A_nbs, full_domain, full_codomain, F );
      graph_size += 10 * fiber . size ();
      if ( d == D ) edge_graph_size += 10 * fiber . size ( D );
    }
  }
  PRINT "RMH: Graph Size = " << graph_size << "\n";
  if ( highest_nontrivial_dim > -1 ) {
  std::cout << "# ";
  std::cout << D << " ";
  std::cout << highest_nontrivial_dim << " ";
  std::cout << full_domain . size (D) << " ";  
  std::cout << full_domain . size () << " ";
  std::cout << full_codomain . size (D) << " ";
  std::cout << full_codomain . size () << " ";
  std::cout << edge_graph_size << " ";
  std::cout << graph_size << " ";
  std::cout << explored_graph_complex << " ";
  std::cout << number_of_analyzed_fibers << " ";
  std::cout << construction_time << " ";
  std::cout << homology_time << " ";
  std::cout << lift_time << " ";
  std::cout << algebra_time << " ";
  std::cout << total_time << " ";
  std::cout << max_chain_memory << " ";

  std::cout << "\n";
  }
#endif
  PRINT "RMH: Returning.\n";
  return 0;
}

int RelativeSelfMapHomology (RelativeMapHomology_t * output,
                             const std::vector < uint64_t > & X_cubes,
                             const std::vector < uint64_t > & A_cubes,
                             const std::vector < uint64_t > & sizes,
                             const std::vector < bool > & periodic,
                             const std::unordered_map < uint64_t, std::vector < uint64_t > > & F_map,
                             bool acyclic_check = true) {
  // Cubical complex dimension
  int D = sizes . size ();

  // std::cout << "X_cubes size: " << X_cubes . size () << std::endl;
  // std::cout << "A_cubes size: " << A_cubes . size () << std::endl;
  // std::cout << "Dimension: " << D << std::endl;

  // Create the Relative Complexes

  // Produce the full complex
  CubicalComplex * full_complex = new CubicalComplex;
  CubicalComplex & X = *full_complex;

  // Cubical complex coordinates from index
  auto coordinates = [&] ( uint64_t index ) {
    std::vector < uint64_t > coords ( D );
    for ( int d = 0; d < D; ++ d ) {
      coords [ d ] = index % sizes [ d ];
      index = index / sizes [ d ];
    }
    return coords;
  };

  // Learn the bounds of the set of cubes
  std::vector < uint64_t > min_cube ( D, -1 );
  std::vector < uint64_t > max_cube ( D, 0 );
  BOOST_FOREACH ( uint64_t index, X_cubes ) {
    std::vector <uint64_t> cube = coordinates ( index );
    for ( int d = 0; d < D; ++ d ) {
      if ( min_cube [ d ] > cube [ d ] ) min_cube [ d ] = cube [ d ];
      if ( max_cube [ d ] < cube [ d ] ) max_cube [ d ] = cube [ d ];
    }
  }

  // Get X sizes and reset periodic flag if not touching domain boundary
  std::vector < uint64_t > dimension_sizes ( D, 1 );
  std::vector < bool > is_periodic = periodic;
  for ( int d = 0; d < D; ++ d ) {
    dimension_sizes [ d ] = max_cube [ d ] - min_cube [ d ] + 1;
    if ( min_cube [ d ] > 0 )               is_periodic [ d ] = false;
    if ( max_cube [ d ] < sizes [ d ] - 1 ) is_periodic [ d ] = false;
  }

  // Initialize cubical complex X
  X . initialize ( dimension_sizes, is_periodic );

  // Add list of cubes X_cubes to complex X
  BOOST_FOREACH ( uint64_t index, X_cubes ) {
    std::vector <uint64_t> cube = coordinates ( index );
    std::vector < uint64_t > offset ( D );
    for ( int d = 0; d < D; ++ d )
      offset [ d ] = cube [ d ] - min_cube [ d ];
    X . addFullCube ( offset );
  }

  // Finalize X construction
  X . finalize ();

  // std::cout << "X dimension: " << X . dimension () << std::endl;
  // std::cout << "X size: " << X . size () << std::endl;

  // Get mapping from cubical complex index to chomp index
  std::unordered_map < uint64_t, uint64_t > index2chomp;
  BOOST_FOREACH ( uint64_t index, X_cubes ) {
    std::vector <uint64_t> cube = coordinates ( index );
    std::vector < uint64_t > offset ( D );
    for ( int d = 0; d < D; ++ d )
      offset [ d ] = cube [ d ] - min_cube [ d ];
    uint64_t chomp_index = X . cubeIndex ( offset );
    index2chomp [ index ] = chomp_index;
  }

  // Get map F in terms of cells indices in X
  std::unordered_map < uint64_t, std::vector < uint64_t > > F;
  for ( const auto& F_pair : F_map ) {
    uint64_t key = F_pair . first;
    uint64_t chomp_key = index2chomp [ key ];
    F [ chomp_key ] . reserve ( F_pair . second . size () );
    for ( uint64_t index : F_pair.second ) {
      uint64_t chomp_index = index2chomp [ index ];
      F [ chomp_key ] . push_back ( chomp_index );
    }
  }

  // Produce the relative complex
  BitmapSubcomplex * pair_complex = new BitmapSubcomplex ( X, true );
  BitmapSubcomplex * rel_complex = new BitmapSubcomplex ( X, false );
  BitmapSubcomplex & XA = * pair_complex;
  BitmapSubcomplex & A = * rel_complex;

  // Add cells from A_cubes to A and remove from X
  BOOST_FOREACH ( uint64_t index, A_cubes ) {
    std::vector <uint64_t> cube = coordinates ( index );
    std::vector < uint64_t > offset ( D );
    for ( int d = 0; d < D; ++ d )
      offset [ d ] = cube [ d ] - min_cube [ d ];
    std::vector < std::vector < uint64_t > > cells = X . fullCubeIndexes ( offset );
    for ( int d = 0; d <= D; ++ d ) {
      BOOST_FOREACH ( uint64_t cell, cells [ d ] ) {
        XA . erase ( cell, d );
        A . insert ( cell, d );
      }
    }
  }

  XA . finalize ();
  A . finalize ();

  // Don't create domain_pair, just use X, XA, and A
  RelativePair domain_pair ( &X, &XA, &A );

  // Set domain and co-domain (same as domain)
  CubicalComplex & full_domain = domain_pair . base ();
  BitmapSubcomplex & domain = domain_pair . pair ();
  BitmapSubcomplex & domain_A = domain_pair . relative ();
  CubicalComplex & full_codomain = domain_pair . base ();
  BitmapSubcomplex & codomain = domain_pair . pair ();

  // std::cout << "Domain size: " << full_domain . size () << std::endl;
  // std::cout << "(X, A) size = (" << domain . size() << ", " << domain_A . size () << ")" << std::endl;
  // std::cout << "(X, A) top size = (" << domain . size(D) << ", " << domain_A . size (D) << ")" << std::endl;

  /// Generate domain_grid_elements_X and domain_grid_elements_A tables.
  /// These are for determining which top cells are involved in constructing
  /// a fiber (and which are relative).
  std::vector < std::vector < boost::unordered_set < Index > > > domain_grid_elements_X ( D + 1 );
  std::vector < std::vector < boost::unordered_set < Index > > > domain_grid_elements_A ( D + 1 );

  for ( int d = 0; d <= D; ++ d ) {
    domain_grid_elements_X [ d ] . resize ( full_domain . size ( d ) );
    domain_grid_elements_A [ d ] . resize ( full_domain . size ( d ) );
  }

  uint64_t N = full_domain . size ( D );
  for ( uint64_t i = 0; i < N; ++ i ) {
    std::vector < boost::unordered_set < uint64_t > > closure_of_i ( D + 1 );
    closure_of_i [ D ] . insert ( i );
    closure ( closure_of_i, full_domain );
    for ( int d = 0; d <= D; ++ d ) {
      BOOST_FOREACH ( uint64_t j, closure_of_i [ d ] ) {
        domain_grid_elements_X [ d ] [ j ] . insert ( i );
      }
    }
  }

  N = domain_A . size ( D );
  for ( uint64_t k = 0; k < N; ++ k ) {
    uint64_t i = domain_A . indexToCell ( k, D );
    std::vector < boost::unordered_set < uint64_t > > closure_of_i ( D + 1 );
    closure_of_i [ D ] . insert ( i );
    closure ( closure_of_i, full_domain );
    for ( int d = 0; d <= D; ++ d ) {
      BOOST_FOREACH ( uint64_t j, closure_of_i [ d ] ) {
        domain_grid_elements_A [ d ] [ j ] . insert ( i );
      }
    }
  }

  // Compute the homology generators
  int cutoff_dimension = full_domain . dimension ();
  // Computing Domain Generators with MorseGenerators
  Generators_t domain_gen = MorseGenerators ( domain, cutoff_dimension );
  // Computing Morse Complex of codomain (codomain_morse)
  MorseComplex codomain_morse ( codomain );
  // Computing Generators of codomain_morse
  Generators_t codomain_morse_gen = SmithGenerators ( codomain_morse, cutoff_dimension );

  // Compute the cycles projected into the codomain through the graph
  std::vector < std::vector < Chain > > codomain_cycles ( D + 1 );

  // Main Loop: Loop through each domain generator
  bool acyclic_map = true;
  for ( int d = 0; d <= D; ++ d ) {
    if ( not acyclic_map ) break;
    int number_of_domain_gen = domain_gen [ d ] . size ();
    codomain_cycles [ d ] . resize ( number_of_domain_gen );
    // std::cout << "Number of generators at dimension " << d << ": " << number_of_domain_gen << std::endl;

    for ( int gi = 0; gi < number_of_domain_gen; ++ gi ) {
      if ( not acyclic_map ) break;
      // std::cout << "Dimension = " << d << ", generator = " << gi << std::endl;

      // Lift the domain generator into the Relative Graph Complex
      Chain answer;
      answer . dimension () = d;
      std::vector < boost::unordered_map < Index, Chain > > graph_boundary ( D + 1 );
      // Begin by Chain Lift
      // std::cout << "Chain Lift" << std::endl;
      // First, get the domain generator and include it into full_domain
      Chain domain_generator = domain . include ( domain_gen [ d ] [ gi ] . first );
      // std::cout << "Chain size = " << domain_generator () . size () << std::endl;

      // Now lift it
      BOOST_FOREACH ( const Term & t, domain_generator () ) {
        if ( not acyclic_map ) break;
        boost::unordered_set < Index > X_nbs, A_nbs;
        X_nbs = domain_grid_elements_X [ d ] [ t . index () ];
        A_nbs = domain_grid_elements_A [ d ] [ t . index () ];
        FiberComplex fiber ( X_nbs, A_nbs, full_domain, full_codomain, F );

        if ( fiber . size () == 0 ) {
          std::cout << "UNEXPECTED: CANT LIFT CHAIN DUE TO EMPTY FIBER." << std::endl;
          acyclic_map = false;
          break;
        }

        if ( acyclic_check && ( not fiber . acyclic () ) ) {
          acyclic_map = false;
          break;
        }

        // Choose a vertex in fiber
        Index choice = fiber . choose ();
        if ( d == 0 ) answer += Term ( choice, t . coef () );
        Chain bd = full_domain . boundary ( t . index (), d );
        BOOST_FOREACH ( const Term & s, bd () ) {
          graph_boundary [d-1] [ s . index () ] . dimension () = 0;
          graph_boundary [d-1] [ s . index () ] += Term ( choice, t . coef () * s . coef () );
        }
      }
      // std::cout << "Cycle Lift -- iterative preboundary loop begins now" << std::endl;

      // Now process every fiber we can see in graph_boundary
      // Loop through every fiber dimension (fd)
      for ( int fd = d - 1; fd >= 0; -- fd ) {
        if ( not acyclic_map ) break;
        // std::cout << "Fiber dimension = " << fd << std::endl;

        // Loop through every non-empty fiber
        typedef std::pair < Index, Chain > IndexChainPair;
        BOOST_FOREACH ( const IndexChainPair & fiberchain, graph_boundary [fd] ) {
          if ( not acyclic_map ) break;
          // Determine fiber
          boost::unordered_set < Index > X_nbs, A_nbs;
          X_nbs = domain_grid_elements_X [ fd ] [ fiberchain . first ];
          A_nbs = domain_grid_elements_A [ fd ] [ fiberchain . first ];
          FiberComplex fiber ( X_nbs, A_nbs, full_domain, full_codomain, F );

          if ( acyclic_check && ( not fiber . acyclic_or_trivial () ) ) {
            acyclic_map = false;
            break;
          }

          // Determine chain in fiber
          Chain projected = fiber . project ( fiberchain . second );
          // Work out the preboundary
          Chain preboundary = fiber . preboundary ( projected );
          // Include the preboundary back into the codomain
          Chain included_preboundary = fiber . include ( preboundary );

          // If this is a zero-dim fiber, add to answer
          if ( fd == 0 ) answer -= included_preboundary;
          // Add preboundary to adjacent fibers
          Chain bd = full_domain . boundary (fiberchain . first, fd );
          Ring grade;
          if ( fd % 2 ) grade = Ring ( -1 ); else grade = Ring ( 1 );
          BOOST_FOREACH ( const Term & s, bd () ) {
            graph_boundary [fd-1] [ s . index () ] . dimension () = d - fd;
            graph_boundary [fd-1] [ s . index () ] -= included_preboundary * (grade * s . coef ());
          }
        }
      }

      // First project into the relative codomain complex
      Chain relative_answer = simplify ( codomain . project ( answer ) );

      // Project the Relative Graph Cycle to the Codomain
      codomain_cycles [ d ] [ gi ] = codomain_morse . lower ( relative_answer );
    }
  }  // Main Loop: for ( int d = 0; d <= D; ++ d )

  if ( not acyclic_map ) {
    std::cout << "Returning no answer due to lack of acyclicity in map." << std::endl;
    return 1;
  }

  // Loop through each Codomain Cycle
  for ( int d = 0; d <= codomain_morse . dimension (); ++ d ) {
    Matrix G = chainsToMatrix ( codomain_morse_gen [ d ], codomain_morse, d );
    Matrix Z = chainsToMatrix ( codomain_cycles [ d ], codomain_morse, d );
    Matrix MapHom = SmithSolve (G, Z);
    // if ( MapHom . number_of_rows () != 0 ) {
    //   std::cout << "Dimension " << d << std::endl;
    //   print_matrix (MapHom);
    // }
    output -> push_back ( MapHom );
  }
  for ( int d = codomain_morse . dimension () + 1; d <= D; ++ d ) {
    output -> push_back ( Matrix (0, 0) );
  }

  return 0;
}

} // namespace chomp

#endif

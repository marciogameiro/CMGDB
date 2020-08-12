#ifndef CMDB_JOIN_H
#define CMDB_JOIN_H

#include <memory>
#include <exception>
#include "TreeGrid.h"
#include "Atlas.h"

template < class GridType, class InputIterator>
void join ( std::shared_ptr<GridType> output, 
	          InputIterator start, 
	          InputIterator stop );

// Note: partial template specialization of a function
//       requires a trick, since the language only provides
//       for partial specialization of a class template

template < class GridType, class InputIterator>
struct joinImpl { 
	static void act( std::shared_ptr<GridType> output, 
	    						 InputIterator start, 
	    						 InputIterator stop ) { 
		// Dynamic Dispatch
		if ( std::shared_ptr<TreeGrid> ptr = 
			   std::dynamic_pointer_cast<TreeGrid> ( output ) ) {
			return joinImpl<TreeGrid,InputIterator>::act ( ptr, start, stop );
		}
		if ( std::shared_ptr<Atlas> ptr = 
			   std::dynamic_pointer_cast<Atlas> ( output ) ) {
			return joinImpl<Atlas,InputIterator>::act ( ptr, start, stop );
		}
		throw std::logic_error ( "Error: joinImpl specialization not "
			                       " written for this Grid class.\n" );
	} 
};

template < class GridType, class InputIterator>
void join ( std::shared_ptr<GridType> output, 
	          InputIterator start, 
	          InputIterator stop ) {
	return joinImpl<GridType,InputIterator>::act ( output, start, stop );
}

// joinImpl specializations:

template < class InputIterator >
struct joinImpl < TreeGrid, InputIterator > { 
	static void act ( std::shared_ptr<TreeGrid> output, 
	    				  		InputIterator start, 
	    					  	InputIterator stop ) { 
		std::shared_ptr<CompressedTreeGrid> joinup 
			( TreeGrid::join ( start, stop ) );
  	output -> assign ( joinup );
	} 
};

template < class InputIterator >
struct joinImpl < Atlas, InputIterator > { 
	static void act ( std::shared_ptr<Atlas> output, 
	    						  InputIterator start, 
	    						  InputIterator stop ) { 

		//std::cout << "Atlas join.\n";
		output -> clear ();
		if ( start == stop ) return;
		std::shared_ptr<Atlas> start_ptr = std::dynamic_pointer_cast<Atlas> ( *start );
		// assert ( start_ptr );
		// Note: It appears we must be joining Atlases with same chart structure
		std::vector < Atlas::size_type > chart_ids;
		for ( Atlas::IdChartPair const& pair : start_ptr -> charts () ) {
			chart_ids . push_back ( pair . first );
		}

		for ( Atlas::size_type chart_id : chart_ids ) {
			//std::cout << "Atlas join, top of loop, chart_id=" << chart_id << ".\n";

			std::vector<std::shared_ptr<TreeGrid> > charts;
			//int atlas_debug_count = 0;
			for ( InputIterator it = start; it != stop; ++ it ) {
				//std::cout << "Looping through Atlas " << atlas_debug_count++ << ".\n";
				std::shared_ptr<Atlas> it_ptr = 
					std::dynamic_pointer_cast<Atlas> ( *it );
				if ( not it_ptr ) {
					throw std::logic_error ( "Atlas::join error: not looping through a container "
			                       			 " of std::shared_ptr<Atlas>.\n" );
				}
				charts . push_back ( it_ptr -> chart ( chart_id ) );
				if ( not std::dynamic_pointer_cast<TreeGrid> ( it_ptr -> chart ( chart_id ) ) ) {
					throw std::logic_error ( "Atlas::join error: just pushed back"
			                       			 " a nonconformant chart.\n" );
				}
			}
			//std::cout << "Atlas join about to reset\n";
			output -> chart ( chart_id ) . reset ( new PointerGrid );
			//std::cout << "Atlas join about to recurse on join\n";
			//std::cout << "Giving it " << charts . size () << " charts to join.\n";
			join ( output -> chart ( chart_id ), charts . begin (), charts . end () );
			//std::cout << "Atlas join returned from join recursion\n";
		}
		output -> finalize ();
	} 
};
#endif

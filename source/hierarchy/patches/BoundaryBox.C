//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-4-4/source/hierarchy/patches/BoundaryBox.C $
// Package:	SAMRAI hierarchy
// Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 1917 $
// Modified:	$LastChangedDate: 2008-01-25 13:28:01 -0800 (Fri, 25 Jan 2008) $
// Description:	BoundaryBox representing a portion of the physical boundary
//

#ifndef included_hier_BoundaryBox_C
#define included_hier_BoundaryBox_C

#include "BoundaryBox.h"

#ifdef DEBUG_CHECK_ASSERTIONS
#include "tbox/Utilities.h"
#endif

#include "BoundaryLookupTable.h"

#ifdef DEBUG_NO_INLINE
#include "BoundaryBox.I"
#endif


namespace SAMRAI {
   namespace hier {


template<int DIM>  BoundaryBox<DIM>::BoundaryBox(const Box<DIM>& box,
                                     const int bdry_type,
                                     const int location_index)
:  d_box(box)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   BoundaryLookupTable<DIM>* blut =
      BoundaryLookupTable<DIM>::getLookupTable();


   const tbox::Array<int>& location_index_max = blut->getMaxLocationIndices();

   TBOX_ASSERT( (bdry_type >= 1) && (bdry_type <= DIM) );

   TBOX_ASSERT(location_index >= 0);
   TBOX_ASSERT(location_index < location_index_max[bdry_type-1]);
#endif

   d_bdry_type = bdry_type;

   d_location_index = location_index;

   d_is_mblk_singularity = false;
}

template<int DIM>  BoundaryBox<DIM>::~BoundaryBox()
{
}

}
}

#endif

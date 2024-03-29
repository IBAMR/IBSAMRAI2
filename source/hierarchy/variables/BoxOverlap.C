//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-4-4/source/hierarchy/variables/BoxOverlap.C $
// Package:	SAMRAI hierarchy
// Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 1917 $
// Modified:	$LastChangedDate: 2008-01-25 13:28:01 -0800 (Fri, 25 Jan 2008) $
// Description:	Base class that describes intersections between AMR boxes
//

#ifndef included_hier_BoxOverlap_C
#define included_hier_BoxOverlap_C

#include "BoxOverlap.h"

#ifdef DEBUG_NO_INLINE
#include "BoxOverlap.I"
#endif
namespace SAMRAI {
    namespace hier {

template<int DIM>  BoxOverlap<DIM>::~BoxOverlap()
{
}

template<int DIM> void BoxOverlap<DIM>::print(std::ostream& os) const
{
   os << "print() method not implemented for this overlap type" << std::endl;
}

}
}
#endif

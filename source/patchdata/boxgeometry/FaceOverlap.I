//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-4-4/source/patchdata/boxgeometry/FaceOverlap.I $
// Package:	SAMRAI patch data geometry
// Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 1917 $
// Modified:	$LastChangedDate: 2008-01-25 13:28:01 -0800 (Fri, 25 Jan 2008) $
// Description:	hier::Box intersection information for face centered objects
//

#ifdef DEBUG_CHECK_ASSERTIONS
#include "tbox/Utilities.h"
#endif

#ifdef DEBUG_NO_INLINE
#define inline
#endif
namespace SAMRAI {
    namespace pdat {

template<int DIM> inline
const hier::BoxList<DIM>&
FaceOverlap<DIM>::getDestinationBoxList(const int axis) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT((axis >= 0) && (axis < DIM));
#endif
   return(d_dst_boxes[axis]);
}

template<int DIM> inline
const hier::IntVector<DIM>& FaceOverlap<DIM>::getSourceOffset() const
{
   return(d_offset);
}

#ifdef DEBUG_NO_INLINE
#undef inline
#endif
}
}

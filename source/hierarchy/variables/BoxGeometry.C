//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-4-4/source/hierarchy/variables/BoxGeometry.C $
// Package:	SAMRAI hierarchy
// Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 1917 $
// Modified:	$LastChangedDate: 2008-01-25 13:28:01 -0800 (Fri, 25 Jan 2008) $
// Description:	Box geometry description for overlap computations
//

#ifndef included_hier_BoxGeometry_C
#define included_hier_BoxGeometry_C

#include "BoxGeometry.h"
#include "tbox/Utilities.h"

#ifdef DEBUG_NO_INLINE
#include "BoxGeometry.I"
#endif

namespace SAMRAI {
   namespace hier {


template<int DIM>  BoxGeometry<DIM>::~BoxGeometry()
{
}

template<int DIM> tbox::Pointer< BoxOverlap<DIM> >
BoxGeometry<DIM>::restrictOverlapToOwnedData(
   const tbox::Pointer< BoxOverlap<DIM> >& overlap,
   const Box<DIM>& src_box,
   const BoxList<DIM>& level_boxes) const
{
   NULL_USE(src_box);
   NULL_USE(level_boxes);
   return(overlap);
}

template<int DIM> tbox::Pointer< BoxOverlap<DIM> >
BoxGeometry<DIM>::removeOverlapOnBoxes(
   const tbox::Pointer< BoxOverlap<DIM> >& overlap,
   const BoxList<DIM>& boxes) const
{
   NULL_USE(boxes);
   return(overlap);
}

template<int DIM> void BoxGeometry<DIM>::computeOwnedDataBoxes(
   BoxList<DIM>& owned_boxes,
   const Box<DIM>& box,
   const BoxList<DIM>& level_boxes,
   const tbox::Array< IntVector<DIM> >& offsets)
{
   /*
    * The datum at a given offset from a cell of the box is owned by the
    * box unless a cell of the level touches it at an earlier offset.
    */
   for (int k = 0; k < offsets.getSize(); k++) {
      BoxList<DIM> boxes(Box<DIM>::shift(box, offsets[k]));
      for (int j = 0; j < k; j++) {
         BoxList<DIM> taken(level_boxes);
         taken.shift(offsets[j]);
         boxes.removeIntersections(taken);
      }
      for (typename BoxList<DIM>::Iterator b(boxes); b; b++) {
         owned_boxes.appendItem(b());
      }
   }
   owned_boxes.coalesceBoxes();
}

template<int DIM> void BoxGeometry<DIM>::intersectOverlapBoxes(
   BoxList<DIM>& result_boxes,
   const BoxList<DIM>& overlap_boxes,
   const BoxList<DIM>& owned_boxes)
{
   for (typename BoxList<DIM>::Iterator b(overlap_boxes); b; b++) {
      BoxList<DIM> boxes(b());
      boxes.intersectBoxes(owned_boxes);
      boxes.coalesceBoxes();
      for (typename BoxList<DIM>::Iterator r(boxes); r; r++) {
         result_boxes.appendItem(r());
      }
   }
}

}
}
#endif

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
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(offsets.getSize() > 0 && offsets[0] == IntVector<DIM>(0));
#endif

   /*
    * The datum at a given offset from a cell of the box is owned by the
    * box unless a cell of the level touches it at an earlier offset.  The
    * first offset is zero, so the box owns every datum that has the index
    * of one of its cells.
    */
   owned_boxes.appendItem(box);

   /*
    * The other data are on the upper border of the box, and the highest
    * cell of the box that touches one of them is the one below it in the
    * directions in which it is on the border.  So the data whose highest
    * such cell is at a given offset form one box, and the box owns those
    * that no level cell outside it touches at an earlier offset.  Only
    * level cells next to the box can do that.
    */
   const Box<DIM> nearby_region(Box<DIM>::grow(box, IntVector<DIM>(1)));
   BoxList<DIM> outside_boxes;
   for (typename BoxList<DIM>::Iterator b(level_boxes); b; b++) {
      const Box<DIM> nearby_box(b() * nearby_region);
      if (nearby_box.empty() || box.contains(nearby_box)) {
         continue;
      }
      if ((nearby_box * box).empty()) {
         outside_boxes.appendItem(nearby_box);
      } else {
         BoxList<DIM> pieces(nearby_box);
         pieces.removeIntersections(box);
         for (typename BoxList<DIM>::Iterator p(pieces); p; p++) {
            outside_boxes.appendItem(p());
         }
      }
   }

   for (int k = 1; k < offsets.getSize(); k++) {
      Box<DIM> border_box(box);
      for (int d = 0; d < DIM; d++) {
#ifdef DEBUG_CHECK_ASSERTIONS
         TBOX_ASSERT(offsets[k](d) == 0 || offsets[k](d) == 1);
#endif
         if (offsets[k](d) != 0) {
            border_box.lower(d) = border_box.upper(d) = box.upper(d) + 1;
         }
      }
      BoxList<DIM> boxes(border_box);
      for (int j = 0; j < k && !boxes.isEmpty(); j++) {
         for (typename BoxList<DIM>::Iterator b(outside_boxes); b; b++) {
            const Box<DIM> taken(Box<DIM>::shift(b(), offsets[j]));
            if (taken.intersects(border_box)) {
               boxes.removeIntersections(taken);
            }
         }
      }
      for (typename BoxList<DIM>::Iterator b(boxes); b; b++) {
         owned_boxes.appendItem(b());
      }
   }
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

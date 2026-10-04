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

template<int DIM> void
BoxGeometry<DIM>::computeOwnedBorderData(
   tbox::Array< BoxList<DIM> >& owned_border_data,
   const BoxList<DIM>& level_boxes,
   const Box<DIM>& owner_box) const
{
   NULL_USE(level_boxes);
   NULL_USE(owner_box);
   owned_border_data.resizeArray(0);
}

template<int DIM> tbox::Pointer< BoxOverlap<DIM> >
BoxGeometry<DIM>::restrictOverlapToOwnedData(
   const tbox::Pointer< BoxOverlap<DIM> >& overlap,
   const Box<DIM>& src_box,
   const tbox::Array< BoxList<DIM> >& owned_border_data) const
{
   NULL_USE(src_box);
   NULL_USE(owned_border_data);
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

template<int DIM> void BoxGeometry<DIM>::computeOwnedBorderBoxes(
   tbox::Array< BoxList<DIM> >& owned_border_boxes,
   const int first,
   const BoxList<DIM>& level_boxes,
   const Box<DIM>& owner_box,
   const IntVector<DIM>* offsets,
   const int num_offsets)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(num_offsets > 0 && offsets[0] == IntVector<DIM>(0));
   TBOX_ASSERT(owned_border_boxes.getSize() >= first + num_offsets - 1);
#endif

   /*
    * The datum at a given offset from a level cell is owned by that cell
    * unless a level cell touches it at an earlier offset.  The first
    * offset is zero, so every level cell owns the datum with its own
    * index, and what is left for the other offsets is on the upper border
    * of the level boxes.
    */
   for (int k = 1; k < num_offsets; k++) {
      BoxList<DIM>& boxes = owned_border_boxes[first + k - 1];
      boxes = level_boxes;
      boxes.intersectBoxes(owner_box);
      boxes.shift(offsets[k]);
      for (int j = 0; j < k && !boxes.isEmpty(); j++) {
         if (j == 0) {
            boxes.removeIntersections(level_boxes);
         } else {
            BoxList<DIM> taken(level_boxes);
            taken.shift(offsets[j]);
            boxes.removeIntersections(taken);
         }
      }
   }
}

template<int DIM> void BoxGeometry<DIM>::intersectOverlapBoxes(
   BoxList<DIM>& result_boxes,
   const BoxList<DIM>& overlap_boxes,
   const Box<DIM>& src_box,
   const tbox::Array< BoxList<DIM> >& owned_border_boxes,
   const int first,
   const IntVector<DIM>* offsets,
   const int num_offsets)
{
   for (typename BoxList<DIM>::Iterator b(overlap_boxes); b; b++) {
      /*
       * Usually only the data that have the index of a source cell are
       * owned, and they form one box.  Otherwise collect the pieces and
       * merge those that came from this overlap box.
       */
      const Box<DIM> own_index_box(b() * src_box);
      BoxList<DIM> border_pieces;
      for (int k = 1; k < num_offsets; k++) {
         const BoxList<DIM>& border_boxes = owned_border_boxes[first + k - 1];
         if (border_boxes.isEmpty()) {
            continue;
         }
         const Box<DIM> border_box(b() * Box<DIM>::shift(src_box, offsets[k]));
         if (border_box.empty()) {
            continue;
         }
         for (typename BoxList<DIM>::Iterator o(border_boxes); o; o++) {
            const Box<DIM> owned_box(border_box * o());
            if (!owned_box.empty()) {
               border_pieces.appendItem(owned_box);
            }
         }
      }
      if (border_pieces.isEmpty()) {
         if (!own_index_box.empty()) {
            result_boxes.appendItem(own_index_box);
         }
      } else {
         if (!own_index_box.empty()) {
            border_pieces.appendItem(own_index_box);
            border_pieces.coalesceBoxes();
         }
         for (typename BoxList<DIM>::Iterator r(border_pieces); r; r++) {
            result_boxes.appendItem(r());
         }
      }
   }
}

}
}
#endif

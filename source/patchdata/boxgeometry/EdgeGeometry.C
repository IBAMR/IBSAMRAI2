//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-4-4/source/patchdata/boxgeometry/EdgeGeometry.C $
// Package:	SAMRAI patch data geometry
// Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 1917 $
// Modified:	$LastChangedDate: 2008-01-25 13:28:01 -0800 (Fri, 25 Jan 2008) $
// Description:	hier::Box geometry information for edge centered objects
//

#ifndef included_pdat_EdgeGeometry_C
#define included_pdat_EdgeGeometry_C

#include "EdgeGeometry.h"
#include "BoxList.h"
#include "EdgeOverlap.h"

#ifdef DEBUG_CHECK_ASSERTIONS
#include "tbox/Utilities.h"
#endif

#ifdef DEBUG_NO_INLINE
#include "EdgeGeometry.I"
#endif
namespace SAMRAI {
    namespace pdat {

/*
*************************************************************************
*									*
* Create a edge geometry object given the box and ghost cell width.	*
*									*
*************************************************************************
*/

template<int DIM>  
EdgeGeometry<DIM>::EdgeGeometry(
   const hier::Box<DIM>& box, 
   const hier::IntVector<DIM>& ghosts)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(ghosts.min() >= 0);
#endif
   d_box    = box;
   d_ghosts = ghosts;
}

template<int DIM>  EdgeGeometry<DIM>::~EdgeGeometry()
{
}

/*
*************************************************************************
*									*
* Attempt to calculate the intersection between two edge centered box	*
* geometries.  The calculateOverlap() checks whether both arguments are	*
* edge geometries; if so, it compuates the intersection.  If not, then	*
* it calls calculateOverlap() on the source object (if retry is true)	*
* to allow the source a chance to calculate the intersection.  See the	*
* hier::BoxGeometry<DIM> base class for more information about the protocol.	*
* A pointer to null is returned if the intersection cannot be computed.	*
* 									*
*************************************************************************
*/

template<int DIM> 
tbox::Pointer< hier::BoxOverlap<DIM> > EdgeGeometry<DIM>::calculateOverlap(
   const hier::BoxGeometry<DIM>& dst_geometry,
   const hier::BoxGeometry<DIM>& src_geometry,
   const hier::Box<DIM>& src_mask,
   const bool overwrite_interior,
   const hier::IntVector<DIM>& src_offset,
   const bool retry) const
{
   const EdgeGeometry<DIM> *t_dst = 
      dynamic_cast<const EdgeGeometry<DIM> *>(&dst_geometry);
   const EdgeGeometry<DIM> *t_src =
      dynamic_cast<const EdgeGeometry<DIM> *>(&src_geometry);

   tbox::Pointer< hier::BoxOverlap<DIM> > over = NULL;

   if ((t_src != NULL) && (t_dst != NULL)) {
      over = doOverlap(*t_dst, *t_src, src_mask, overwrite_interior, 
		       src_offset);
   } else if (retry) {
      over = src_geometry.calculateOverlap(dst_geometry, src_geometry,
                                           src_mask, overwrite_interior,
                                           src_offset, false);
   }
   return(over);
}

/*
*************************************************************************
*									*
* Convert an AMR-index space hier::Box into a edge-index space box by a	*
* cyclic shift of indices.						*
*									*
*************************************************************************
*/

template<int DIM> 
hier::Box<DIM> EdgeGeometry<DIM>::toEdgeBox(
   const hier::Box<DIM>& box, 
   int axis)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(0 <= axis && axis < DIM);
#endif

   hier::Box<DIM> edge_box;

   if (!box.empty()) {
      edge_box = box;
      for (int i = 0; i < DIM; i++) {
         if (axis != i) {
            edge_box.upper(i) += 1;
         }
      }
   }

   return(edge_box);
}

/*
*************************************************************************
*									*
* Compute the overlap between two edge centered boxes.  The algorithm	*
* is fairly straight-forward.  First, we perform a quick-and-dirty	*
* intersection to see if the boxes might overlap.  If that intersection	*
* is not empty, then we need to do a better job calculating the overlap	*
* for each dimension.  Note that the AMR index space boxes must be	*
* shifted into the edge centered space before we calculate the proper	*
* intersections.							*
*									*
*************************************************************************
*/

template<int DIM> 
tbox::Pointer< hier::BoxOverlap<DIM> > EdgeGeometry<DIM>::doOverlap(
   const EdgeGeometry<DIM>& dst_geometry,
   const EdgeGeometry<DIM>& src_geometry,
   const hier::Box<DIM>& src_mask,
   const bool overwrite_interior,
   const hier::IntVector<DIM>& src_offset)
{
   hier::BoxList<DIM> dst_boxes[DIM];

   // Perform a quick-and-dirty intersection to see if the boxes might overlap

   const hier::Box<DIM> src_box =
      hier::Box<DIM>::grow(src_geometry.d_box, src_geometry.d_ghosts) * src_mask;
   const hier::Box<DIM> src_shift =
      hier::Box<DIM>::shift(src_box, src_offset);
   const hier::Box<DIM> dst_ghost =
      hier::Box<DIM>::grow(dst_geometry.d_box, dst_geometry.d_ghosts);

   // Compute the intersection (if any) for each of the edge directions

   const hier::Box<DIM> quick_check =
      hier::Box<DIM>::grow(src_shift, 1) * hier::Box<DIM>::grow(dst_ghost, 1);

   if (!quick_check.empty()) {

      for (int d = 0; d < DIM; d++) {

         const hier::Box<DIM> dst_edge = toEdgeBox(dst_ghost, d);
         const hier::Box<DIM> src_edge = toEdgeBox(src_shift, d);
         const hier::Box<DIM> together = dst_edge * src_edge;

         if (!together.empty()) {

            if (!overwrite_interior) {
               const hier::Box<DIM> int_edge = toEdgeBox(dst_geometry.d_box, d);
               dst_boxes[d].removeIntersections(together,int_edge);
            } else {
               dst_boxes[d].appendItem(together);
            }

         }  // if (!together.empty())

      }  // loop over dim

   }  // if (!quick_check.empty())

   // Create the edge overlap data object using the boxes and source shift

   hier::BoxOverlap<DIM> *overlap = new EdgeOverlap<DIM>(dst_boxes, src_offset);
   return(tbox::Pointer< hier::BoxOverlap<DIM> >(overlap));
}

/*
*************************************************************************
*                                                                       *
* Restrict an overlap to the data owned by the source box.  An edge     *
* in direction d is touched by the cell with the same index and by the  *
* cells below it in any combination of the other directions.            *
*                                                                       *
*************************************************************************
*/

template<int DIM> void
EdgeGeometry<DIM>::computeOwnedBorderData(
   tbox::Array< hier::BoxList<DIM> >& owned_border_data,
   const hier::BoxList<DIM>& level_boxes,
   const hier::Box<DIM>& owner_box) const
{
   hier::BoxList<DIM> boxes(level_boxes);
   boxes.coalesceBoxes();

   const int num_offsets = 1 << (DIM - 1);
   owned_border_data.resizeArray(DIM * (num_offsets - 1));
   hier::IntVector<DIM> offsets[1 << (DIM - 1)];
   for (int d = 0; d < DIM; d++) {
      /*
       * Bit i of k is the offset in direction i, so that the offsets are
       * sorted with the last coordinate compared first.
       */
      int n = 0;
      for (int k = 0; k < (1 << DIM); k++) {
         if ((k >> d) & 1) continue;
         for (int i = 0; i < DIM; i++) {
            offsets[n](i) = (k >> i) & 1;
         }
         n++;
      }
      hier::BoxGeometry<DIM>::computeOwnedBorderBoxes(
         owned_border_data, d * (num_offsets - 1), boxes, owner_box, offsets, num_offsets);
   }
}

template<int DIM> tbox::Pointer< hier::BoxOverlap<DIM> >
EdgeGeometry<DIM>::restrictOverlapToOwnedData(
   const tbox::Pointer< hier::BoxOverlap<DIM> >& overlap,
   const hier::Box<DIM>& src_box,
   const tbox::Array< hier::BoxList<DIM> >& owned_border_data) const
{
   const EdgeOverlap<DIM>* t_overlap =
      dynamic_cast<const EdgeOverlap<DIM>*>(overlap.getPointer());
   const int num_offsets = 1 << (DIM - 1);
   if (t_overlap == NULL ||
       owned_border_data.getSize() != DIM * (num_offsets - 1)) {
      return(overlap);
   }

   hier::BoxList<DIM> dst_boxes[DIM];
   hier::IntVector<DIM> offsets[1 << (DIM - 1)];
   for (int d = 0; d < DIM; d++) {
      /*
       * Bit i of k is the offset in direction i, so that the offsets are
       * sorted with the last coordinate compared first.
       */
      int n = 0;
      for (int k = 0; k < (1 << DIM); k++) {
         if ((k >> d) & 1) continue;
         for (int i = 0; i < DIM; i++) {
            offsets[n](i) = (k >> i) & 1;
         }
         n++;
      }
      hier::BoxGeometry<DIM>::intersectOverlapBoxes(
         dst_boxes[d], t_overlap->getDestinationBoxList(d), src_box,
         owned_border_data, d * (num_offsets - 1), offsets, num_offsets);
   }

   return(new EdgeOverlap<DIM>(dst_boxes, t_overlap->getSourceOffset()));
}

/*
*************************************************************************
*                                                                       *
* Compute the overlap between two edge centered boxes as doOverlap()    *
* does, but keep only the data owned by the source box, as              *
* restrictOverlapToOwnedData() does.                                    *
*                                                                       *
*************************************************************************
*/

template<int DIM> tbox::Pointer< hier::BoxOverlap<DIM> >
EdgeGeometry<DIM>::calculateOwnedOverlap(
   const hier::BoxGeometry<DIM>& src_geometry,
   const hier::Box<DIM>& src_mask,
   const bool overwrite_interior,
   const hier::IntVector<DIM>& src_offset,
   const hier::Box<DIM>& src_box,
   const tbox::Array< hier::BoxList<DIM> >& owned_border_data) const
{
   const EdgeGeometry<DIM> *t_src =
      dynamic_cast<const EdgeGeometry<DIM> *>(&src_geometry);
   const int num_offsets = 1 << (DIM - 1);
   if (t_src == NULL ||
       owned_border_data.getSize() != DIM * (num_offsets - 1)) {
      return(hier::BoxGeometry<DIM>::calculateOwnedOverlap(
                src_geometry, src_mask, overwrite_interior, src_offset,
                src_box, owned_border_data));
   }

   hier::BoxList<DIM> dst_boxes[DIM];

   const hier::Box<DIM> src_ghost =
      hier::Box<DIM>::grow(t_src->d_box, t_src->d_ghosts) * src_mask;
   const hier::Box<DIM> src_shift =
      hier::Box<DIM>::shift(src_ghost, src_offset);
   const hier::Box<DIM> dst_ghost =
      hier::Box<DIM>::grow(d_box, d_ghosts);

   const hier::Box<DIM> quick_check =
      hier::Box<DIM>::grow(src_shift, 1) * hier::Box<DIM>::grow(dst_ghost, 1);

   if (!quick_check.empty()) {
      hier::IntVector<DIM> offsets[1 << (DIM - 1)];
      for (int d = 0; d < DIM; d++) {
         /*
          * Bit i of k is the offset in direction i, so that the offsets
          * are sorted with the last coordinate compared first.
          */
         int n = 0;
         for (int k = 0; k < (1 << DIM); k++) {
            if ((k >> d) & 1) continue;
            for (int i = 0; i < DIM; i++) {
               offsets[n](i) = (k >> i) & 1;
            }
            n++;
         }
         const hier::Box<DIM> together =
            toEdgeBox(dst_ghost, d) * toEdgeBox(src_shift, d);
         if (!together.empty()) {
            if (!overwrite_interior) {
               hier::BoxList<DIM> boxes;
               boxes.removeIntersections(together, toEdgeBox(d_box, d));
               hier::BoxGeometry<DIM>::intersectOverlapBoxes(
                  dst_boxes[d], boxes, src_box,
                  owned_border_data, d * (num_offsets - 1), offsets, num_offsets);
            } else {
               hier::BoxGeometry<DIM>::intersectOverlapBox(
                  dst_boxes[d], together, src_box,
                  owned_border_data, d * (num_offsets - 1), offsets, num_offsets);
            }
         }
      }
   }

   return(new EdgeOverlap<DIM>(dst_boxes, src_offset));
}

/*
*************************************************************************
*                                                                       *
* Remove from an overlap the data touched by the cells of some boxes.   *
*                                                                       *
*************************************************************************
*/

template<int DIM> tbox::Pointer< hier::BoxOverlap<DIM> >
EdgeGeometry<DIM>::removeOverlapOnBoxes(
   const tbox::Pointer< hier::BoxOverlap<DIM> >& overlap,
   const hier::BoxList<DIM>& boxes) const
{
   const EdgeOverlap<DIM>* t_overlap =
      dynamic_cast<const EdgeOverlap<DIM>*>(overlap.getPointer());
   if (t_overlap == NULL) {
      return(overlap);
   }

   hier::BoxList<DIM> dst_boxes[DIM];
   for (int d = 0; d < DIM; d++) {
      dst_boxes[d] = t_overlap->getDestinationBoxList(d);
      for (typename hier::BoxList<DIM>::Iterator b(boxes); b; b++) {
         dst_boxes[d].removeIntersections(toEdgeBox(b(), d));
      }
   }

   return(new EdgeOverlap<DIM>(dst_boxes, t_overlap->getSourceOffset()));
}

}
}
#endif

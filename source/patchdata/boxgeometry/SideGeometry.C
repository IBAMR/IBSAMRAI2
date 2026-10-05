//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-4-4/source/patchdata/boxgeometry/SideGeometry.C $
// Package:	SAMRAI patch data geometry
// Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 2856 $
// Modified:	$LastChangedDate: 2009-01-30 13:58:39 -0800 (Fri, 30 Jan 2009) $
// Description:	hier::Box geometry information for side centered objects
//

#ifndef included_pdat_SideGeometry_C
#define included_pdat_SideGeometry_C

#include "SideGeometry.h"
#include "BoxList.h"
#include "SideOverlap.h"

#ifdef DEBUG_CHECK_ASSERTIONS
#include "tbox/Utilities.h"
#endif

#ifdef DEBUG_NO_INLINE
#include "SideGeometry.I"
#endif
namespace SAMRAI {
    namespace pdat {

/*
*************************************************************************
*									*
* Create a side geometry object given the box, ghost cell width, and    *
* direction information.                                                *
*									*
*************************************************************************
*/

template<int DIM>  SideGeometry<DIM>::SideGeometry(
   const hier::Box<DIM>& box,
   const hier::IntVector<DIM>& ghosts,
   const hier::IntVector<DIM>& directions)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(ghosts.min() >= 0);
   TBOX_ASSERT(directions.min() >= 0);
#endif
   d_box    = box;
   d_ghosts = ghosts;
   d_directions = directions;
}

template<int DIM>  SideGeometry<DIM>::~SideGeometry()
{
}

/*
*************************************************************************
*									*
* Attempt to calculate the intersection between two side centered box	*
* geometries.  The calculateOverlap() checks whether both arguments are	*
* side geometries; if so, it compuates the intersection.  If not, then	*
* it calls calculateOverlap() on the source object (if retry is true)	*
* to allow the source a chance to calculate the intersection.  See the	*
* hier::BoxGeometry<DIM> base class for more information about the protocol.	*
* A pointer to null is returned if the intersection cannot be computed.	*
* 									*
*************************************************************************
*/

template<int DIM> tbox::Pointer< hier::BoxOverlap<DIM> > SideGeometry<DIM>::calculateOverlap(
   const hier::BoxGeometry<DIM>& dst_geometry,
   const hier::BoxGeometry<DIM>& src_geometry,
   const hier::Box<DIM>& src_mask,
   const bool overwrite_interior,
   const hier::IntVector<DIM>& src_offset,
   const bool retry) const
{
   const SideGeometry<DIM> *t_dst =
      dynamic_cast<const SideGeometry<DIM> *>(&dst_geometry);
   const SideGeometry<DIM> *t_src =
      dynamic_cast<const SideGeometry<DIM> *>(&src_geometry);

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
* Convert an AMR-index space hier::Box into a side-index space box by a	*
* increasing the index size by one in the axis direction.		*
*									*
*************************************************************************
*/

template<int DIM> hier::Box<DIM> SideGeometry<DIM>::toSideBox(
   const hier::Box<DIM>& box,
   int side_normal)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT( (side_normal >= 0) && (side_normal < DIM) );
#endif
   hier::Box<DIM> side_box;

   if (!box.empty()) {
      side_box = box;
      side_box.upper(side_normal) += 1;
   }

   return(side_box);
}

/*
*************************************************************************
*									*
* Compute the overlap between two side centered boxes.  The algorithm	*
* is fairly straight-forward.  First, we perform a quick-and-dirty	*
* intersection to see if the boxes might overlap.  If that intersection	*
* is not empty, then we need to do a better job calculating the overlap	*
* for each dimension.  Note that the AMR index space boxes must be	*
* shifted into the side centered space before we calculate the proper	*
* intersections.							*
*									*
*************************************************************************
*/

template<int DIM> tbox::Pointer< hier::BoxOverlap<DIM> > SideGeometry<DIM>::doOverlap(
   const SideGeometry<DIM>& dst_geometry,
   const SideGeometry<DIM>& src_geometry,
   const hier::Box<DIM>& src_mask,
   const bool overwrite_interior,
   const hier::IntVector<DIM>& src_offset)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(dst_geometry.getDirectionVector() 
           == src_geometry.getDirectionVector());
#endif

   hier::BoxList<DIM> dst_boxes[DIM];

   // Perform a quick-and-dirty intersection to see if the boxes might overlap

   const hier::Box<DIM> src_box =
      hier::Box<DIM>::grow(src_geometry.d_box, src_geometry.d_ghosts) * src_mask;
   const hier::Box<DIM> src_shift =
      hier::Box<DIM>::shift(src_box, src_offset);
   const hier::Box<DIM> dst_ghost =
      hier::Box<DIM>::grow(dst_geometry.d_box, dst_geometry.d_ghosts);

   // Compute the intersection (if any) for each of the side directions

   const hier::Box<DIM> quick_check =
      hier::Box<DIM>::grow(src_shift, 1) * hier::Box<DIM>::grow(dst_ghost, 1);

   if (!quick_check.empty()) {

      const hier::IntVector<DIM>& dirs = src_geometry.getDirectionVector();
      for (int d = 0; d < DIM; d++) {
         if ( dirs(d) ) {
            const hier::Box<DIM> dst_side = toSideBox(dst_ghost, d);
            const hier::Box<DIM> src_side = toSideBox(src_shift, d);
            const hier::Box<DIM> together = dst_side * src_side;
            if (!together.empty()) {
               if (!overwrite_interior) {
                  const hier::Box<DIM> int_side = toSideBox(dst_geometry.d_box, d);
                  dst_boxes[d].removeIntersections(together,int_side);
               } else {
                  dst_boxes[d].appendItem(together);
               }
            }  // if (!together.empty())
         } // if (dirs(d))
      }  // loop over dim && dirs(d)

   }  // if (!quick_check.empty())

   // Create the side overlap data object using the boxes and source shift

   hier::BoxOverlap<DIM> *overlap = new SideOverlap<DIM>(dst_boxes, src_offset);
   return(tbox::Pointer< hier::BoxOverlap<DIM> >(overlap));
}

/*
*************************************************************************
*                                                                       *
* Restrict an overlap to the data owned by the source box.  A side      *
* with normal direction d is touched by the cell with the same index    *
* and by the cell below it in direction d.                              *
*                                                                       *
*************************************************************************
*/

template<int DIM> void
SideGeometry<DIM>::computeOwnedBorderData(
   tbox::Array< hier::BoxList<DIM> >& owned_border_data,
   const hier::BoxList<DIM>& level_boxes,
   const hier::Box<DIM>& owner_box) const
{
   hier::BoxList<DIM> boxes(level_boxes);
   boxes.coalesceBoxes();

   owned_border_data.resizeArray(DIM);
   hier::IntVector<DIM> offsets[2];
   for (int d = 0; d < DIM; d++) {
      offsets[0] = hier::IntVector<DIM>(0);
      offsets[1] = hier::IntVector<DIM>(0);
      offsets[1](d) = 1;
      hier::BoxGeometry<DIM>::computeOwnedBorderBoxes(
         owned_border_data, d, boxes, owner_box, offsets, 2);
   }
}

template<int DIM> tbox::Pointer< hier::BoxOverlap<DIM> >
SideGeometry<DIM>::restrictOverlapToOwnedData(
   const tbox::Pointer< hier::BoxOverlap<DIM> >& overlap,
   const hier::Box<DIM>& src_box,
   const tbox::Array< hier::BoxList<DIM> >& owned_border_data) const
{
   const SideOverlap<DIM>* t_overlap =
      dynamic_cast<const SideOverlap<DIM>*>(overlap.getPointer());
   if (t_overlap == NULL || owned_border_data.getSize() != DIM) {
      return(overlap);
   }

   hier::BoxList<DIM> dst_boxes[DIM];
   hier::IntVector<DIM> offsets[2];
   for (int d = 0; d < DIM; d++) {
      if (d_directions(d)) {
         offsets[0] = hier::IntVector<DIM>(0);
         offsets[1] = hier::IntVector<DIM>(0);
         offsets[1](d) = 1;
         hier::BoxGeometry<DIM>::intersectOverlapBoxes(
            dst_boxes[d], t_overlap->getDestinationBoxList(d), src_box,
            owned_border_data, d, offsets, 2);
      }
   }

   return(new SideOverlap<DIM>(dst_boxes, t_overlap->getSourceOffset()));
}

/*
*************************************************************************
*                                                                       *
* Compute the overlap between two side centered boxes as doOverlap()    *
* does, but keep only the data owned by the source box, as              *
* restrictOverlapToOwnedData() does.                                    *
*                                                                       *
*************************************************************************
*/

template<int DIM> tbox::Pointer< hier::BoxOverlap<DIM> >
SideGeometry<DIM>::calculateOwnedOverlap(
   const hier::BoxGeometry<DIM>& src_geometry,
   const hier::Box<DIM>& src_mask,
   const bool overwrite_interior,
   const hier::IntVector<DIM>& src_offset,
   const hier::Box<DIM>& src_box,
   const tbox::Array< hier::BoxList<DIM> >& owned_border_data) const
{
   const SideGeometry<DIM> *t_src =
      dynamic_cast<const SideGeometry<DIM> *>(&src_geometry);
   if (t_src == NULL || owned_border_data.getSize() != DIM) {
      return(hier::BoxGeometry<DIM>::calculateOwnedOverlap(
                src_geometry, src_mask, overwrite_interior, src_offset,
                src_box, owned_border_data));
   }
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(d_directions == t_src->d_directions);
#endif

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
      hier::IntVector<DIM> offsets[2];
      for (int d = 0; d < DIM; d++) {
         if (d_directions(d)) {
            offsets[0] = hier::IntVector<DIM>(0);
            offsets[1] = hier::IntVector<DIM>(0);
            offsets[1](d) = 1;
            const hier::Box<DIM> together =
               toSideBox(dst_ghost, d) * toSideBox(src_shift, d);
            if (!together.empty()) {
               if (!overwrite_interior) {
                  hier::BoxList<DIM> boxes;
                  boxes.removeIntersections(together, toSideBox(d_box, d));
                  hier::BoxGeometry<DIM>::intersectOverlapBoxes(
                     dst_boxes[d], boxes, src_box,
                     owned_border_data, d, offsets, 2);
               } else {
                  hier::BoxGeometry<DIM>::intersectOverlapBox(
                     dst_boxes[d], together, src_box,
                     owned_border_data, d, offsets, 2);
               }
            }
         }
      }
   }

   return(new SideOverlap<DIM>(dst_boxes, src_offset));
}

/*
*************************************************************************
*                                                                       *
* Remove from an overlap the data touched by the cells of some boxes.   *
*                                                                       *
*************************************************************************
*/

template<int DIM> tbox::Pointer< hier::BoxOverlap<DIM> >
SideGeometry<DIM>::removeOverlapOnBoxes(
   const tbox::Pointer< hier::BoxOverlap<DIM> >& overlap,
   const hier::BoxList<DIM>& boxes) const
{
   const SideOverlap<DIM>* t_overlap =
      dynamic_cast<const SideOverlap<DIM>*>(overlap.getPointer());
   if (t_overlap == NULL) {
      return(overlap);
   }

   hier::BoxList<DIM> dst_boxes[DIM];
   for (int d = 0; d < DIM; d++) {
      dst_boxes[d] = t_overlap->getDestinationBoxList(d);
      for (typename hier::BoxList<DIM>::Iterator b(boxes); b; b++) {
         dst_boxes[d].removeIntersections(toSideBox(b(), d));
      }
   }

   return(new SideOverlap<DIM>(dst_boxes, t_overlap->getSourceOffset()));
}

}
}
#endif

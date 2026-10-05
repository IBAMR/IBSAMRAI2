//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-4-4/source/patchdata/boxgeometry/FaceGeometry.C $
// Package:	SAMRAI patch data geometry
// Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 1917 $
// Modified:	$LastChangedDate: 2008-01-25 13:28:01 -0800 (Fri, 25 Jan 2008) $
// Description:	hier::Box geometry information for face centered objects
//

#ifndef included_pdat_FaceGeometry_C
#define included_pdat_FaceGeometry_C

#include "FaceGeometry.h"
#include "BoxList.h"
#include "FaceOverlap.h"

#ifdef DEBUG_CHECK_ASSERTIONS
#include "tbox/Utilities.h"
#endif

#ifdef DEBUG_NO_INLINE
#include "FaceGeometry.I"
#endif
namespace SAMRAI {
    namespace pdat {

/*
*************************************************************************
*									*
* Create a face geometry object given the box and ghost cell width.	*
*									*
*************************************************************************
*/

template<int DIM>  FaceGeometry<DIM>::FaceGeometry(
   const hier::Box<DIM>& box, const hier::IntVector<DIM>& ghosts)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(ghosts.min() >= 0);
#endif
   d_box    = box;
   d_ghosts = ghosts;
}

template<int DIM>  FaceGeometry<DIM>::~FaceGeometry()
{
}

/*
*************************************************************************
*									*
* Attempt to calculate the intersection between two face centered box	*
* geometries.  The calculateOverlap() checks whether both arguments are	*
* face geometries; if so, it compuates the intersection.  If not, then	*
* it calls calculateOverlap() on the source object (if retry is true)	*
* to allow the source a chance to calculate the intersection.  See the	*
* hier::BoxGeometry<DIM> base class for more information about the protocol.	*
* A pointer to null is returned if the intersection cannot be computed.	*
* 									*
*************************************************************************
*/

template<int DIM> 
tbox::Pointer< hier::BoxOverlap<DIM> > FaceGeometry<DIM>::calculateOverlap(
   const hier::BoxGeometry<DIM>& dst_geometry,
   const hier::BoxGeometry<DIM>& src_geometry,
   const hier::Box<DIM>& src_mask,
   const bool overwrite_interior,
   const hier::IntVector<DIM>& src_offset,
   const bool retry) const
{
   const FaceGeometry<DIM> *t_dst = 
      dynamic_cast<const FaceGeometry<DIM> *>(&dst_geometry);
   const FaceGeometry<DIM> *t_src =
      dynamic_cast<const FaceGeometry<DIM> *>(&src_geometry);

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
* Convert an AMR-index space hier::Box into a face-index space box by a	*
* cyclic shift of indices.						*
*									*
*************************************************************************
*/

template<int DIM> hier::Box<DIM> 
FaceGeometry<DIM>::toFaceBox(
   const hier::Box<DIM>& box, 
   int face_normal)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT( (face_normal >= 0) && (face_normal < DIM) );
#endif

   hier::Box<DIM> face_box;

   if (!box.empty()) {
      const int x = face_normal;
      face_box.lower(0) = box.lower(x);
      face_box.upper(0) = box.upper(x)+1;
      for (int i = 1; i < DIM; i++) {
         const int y = (face_normal + i) % DIM;
         face_box.lower(i) = box.lower(y);
         face_box.upper(i) = box.upper(y);
      }
   }

   return(face_box);
}

/*
*************************************************************************
*									*
* Compute the overlap between two face centered boxes.  The algorithm	*
* is fairly straight-forward.  First, we perform a quick-and-dirty	*
* intersection to see if the boxes might overlap.  If that intersection	*
* is not empty, then we need to do a better job calculating the overlap	*
* for each dimension.  Note that the AMR index space boxes must be	*
* shifted into the face centered space before we calculate the proper	*
* intersections.							*
*									*
*************************************************************************
*/

template<int DIM> 
tbox::Pointer< hier::BoxOverlap<DIM> > FaceGeometry<DIM>::doOverlap(
   const FaceGeometry<DIM>& dst_geometry,
   const FaceGeometry<DIM>& src_geometry,
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

   // Compute the intersection (if any) for each of the face directions

   const hier::Box<DIM> quick_check =
      hier::Box<DIM>::grow(src_shift, 1) * hier::Box<DIM>::grow(dst_ghost, 1);

   if (!quick_check.empty()) {
      for (int d = 0; d < DIM; d++) {
         const hier::Box<DIM> dst_face = toFaceBox(dst_ghost, d);
         const hier::Box<DIM> src_face = toFaceBox(src_shift, d);
         const hier::Box<DIM> together = dst_face * src_face;
         if (!together.empty()) {
            if (!overwrite_interior) {
               const hier::Box<DIM> int_face = toFaceBox(dst_geometry.d_box, d);
               dst_boxes[d].removeIntersections(together,int_face);
            } else {
               dst_boxes[d].appendItem(together);
            }
         }  // if (!together.empty())
      }  // loop over dim
   }  // !quick_check.empty()

   // Create the face overlap data object using the boxes and source shift

   hier::BoxOverlap<DIM> *overlap = new FaceOverlap<DIM>(dst_boxes, src_offset);
   return(tbox::Pointer< hier::BoxOverlap<DIM> >(overlap));
}

/*
*************************************************************************
*                                                                       *
* Restrict an overlap to the data owned by the source box.  A face      *
* with normal direction d is touched by the cell with the same index    *
* and by the cell below it in direction d.  Face indices are permuted   *
* so that the normal direction comes first.                             *
*                                                                       *
*************************************************************************
*/

template<int DIM> void
FaceGeometry<DIM>::computeOwnedBorderData(
   tbox::Array< hier::BoxList<DIM> >& owned_border_data,
   const hier::BoxList<DIM>& level_boxes,
   const hier::Box<DIM>& owner_box) const
{
   hier::BoxList<DIM> boxes(level_boxes);
   boxes.coalesceBoxes();

   /*
    * The index of a face with normal direction d starts with its
    * component in direction d, so store the boxes in that order.
    */
   owned_border_data.resizeArray(DIM);
   hier::IntVector<DIM> offsets[2];
   for (int d = 0; d < DIM; d++) {
      offsets[0] = hier::IntVector<DIM>(0);
      offsets[1] = hier::IntVector<DIM>(0);
      offsets[1](d) = 1;
      hier::BoxGeometry<DIM>::computeOwnedBorderBoxes(
         owned_border_data, d, boxes, owner_box, offsets, 2);

      hier::BoxList<DIM> face_boxes;
      for (typename hier::BoxList<DIM>::Iterator b(owned_border_data[d]); b; b++) {
         hier::Box<DIM> face_box;
         for (int i = 0; i < DIM; i++) {
            face_box.lower(i) = b().lower((d + i) % DIM);
            face_box.upper(i) = b().upper((d + i) % DIM);
         }
         face_boxes.appendItem(face_box);
      }
      owned_border_data[d] = face_boxes;
   }
}

template<int DIM> tbox::Pointer< hier::BoxOverlap<DIM> >
FaceGeometry<DIM>::restrictOverlapToOwnedData(
   const tbox::Pointer< hier::BoxOverlap<DIM> >& overlap,
   const hier::Box<DIM>& src_box,
   const tbox::Array< hier::BoxList<DIM> >& owned_border_data) const
{
   const FaceOverlap<DIM>* t_overlap =
      dynamic_cast<const FaceOverlap<DIM>*>(overlap.getPointer());
   if (t_overlap == NULL || owned_border_data.getSize() != DIM) {
      return(overlap);
   }

   hier::BoxList<DIM> dst_boxes[DIM];
   hier::IntVector<DIM> offsets[2];
   offsets[0] = hier::IntVector<DIM>(0);
   offsets[1] = hier::IntVector<DIM>(0);
   offsets[1](0) = 1;
   for (int d = 0; d < DIM; d++) {
      hier::Box<DIM> src_face_box;
      for (int i = 0; i < DIM; i++) {
         src_face_box.lower(i) = src_box.lower((d + i) % DIM);
         src_face_box.upper(i) = src_box.upper((d + i) % DIM);
      }
      hier::BoxGeometry<DIM>::intersectOverlapBoxes(
         dst_boxes[d], t_overlap->getDestinationBoxList(d), src_face_box,
         owned_border_data, d, offsets, 2);
   }

   return(new FaceOverlap<DIM>(dst_boxes, t_overlap->getSourceOffset()));
}

/*
*************************************************************************
*                                                                       *
* Compute the overlap between two face centered boxes as doOverlap()    *
* does, but keep only the data owned by the source box, as              *
* restrictOverlapToOwnedData() does.                                    *
*                                                                       *
*************************************************************************
*/

template<int DIM> tbox::Pointer< hier::BoxOverlap<DIM> >
FaceGeometry<DIM>::calculateOwnedOverlap(
   const hier::BoxGeometry<DIM>& src_geometry,
   const hier::Box<DIM>& src_mask,
   const bool overwrite_interior,
   const hier::IntVector<DIM>& src_offset,
   const hier::Box<DIM>& src_box,
   const tbox::Array< hier::BoxList<DIM> >& owned_border_data) const
{
   const FaceGeometry<DIM> *t_src =
      dynamic_cast<const FaceGeometry<DIM> *>(&src_geometry);
   if (t_src == NULL || owned_border_data.getSize() != DIM) {
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
      hier::IntVector<DIM> offsets[2];
      offsets[0] = hier::IntVector<DIM>(0);
      offsets[1] = hier::IntVector<DIM>(0);
      offsets[1](0) = 1;
      for (int d = 0; d < DIM; d++) {
         hier::Box<DIM> src_face_box;
         for (int i = 0; i < DIM; i++) {
            src_face_box.lower(i) = src_box.lower((d + i) % DIM);
            src_face_box.upper(i) = src_box.upper((d + i) % DIM);
         }
         const hier::Box<DIM> together =
            toFaceBox(dst_ghost, d) * toFaceBox(src_shift, d);
         if (!together.empty()) {
            if (!overwrite_interior) {
               hier::BoxList<DIM> boxes;
               boxes.removeIntersections(together, toFaceBox(d_box, d));
               hier::BoxGeometry<DIM>::intersectOverlapBoxes(
                  dst_boxes[d], boxes, src_face_box,
                  owned_border_data, d, offsets, 2);
            } else {
               hier::BoxGeometry<DIM>::intersectOverlapBox(
                  dst_boxes[d], together, src_face_box,
                  owned_border_data, d, offsets, 2);
            }
         }
      }
   }

   return(new FaceOverlap<DIM>(dst_boxes, src_offset));
}

/*
*************************************************************************
*                                                                       *
* Remove from an overlap the data touched by the cells of some boxes.   *
*                                                                       *
*************************************************************************
*/

template<int DIM> tbox::Pointer< hier::BoxOverlap<DIM> >
FaceGeometry<DIM>::removeOverlapOnBoxes(
   const tbox::Pointer< hier::BoxOverlap<DIM> >& overlap,
   const hier::BoxList<DIM>& boxes) const
{
   const FaceOverlap<DIM>* t_overlap =
      dynamic_cast<const FaceOverlap<DIM>*>(overlap.getPointer());
   if (t_overlap == NULL) {
      return(overlap);
   }

   hier::BoxList<DIM> dst_boxes[DIM];
   for (int d = 0; d < DIM; d++) {
      dst_boxes[d] = t_overlap->getDestinationBoxList(d);
      for (typename hier::BoxList<DIM>::Iterator b(boxes); b; b++) {
         dst_boxes[d].removeIntersections(toFaceBox(b(), d));
      }
   }

   return(new FaceOverlap<DIM>(dst_boxes, t_overlap->getSourceOffset()));
}

}
}
#endif

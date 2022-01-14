//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-4-4/source/patchdata/boxgeometry/FaceGeometry.h $
// Package:	SAMRAI patch data geometry
// Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 2224 $
// Modified:	$LastChangedDate: 2008-06-20 17:51:16 -0700 (Fri, 20 Jun 2008) $
// Description:	hier::Box geometry information for face centered objects
//

#ifndef included_pdat_FaceGeometry
#define included_pdat_FaceGeometry

#include "SAMRAI_config.h"
#include "Box.h"
#include "BoxGeometry.h"
#include "BoxOverlap.h"
#include "FaceOverlap.h"
#include "IntVector.h"
#include "tbox/Pointer.h"

namespace SAMRAI {
    namespace pdat {

/*!
 * Class FaceGeometry<DIM> manages the mapping between the AMR index space
 * and the face-centered geometry index space.  It is a subclass of
 * hier::BoxGeometry<DIM> and it computes intersections between face-
 * centered box geometries for communication operations.
 *
 * See header file for FaceData<DIM> class for a more detailed
 * description of the data layout.
 *
 * @see hier::BoxGeometry
 * @see pdat::FaceOverlap
 */

template<int DIM> class FaceGeometry : public hier::BoxGeometry<DIM>
{
public:

   /*!
    * The BoxOverlap implemenation for this geometry.
    */
   typedef FaceOverlap<DIM> Overlap;

  /*!
    * @brief Convert an AMR index box space box into an face geometry box.
    * An face geometry box extends the given AMR index box space box
    * by one in upper dimension for the face normal coordinate direction.
    * 
    * Recall that box indices are cyclically shifted such that the face normal
    * direction is the first coordinate index.  See SideData header file.
    */
   static hier::Box<DIM> toFaceBox(const hier::Box<DIM>& box,
                                   int face_normal);

   /*!
    * @brief Construct the face geometry object given an AMR index
    * space box and ghost cell width.
    */
   FaceGeometry(const hier::Box<DIM>& box,
                const hier::IntVector<DIM>& ghosts);

   /*!
    * @brief The virtual destructor does nothing interesting.
    */
   virtual ~FaceGeometry<DIM>();

   /*!
    * @brief Compute the overlap in face-centered index space between
    * the source box geometry and the destination box geometry.
    */
   virtual tbox::Pointer< hier::BoxOverlap<DIM> > calculateOverlap(
      const hier::BoxGeometry<DIM>& dst_geometry,
      const hier::BoxGeometry<DIM>& src_geometry,
      const hier::Box<DIM>& src_mask,
      const bool overwrite_interior,
      const hier::IntVector<DIM>& src_offset,
      const bool retry) const;

   /*!
    * @brief Return the box for this face centered box geometry
    * object.
    */
   const hier::Box<DIM>& getBox() const;

   /*!
    * @brief Return the ghost cell width for this face centered box
    * geometry object.
    */
   const hier::IntVector<DIM>& getGhosts() const;

private:
   /**
    * Function doOverlap() is the function that computes the overlap
    * between the source and destination objects, where both box geometry
    * objects are guaranteed to have face centered geometry.
    */
   static tbox::Pointer< hier::BoxOverlap<DIM> > doOverlap(
      const FaceGeometry<DIM>& dst_geometry,
      const FaceGeometry<DIM>& src_geometry,
      const hier::Box<DIM>& src_mask,
      const bool overwrite_interior,
      const hier::IntVector<DIM>& src_offset);

   FaceGeometry(const FaceGeometry<DIM>&);	// not implemented
   void operator=(const FaceGeometry<DIM>&);		// not implemented

   hier::Box<DIM> d_box;
   hier::IntVector<DIM> d_ghosts;

};

}
}
#ifndef DEBUG_NO_INLINE
#include "FaceGeometry.I"
#endif
#endif

#ifdef INCLUDE_TEMPLATE_IMPLEMENTATION
#include "FaceGeometry.C"
#endif

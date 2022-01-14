//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-4-4/source/patchdata/boxgeometry/NodeGeometry.h $
// Package:	SAMRAI patch data geometry
// Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 3061 $
// Modified:	$LastChangedDate: 2009-03-19 16:03:30 -0700 (Thu, 19 Mar 2009) $
// Description:	hier::Box geometry information for node centered objects
//

#ifndef included_pdat_NodeGeometry
#define included_pdat_NodeGeometry

#include "SAMRAI_config.h"
#include "Box.h"
#include "BoxGeometry.h"
#include "BoxOverlap.h"
#include "NodeOverlap.h"
#include "IntVector.h"
#include "tbox/Pointer.h"

namespace SAMRAI {
    namespace pdat {

/*!
 * Class NodeGeometry<DIM> manages the mapping between the AMR index space
 * and the node-centered geometry index space.  It is a subclass of
 * hier::BoxGeometry<DIM> and it computes intersections between node-
 * centered box geometries for communication operations.
 *
 * See header file for NodeData<DIM> class for a more detailed
 * description of the data layout.
 *
 * @see hier::BoxGeometry
 * @see pdat::NodeOverlap
 */

template<int DIM> class NodeGeometry : public hier::BoxGeometry<DIM>
{
public:

   /*!
    * The BoxOverlap implemenation for this geometry.
    */
   typedef NodeOverlap<DIM> Overlap;

   /*!
    * @brief Convert an AMR index box space box into a node geometry box.
    * A node geometry box is extends the given AMR index box space box
    * by one in upper dimension for each coordinate direction.
    */
   static hier::Box<DIM> toNodeBox(const hier::Box<DIM>& box);

   /*!
    * @brief Construct the node geometry object given an AMR index
    * space box and ghost cell width. 
    */
   NodeGeometry(const hier::Box<DIM>& box, 
                const hier::IntVector<DIM>& ghosts);

   /*!
    * @brief The virtual destructor does nothing interesting.
    */
   virtual ~NodeGeometry<DIM>();

   /*!
    * @brief Compute the overlap in node-centered index space between
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
    * @brief Compute the node-centered destination boxes that represent
    * the overlap between the source box geometry and the destination
    * box geometry.
    */
   void computeDestinationBoxes(
      hier::BoxList<DIM>& dst_boxes,
      const NodeGeometry<DIM>& src_geometry,
      const hier::Box<DIM>& src_mask,
      const bool overwrite_interior,
      const hier::IntVector<DIM>& src_offset) const;

   /*!
    * @brief Return the box for this node centered box geometry
    * object.
    */
   const hier::Box<DIM>& getBox() const;
 
   /*!
    * @brief Return the ghost cell width for this node centered box
    * geometry object.
    */
   const hier::IntVector<DIM>& getGhosts() const;

private:
   /**
    * Function doOverlap() is the function that computes the overlap
    * between the source and destination objects, where both box geometry
    * objects are guaranteed to have node centered geometry.
    */
   static tbox::Pointer< hier::BoxOverlap<DIM> > doOverlap(
      const NodeGeometry<DIM>& dst_geometry,
      const NodeGeometry<DIM>& src_geometry,
      const hier::Box<DIM>& src_mask,
      const bool overwrite_interior,
      const hier::IntVector<DIM>& src_offset);

   NodeGeometry(const NodeGeometry<DIM>&);	// not implemented
   void operator=(const NodeGeometry<DIM>&);		// not implemented

   hier::Box<DIM> d_box;
   hier::IntVector<DIM> d_ghosts;

};

}
}
#ifndef DEBUG_NO_INLINE
#include "NodeGeometry.I"
#endif
#endif

#ifdef INCLUDE_TEMPLATE_IMPLEMENTATION
#include "NodeGeometry.C"
#endif

//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-4-4/source/hierarchy/variables/BoxGeometry.h $
// Package:	SAMRAI hierarchy
// Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 2132 $
// Modified:	$LastChangedDate: 2008-04-14 14:51:47 -0700 (Mon, 14 Apr 2008) $
// Description:	Box geometry description for overlap computations
//

#ifndef included_hier_BoxGeometry
#define included_hier_BoxGeometry

#include "SAMRAI_config.h"
#include "Box.h"
#include "BoxList.h"
#include "BoxOverlap.h"
#include "IntVector.h"
#include "tbox/Array.h"
#include "tbox/Pointer.h"
#include "tbox/DescribedClass.h"

namespace SAMRAI {
   namespace hier {

/**
 * Class BoxGeometry<DIM> encapsulates the geometry information associated
 * with data defined over a box region, such as its ghost cell width and the 
 * centering or geometry of the data index space.  The intersection (or overlap) 
 * of two box geometries is generated via member function calculateOverlap().  The
 * form of this overlap depends on the particular geometry represented by the
 * subclass.
 *
 * Box geometry objects are created by the patch data factories since patch
 * data objects may not be available for patches that are distributed across
 * processor memories (patch data factories are always replicated).
 *
 * The concept of ``overlap'' or data dependency is more complex for generic
 * box geometry objects than for just cell-centered box indices in the abstract
 * AMR index space.  Problems arise in cases where data lies on the outside
 * corners, faces, or edges of a box.  For these data types, it is likely 
 * that there will exist duplicate data values on different patches.
 *
 * The solution implemented here introduces the concept of ``priority''
 * between patches.  Data of patches with higher priority can overwrite
 * the interiors (face, node, or edge values associated with cells that
 * constitute the interior of the patch) of patches with lower priorities,
 * but lower priority patches can never overwrite the interiors of higher
 * priority patches.  This scheme introduces a total ordering of data and
 * therefore eliminates the duplicate information problem.
 *
 * In practice, this protocol means two things: (1) the communication
 * routines must always process copies from low priority sources to high
 * priority sources, and (2) patches must be given special permission to
 * overwrite their interior values during a write.  All destinations are
 * therefore represented by three quantities: (1) the box geometry of the
 * destination (which encodes the box, ghost cells, and geometry), (2) the
 * box geometry of the source, and (3) a flag indicating whether the source
 * has a higher priority than the destination (that is, whether the source
 * can overwrite the interior of the destination).  If the overwrite flag is
 * set, then data will be copied over the specified box domain and may write
 * into the interior of the destination.  If the overwrite flag is not set,
 * then data will be copied only into the ghost cell values and not the
 * interior values of the patch.
 *
 * @see hier::BoxOverlap
 * @see hier::PatchDataFactory
 * @see hier::PatchData
 */

template<int DIM> class BoxGeometry : public tbox::DescribedClass
{
public:
   /**
    * The default constructor for BoxGeometry<DIM> does nothing interesting.
    */
   BoxGeometry();

   /**
    * The virtual destructor does nothing interesting.
    */
   virtual ~BoxGeometry();

   /**
    * Calculate the overlap between two box geometry objects given the
    * source and destination (given by this) geometries, a source mask,
    * the priority overwrite flag, and an offset between source and
    * destination index spaces.  The box overlap description returned by
    * this function will be used in later copy and pack/unpack calls on
    * the patch data object.  The offset is from the source space 
    * into the destination index space.  That is, if p is in the source
    * index space, then p + sourceOffset is the corresponding point in the
    * destination index space.  The overwrite flag is used to represent
    * priority between patches.  If it is set, then the copy is allowed
    * to modify the interior of the destination region.  Note that the
    * source and destination box geometries encode the geometry of the box
    * that they represent; thus, it is possible to calculate intersections
    * between different geometries.  This will be necessary when copying
    * data from flux sum counters into a face centered array in the AMR
    * flux synchronization algorithm.
    */
   tbox::Pointer< BoxOverlap<DIM> >
   calculateOverlap(const BoxGeometry<DIM>& src_geometry,
                    const Box<DIM>& src_mask,
                    const bool overwrite_interior,
                    const IntVector<DIM>& src_offset) const;
 
   /**
    * Calculate the overlap between two box geometry objects given the
    * source and destination geometries.  This form calculateOverlap() is
    * redefined by the subclasses of BoxGeometry<DIM> for the appropriate
    * intersection algorithms.  If calculateOverlap() cannot compute the
    * intersection between the two given geometries and retry is true, then
    * calculateOverlap() is called on the destination geometry object with
    * retry set to false (to avoid infinite recursion).  This protocol
    * makes it possible to add new box geometry types and still calculate
    * intersections with existing box geometry types.
    */
   virtual tbox::Pointer< BoxOverlap<DIM> >
   calculateOverlap(const BoxGeometry<DIM>& dst_geometry,
                    const BoxGeometry<DIM>& src_geometry,
                    const Box<DIM>& src_mask,
                    const bool overwrite_interior,
                    const IntVector<DIM>& src_offset,
                    const bool retry) const = 0;

   /**
    * Return the part of the given overlap whose destination data are owned
    * by the given source box.  Data on the border of a box are shared by
    * every box that touches them.  A box owns a datum if, among the cells
    * of the given level boxes that touch the datum, the one with the
    * largest index (comparing the last coordinate first) is in the box.
    * In particular, a box owns the data in its interior and on its lower
    * border.  The source box and the level boxes are given in the
    * destination index space, and the level boxes must cover the source
    * box.
    *
    * The default implementation returns the given overlap, which is
    * appropriate for data that are not shared between boxes.
    */
   virtual tbox::Pointer< BoxOverlap<DIM> >
   restrictOverlapToOwnedData(const tbox::Pointer< BoxOverlap<DIM> >& overlap,
                              const Box<DIM>& src_box,
                              const BoxList<DIM>& level_boxes) const;

   /**
    * Return the part of the given overlap whose destination data are not
    * touched by any cell of the given boxes, which are given in the
    * destination index space.
    *
    * The default implementation returns the given overlap, which is
    * appropriate for data that are not shared between boxes.
    */
   virtual tbox::Pointer< BoxOverlap<DIM> >
   removeOverlapOnBoxes(const tbox::Pointer< BoxOverlap<DIM> >& overlap,
                        const BoxList<DIM>& boxes) const;

protected:
   /**
    * Compute the indices of the data owned by the given box, as described
    * in restrictOverlapToOwnedData(), for data located at the given offsets
    * from the index of each cell that touches them.  The offsets must be
    * sorted so that an earlier offset corresponds to a cell with a larger
    * index; the first is therefore zero.
    */
   static void computeOwnedDataBoxes(BoxList<DIM>& owned_boxes,
                                     const Box<DIM>& box,
                                     const BoxList<DIM>& level_boxes,
                                     const tbox::Array< IntVector<DIM> >& offsets);

   /**
    * Intersect each of the given overlap boxes with the given owned boxes
    * and append the results to a list.  Each resulting box is contained in
    * one of the overlap boxes: boxes that come from different overlap boxes
    * are never merged, since an overlap box may correspond to a separate
    * piece of the source data.
    */
   static void intersectOverlapBoxes(BoxList<DIM>& result_boxes,
                                     const BoxList<DIM>& overlap_boxes,
                                     const BoxList<DIM>& owned_boxes);

private:
   BoxGeometry(const BoxGeometry<DIM>&);	// not implemented
   void operator=(const BoxGeometry<DIM>&);	// not implemented

};

}
}

#ifndef DEBUG_NO_INLINE
#include "BoxGeometry.I"
#endif
#endif

#ifdef INCLUDE_TEMPLATE_IMPLEMENTATION
#include "BoxGeometry.C"
#endif

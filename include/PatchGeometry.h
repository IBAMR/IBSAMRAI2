//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-4-4/source/hierarchy/patches/PatchGeometry.h $
// Package:	SAMRAI hierarchy package
// Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 2147 $
// Modified:	$LastChangedDate: 2008-04-23 16:48:12 -0700 (Wed, 23 Apr 2008) $
// Description: Base class for geometry management on patches
//

#ifndef included_hier_PatchGeometry
#define included_hier_PatchGeometry

#include "SAMRAI_config.h"

#ifndef included_iostream
#define included_iostream
#include <iostream>
#endif

#include "tbox/Array.h"
#include "BoundaryBox.h"
#include "IntVector.h"
#include "tbox/DescribedClass.h"
#include "tbox/List.h"


namespace SAMRAI {
    namespace hier {

/**
 * Class PatchGeometry is the base class for geometry classes that
 * manage index space and mesh increment information on individual patches.
 * Patch geometry information is used for setting boundary conditions in
 * ghost cells and is used in the inter-level transfer operators for refining
 * or coarsening data between two patches associated with different index 
 * spaces.  The boundary information for patches is actually computed by 
 * the GridGeometry class.
 * 
 * @see hier::BoundaryBox
 * @see hier::GridGeometry
 */

template<int DIM> class PatchGeometry  : public tbox::DescribedClass
{
public:
   /**
    * The default constructor for the patch geometry base class.
    */
   PatchGeometry(
      const IntVector<DIM>& ratio_to_level_zero,
      const tbox::Array< tbox::Array<bool> >& touches_regular_bdry,
      const tbox::Array< tbox::Array<bool> >& touches_periodic_bdry);

   /**
    * The virtual destructor for the patch geometry base class.
    */
   virtual ~PatchGeometry();

   /**
    * Return const reference to patch boundary information.
    */ 
   const tbox::Array< BoundaryBox<DIM> >* getPatchBoundaries() const; 

   /*!
    * @brief Set the boundary box arrays for this patch geometry.
    *
    * An array of length DIM of tbox::Array< BoundaryBox<DIM> > is passed
    * in to be stored as the boundary boxes for this patch geometry.
    *
    * @param bdry The array of BoundaryBox arrays.
    */
   void setBoundaryBoxesOnPatch(
      const tbox::Array< BoundaryBox<DIM> > bdry[DIM]);

   /**
    * Return const reference to ratio to level zero index space.
    */
   const IntVector<DIM>& getRatio() const;

   /**
    * Return a boolean value indicating whether the patch boundary 
    * intersects the physical domain boundary in a non-periodic
    * direction.  In other words, the return value is true when the
    * patch has non-empty boundary boxes that lie outside the physical
    * domain.  Otherwise, the return value is false.  Note that when
    * a patch touches the "boundary" of the physical domain in a periodic 
    * direction, there are no boundary boxes to fill; the data is filled
    * from the proper region of the domain interior in the periodic direction.
    */
   bool intersectsPhysicalBoundary() const;

   /**
    * Return array of boundary box components for patch each of which
    * intersects the patch at a single point (i.e., 0-dim intersection
    * between cells in patch and cells in boundary box).
    */
   const tbox::Array< BoundaryBox<DIM> >& getNodeBoundaries() const; 

   /**
    * Return array of boundary box components for patch each of which
    * intersects the patch along a 1-dim edge (i.e., 1-dim intersection
    * between cells in patch and cells in boundary box).
    *
    * When assertion checking is active, this routine throws an assertion
    * when DIM < 2.
    */
   const tbox::Array< BoundaryBox<DIM> >& getEdgeBoundaries() const;

   /**
    * Return array of boundary box components for patch each of which
    * intersects the patch along a 2-dim face (i.e., 2-dim intersection
    * between cells in patch and cells in boundary box).
    *
    * When assertion checking is active, this routine throws an assertion
    * when DIM < 3.
    */
   const tbox::Array< BoundaryBox<DIM> >& getFaceBoundaries() const;

   /**
    * Return array of boundary box components for patch each of which
    * intersects the patch as a (DIM - codim)-dimensional object.
    * That is,
    *
    * if DIM == 1: (codim == 1) => same components as getNodeBoundaries.
    *  
    * if DIM == 2, (codim == 1) => same components as getEdgeBoundaries.  
    *              (codim == 2) => same components as getNodeBoundaries.  
    *  
    * if DIM == 3, (codim == 1) => same components as getFaceBoundaries.  
    *              (codim == 2) => same components as getEdgeBoundaries.  
    *              (codim == 3) => same components as getNodeBoundaries.  
    *
    * When assertion checking is active, this routine throws an assertion
    * when codim < 0 or codim > DIM.
    */
   const tbox::Array< BoundaryBox<DIM> >& getCodimensionBoundaries(
      const int codim) const;

   /**
    * Set the array of boundary box components of the given codimension
    * for a patch.
    */
   void setCodimensionBoundaries(
      const tbox::Array< BoundaryBox<DIM> >& bdry_boxes,
      const int codim);

#if (INCLUDE_DEPRECATED <= 2)
   /*
    * Deprecated methods.  These are identical to the previous calls
    * and should not be used.
    */
   const tbox::Array< BoundaryBox<DIM> >* getPatchBoundary() const; 
   const tbox::Array< BoundaryBox<DIM> >& getNodeBoundary() const; 
   const tbox::Array< BoundaryBox<DIM> >& getEdgeBoundary() const;
   const tbox::Array< BoundaryBox<DIM> >& getFaceBoundary() const;
   const tbox::Array< BoundaryBox<DIM> >& getCodimensionBoundary(
      const int codim) const;
   void setCodimensionBoundary(
      const tbox::Array< BoundaryBox<DIM> >& bdry_boxes,
      const int codim);
#endif



   /*!
    * @brief Compute a box outside a physical domain that needs to be filled.
    *
    * The patch box will be grown by the given ghost cell width and
    * then intersected with the boundary box.  The resulting intersection
    * will be grown to the needed ghost cell width in the direction
    * normal to the boundary.
    *
    * @param bbox BoundaryBox representing location and type of boundary
    * @param patch_box The box for the patch where data is being filled
    * @param gcw ghost cell width to fill
    */ 
   Box<DIM> getBoundaryFillBox(const BoundaryBox<DIM>& bbox,
                                const Box<DIM>& patch_box,
                                const IntVector<DIM>& gcw) const;

   /*!
    * @brief Query whether patch touches a regular boundary
    *
    * Returns true if the Patch touches any non-periodic physical boundary
    */
   bool getTouchesRegularBoundary() const;

   /*!
    * @brief Query whether patch touches a regular boundary
    *
    * Returns true if the Patch touches any periodic boundary
    */
   bool getTouchesPeriodicBoundary() const;

   /*!
    * @brief Query whether patch touches a specific regular boundary
    *
    * Returns true if the Patch touches a non-periodic physical boundary
    * on the side of the Patch specified in the argument list.  The side
    * is specified by an axis direction and a flag specified the upper or
    * lower side.
    *
    * @param axis       Axis direction normal to the side being checked
    * @param upperlower Flag should be 0 if checking the lower side in the
    *                   axis direction, or 1 if checking the upper side.
    */
   bool getTouchesRegularBoundary(int axis, int upperlower) const;

   /*!
    * @brief Query whether patch touches a specific regular boundary
    *
    * Returns true if the Patch touches a periodic boundary
    * on the side of the Patch specified in the argument list.  The side
    * is specified by an axis direction and a flag specified the upper or
    * lower side.
    *
    * @param axis       Axis direction normal to the side being checked
    * @param upperlower Flag should be 0 if checking the lower side in the
    *                   axis direction, or 1 if checking the upper side.
    */
   bool getTouchesPeriodicBoundary(int axis, int upperlower) const;

   /**
    * Print object data to the specified output stream.
    */
   virtual void printClassData(std::ostream& stream) const;

private:
   bool d_has_regular_boundary;
   bool d_has_periodic_boundary;
   IntVector<DIM> d_ratio_to_level_zero;
   tbox::Array< BoundaryBox<DIM> > d_patch_boundaries[DIM];

   bool d_touches_regular_bdry[DIM][2];
   bool d_touches_periodic_bdry[DIM][2];
};

}
}
#ifndef DEBUG_NO_INLINE
#include "PatchGeometry.I"
#endif
#endif

#ifdef INCLUDE_TEMPLATE_IMPLEMENTATION
#include "PatchGeometry.C"
#endif

//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-4-4/source/patchdata/outerface/OuterfaceDataFactory.h $
// Package:	SAMRAI patch data
// Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 2132 $
// Modified:	$LastChangedDate: 2008-04-14 14:51:47 -0700 (Mon, 14 Apr 2008) $
// Description: Factory class for creating outerface data objects
//

#ifndef included_pdat_OuterfaceDataFactory
#define included_pdat_OuterfaceDataFactory

#include "SAMRAI_config.h"
#include "Box.h"
#include "BoxGeometry.h"
#include "IntVector.h"
#include "PatchDataFactory.h"
#include "tbox/Arena.h"
#include "tbox/Complex.h"
#include "tbox/Pointer.h"

#ifndef NULL
#define NULL (0)
#endif

namespace SAMRAI {
    namespace pdat {

/**
 * Class OuterfaceDataFactory is a factory class used to allocate new
 * instances of OuterfaceData objects.  It is a subclass of the patch
 * data factory class and outerface data is a subclass of patch data.  Both
 * the factory and data classes are templated on the type of the contained
 * object (e.g., double or int).
 *
 * @see pdat::OuterfaceData
 * @see pdat::PatchDataFactory
 */

template<int DIM, class TYPE>
class OuterfaceDataFactory : public hier::PatchDataFactory<DIM>
{
public:
   /**
    * The default constructor for the outerface data factory class.
    * The depth (number of components) gives the default for all of
    * the outerface data objects created with this factory.
    */
   OuterfaceDataFactory(int depth);

   /**
    * Virtual destructor for the outerface data factory class.
    */
   virtual ~OuterfaceDataFactory();

   /**
    * @brief Abstract virtual function to clone a patch data factory.

    * This will return a new instantiation of the abstract factory
    * with the same properties.  The properties of the cloned factory
    * can then be changed without modifying the original.
    *
    * @param ghosts default ghost cell width for concrete classes created from
    * the factory.
    */
   virtual tbox::Pointer< hier::PatchDataFactory<DIM> > cloneFactory(const hier::IntVector<DIM>& ghosts);

   /**
    * Virtual factory function to allocate a concrete outerface data object.
    * The default information about the object (e.g., depth) is taken from
    * the factory.  If no memory pool is provided, then the allocation routine
    * assumes some default memory pool.
    */
   virtual tbox::Pointer< hier::PatchData<DIM> > allocate(
      const hier::Box<DIM>& box,
      tbox::Pointer<tbox::Arena> pool = tbox::Pointer<tbox::Arena>(NULL)) const;

   /**
    * Virtual factory function to allocate a concrete cell data object.
    * Same as above function, except passes in a patch instead of a box.
    */
   virtual tbox::Pointer< hier::PatchData<DIM> > allocate(
      const hier::Patch<DIM>& patch,
      tbox::Pointer<tbox::Arena> pool = tbox::Pointer<tbox::Arena>(NULL)) const;

   /**
    * Allocate the box geometry object associated with the patch data.
    * This information will be used in the computation of intersections
    * and data dependencies between objects.
    */
   virtual tbox::Pointer< hier::BoxGeometry<DIM> >
   getBoxGeometry(const hier::Box<DIM>& box) const;

   /**
    * Get the default depth (number of components).  This is the default
    * depth that will be used in the instantiation of outerface data objects.
    */
   int getDefaultDepth() const;

   /**
    * Set the default depth (number of components).  This is the default
    * depth that will be used in the instantiation of outerface data objects.
    */
   void setDefaultDepth(const int depth);

   /**
    * Calculate the amount of memory needed to store the outerface data
    * object, including object data and dynamically allocated data.
    */
   virtual size_t getSizeOfMemory(const hier::Box<DIM>& box) const;

   /**
    * Return a boolean true value indicating that fine data for the outerface quantity will 
    * take precedence on coarse-fine interfaces.  See the OuterfaceVariable<DIM> class 
    * header file for more information.
    */
   bool fineBoundaryRepresentsVariable() const {return true;}

   /**
    * Return true since the outerface data index space extends beyond the interior of
    * patches.  That is, outerface data lives on patch borders.
    */
   bool dataLivesOnPatchBorder() const {return true;}

   /**
    * Return whether it is valid to copy this OuterfaceDataFactory to the 
    * supplied destination patch data factory.  It will return true if 
    * dst_pdf is FaceDataFactory or OuterfaceDataFactory, false otherwise.
    */
   bool validCopyTo(
      const tbox::Pointer< hier::PatchDataFactory<DIM> >& dst_pdf) const;

private:
   int d_depth;

   hier::IntVector<DIM> d_no_ghosts;
};

}
}
#ifndef DEBUG_NO_INLINE
#include "OuterfaceDataFactory.I"
#endif
#endif

#ifdef INCLUDE_TEMPLATE_IMPLEMENTATION
#include "OuterfaceDataFactory.C"
#endif

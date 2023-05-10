//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-4-4/source/patchdata/side/SideIterator.h $
// Package:	SAMRAI patch data
// Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 2132 $
// Modified:	$LastChangedDate: 2008-04-14 14:51:47 -0700 (Mon, 14 Apr 2008) $
// Description:	Iterator for side centered patch data types
//

#ifndef included_pdat_SideIterator
#define included_pdat_SideIterator

#include "SAMRAI_config.h"
#include "Box.h"
#include "SideGeometry.h"
#include "SideIndex.h"

namespace SAMRAI {
    namespace pdat {

/**
 * Class SideIterator<DIM> is an iterator that provides methods for
 * stepping through the index space associated with a side centered box.
 * The indices are enumerated in column-major (e.g., Fortran) order.
 * The iterator should be used as follows:
   \verbatim
   hier::Box<DIM> box;
   ...
   for (SideIterator<DIM> c(box, axis); c; c++) {
      // use index c of the box
   }
   \endverbatim
 * Note that the side iterator may not compile to efficient code, depending
 * on your compiler.  Many compilers are not smart enough to optimize the
 * looping constructs and indexing operations.
 *
 * @see pdat::SideData
 * @see pdat::SideGeometry
 * @see pdat::SideIndex
 */

template<int DIM> class SideIterator
{
public:
   /**
    * Default constructor for the side iterator.  The iterator must
    * be initialized before it can be used to iterate over a box.
    */
   SideIterator();

   /**
    * Constructor for the side iterator.  The iterator will enumerate
    * the indices in the argument box.
    */
   SideIterator(const hier::Box<DIM>& box, const int axis);

   /**
    * Copy constructor for the side iterator
    */
   SideIterator(const SideIterator<DIM>& iterator);

   /**
    * Assignment operator for the side iterator.
    */
   SideIterator<DIM>& operator=(const SideIterator<DIM>& iterator);

   /**
    * Destructor for the side iterator.
    */
   ~SideIterator();

   /**
    * Extract the side index corresponding to the iterator position in the box.
    */
   const SideIndex<DIM>& operator*() const;

   /**
    * Extract the side index corresponding to the iterator position in the box.
    */
   const SideIndex<DIM>& operator()() const;

   /**
    * Return true if the iterator points to a valid index within the box.
    */
   operator bool() const;

#ifndef LACKS_BOOL_VOID_RESOLUTION
   /**
    * Return a non-NULL if the iterator points to a valid index within the box.
    */
   operator const void*() const;
#endif

   /**
    * Return whether the iterator points to a valid index within the box.
    * This operator mimics the !p operation applied to a pointer p.
    */
   bool operator!() const;

   /**
    * Increment the iterator to point to the next index in the box.
    */
   void operator++(int);

   /**
    * Test two iterators for equality (same index value).
    */
   bool operator==(const SideIterator<DIM>& iterator) const;

   /**
    * Test two iterators for inequality (different index values).
    */
   bool operator!=(const SideIterator<DIM>& iterator) const;

private:
   SideIndex<DIM> d_index;
   hier::Box<DIM> d_box;
};

}
}
#ifndef DEBUG_NO_INLINE
#include "SideIterator.I"
#endif
#endif

#ifdef INCLUDE_TEMPLATE_IMPLEMENTATION
#include "SideIterator.C"
#endif

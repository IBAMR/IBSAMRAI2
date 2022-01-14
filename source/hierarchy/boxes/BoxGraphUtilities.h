//
// File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-4-4/source/hierarchy/boxes/BoxGraphUtilities.h $
// Package:     SAMRAI hierarchy
// Copyright:   (c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:    $LastChangedRevision: 2132 $
// Modified:    $LastChangedDate: 2008-04-14 14:51:47 -0700 (Mon, 14 Apr 2008) $
// Description: Utility class for operations that reduce complexity of box calculus
//

#ifndef included_hier_BoxGraphUtilities
#define included_hier_BoxGraphUtilities

#include "SAMRAI_config.h"

#include "BoxArray.h"
#include "IntVector.h"
#include "tbox/Array.h"
#include "tbox/List.h"

#ifndef included_iostream
#define included_iostream
#include <iostream>
#endif


namespace SAMRAI {
   namespace hier {

/*!
 * Class BoxGraphUtilities is a utility class that provides
 * methods for "expanding" an array of boxes when some of the
 * boxes touch periodic boundaries.  If a box touches a periodic
 * boundary, then there are one or more "virtual" boxes that
 * "wrap around" to the opposite side(s) of the domain:
 * 
 * 
 * @verbatim
 * 
 *         ---------------------------
 *         |                         |
 *         |                         |
 *         |                         |
 *         |-------------            |-------------
 *         | this       |            | the virtual|
 *         | box touches|            | wrap around|
 *         | a periodic |            | box        |
 *         | boundary   |            |            |
 *         |-------------            |-------------
 *         |                         |
 *         |       domain            |
 *         |        box              |
 *         ---------------------------
 * 
 * @endverbatim
 * 
 * 
 * The BoxTop, BoxGraph, and possibly other classes require
 * as input arrays of boxes in which the "virtual" boxes are
 * explicitly represented.
 */

template<int DIM> struct BoxGraphUtilities
{
   /*!
    * @brief Returns an array of boxes that includes entries for each
    * box that touches a periodic boundary.
    * 
    * The shift array must either have zero length, or have the same length
    * as \b in_boxes, otherwise an unrecoverable error will be thrown.
    * 
    * @param out_boxes contains all items in the \b in_boxes array, and contains
    *                  one or more additional items for any box that touches
    *                  a periodic boundary.
    * @param in_boxes  array of input boxes.
    * @param shifts    shift information for each of the input boxes.
    */
   static void makeBoxesPlusPeriodicBoxes(
      BoxArray<DIM>& out_boxes,
      const BoxArray<DIM>& in_boxes,
      const tbox::Array< tbox::List< IntVector<DIM> > >& shifts);

   /*!
    * @brief Returns an array of boxes that includes entries for each
    * box that touches a periodic boundary.
    * 
    * The shift array must either have zero length, or have the same length
    * as \b in_boxes, otherwise an unrecoverable error will be thrown.
    * 
    * 
    * @param out_boxes contains all items in the \b in_boxes array, and contains
    *                  one or more additional items for any box that touches
    *                  a periodic boundary.
    * @param out_indices contains an entry for each box in \b out_boxes;
    *                    the entry indicates the box, w.r.t \b in_boxes,
    *                    from which the box was derived.
    * @param in_boxes  array of input boxes.
    * @param shifts    shift information for each of the input boxes.
    */
   static void makeBoxesPlusPeriodicBoxes(
      BoxArray<DIM>& out_boxes,
      tbox::Array<int>& out_indices,
      const BoxArray<DIM>& in_boxes,
      const tbox::Array< tbox::List< IntVector<DIM> > >& shifts);

   /*!
    * @brief Returns the sum of shifts[j].getNumberOfItems().
    * 
    * This function is called by makeBoxesPlusPeriodicBoxes().
    * 
    * @param shifts periodic shift information for each box.
    */
   static int countPeriodicBoxes(
      const tbox::Array< tbox::List< IntVector<DIM> > >& shifts);

   /*!
    * @brief Compare function for use with qsort when sorting integers
    * in ascending order.
    * 
    * Sample usage:
    * 
    * @verbatim
    *   intarray[len];
    *   ...
    *   qsort(array, len, sizeof(int), qsortIntCompare);
    * @endverbatim
    * 
    */
   static int qsortIntCompare(const void *v, const void *w);
};


}
}

#ifdef INCLUDE_TEMPLATE_IMPLEMENTATION
#include "BoxGraphUtilities.C"
#endif

#endif


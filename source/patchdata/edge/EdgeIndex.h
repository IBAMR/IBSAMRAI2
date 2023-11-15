//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-4-4/source/patchdata/edge/EdgeIndex.h $
// Package:	SAMRAI patch data
// Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 2132 $
// Modified:	$LastChangedDate: 2008-04-14 14:51:47 -0700 (Mon, 14 Apr 2008) $
// Description:	hier::Index for edge centered patch data types
//

#ifndef included_pdat_EdgeIndex
#define included_pdat_EdgeIndex

#include "SAMRAI_config.h"
#include "IntVector.h"
#include "Index.h"

namespace SAMRAI {
    namespace pdat {

/**
 * Class EdgeIndex<DIM> implements a simple n-dimensional integer
 * vector for edge centered variables.  Edge indices contain an integer
 * index location in AMR index space along with the designated edge axis
 * (X=0, Y=1, or Z=2).  See the edge box geometry class for more information
 * about the mapping between the AMR index space and the edge indices.
 *
 * @see hier::Index
 * @see pdat::EdgeData
 * @see pdat::EdgeGeometry
 * @see pdat::EdgeIterator
 */

template<int DIM> class EdgeIndex : public hier::Index<DIM>
{
public:
   /**
    * The default constructor for a edge index creates an uninitialized index.
    */
   EdgeIndex();

   /**
    * Construct a edge index from a regular index, axis, and edge.  The axis
    * can be one of EdgeIndex<DIM>::X (0), EdgeIndex<DIM>::Y (1), or
    * EdgeIndex<DIM>::Z (2). The edge argument can be one of the constants
    * EdgeIndex<DIM>::LowerLeft (0), EdgeIndex<DIM>::LowerRight(1),
    * EdgeIndex<DIM>::UpperLeft (2) or EdgeIndex<DIM>::UpperRight(3).
    */
   EdgeIndex(const hier::Index<DIM>& rhs, const int axis, const int edge);

   /**
    * The copy constructor creates a edge index equal to the argument.
    */
   EdgeIndex(const EdgeIndex<DIM>& rhs);

   /**
    * The assignment operator sets the edge index equal to the argument.
    */
   EdgeIndex<DIM>& operator=(const EdgeIndex<DIM>& rhs);

   /**
    * Get the axis for which this edge index is defined (X=0, Y=1, Z=2).
    */
   int getAxis() const;

   /**
    * Set the edge axis (X=0, Y=1, Z=2).
    */
   void setAxis(const int axis);

   /**
    * For dimension 2, converts the edge index into the index on the
    * lower side (argument edge = 0) or the upper side (argument
    * edge = 1) of the edge. 
    * For dimension 3, converts the edge index into the index on the lower 
    * left (argument edge == 0), the lower right (argument edge == 1),
    * the upper left (argument edge == 2), or the upper right 
    * (argument edge == 3) of thw edge.
    */
   hier::Index<DIM> toCell(const int edge) const;

   /**
    * Plus-equals operator for a edge index and an integer vector.
    */
   EdgeIndex<DIM>& operator+=(const hier::IntVector<DIM>& rhs);

   /**
    * Plus operator for a edge index and an integer vector.
    */
   EdgeIndex<DIM> operator+(const hier::IntVector<DIM>& rhs) const;

   /**
    * Plus-equals operator for a edge index and an integer.
    */
   EdgeIndex<DIM>& operator+=(const int rhs);

   /**
    * Plus operator for a edge index and an integer.
    */
   EdgeIndex<DIM> operator+(const int rhs) const;

   /**
    * Minus-equals operator for a edge index and an integer vector.
    */
   EdgeIndex<DIM>& operator-=(const hier::IntVector<DIM>& rhs);

   /**
    * Minus operator for a edge index and an integer vector.
    */
   EdgeIndex<DIM> operator-(const hier::IntVector<DIM>& rhs) const;

   /**
    * Minus-equals operator for a edge index and an integer.
    */
   EdgeIndex<DIM>& operator-=(const int rhs);

   /**
    * Minus operator for a edge index and an integer.
    */
   EdgeIndex<DIM> operator-(const int rhs) const;

   /**
    * Times-equals operator for a edge index and an integer vector.
    */
   EdgeIndex<DIM>& operator*=(const hier::IntVector<DIM>& rhs);

   /**
    * Times operator for a edge index and an integer vector.
    */
   EdgeIndex<DIM> operator*(const hier::IntVector<DIM>& rhs) const;

   /**
    * Times-equals operator for a edge index and an integer.
    */
   EdgeIndex<DIM>& operator*=(const int rhs);

   /**
    * Times operator for a edge index and an integer.
    */
   EdgeIndex<DIM> operator*(const int rhs) const;

   /**
    * Returns true if two edge index objects are equal.  All components
    * and the corresponding edge axes must be the same for equality.
    */
   bool operator==(const EdgeIndex<DIM>& rhs) const;

   /**
    * Returns true if two edge index objects are not equal.  Any of
    * the components or axes may be different for inequality.
    */
   bool operator!=(const EdgeIndex<DIM>& rhs) const;

   enum {

      X = 0,
      Y = 1,
      Z = 2,
      Lower = 0,
      Upper = 1,
      LowerLeft = 0,
      LowerRight = 1,
      UpperLeft = 2,
      UpperRight = 3

   };

private:
   int d_axis;
};

}
}
#ifndef DEBUG_NO_INLINE
#include "EdgeIndex.I"
#endif
#endif

#ifdef INCLUDE_TEMPLATE_IMPLEMENTATION
#include "EdgeIndex.C"
#endif

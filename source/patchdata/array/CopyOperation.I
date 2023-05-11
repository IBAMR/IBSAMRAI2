//
// File:	$URL$
// Package:	SAMRAI patch data
// Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 1917 $
// Modified:	$LastChangedDate: 2008-01-25 13:28:01 -0800 (Fri, 25 Jan 2008) $
// Description:	Copy operation on single array data elements templated on data type
//

#include <algorithm>

#ifdef DEBUG_NO_INLINE
#define inline
#endif

namespace SAMRAI {
    namespace pdat {

/*
 * Member functions for CopyOperation
 */

template <class TYPE> inline
CopyOperation<TYPE>::CopyOperation()
{
}

template <class TYPE> inline
CopyOperation<TYPE>::~CopyOperation()
{
}

template <class TYPE> inline
void CopyOperation<TYPE>::operator()(TYPE* vdst, const TYPE* vsrc, const int size) const
{
   // This assumes vdst and vsrc do not overlap, which is reasonable for SAMRAI: there's
   // no situation in which we would copy data from one part of a patch into another
   // overlapping part of the same patch.
   if (vsrc != vdst)
      std::copy(vsrc, vsrc + size, vdst);
}

#ifdef DEBUG_NO_INLINE
#undef inline
#endif

}
}

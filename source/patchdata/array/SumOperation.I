//
// File:	$URL$
// Package:	SAMRAI patch data
// Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 1917 $
// Modified:	$LastChangedDate: 2008-01-25 13:28:01 -0800 (Fri, 25 Jan 2008) $
// Description:	Sum operation on single array data elements templated on data type
//

#ifdef DEBUG_NO_INLINE
#define inline
#endif

namespace SAMRAI {
    namespace pdat {

/*
 * Member functions for SumOperation
 */

template <class TYPE> inline
SumOperation<TYPE>::SumOperation()
{
}

template <class TYPE> inline
SumOperation<TYPE>::~SumOperation()
{
}

template <class TYPE> inline
void SumOperation<TYPE>::operator()(TYPE* vdst, const TYPE* vsrc, const int size) const
{
   std::transform(vdst, vdst+size, vsrc, vdst, std::plus<TYPE>());
}

#ifdef DEBUG_NO_INLINE
#undef inline
#endif

}
}

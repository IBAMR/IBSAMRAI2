//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-4-4/source/toolbox/memory/Array.C $
// Package:	SAMRAI toolbox for memory management
// Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 2195 $
// Modified:	$LastChangedDate: 2008-05-14 11:33:30 -0700 (Wed, 14 May 2008) $
// Description:	A simple array template class
//

#ifndef included_tbox_Array_C
#define included_tbox_Array_C

#include "tbox/Array.h"
#include "tbox/Arena.h"
#include "tbox/Pointer.h"

#ifdef DEBUG_NO_INLINE
#include "tbox/Array.I"
#endif

namespace SAMRAI {
   namespace tbox {

template <class TYPE>
bool Array<TYPE>::Allocator::s_is_available = false;

template <class TYPE>
std::size_t Array<TYPE>::Allocator::s_number_of_allocations = 0;

template <class TYPE>
std::vector<std::vector<TYPE *>> Array<TYPE>::Allocator::s_block_stacks;

template <class TYPE>
Array<TYPE>::Array(const int n)
{
   d_objects  = Allocator::getAllocator().allocate(n);
   d_counter  = new ReferenceCounter;
   d_elements = n;
}

template <class TYPE>
Array<TYPE>::Array(const int n, const Pointer<Arena>& /*pool*/)
{
   d_objects  = Allocator::getAllocator().allocate(n);
   d_counter  = new ReferenceCounter;
   d_elements = n;
}

template <class TYPE>
Array<TYPE>& Array<TYPE>::operator=(const Array<TYPE>& rhs)
{
   if (this != &rhs) {
      if (d_counter && d_counter->deleteReference()) deleteObjects();
      d_objects  = rhs.d_objects;
      d_counter  = rhs.d_counter;
      d_elements = rhs.d_elements;
      if (d_counter) d_counter->addReference();
   }
   return(*this);
}

template <class TYPE>
Array<TYPE>& Array<TYPE>::operator=(Array<TYPE>&& rhs)
{
   if (this != &rhs) {
      Array<TYPE> temp(std::move(rhs));

      std::swap(this->d_objects, temp.d_objects);
      std::swap(this->d_counter, temp.d_counter);
      std::swap(this->d_elements, temp.d_elements);
   }
   return(*this);
}

template <class TYPE>
TYPE *Array<TYPE>::allocateObjects(const int n, Arena */*arena*/)
{
   return Allocator::getAllocator().allocate(n);
}

template <class TYPE>
void Array<TYPE>::deleteObjects()
{
   Allocator::deallocate(d_objects, d_elements);
   delete d_counter;

   d_objects  = (TYPE *) NULL;
   d_counter  = (ReferenceCounter *) NULL;
   d_elements = 0;
}

template <class TYPE>
void Array<TYPE>::resizeArray(const int n)
{
   if (n != d_elements) {
      Array<TYPE> array(n);
      const int s = (d_elements < n ? d_elements : n);
      for (int i = 0; i < s; i++) {
         array.d_objects[i] = d_objects[i];
      }
      this->operator=(array);
   }
}

template <class TYPE>
void Array<TYPE>::resizeArray(
   const int n, const Pointer<Arena>& /*pool*/)
{
   Array<TYPE>::resizeArray(n);
}

template <class TYPE>
std::size_t Array<TYPE>::getNumberOfAllocations()
{
   return Allocator::getNumberOfAllocations();
}


}
}

#endif

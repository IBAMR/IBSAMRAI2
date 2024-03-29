//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-4-4/source/mathops/cell/PatchCellDataOpsReal.C $
// Package:	SAMRAI mathops
// Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 2249 $
// Modified:	$LastChangedDate: 2008-07-03 08:17:20 -0700 (Thu, 03 Jul 2008) $
// Description:	Templated operations for real cell-centered patch data.
//

#ifndef included_math_PatchCellDataOpsReal_C
#define included_math_PatchCellDataOpsReal_C

#include "PatchCellDataOpsReal.h"
#include "tbox/Utilities.h"

namespace SAMRAI {
    namespace math {

template<int DIM, class TYPE>
PatchCellDataOpsReal<DIM,TYPE>::PatchCellDataOpsReal()
{
}

template<int DIM, class TYPE>
PatchCellDataOpsReal<DIM,TYPE>::~PatchCellDataOpsReal()
{
}

/*
*************************************************************************
*                                                                       *
* The const constructor and assignment operator are not actually used   *
* but are defined here for compilers that require an implementation for *
* every declaration.                                                    *
*                                                                       *
*************************************************************************
*/

template<int DIM, class TYPE>
PatchCellDataOpsReal<DIM,TYPE>::PatchCellDataOpsReal(
   const PatchCellDataOpsReal<DIM,TYPE>& foo)
{  
   NULL_USE(foo); 
}

template<int DIM, class TYPE>
void PatchCellDataOpsReal<DIM,TYPE>::operator=(
   const PatchCellDataOpsReal<DIM,TYPE>& foo)
{
   NULL_USE(foo); 
}

/*
*************************************************************************
*                                                                       *
* General templated operations for real cell-centered patch data.       *
*                                                                       *
*************************************************************************
*/

template<int DIM, class TYPE>
void PatchCellDataOpsReal<DIM,TYPE>::swapData(
   tbox::Pointer< hier::Patch<DIM> > patch,
   const int data1_id,
   const int data2_id) const
{
   tbox::Pointer< pdat::CellData<DIM,TYPE> > d1 = patch->getPatchData(data1_id);
   tbox::Pointer< pdat::CellData<DIM,TYPE> > d2 = patch->getPatchData(data2_id);
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!d1.isNull() && !d2.isNull());
   TBOX_ASSERT(d1->getDepth() && d2->getDepth());
   TBOX_ASSERT(d1->getBox() == d2->getBox());
   TBOX_ASSERT(d1->getGhostBox() == d2->getGhostBox());
#endif
   patch->setPatchData( data1_id, d2 );
   patch->setPatchData( data2_id, d1 );
}

template<int DIM, class TYPE>
void PatchCellDataOpsReal<DIM,TYPE>::printData(
   const tbox::Pointer< pdat::CellData<DIM,TYPE> >& data,
   const hier::Box<DIM>& box,
   std::ostream& s) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!data.isNull());
#endif
   s << "Data box = " << box << std::endl;
   data->print(box, s);
   s << "\n";
}

template<int DIM, class TYPE>
void PatchCellDataOpsReal<DIM,TYPE>::copyData(
   tbox::Pointer< pdat::CellData<DIM,TYPE> >& dst,
   const tbox::Pointer< pdat::CellData<DIM,TYPE> >& src,
   const hier::Box<DIM>& box) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!dst.isNull() && !src.isNull());
#endif
   (dst->getArrayData()).copy(src->getArrayData(), box);
}

template<int DIM, class TYPE>
void PatchCellDataOpsReal<DIM,TYPE>::setToScalar(
   tbox::Pointer< pdat::CellData<DIM,TYPE> >& dst,
   const TYPE& alpha,
   const hier::Box<DIM>& box) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!dst.isNull());
#endif
   dst->fillAll(alpha, box);
}

}
}
#endif

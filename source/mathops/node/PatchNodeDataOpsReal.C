//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-4-4/source/mathops/node/PatchNodeDataOpsReal.C $
// Package:	SAMRAI mathops
// Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 2249 $
// Modified:	$LastChangedDate: 2008-07-03 08:17:20 -0700 (Thu, 03 Jul 2008) $
// Description:	Templated operations for real node-centered patch data.
//

#ifndef included_math_PatchNodeDataOpsReal_C
#define included_math_PatchNodeDataOpsReal_C

#include "PatchNodeDataOpsReal.h"
#include "tbox/Utilities.h"
#include "NodeGeometry.h"

namespace SAMRAI {
    namespace math {

template<int DIM, class TYPE>
PatchNodeDataOpsReal<DIM,TYPE>::PatchNodeDataOpsReal()
{
}

template<int DIM, class TYPE>
PatchNodeDataOpsReal<DIM,TYPE>::~PatchNodeDataOpsReal()
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
PatchNodeDataOpsReal<DIM,TYPE>::PatchNodeDataOpsReal(
   const PatchNodeDataOpsReal<DIM,TYPE>& foo)
{
   NULL_USE(foo);
}

template<int DIM, class TYPE>
void PatchNodeDataOpsReal<DIM,TYPE>::operator=(
   const PatchNodeDataOpsReal<DIM,TYPE>& foo)
{
   NULL_USE(foo);
}

/*
*************************************************************************
*                                                                       *
* General templated operations for real node-centered patch data.       *
*                                                                       *
*************************************************************************
*/

template<int DIM, class TYPE>
void PatchNodeDataOpsReal<DIM,TYPE>::swapData(
   tbox::Pointer< hier::Patch<DIM> > patch,
   const int data1_id,
   const int data2_id) const
{
   tbox::Pointer< pdat::NodeData<DIM,TYPE> > d1 = patch->getPatchData(data1_id);
   tbox::Pointer< pdat::NodeData<DIM,TYPE> > d2 = patch->getPatchData(data2_id);
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
void PatchNodeDataOpsReal<DIM,TYPE>::printData(
   const tbox::Pointer< pdat::NodeData<DIM,TYPE> >& data,
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
void PatchNodeDataOpsReal<DIM,TYPE>::copyData(
   tbox::Pointer< pdat::NodeData<DIM,TYPE> >& dst,
   const tbox::Pointer< pdat::NodeData<DIM,TYPE> >& src,
   const hier::Box<DIM>& box) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!dst.isNull() && !src.isNull());
#endif
   const hier::Box<DIM> node_box = pdat::NodeGeometry<DIM>::toNodeBox(box);
   (dst->getArrayData()).copy(src->getArrayData(), node_box);
}

template<int DIM, class TYPE>
void PatchNodeDataOpsReal<DIM,TYPE>::setToScalar(
   tbox::Pointer< pdat::NodeData<DIM,TYPE> >& dst,
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

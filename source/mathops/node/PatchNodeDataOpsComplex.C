//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-4-4/source/mathops/node/PatchNodeDataOpsComplex.C $
// Package:	SAMRAI mathops
// Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 1917 $
// Modified:	$LastChangedDate: 2008-01-25 13:28:01 -0800 (Fri, 25 Jan 2008) $
// Description:	Operations for complex node-centered patch data.
//

#ifndef included_math_PatchNodeDataOpsComplex_C
#define included_math_PatchNodeDataOpsComplex_C

#include "PatchNodeDataOpsComplex.h"
#include "NodeGeometry.h"
#ifdef DEBUG_CHECK_ASSERTIONS
#include "tbox/Utilities.h"
#endif

namespace SAMRAI {
    namespace math {

template<int DIM>  PatchNodeDataOpsComplex<DIM>::PatchNodeDataOpsComplex()
{
}

template<int DIM>  PatchNodeDataOpsComplex<DIM>::~PatchNodeDataOpsComplex()
{
}

/*
*************************************************************************
*                                                                       *
* General operations for complex node-centered patch data.              *
*                                                                       *
*************************************************************************
*/

template<int DIM> void PatchNodeDataOpsComplex<DIM>::swapData(
   tbox::Pointer< hier::Patch<DIM> > patch,
   const int data1_id,
   const int data2_id) const
{
   tbox::Pointer< pdat::NodeData<DIM,dcomplex> > d1 = patch->getPatchData(data1_id);
   tbox::Pointer< pdat::NodeData<DIM,dcomplex> > d2 = patch->getPatchData(data2_id);
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!d1.isNull() && !d2.isNull());
   TBOX_ASSERT(d1->getDepth() && d2->getDepth());
   TBOX_ASSERT(d1->getBox() == d2->getBox());
   TBOX_ASSERT(d1->getGhostBox() == d2->getGhostBox());
#endif
   patch->setPatchData( data1_id, d2 );
   patch->setPatchData( data2_id, d1 );
}

template<int DIM> void PatchNodeDataOpsComplex<DIM>::printData(
   const tbox::Pointer< pdat::NodeData<DIM,dcomplex> >& data,
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

template<int DIM> void PatchNodeDataOpsComplex<DIM>::copyData(
   tbox::Pointer< pdat::NodeData<DIM,dcomplex> >& dst,
   const tbox::Pointer< pdat::NodeData<DIM,dcomplex> >& src,
   const hier::Box<DIM>& box) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!dst.isNull() && !src.isNull());
#endif
   const hier::Box<DIM> node_box = pdat::NodeGeometry<DIM>::toNodeBox(box);
   (dst->getArrayData()).copy(src->getArrayData(), node_box);
}

template<int DIM> void PatchNodeDataOpsComplex<DIM>::setToScalar(
   tbox::Pointer< pdat::NodeData<DIM,dcomplex> >& dst,
   const dcomplex& alpha,
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

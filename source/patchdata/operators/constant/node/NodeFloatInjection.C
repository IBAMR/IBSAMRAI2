//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-4-4/source/patchdata/operators/constant/node/NodeFloatInjection.C $
// Package:	SAMRAI patchdata
// Copyright:   (c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 1917 $
// Modified:	$LastChangedDate: 2008-01-25 13:28:01 -0800 (Fri, 25 Jan 2008) $
// Description: Constant averaging operator for node-centered float data on
//              a  mesh.
//

#ifndef included_pdat_NodeFloatInjection_C
#define included_pdat_NodeFloatInjection_C

#include "NodeFloatInjection.h"

#include<float.h>
#include<math.h>
#include "tbox/Utilities.h"
#include "Index.h"
#include "NodeData.h"
#include "NodeVariable.h"

/*
*************************************************************************
*                                                                       *
* External declarations for FORTRAN  routines.                          *
*                                                                       *
*************************************************************************
*/
extern "C" {
// in concoarsen1d.f:
   void conavgnodeflot1d_( const int&, const int&,
                           const int&, const int&,
                           const int&, const int&,
                           const int*,
                           const float*, float* );
// in concoarsen2d.f:
   void conavgnodeflot2d_( const int&, const int&, const int&, const int&,
                           const int&, const int&, const int&, const int&,
                           const int&, const int&, const int&, const int&,
                           const int*,
                           const float*, float* );
// in concoarsen3d.f:
   void conavgnodeflot3d_( const int&, const int&, const int&,
                           const int&, const int&, const int&,
                           const int&, const int&, const int&,
                           const int&, const int&, const int&,
                           const int&, const int&, const int&,
                           const int&, const int&, const int&,
                           const int*,
                           const float*, float* );
}

namespace SAMRAI {
    namespace pdat {

template<int DIM> NodeFloatInjection<DIM>::NodeFloatInjection()
: xfer::CoarsenOperator<DIM>()
{
   d_name_id = "CONSTANT_COARSEN";
}

template<int DIM> NodeFloatInjection<DIM>::~NodeFloatInjection()
{
}

template<int DIM> bool NodeFloatInjection<DIM>::findCoarsenOperator(
   const tbox::Pointer< hier::Variable<DIM> >& var,
   const std::string &op_name) const
{
   const tbox::Pointer< NodeVariable<DIM,float> > cast_var(var);
   if ( !cast_var.isNull() && (op_name == d_name_id) ) {
      return(true);
   } else {
      return(false);
   }
}

template<int DIM> const std::string&
NodeFloatInjection<DIM>::getOperatorName() const
{
   return(d_name_id);
}

template<int DIM> int NodeFloatInjection<DIM>::getOperatorPriority() const
{
   return(0);
}

template<int DIM> hier::IntVector<DIM>
NodeFloatInjection<DIM>::getStencilWidth() const {
   return(hier::IntVector<DIM>(0));
}

template<int DIM> void NodeFloatInjection<DIM>::coarsen(
   hier::Patch<DIM>& coarse,
   const hier::Patch<DIM>& fine,
   const int dst_component,
   const int src_component,
   const hier::Box<DIM>& coarse_box,
   const hier::IntVector<DIM>& ratio) const
{
   tbox::Pointer< NodeData<DIM,float> >
      fdata = fine.getPatchData(src_component);
   tbox::Pointer< NodeData<DIM,float> >
      cdata = coarse.getPatchData(dst_component);
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!fdata.isNull());
   TBOX_ASSERT(!cdata.isNull());
   TBOX_ASSERT(cdata->getDepth() == fdata->getDepth());
#endif

   const hier::Index<DIM> filo = fdata->getGhostBox().lower();
   const hier::Index<DIM> fihi = fdata->getGhostBox().upper();
   const hier::Index<DIM> cilo = cdata->getGhostBox().lower();
   const hier::Index<DIM> cihi = cdata->getGhostBox().upper();

   const hier::Index<DIM> ifirstc = coarse_box.lower();
   const hier::Index<DIM> ilastc = coarse_box.upper();

   for (int d = 0; d < cdata->getDepth(); d++) {
      if (DIM == 1) {
	 conavgnodeflot1d_(ifirstc(0),ilastc(0),
			   filo(0),fihi(0),
			   cilo(0),cihi(0),
			   ratio,
			   fdata->getPointer(d),
			   cdata->getPointer(d));
      } else if (DIM == 2) {
	 conavgnodeflot2d_(ifirstc(0),ifirstc(1),ilastc(0),ilastc(1),
			   filo(0),filo(1),fihi(0),fihi(1),
			   cilo(0),cilo(1),cihi(0),cihi(1),
			   ratio,
			   fdata->getPointer(d),
			   cdata->getPointer(d));
      } else if (DIM == 3) {
	 conavgnodeflot3d_(ifirstc(0),ifirstc(1),ifirstc(2),
			   ilastc(0),ilastc(1),ilastc(2),
			   filo(0),filo(1),filo(2),
			   fihi(0),fihi(1),fihi(2),
			   cilo(0),cilo(1),cilo(2),
			   cihi(0),cihi(1),cihi(2),
			   ratio,
			   fdata->getPointer(d),
			   cdata->getPointer(d));
      } else {
	 TBOX_ERROR("NodeFloatConstantRefine::coarsen DIM > 3 not supported" << std::endl);
      }
   }
}

}
}
#endif

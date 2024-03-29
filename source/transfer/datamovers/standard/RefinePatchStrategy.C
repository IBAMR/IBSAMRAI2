//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-4-4/source/transfer/datamovers/standard/RefinePatchStrategy.C $
// Package:	SAMRAI data transfer
// Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 1917 $
// Modified:	$LastChangedDate: 2008-01-25 13:28:01 -0800 (Fri, 25 Jan 2008) $
// Description:	Strategy interface to user routines for refining AMR data.
//

#ifndef included_xfer_RefinePatchStrategy_C
#define included_xfer_RefinePatchStrategy_C

#include "RefinePatchStrategy.h"

namespace SAMRAI {
    namespace xfer {

/*
*************************************************************************
*									*
* The default constructor and virtual destructor do nothing             *
* particularly interesting.		                                *
*									*
*************************************************************************
*/

template<int DIM>  RefinePatchStrategy<DIM>::RefinePatchStrategy()
{
}

template<int DIM>  RefinePatchStrategy<DIM>::~RefinePatchStrategy()
{
}

/*
*************************************************************************
*									*
* Loop over all fill boxes and call the user-defined preprocesses.	*
*									*
*************************************************************************
*/

template<int DIM> void RefinePatchStrategy<DIM>::preprocessRefineBoxes(
   hier::Patch<DIM>& fine,
   const hier::Patch<DIM>& coarse,
   const hier::BoxList<DIM>& fine_boxes,
   const hier::IntVector<DIM>& ratio)
{
   for (typename hier::BoxList<DIM>::Iterator b(fine_boxes); b; b++) {
      this->preprocessRefine(fine, coarse, b(), ratio);
   }
}

/*
*************************************************************************
*									*
* Loop over all fill boxes and call the user-defined postprocesses.	*
*									*
*************************************************************************
*/

template<int DIM> void RefinePatchStrategy<DIM>::postprocessRefineBoxes(
   hier::Patch<DIM>& fine,
   const hier::Patch<DIM>& coarse,
   const hier::BoxList<DIM>& fine_boxes,
   const hier::IntVector<DIM>& ratio)
{
   for (typename hier::BoxList<DIM>::Iterator b(fine_boxes); b; b++) {
      this->postprocessRefine(fine, coarse, b(), ratio);
   }
}

}
}
#endif

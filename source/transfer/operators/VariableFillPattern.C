//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/trunk/source/SAMRAI/xfer/VariableFillPattern.C $
// Package:	SAMRAI data transfer
// Copyright:	(c) 1997-2009 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 2861 $
// Modified:	$LastChangedDate: 2009-02-02 15:22:36 -0800 (Mon, 02 Feb 2009) $
// Description:	Abstract fill pattern class to provide interface for stencils
//

#ifndef included_xfer_VariableFillPattern_C
#define included_xfer_VariableFillPattern_C

#include "VariableFillPattern.h"

namespace SAMRAI {
    namespace xfer {

/*
*************************************************************************
*                                                                       *
* Default constructor                                                   *
*                                                                       *
*************************************************************************
*/
template<int DIM>
VariableFillPattern<DIM>::VariableFillPattern()
{
}

/*
*************************************************************************
*									*
* Destructor                                                            *
*									*
*************************************************************************
*/
template<int DIM>
VariableFillPattern<DIM>::~VariableFillPattern()
{
}

template<int DIM>
tbox::Pointer< hier::BoxOverlap<DIM> >
VariableFillPattern<DIM>::calculateOverlapOnLevel(const hier::BoxGeometry<DIM>& dst_geometry,
                                                  const hier::BoxGeometry<DIM>& src_geometry,
                                                  const hier::Box<DIM>& dst_patch_box,
                                                  const hier::Box<DIM>& src_mask,
                                                  const bool overwrite_interior,
                                                  const hier::IntVector<DIM>& src_offset,
                                                  const int dst_level_num,
                                                  const int src_level_num) const
{
   NULL_USE(dst_level_num);
   NULL_USE(src_level_num);
   return calculateOverlap(dst_geometry,
                           src_geometry,
                           dst_patch_box,
                           src_mask,
                           overwrite_interior,
                           src_offset);
}

template<int DIM>
void VariableFillPattern<DIM>::setTargetPatchLevelNumber(const int level_num)
{
   NULL_USE(level_num);
   return;
}

}
}
#endif

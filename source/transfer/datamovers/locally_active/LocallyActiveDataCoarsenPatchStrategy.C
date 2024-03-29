//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-4-4/source/transfer/datamovers/locally_active/LocallyActiveDataCoarsenPatchStrategy.C $
// Package:	SAMRAI data transfer
// Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 1917 $
// Modified:	$LastChangedDate: 2008-01-25 13:28:01 -0800 (Fri, 25 Jan 2008) $
// Description:	Strategy interface to user routines for coarsening locally-active AMR data.
//

#ifndef included_xfer_LocallyActiveCoarsenPatchStrategy_C
#define included_xfer_LocallyActiveCoarsenPatchStrategy_C
 
#include "LocallyActiveDataCoarsenPatchStrategy.h"

namespace SAMRAI {
    namespace xfer {

template<int DIM>
LocallyActiveDataCoarsenPatchStrategy<DIM>::LocallyActiveDataCoarsenPatchStrategy()
{
}

template<int DIM>
LocallyActiveDataCoarsenPatchStrategy<DIM>::~LocallyActiveDataCoarsenPatchStrategy()
{
}

}
}
#endif

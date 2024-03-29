/*
 * File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-4-4/source/hierarchy/dlbg/LayerNode.C $
 * Copyright:   (c) 1997-2003 Lawrence Livermore National Security, LLC
 * Revision:    $LastChangedRevision: 2155 $
 * Modified:    $LastChangedDate: 2008-04-28 09:43:00 -0700 (Mon, 28 Apr 2008) $
 * Description: Node in the distribued box graph.
 */

#ifndef included_hier_LayerNode_C
#define included_hier_LayerNode_C


#include "LayerNode.h"

#include <iostream>
#include <iomanip>

#ifdef DEBUG_NO_INLINE
#include "LayerNode.I"
#endif

namespace SAMRAI {
   namespace hier {

template<int DIM>
LayerNode<DIM>::~LayerNode()
{
   return;
}

template<int DIM>
std::ostream &operator<<( std::ostream &co, const LayerNode<DIM> &r )
{
   co << r.d_owner_rank << '#' << r.d_local_index << ':' << r.getBox();
   return co;
}




}
}
#endif

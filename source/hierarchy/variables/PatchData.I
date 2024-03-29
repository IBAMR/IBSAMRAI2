//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-4-4/source/hierarchy/variables/PatchData.I $
// Package:	SAMRAI hierarchy
// Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 1917 $
// Modified:	$LastChangedDate: 2008-01-25 13:28:01 -0800 (Fri, 25 Jan 2008) $
// Description:	Abstract base class for patch data objects
//


#ifdef DEBUG_NO_INLINE
#define inline
#endif
namespace SAMRAI {
    namespace hier {

template<int DIM> inline
const Box<DIM>& PatchData<DIM>::getBox() const
{
   return(d_box);
}

template<int DIM> inline
const Box<DIM>& PatchData<DIM>::getGhostBox() const
{
   return(d_ghost_box);
}

template<int DIM> inline
void PatchData<DIM>::setGhostBox(const Box<DIM>& ghost_box)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(ghost_box * d_box == d_box);
#endif
   d_ghost_box = ghost_box;
}

template<int DIM> inline
const IntVector<DIM>& PatchData<DIM>::getGhostCellWidth() const
{
   return(d_ghosts);
}

template<int DIM> inline
void PatchData<DIM>::setTime(const double timestamp)
{
   d_timestamp = timestamp;
}

template<int DIM> inline
double PatchData<DIM>::getTime() const
{
   return(d_timestamp);
}

#ifdef DEBUG_NO_INLINE
#undef inline
#endif
}
}

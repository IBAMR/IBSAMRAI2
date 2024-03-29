//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-4-4/source/patchdata/side/SideDataFactory.I $
// Package:	SAMRAI patch data
// Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 1917 $
// Modified:	$LastChangedDate: 2008-01-25 13:28:01 -0800 (Fri, 25 Jan 2008) $
// Description: Factory class for creating side data objects
//


#ifdef DEBUG_NO_INLINE
#define inline
#endif
namespace SAMRAI {
    namespace pdat {

template<int DIM, class TYPE> inline
int SideDataFactory<DIM,TYPE>::getDefaultDepth() const
{
   return(d_depth);
}

template<int DIM, class TYPE> inline
void SideDataFactory<DIM,TYPE>::setDefaultDepth(const int depth)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(depth > 0);
#endif
   d_depth = depth;
}

template<int DIM, class TYPE> inline
const hier::IntVector<DIM>& SideDataFactory<DIM,TYPE>::getDirectionVector() const
{
   return(d_directions);
}

template<int DIM, class TYPE> inline
hier::MultiblockDataTranslator<DIM>*
SideDataFactory<DIM,TYPE>::getMultiblockDataTranslator()
{
   if (d_mb_trans == NULL) {
      d_mb_trans = new MultiblockSideDataTranslator<DIM,TYPE>();
   }
   return d_mb_trans;
}


#ifdef DEBUG_NO_INLINE
#undef inline
#endif
}
}

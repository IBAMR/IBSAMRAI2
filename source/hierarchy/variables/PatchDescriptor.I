//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-4-4/source/hierarchy/variables/PatchDescriptor.I $
// Package:	SAMRAI hierarchy
// Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 1917 $
// Modified:	$LastChangedDate: 2008-01-25 13:28:01 -0800 (Fri, 25 Jan 2008) $
// Description:	Factory class for patch data objects that live on a patch
//


#ifdef DEBUG_NO_INLINE
#define inline
#endif
namespace SAMRAI {
    namespace hier {

template<int DIM> inline
tbox::Pointer< PatchDataFactory<DIM> >
PatchDescriptor<DIM>::getPatchDataFactory(const int id) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT( (id >= 0) && (id < d_max_number_registered_components) );
#endif
   return(d_factories[id]);
}

template<int DIM> inline
int PatchDescriptor<DIM>::getMaxNumberRegisteredComponents() const
{
   return(d_max_number_registered_components);
}

template<int DIM> inline
const std::string& PatchDescriptor<DIM>::mapIndexToName(const int id) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT( (id >= 0) && (id < d_max_number_registered_components) );
#endif
   return(d_names[id]);
}

#ifdef DEBUG_NO_INLINE
#undef inline
#endif
}
}

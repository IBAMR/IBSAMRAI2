//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-4-4/source/patchdata/outeredge/OuteredgeVariable.C $
// Package:	SAMRAI patch data
// Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 1973 $
// Modified:	$LastChangedDate: 2008-02-11 16:39:15 -0800 (Mon, 11 Feb 2008) $
// Description:	Variable class for defining outeredge centered variables
//

#ifndef included_pdat_OuteredgeVariable_C
#define included_pdat_OuteredgeVariable_C

#include "OuteredgeVariable.h"
#include "OuteredgeDataFactory.h"
#include "tbox/Utilities.h"

#ifndef NULL
#define NULL (0)
#endif

namespace SAMRAI {
    namespace pdat {

/*
*************************************************************************
*                                                                       *
* Constructor and destructor for side variable objects                  *
*                                                                       *
*************************************************************************
*/

template <int DIM, class TYPE>
OuteredgeVariable<DIM,TYPE>::OuteredgeVariable(
   const std::string &name, int depth)
:  hier::Variable<DIM>(name, new OuteredgeDataFactory<DIM,TYPE>(depth))
{
}

template <int DIM, class TYPE>
OuteredgeVariable<DIM,TYPE>::~OuteredgeVariable()
{
}

/*
*************************************************************************
*                                                                       *
* These are private and should not be used.  They are defined here      *
* because some template instantiation methods fail if some member       *
* functions are left undefined.                                         *
*                                                                       *
*************************************************************************
*/

template <int DIM, class TYPE>
OuteredgeVariable<DIM,TYPE>::OuteredgeVariable(
   const OuteredgeVariable<DIM,TYPE>& foo)
:  hier::Variable<DIM>(NULL, NULL)
{
   NULL_USE(foo);
}

template <int DIM, class TYPE>
void OuteredgeVariable<DIM,TYPE>::operator=(
   const OuteredgeVariable<DIM,TYPE>& foo)
{
   NULL_USE(foo);
}

}
}
#endif


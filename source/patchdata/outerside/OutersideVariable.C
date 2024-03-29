//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-4-4/source/patchdata/outerside/OutersideVariable.C $
// Package:	SAMRAI patch data
// Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 1917 $
// Modified:	$LastChangedDate: 2008-01-25 13:28:01 -0800 (Fri, 25 Jan 2008) $
// Description:	hier::Variable class for defining outerside centered variables
//

#ifndef included_pdat_OutersideVariable_C
#define included_pdat_OutersideVariable_C

#include "OutersideVariable.h"
#include "OutersideDataFactory.h"
#include "tbox/Utilities.h"

#ifndef NULL
#define NULL (0)
#endif

namespace SAMRAI {
    namespace pdat {

/*
*************************************************************************
*									*
* Constructor and destructor for side variable objects			*
*									*
*************************************************************************
*/

template<int DIM, class TYPE>
OutersideVariable<DIM,TYPE>::OutersideVariable(
   const std::string &name, int depth)
:  hier::Variable<DIM>(name, new OutersideDataFactory<DIM,TYPE>(depth))
{
}

template<int DIM, class TYPE>
OutersideVariable<DIM,TYPE>::~OutersideVariable()
{
}

/*
*************************************************************************
*									*
* These are private and should not be used.  They are defined here	*
* because some template instantiation methods fail if some member	*
* functions are left undefined.						*
*									*
*************************************************************************
*/

template<int DIM, class TYPE>
OutersideVariable<DIM,TYPE>::OutersideVariable(
   const OutersideVariable<DIM,TYPE>& foo)
:  hier::Variable<DIM>(NULL, NULL)
{
   NULL_USE(foo);
}

template<int DIM, class TYPE>
void OutersideVariable<DIM,TYPE>::operator=(
   const OutersideVariable<DIM,TYPE>& foo)
{
   NULL_USE(foo);
}

}
}
#endif

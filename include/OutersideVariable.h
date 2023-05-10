//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-4-4/source/patchdata/outerside/OutersideVariable.h $
// Package:	SAMRAI patch data
// Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 2132 $
// Modified:	$LastChangedDate: 2008-04-14 14:51:47 -0700 (Mon, 14 Apr 2008) $
// Description:	hier::Variable class for defining outerside centered variables
//

#ifndef included_pdat_OutersideVariable
#define included_pdat_OutersideVariable

#include "SAMRAI_config.h"
#include "tbox/Complex.h"
#ifndef included_String
#include <string>
#define included_String
#endif
#include "Variable.h"

namespace SAMRAI {
    namespace pdat {

/*!
 * @brief Class OutersideVariable<DIM> is a templated variable class
 * used to define side-centered data quantities only on patch boundaries.
 * It is a subclass of hier::Variable and is templated on the type
 * of the underlying data (e.g., double, int, bool, etc.).
 *
 * Note that the data layout in the outerside data arrays matches the corresponding 
 * array sections provided by the side data implementation. See header file for 
 * the OutersideData<DIM> class for a more detailed description of the data layout.
 *
 * IMPORTANT: The class OuterfaceVariable<DIM> and associated "outerface
 * data" classes define the same storage as this outerside variable class,
 * except that the individual array indices are permuted in the outerface
 * data type.
 *
 * @see pdat::SideData
 * @see pdat::OutersideData
 * @see pdat::OutersideDataFactory
 * @see hier::Variable
 */

template<int DIM, class TYPE>
class OutersideVariable : public hier::Variable<DIM>
{
public:
   /*!
    * @brief Create an outerside variable object having properties
    * specified by the name and depth (i.e., number of data values
    * at each index location).  The default depth is one.
    *
    * Note that The ghost cell width for all outerside data is currently
    * fixed at zero; this may be changed in the future if needed.
    */
   OutersideVariable(const std::string &name,
                     int depth = 1);

   /*!
    * @brief Virtual destructor for outerside variable objects.
    */
   virtual ~OutersideVariable();

   /*!
    * @brief Return a boolean true value indicating that fine patch
    * values take precedence on coarse-fine interfaces.
    */
   bool fineBoundaryRepresentsVariable() const {return true;}

   /*!
    * @brief Return true indicating that outerside data
    * exists on the patch boundary.
    */
   bool dataLivesOnPatchBorder() const {return true;}

private:
   // neither of the following functions are implemented
   OutersideVariable(const OutersideVariable<DIM,TYPE>&);
   void operator=(const OutersideVariable<DIM,TYPE>&);

};

}
}
#endif

#ifdef INCLUDE_TEMPLATE_IMPLEMENTATION
#include "OutersideVariable.C"
#endif

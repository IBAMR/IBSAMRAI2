//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-4-4/source/patchdata/edge/EdgeVariable.h $
// Package:	SAMRAI patch data
// Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 2132 $
// Modified:	$LastChangedDate: 2008-04-14 14:51:47 -0700 (Mon, 14 Apr 2008) $
// Description:	hier::Variable class for defining edge centered variables
//

#ifndef included_pdat_EdgeVariable
#define included_pdat_EdgeVariable

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
 * Class EdgeVariable<DIM> is a templated variable class used to define
 * edge-centered quantities on an AMR mesh.   It is a subclass of
 * hier::Variable and is templated on the type of the underlying data
 * (e.g., double, int, bool, etc.).
 *
 * See header file for EdgeData<DIM> class for a more detailed
 * description of the data layout.
 *
 * @see pdat::EdgeData
 * @see pdat::EdgeDataFactory
 * @see hier::Variable
 */

template<int DIM, class TYPE>
class EdgeVariable : public hier::Variable<DIM>
{
public:
   /*!
    * @brief Create an edge-centered variable object with the given name and
    * depth (i.e., number of data values at each edge index location).
    * A default depth of one is provided.   The fine boundary representation
    * boolean argument indicates which values (either coarse or fine) take
    * precedence at coarse-fine mesh boundaries during coarsen and refine  
    * operations.  The default is that fine data values take precedence  
    * on coarse-fine interfaces.
    */
   EdgeVariable(const std::string &name,
                int depth = 1,
                bool fine_boundary_represents_var = true);

   /*!
    * @brief Virtual destructor for edge variable objects.
    */
   virtual ~EdgeVariable<DIM,TYPE>();

   /*!
    * @brief Return boolean indicating which edge data values (coarse
    * or fine) take precedence at coarse-fine mesh interfaces.  The
    * value is set in the constructor.
    */
   bool fineBoundaryRepresentsVariable() const
       {return d_fine_boundary_represents_var;}
 
   /*!
    * @brief Return true indicating that edge data on a patch interior
    * exists on the patch boundary.
    */
   bool dataLivesOnPatchBorder() const {return true;}

private:
   bool d_fine_boundary_represents_var; 
 
   EdgeVariable(const EdgeVariable<DIM,TYPE>&);// not implemented
   void operator=(const EdgeVariable<DIM,TYPE>&);	// not implemented

};

}
}
#endif

#ifdef INCLUDE_TEMPLATE_IMPLEMENTATION
#include "EdgeVariable.C"
#endif

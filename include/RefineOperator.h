//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-4-4/source/transfer/operators/RefineOperator.h $
// Package:	SAMRAI transfer 
// Copyright:   (c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 2132 $
// Modified:	$LastChangedDate: 2008-04-14 14:51:47 -0700 (Mon, 14 Apr 2008) $
// Description: Abstract base class for spatial refinement operators.
//

#ifndef included_xfer_RefineOperator
#define included_xfer_RefineOperator

#include "SAMRAI_config.h"
#include "Box.h"
#include "IntVector.h"
#include "Patch.h"
#include "Variable.h"
#include "tbox/DescribedClass.h"
#ifndef included_String
#include <string>
#define included_String
#endif

namespace SAMRAI {
    namespace xfer {

/**
 * Class RefineOperator<DIM> is an abstract base class for each
 * spatial refinement operator used in the SAMRAI framework.  This class
 * defines the interface between numerical refinement routines and the
 * rest of the framework.  Each concrete refinement operator subclass
 * must provide three four operations:
 * 


 * - \b (1) {an implementation of the refinement operation
 *            appropriate for its corresponding patch data type.}
 * - \b (2) {a function that determines whether or not the operator
 *             matches an arbitrary request for a refinement operator.}
 * - \b (3) {a function that returns the stencil width of the operator
 *             (i.e., the number of ghost cells needed by the operator).}
 * - \b (4) {a function that returns an integer stating the priority of the
 *             operator with respect to other refinement operators.}
 * 


 * To add a new refinement operator (either for a new patch data type
 * or for a new time refinement routine on an existing type), define
 * the operator by inheriting from this abstract base class.  The operator
 * subclass must implement the refinement operation in the refine()
 * fnction, and provide a response to a general operator request in the
 * findRefineOperator() function.  The stencil width and operator priority
 * must be returned from the getStencilWidth() and getOperatorPriority()
 * functions, respectively.  Then, the new operator must be added to the
 * operator list for the appropriate transfer geometry object using the
 * Geometry<DIM>::addSpatialRefineOperator() function.
 *
 * Since spatial refinement operators usually depend on patch data centering
 * and data type as well as the mesh coordinate system, they are defined
 * in the <EM>geometry</EM> package.
 *
 * @see xfer::Geometry
 */

template<int DIM> class RefineOperator : public tbox::DescribedClass
{
public:
   /**
    * The default constructor for the refinement operator does
    * nothing interesting.
    */
   RefineOperator();

   /**
    * The virtual destructor for the refinement operator does
    * nothing interesting.
    */
   virtual ~RefineOperator();

   /**
    * Return true if the refinement operation matches the variable and
    * name string identifier request; false, otherwise.
    */
   virtual bool findRefineOperator(
      const tbox::Pointer< hier::Variable<DIM> >& var,
      const std::string &op_name) const = 0;

   /**
    * Return name string identifier of the refinement operation.
    */
   virtual const std::string& getOperatorName() const = 0;

   /**
    * Return the priority of this operator relative to other refinement
    * operators.  The SAMRAI transfer routines guarantee that refinement
    * using operators with lower priority will be performed before those
    * with higher priority.
    */
   virtual int getOperatorPriority() const = 0;

   /**
    * Return the stencil width associated with the refinement operator.
    * The SAMRAI transfer routines guarantee that the source patch will
    * contain sufficient ghost cell data surrounding the interior to
    * satisfy the stencil width requirements for each refinement operator.
    */
   virtual hier::IntVector<DIM> getStencilWidth() const = 0;

   /**
    * Refine the source component on the coarse patch to the destination
    * component on the fine patch. The refinement operation is performed
    * on the intersection of the destination patch and the fine box.
    * The coarse patch is guaranteed to contain sufficient data for the
    * stencil width of the refinement operator.
    */
   virtual void refine(
      hier::Patch<DIM>& fine,
      const hier::Patch<DIM>& coarse,
      const int dst_component,
      const int src_component,
      const hier::Box<DIM>& fine_box,
      const hier::IntVector<DIM>& ratio) const = 0;

private:
   RefineOperator(const RefineOperator<DIM>&);	// not implemented
   void operator=(const RefineOperator<DIM>&);		// not implemented

};

}
}
#endif

#ifdef INCLUDE_TEMPLATE_IMPLEMENTATION
#include "RefineOperator.C"
#endif

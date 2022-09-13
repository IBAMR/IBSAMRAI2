//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-4-4/source/patchdata/operators/constant/node/NodeFloatInjection.h $
// Package:	SAMRAI patchdata
// Copyright:   (c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 2132 $
// Modified:	$LastChangedDate: 2008-04-14 14:51:47 -0700 (Mon, 14 Apr 2008) $
// Description: Constant averaging operator for node-centered float data on
//              a  mesh.
//

#ifndef included_pdat_NodeFloatInjection
#define included_pdat_NodeFloatInjection

#include "SAMRAI_config.h"
#include "Box.h"
#include "IntVector.h"
#include "Patch.h"
#include "tbox/Pointer.h"
#ifndef included_String
#include <string>
#define included_String
#endif
#include "CoarsenOperator.h"

namespace SAMRAI {
    namespace pdat {

/**
 * Class NodeFloatInjection<DIM> implements constant
 * averaging (i.e., injection) for node-centered float patch data defined
 * over a  mesh.  It is derived from the xfer::CoarsenOperator<DIM> base
 * class.  The numerical operations for theaveraging use FORTRAN numerical
 * routines.
 *
 * The findCoarsenOperator() operator function returns true if the input
 * variable is node-centered float, and the string is "CONSTANT_COARSEN".
 *
 * @see xfer::CoarsenOperator
 */

template<int DIM> class NodeFloatInjection
: public xfer::CoarsenOperator<DIM>
{
public:
   /**
    * Uninteresting default constructor.
    */
   NodeFloatInjection();

   /**
    * Uninteresting virtual destructor.
    */
   virtual ~NodeFloatInjection();

   /**
    * Return true if the variable and name string match the node-centered
    * constant averaging; otherwise, return false.
    */
   bool findCoarsenOperator(const tbox::Pointer< hier::Variable<DIM> >& var,
                            const std::string &op_name) const;

   /**
    * Return name string identifier of this coarsening operator.
    */
   const std::string& getOperatorName() const;

   /**
    * The priority of node-centered constant averaging is 0.
    * It will be performed before any user-defined coarsen operations.
    */
   int getOperatorPriority() const;

   /**
    * The stencil width of the constant averaging operator is the vector of
    * zeros.  That is, its stencil does not extend outside the fine box.
    */
   hier::IntVector<DIM> getStencilWidth() const;

   /**
    * Coarsen the source component on the fine patch to the destination
    * component on the coarse patch using the node-centered float constant
    * averaging operator.  Coarsening is performed on the intersection of
    * the destination patch and the coarse box.  It is assumed that the
    * fine patch contains sufficient data for the stencil width of the
    * coarsening operator.
    */
   void coarsen(hier::Patch<DIM>& coarse,
                const hier::Patch<DIM>& fine,
                const int dst_component,
                const int src_component,
                const hier::Box<DIM>& coarse_box,
                const hier::IntVector<DIM>& ratio) const;

private:
   std::string d_name_id;

};

}
}
#endif

#ifdef INCLUDE_TEMPLATE_IMPLEMENTATION
#include "NodeFloatInjection.C"
#endif

//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-4-4/source/mathops/cell/PatchCellDataOpsReal.h $
// Package:	SAMRAI mathops
// Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 2132 $
// Modified:	$LastChangedDate: 2008-04-14 14:51:47 -0700 (Mon, 14 Apr 2008) $
// Description:	Templated operations for real cell-centered patch data.
//

#ifndef included_math_PatchCellDataOpsReal
#define included_math_PatchCellDataOpsReal

#include "SAMRAI_config.h"
#ifndef included_iostream
#define included_iostream
#include <iostream>
#endif
#include "PatchCellDataBasicOps.h"
#include "PatchCellDataMiscellaneousOpsReal.h"
#include "PatchCellDataNormOpsReal.h"
#include "Box.h"
#include "Patch.h"
#include "CellData.h"
#include "tbox/PIO.h"
#include "tbox/Pointer.h"

namespace SAMRAI {
    namespace math {

/**
 * Class PatchCellDataOpsReal<DIM> provides a collection of operations
 * to manipulate float and double numerical cell-centered patch data.  The
 * operations include basic arithmetic, norms and ordering, and assorted 
 * miscellaneous operations.  With the exception of a few basic routines, 
 * this class inherits its interface (and thus its functionality) from the 
 * base classes PatchCellDataBasicOps<DIM>, PatchCellDataNormOpsReal<DIM>,
 * and PatchCellDataMiscellaneousOpsReal<DIM> from which it is derived.  The 
 * name of each of these base classes is indicative of the set of 
 * cell-centered patch data operations that it provides.  
 *
 * Note that this templated class should only be used to instantiate 
 * objects with double or float as the template parameter.  A similar set of 
 * operations is implemented for complex and integer patch data in the classes 
 * PatchCellDataOpsComplex<DIM> and PatchCellDataOpsInteger<DIM>, 
 * respectively. 
 *
 * @see math::PatchCellDataBasicOps
 * @see math::PatchCellDataMiscellaneousOpsReal
 * @see math::PatchCellDataNormOpsReal
 */

template<int DIM, class TYPE>
class PatchCellDataOpsReal : 
   public tbox::DescribedClass,
   public PatchCellDataBasicOps<DIM,TYPE>,
   public PatchCellDataMiscellaneousOpsReal<DIM,TYPE>,
   public PatchCellDataNormOpsReal<DIM,TYPE>
{
public:
   /** 
    * Empty constructor and destructor.
    */
   PatchCellDataOpsReal();

   virtual ~PatchCellDataOpsReal();

   /**
    * Copy dst data to src data over given box.
    */
   void copyData(tbox::Pointer< pdat::CellData<DIM,TYPE> >& dst,
                 const tbox::Pointer< pdat::CellData<DIM,TYPE> >& src,
                 const hier::Box<DIM>& box) const;

   /**
    * Swap pointers for patch data objects.  Objects are checked for 
    * consistency of depth, box, and ghost box.
    */
   void swapData(tbox::Pointer< hier::Patch<DIM> > patch,
                 const int data1_id,
                 const int data2_id) const;

   /**
    * Print data entries over given box to given output stream.
    */
   void printData(const tbox::Pointer< pdat::CellData<DIM,TYPE> >& data,
                  const hier::Box<DIM>& box,
                  std::ostream& s = tbox::plog) const;

   /**
    * Initialize data to given scalar over given box.
    */
   void setToScalar(tbox::Pointer< pdat::CellData<DIM,TYPE> >& dst,
                    const TYPE& alpha,
                    const hier::Box<DIM>& box) const;

private:
   // The following are not implemented:
   PatchCellDataOpsReal(const PatchCellDataOpsReal<DIM,TYPE>&);
   void operator=(const PatchCellDataOpsReal<DIM,TYPE>&);

};

}
}
#endif

#ifdef INCLUDE_TEMPLATE_IMPLEMENTATION
#include "PatchCellDataOpsReal.C"
#endif

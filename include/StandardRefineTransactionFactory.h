//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-4-4/source/transfer/datamovers/standard/StandardRefineTransactionFactory.h $
// Package:	SAMRAI transfer
// Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 2132 $
// Modified:	$LastChangedDate: 2008-04-14 14:51:47 -0700 (Mon, 14 Apr 2008) $
// Description:	Concrete factory to create standard copy and time transactions 
//              for refine schedules. 
//

#ifndef included_xfer_StandardRefineTransactionFactory
#define included_xfer_StandardRefineTransactionFactory

#include "SAMRAI_config.h"
#include "Box.h"
#include "BoxOverlap.h"
#include "PatchLevel.h"
#include "tbox/Arena.h"
#include "tbox/Pointer.h"
#include "tbox/Transaction.h"
#include "RefineClasses.h"
#include "RefineTransactionFactory.h"


namespace SAMRAI {
    namespace xfer {

/*!
 * @brief Concrete subclass of RefineTransactionFactory<DIM> base class that allocates
 * RefineCopyTransaction<DIM> and RefineTimeTransaction<DIM> objects for a
 * RefineSchedule<DIM> object.
 * 
 * @see xfer::RefineCopyTransaction
 * @see xfer::RefineTimeTransaction
 * @see xfer::RefineTransactionFactory
 */

template<int DIM> 
class StandardRefineTransactionFactory : public RefineTransactionFactory<DIM>
{
public:
   /*!
    * @brief Default constructor.
    */
   StandardRefineTransactionFactory();

   /*!
    * @brief Virtual destructor.
    */
   virtual ~StandardRefineTransactionFactory();

   /*!
    * @brief Set the array of RefineClass::Data items used by the transactions.
    */
   void setRefineItems(const typename RefineClasses<DIM>::Data** refine_items,
                       int num_refine_items);

   /*!
    * @brief Clear the array of RefineClass::Data items used by the transactions.
    */
   void unsetRefineItems();

   /*!
    * @brief Set simulation time used by the refine time transaction objects.
    */
   void setTransactionTime(double fill_time);

   /*!
    * @brief Allocate an appropriate refine copy or time transaction object.
    * When time interpolation flag is passed as true a RefineTimeTransaction<DIM>
    * object will be created.  Otherwise, a RefineCopyTransaction<DIM> aill be
    * created.
    *
    * @param dst_level      tbox::Pointer to destination patch level.
    * @param src_level      tbox::Pointer to source patch level.
    * @param overlap        tbox::Pointer to overlap region between patches.
    * @param dst_patch_id   Integer index of destination patch in destination
    *                       patch level.
    * @param src_patch_id   Integer index of source patch in source patch level.
    * @param ritem_id       Integer index of RefineClass::Data item associated
    *                       with transaction.
    * @param box            Optional const reference to box defining region of 
    *                       refine transaction.  Default is an empty box.
    * @param use_time_interpolation  Optional boolean flag indicating whether the
    *                       refine transaction involves time interpolation. 
    *                       Default is false.
    * @param pool           Optional pointer to memory pool from which the
    *                       refine transaction may be allocated. Default is null.
    */
   tbox::Pointer<tbox::Transaction>
   allocate(tbox::Pointer< hier::PatchLevel<DIM> > dst_level,
            tbox::Pointer< hier::PatchLevel<DIM> > src_level,
            tbox::Pointer< hier::BoxOverlap<DIM> > overlap,
            int dst_patch_id,
            int src_patch_id,
            int ritem_id,
            const hier::Box<DIM>& box = hier::Box<DIM>(),
            bool use_time_interpolation = false,
            tbox::Pointer<tbox::Arena> pool = (tbox::Arena *) NULL) const;

private:
   // The following two functions are not implemented
   StandardRefineTransactionFactory(const StandardRefineTransactionFactory<DIM>&);
   void operator=(const StandardRefineTransactionFactory<DIM>&);

   const typename RefineClasses<DIM>::Data** d_refine_items;
   int d_num_refine_items;

};

}
}
#endif

#ifdef INCLUDE_TEMPLATE_IMPLEMENTATION
#include "StandardRefineTransactionFactory.C"
#endif

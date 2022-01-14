//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-4-4/source/transfer/datamovers/locally_active/StandardLocallyActiveDataCoarsenTransactionFactory.h $
// Package:	SAMRAI transfer
// Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 2195 $
// Modified:	$LastChangedDate: 2008-05-14 11:33:30 -0700 (Wed, 14 May 2008) $
// Description:	Concrete factory for create standard copy transactions 
//              for locally-active data coarsen schedules. 
//

#ifndef included_xfer_StandardLocallyActiveDataCoarsenTransactionFactory
#define included_xfer_StandardLocallyActiveDataCoarsenTransactionFactory

#include "SAMRAI_config.h"
#include "Box.h"
#include "BoxOverlap.h"
#include "PatchLevel.h"
#include "tbox/Arena.h"
#include "tbox/Pointer.h"
#include "tbox/Transaction.h"
#include "CoarsenClasses.h"
#include "LocallyActiveDataCoarsenTransactionFactory.h"


namespace SAMRAI {
    namespace xfer {

/*!
 * @brief Concrete subclass of LocallyActiveDataCoarsenTransactionFactory<DIM> base 
 * class that allocates CoarsenCopyTransaction<DIM> objects for a 
 * LocallyActiveDataCoarsenSchedule<DIM> object.
 * 
 * @see xfer::CoarsenCopyTransaction
 * @see xfer::LocallyActiveDataCoarsenTransactionFactory
 */

template<int DIM> 
class StandardLocallyActiveDataCoarsenTransactionFactory : 
   public LocallyActiveDataCoarsenTransactionFactory<DIM>
{
public:
   /*!
    * @brief Default constructor.
    */
   StandardLocallyActiveDataCoarsenTransactionFactory();

   /*!
    * @brief Virtual destructor.
    */
   virtual ~StandardLocallyActiveDataCoarsenTransactionFactory<DIM>();

   /*!
    * @brief Set the array of CoarsenClass::Data items used by the transactions.
    */
   void setCoarsenItems(const typename CoarsenClasses<DIM>::Data** coarsen_items,
                        int num_coarsen_items);

   /*!
    * @brief Clear the array of CoarsenClass::Data items used by the transactions.
    */
   void unsetCoarsenItems();

   /*!
    * @brief Allocate a CoarsenCopyTransaction<DIM> object.
    *
    * @param dst_level      tbox::Pointer to destination patch level.
    * @param src_level      tbox::Pointer to source patch level.
    * @param overlap        tbox::Pointer to overlap region between patches.
    * @param dst_patch_id   Integer index of destination patch in destination
    *                       patch level.
    * @param src_patch_id   Integer index of source patch in source patch level.
    * @param citem_id       Integer index of CoarsenClass::Data item associated
    *                       with transaction.
    * @param pool           Optional pointer to memory pool from which the
    *                       coarsen transaction may be allocated. Default is null.
    */
   tbox::Pointer<tbox::Transaction>
   allocate(tbox::Pointer< hier::PatchLevel<DIM> > dst_level,
            tbox::Pointer< hier::PatchLevel<DIM> > src_level,
            tbox::Pointer< hier::BoxOverlap<DIM> > overlap,
            int dst_patch_id,
            int src_patch_id,
            int citem_id,
            tbox::Pointer<tbox::Arena> pool = (tbox::Arena *) NULL) const;

private:
   // The following two functions are not implemented
   StandardLocallyActiveDataCoarsenTransactionFactory(
      const StandardLocallyActiveDataCoarsenTransactionFactory<DIM>&);
   void operator=(const StandardLocallyActiveDataCoarsenTransactionFactory<DIM>&);

   const typename xfer::CoarsenClasses<DIM>::Data** d_coarsen_items;
   int d_num_coarsen_items;

};

}
}

#ifdef INCLUDE_TEMPLATE_IMPLEMENTATION
#include "StandardLocallyActiveDataCoarsenTransactionFactory.C"
#endif

#endif


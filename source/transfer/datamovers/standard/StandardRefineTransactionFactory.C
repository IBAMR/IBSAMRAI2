//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-4-4/source/transfer/datamovers/standard/StandardRefineTransactionFactory.C $
// Package:	SAMRAI transfer
// Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 1917 $
// Modified:	$LastChangedDate: 2008-01-25 13:28:01 -0800 (Fri, 25 Jan 2008) $
// Description:	Concrete factory to create standard copy and time transactions
//              for refine schedules.
//

#ifndef included_xfer_StandardRefineTransactionFactory_C
#define included_xfer_StandardRefineTransactionFactory_C

#include "StandardRefineTransactionFactory.h"

#include "RefineCopyTransaction.h"
#include "RefineTimeTransaction.h"
#include "tbox/ArenaManager.h"

namespace SAMRAI {
    namespace xfer {

/*
*************************************************************************
*                                                                       *
* Default constructor and destructor.                                   * 
*                                                                       *
*************************************************************************
*/

template<int DIM> 
StandardRefineTransactionFactory<DIM>::StandardRefineTransactionFactory()
{
}

template<int DIM> 
StandardRefineTransactionFactory<DIM>::~StandardRefineTransactionFactory()
{
}

/*
*************************************************************************
*                                                                       *
* Set/unset information for transactions managed by this factory class. *
*                                                                       *
*************************************************************************
*/

template<int DIM>
void StandardRefineTransactionFactory<DIM>::setRefineItems(
   const typename RefineClasses<DIM>::Data** refine_items,
   int num_refine_items)
{
   xfer::RefineCopyTransaction<DIM>::setRefineItems(refine_items, 
                                                    num_refine_items);
   xfer::RefineTimeTransaction<DIM>::setRefineItems(refine_items, 
                                                    num_refine_items);
   d_refine_items = refine_items;
   d_num_refine_items = num_refine_items;
}

template<int DIM>
void StandardRefineTransactionFactory<DIM>::unsetRefineItems()
{
   xfer::RefineCopyTransaction<DIM>::unsetRefineItems();
   xfer::RefineTimeTransaction<DIM>::unsetRefineItems();
   d_refine_items = (const typename xfer::RefineClasses<DIM>::Data**)NULL;
   d_num_refine_items = 0;
}

template<int DIM>
void StandardRefineTransactionFactory<DIM>::setTransactionTime(
   double fill_time)
{
   xfer::RefineTimeTransaction<DIM>::setTransactionTime(fill_time);
}

/*
*************************************************************************
*                                                                       *
* Allocate appropriate transaction object.                              *
*                                                                       *
*************************************************************************
*/

template<int DIM>
tbox::Pointer<tbox::Transaction>
StandardRefineTransactionFactory<DIM>::allocate(
   tbox::Pointer< hier::PatchLevel<DIM> > dst_level,
   tbox::Pointer< hier::PatchLevel<DIM> > src_level,
   tbox::Pointer< hier::BoxOverlap<DIM> > overlap,
   int dst_patch_id,
   int src_patch_id,
   int ritem_id,
   const hier::Box<DIM>& box,
   bool use_time_interpolation,
   tbox::Pointer<tbox::Arena> pool) const
{
   if (pool.isNull()) {
      pool = tbox::ArenaManager::getManager()->getStandardAllocator();
   }

   if (use_time_interpolation) {

      RefineTimeTransaction<DIM>* transaction =
         new (pool) RefineTimeTransaction<DIM>(dst_level, src_level,
                                               overlap,
                                               dst_patch_id, 
                                               src_patch_id,
                                               box,
                                               ritem_id);
      return(tbox::Pointer<tbox::Transaction>(transaction, pool));

   } else {
  
      RefineCopyTransaction<DIM>* transaction =
         new (pool) RefineCopyTransaction<DIM>(dst_level, src_level,
                                               overlap,
                                               dst_patch_id,
                                               src_patch_id,
                                               ritem_id);
      return(tbox::Pointer<tbox::Transaction>(transaction, pool));

   }
}

}
}
#endif

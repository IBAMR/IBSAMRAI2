//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-4-4/source/transfer/datamovers/standard/CoarsenCopyTransaction.h $
// Package:	SAMRAI data transfer
// Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 2132 $
// Modified:	$LastChangedDate: 2008-04-14 14:51:47 -0700 (Mon, 14 Apr 2008) $
// Description:	Communication transaction for data copies during data coarsening
//
 
#ifndef included_xfer_CoarsenCopyTransaction
#define included_xfer_CoarsenCopyTransaction

#include "SAMRAI_config.h"
#ifndef included_iostream
#define included_iostream
#include <iostream>
#endif
#include "BoxOverlap.h"
#include "PatchLevel.h"
#include "tbox/Pointer.h"
#include "tbox/Transaction.h"
#include "CoarsenClasses.h"

namespace SAMRAI {
    namespace xfer {

/*!
 * @brief Class CoarsenCopyTransaction<DIM> represents a single copy communication 
 * transaction between two processors or a local data copy for coaren schedules.  
 * Note that to there is an implicit hand-shaking between objects of this class 
 * and the CoarsenSchedule<DIM> object that constructs them.  Following the coarsen
 * schedule implementation, the source patch data index for a copy transaction
 * always refers to the source data, and the destination patch data index for a copy 
 * transaction is always the destination data, all as defined in the 
 * CoarsenClasses<DIM> class.
 *
 * @see xfer::CoarsenSchedule
 * @see xfer::CoarsenClasses
 * @see tbox::Schedule
 * @see tbox::Transaction
 */

template<int DIM> class CoarsenCopyTransaction : public tbox::Transaction
{
public:
   /*!
    * Static member function to set the array of coarsen class data items that 
    * is shared by all object instances of this copy transaction class during 
    * data transfers.  The array must be set before any transactions are executed.  
    * The array is set in the CoarsenSchedule<DIM> class.
    */
   static void setCoarsenItems(const typename CoarsenClasses<DIM>::Data** coarsen_items,
                               int num_coarsen_items);

   /*!
    * Static member function to unset the array of coarsen class data items that
    * is shared by all object instances of this copy transaction class during 
    * data transfers.  The unset function is used to prevent erroneous execution 
    * of different schedules.  The array is unset in the CoarsenSchedule<DIM> class.
    */
   static void unsetCoarsenItems();

   /*!
    * Construct a transaction with the specified source and destination
    * levels, patches, and patch data components found in the coarsen class
    * item with the given id owned by the calling coarsen schedule.  In general, 
    * this constructor is called by a CoarsenSchedule<DIM> object for each data 
    * transaction (not involving time interpolation) that must occur.  This 
    * transaction will be responsible for one of the following: (1) a local data 
    * copy, (2) packing a message stream with source patch data, or (3) unpacking 
    * destination patch data from a message stream.
    * 
    * @param dst_level        tbox::Pointer to destination patch level.
    * @param src_level        tbox::Pointer to source patch level.
    * @param overlap          tbox::Pointer to overlap region between patches.
    * @param dst_patch        Integer index of destination patch in destination 
    *                         patch level. 
    * @param src_patch        Integer index of source patch in source patch level. 
    * @param coarsen_item_id  Integer id of coarsen data item owned by coarsen schedule.
    *
    * When assertion checking is active, an assertion will result if any of the pointer
    * arguments is null, or if any of the integer arguments are invalid (i.e., < 0);
    */
   CoarsenCopyTransaction(
      tbox::Pointer< hier::PatchLevel<DIM> > dst_level,
      tbox::Pointer< hier::PatchLevel<DIM> > src_level,
      tbox::Pointer< hier::BoxOverlap<DIM> > overlap,
      int dst_patch,
      int src_patch,
      int coarsen_item_id);

   /*!
    * The virtual destructor for the copy transaction releases all
    * memory associated with the transaction.
    */
   virtual ~CoarsenCopyTransaction<DIM>();

   /*!
    * Return a boolean indicating whether this transaction can estimate
    * the size of an incoming message.  If this is false, then a different
    * communication protocol kicks in and the message size is transmitted
    * between nodes.
    */
   virtual bool canEstimateIncomingMessageSize();

   /*!
    * Return the integer buffer space (in bytes) needed for the incoming message.
    * This routine is only called if the transaction can estimate the
    * size of the incoming message.  See canEstimateIncomingMessageSize().
    */
   virtual int computeIncomingMessageSize();

   /*!
    * Return the integer buffer space (in bytes) needed for the outgoing message.
    */
   virtual int computeOutgoingMessageSize();

   /*!
    * Return the sending processor number for the communications transaction.
    */
   virtual int getSourceProcessor();

   /*!
    * Return the receiving processor number for the communications transaction.
    */
   virtual int getDestinationProcessor();

   /*!
    * Pack the transaction data into the message stream.
    */
   virtual void packStream(tbox::AbstractStream& stream);

   /*!
    * Unpack the transaction data from the message stream.
    */
   virtual void unpackStream(tbox::AbstractStream& stream);

   /*!
    * Perform the local data copy for the transaction.
    */
   virtual void copyLocalData();

   /*!
    * Print out transaction information.
    */
   virtual void printClassData(std::ostream& stream) const;

private:
   CoarsenCopyTransaction(const CoarsenCopyTransaction<DIM>&); // not implemented
   void operator=(const CoarsenCopyTransaction<DIM>&);	// not implemented

   static const typename CoarsenClasses<DIM>::Data** s_coarsen_items;
   static int s_num_coarsen_items;

   tbox::Pointer< hier::PatchLevel<DIM> >    d_dst_level;
   tbox::Pointer< hier::PatchLevel<DIM> >    d_src_level;
   tbox::Pointer< hier::BoxOverlap<DIM> >    d_overlap;
   int                               d_dst_patch;
   int                               d_src_patch;
   int                               d_coarsen_item_id;
   int                               d_incoming_bytes;
   int                               d_outgoing_bytes;

};

}
}
#endif

#ifdef INCLUDE_TEMPLATE_IMPLEMENTATION
#include "CoarsenCopyTransaction.C"
#endif

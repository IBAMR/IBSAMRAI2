//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-4-4/source/toolbox/parallel/Schedule.C $
// Package:	SAMRAI toolbox
// Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 2410 $
// Modified:	$LastChangedDate: 2008-10-08 14:30:40 -0700 (Wed, 08 Oct 2008) $
// Description:	Schedule of communication transactions between processors
//

#include "tbox/Schedule.h"
#include "tbox/ShutdownRegistry.h"
#include "tbox/TimerManager.h"

#include <vector>

namespace SAMRAI {
   namespace tbox {

#define SCHEDULE_SIZE_TAG (1)
#define SCHEDULE_DATA_TAG (2)

typedef List< Pointer< Transaction > >::Iterator ITERATOR;

/*
*************************************************************************
*									*
* Initialize the communication schedule and allocate the data arrays.	*
*									*
*************************************************************************
*/

Schedule::Schedule()
{
   d_nnodes = SAMRAI_MPI::getNodes();

   d_incoming.resizeArray(d_nnodes);
   d_outgoing.resizeArray(d_nnodes);

   d_send_set.resizeArray(d_nnodes);
   d_recv_set.resizeArray(d_nnodes);
}

/*
*************************************************************************
*									*
* All of the data arrays are automatically deallocated.  Note that the	*
* destructor should not be called during a communication phase.		*
*									*
*************************************************************************
*/

Schedule::~Schedule()
{
}

/*
*************************************************************************
*									*
* Add a transaction to the head of a list of data transactions in       *
* this schedule. The assignment of the transaction to a list depends    *
* on the source and destination processors of the transaction.	        *
*									*
*************************************************************************
*/

void Schedule::addTransaction(
   const Pointer<Transaction>& transaction)
{
   const int my_id  = SAMRAI_MPI::getRank();
   const int src_id = transaction->getSourceProcessor();
   const int dst_id = transaction->getDestinationProcessor();

   if ((my_id == src_id) && (my_id == dst_id)) {
      d_local_set.addItem(transaction);
   } else {
      if (my_id == dst_id) {
         d_recv_set[src_id].addItem(transaction);
      } else if (my_id == src_id) {
         d_send_set[dst_id].addItem(transaction);
      }
   }
}

/*
*************************************************************************
*                                                                       *
* Append a transaction to the tail of a list of data transactions in    *
* this schedule.  The assignment of the transaction to a list depends   *
* on the source and destination processors of the transaction.          *
*                                                                       *
*************************************************************************
*/

void Schedule::appendTransaction(
   const Pointer<Transaction>& transaction)
{
   const int my_id  = SAMRAI_MPI::getRank();
   const int src_id = transaction->getSourceProcessor();
   const int dst_id = transaction->getDestinationProcessor();

   if ((my_id == src_id) && (my_id == dst_id)) {
      d_local_set.appendItem(transaction);
   } else {
      if (my_id == dst_id) {
         d_recv_set[src_id].appendItem(transaction);
      } else if (my_id == src_id) {
         d_send_set[dst_id].appendItem(transaction);
      }
   }
}


/*
*************************************************************************
*									*
* Perform the communication described by the schedule.			*
*									*
*************************************************************************
*/

void Schedule::communicate()
{
   SAMRAI_SETUP_TIMER_AND_SCOPE("tbox::Schedule::communicate()");
   beginCommunication();
   finalizeCommunication();
}

/*
*************************************************************************
*									*
* Begin communication but do not deliver data to the destinations.	*
* This routine calculates message sizes, posts receives, and sends	*
* outgoing messages.  Since we does not wait for message completion,	*
* use finalizeCommunication() to ensure that communication has		*
* finished.								*
*									*
*************************************************************************
*/

void Schedule::beginCommunication()
{
   SAMRAI_SETUP_TIMER_AND_SCOPE("tbox::Schedule::beginCommunication()");
   calculateSendSizes();
   calculateReceiveSizes();
   postMessageReceives();
   sendMessages();
}

/*
*************************************************************************
*									*
* Perform the local data copies, deliver the messages into their	*
* destinations, and deallocate the send buffers.			*
*									*
*************************************************************************
*/

void Schedule::finalizeCommunication()
{
   SAMRAI_SETUP_TIMER_AND_SCOPE("tbox::Schedule::finalizeCommunication()");

#ifdef HAVE_MPI
   std::vector<MPI_Request> requests;
   // Wait for both sends and recvs to complete
   for (int p = 0; p < d_nnodes; p++) {
      if (d_incoming[p].d_stream_in_use) {
         requests.push_back(d_incoming[p].d_request_id);
         d_incoming[p].d_request_id = MPI_REQUEST_NULL;
      }
      if (d_outgoing[p].d_stream_in_use) {
         requests.push_back(d_outgoing[p].d_request_id);
         d_outgoing[p].d_request_id = MPI_REQUEST_NULL;
      }
   }

   {
      SAMRAI_SETUP_TIMER_AND_SCOPE("tbox::Schedule::finalizeCommunication()[waitall]");
      int ierr = MPI_Waitall(requests.size(), requests.data(), MPI_STATUSES_IGNORE);
      TBOX_ASSERT(ierr == 0);
   }
#endif

   performLocalCopies();

   processIncomingMessages();
   deallocateSendBuffers();
}

/*
*************************************************************************
*									*
* Calculate the number of bytes to send to every other processor.	*
*									*
*************************************************************************
*/

void Schedule::calculateSendSizes()
{
   SAMRAI_SETUP_TIMER_AND_SCOPE("tbox::Schedule::calculateSendSizes()");
   for (int p = 0; p < d_nnodes; p++) {
      int bytes = 0;
      d_outgoing[p].d_must_communicate_byte_size = false;
      for (ITERATOR s(d_send_set[p]); s; s++) {
         if (!s()->canEstimateIncomingMessageSize()) {
            d_outgoing[p].d_must_communicate_byte_size = true;
         }
         bytes += s()->computeOutgoingMessageSize();
      }
      d_outgoing[p].d_bytes_in_stream = bytes;
   }
}

/*
*************************************************************************
*									*
* Estimate the number of bytes to be received from every other		*
* processor.  If we cannot estimate the bytes from available local	*
* information, then receive the number of bytes from the other		*
* processor.  We need to be careful here since a global communication	*
* may not work, since not all processors may need size information,	*
* so not all processors may enter the global exchange.  Since each	*
* processor will usually communicate with a small number of other	*
* processors, local communication is also probably more efficient	*
* than a global synchronization.					*
*									*
* The general outline of the algorithm is as follows:			*
*									*
*	(1) Walk the receive list and compute incoming byte sizes	*
*	(2) If we cannot compute incoming bytes, post a message receive	*
*	(3) Walk the send list and communicate byte sizes if needed	*
*	(4) Wait for message receives to finish				*
*									*
*************************************************************************
*/

void Schedule::calculateReceiveSizes()
{
   SAMRAI_SETUP_TIMER_AND_SCOPE("tbox::Schedule::calculateReceiveSizes()");
   /*
    * Walk the receive set list and compute the bytes to be received.
    * If the size of the receive message cannot be determined from local
    * information, then post a message receive from the sending processor.
    */
#ifdef HAVE_MPI
   int ierr = 0;
   std::vector<MPI_Request> requests;
#endif
   for (int p = 0; p < d_nnodes; p++) {
      int bytes = 0;
      bool needs_bytes = false;

      for (ITERATOR r(d_recv_set[p]); r && !needs_bytes; r++) {
         if (!r()->canEstimateIncomingMessageSize()) {
            needs_bytes = true;
         }
         bytes += r()->computeIncomingMessageSize();
      }

      if (!needs_bytes) {
         d_incoming[p].d_must_communicate_byte_size = false;
         d_incoming[p].d_bytes_in_stream            = bytes;
         d_incoming[p].d_stream_in_use              = false;
      } else {
         d_incoming[p].d_must_communicate_byte_size = true;
         d_incoming[p].d_stream_in_use              = true;
#ifdef HAVE_MPI
         SAMRAI_MPI::updateIncomingStatistics(1, sizeof(int));
         ierr = MPI_Irecv(&d_incoming[p].d_bytes_in_stream,
                          1,
                          MPI_INT,
                          p,
                          SCHEDULE_SIZE_TAG,
                          SAMRAI_MPI::getCommunicator(),
                          &d_incoming[p].d_request_id);
         TBOX_ASSERT(ierr == 0);
         requests.push_back(d_incoming[p].d_request_id);
#endif
      }
   }

   /*
    * Walk the send list and send out size messages if size data
    * is needed
    */

   for (int q = 0; q < d_nnodes; q++) {
      if (d_outgoing[q].d_must_communicate_byte_size) {
#ifdef HAVE_MPI
         SAMRAI_MPI::updateOutgoingStatistics(1, sizeof(int));
         ierr = MPI_Isend(&d_outgoing[q].d_bytes_in_stream,
                          1,
                          MPI_INT,
                          q,
                          SCHEDULE_SIZE_TAG,
                          SAMRAI_MPI::getCommunicator(),
                          &d_outgoing[q].d_request_id);
         TBOX_ASSERT(ierr == 0);
         requests.push_back(d_outgoing[q].d_request_id);
#endif
      }
   }

#ifdef HAVE_MPI
   ierr = MPI_Waitall(static_cast<int>(requests.size()),
                      requests.data(),
                      MPI_STATUSES_IGNORE);
   TBOX_ASSERT(ierr == 0);
#endif

   for (int w = 0; w < d_nnodes; w++) {
      if (d_incoming[w].d_stream_in_use) {
         d_incoming[w].d_stream_in_use = false;
      }
      if (d_outgoing[w].d_stream_in_use) {
         d_outgoing[w].d_stream_in_use = false;
      }
   }
}

/*
*************************************************************************
*									*
* Allocate receive buffers and post receives for incoming messages.	*
*									*
*************************************************************************
*/

void Schedule::postMessageReceives()
{
   SAMRAI_SETUP_TIMER_AND_SCOPE("tbox::Schedule::postMessageReceives()");
   for (int p = 0; p < d_nnodes; p++) {
      const int bytes = d_incoming[p].d_bytes_in_stream;
      if (bytes == 0) {
         d_incoming[p].d_stream_in_use = false;
      } else {
         d_incoming[p].d_stream_in_use = true;
         d_incoming[p].d_stream =
            new MessageStream(bytes, MessageStream::Read);
#ifdef HAVE_MPI
         SAMRAI_MPI::updateIncomingStatistics(1, bytes);
         int ierr = MPI_Irecv(d_incoming[p].d_stream->getBufferStart(),
                              bytes,
                              MPI_BYTE,
                              p,
                              SCHEDULE_DATA_TAG,
                              SAMRAI_MPI::getCommunicator(),
                              &d_incoming[p].d_request_id);
         TBOX_ASSERT(ierr == 0);
#endif
      }
   }
}

/*
*************************************************************************
*									*
* Allocate the send buffer, pack the data, and initiate the message	*
* sends.								*
*									*
*************************************************************************
*/

void Schedule::sendMessages()
{
   SAMRAI_SETUP_TIMER_AND_SCOPE("tbox::Schedule::sendMessages()");
   for (int p = 0; p < d_nnodes; p++) {
      const int bytes = d_outgoing[p].d_bytes_in_stream;
      if (bytes == 0) {
         d_outgoing[p].d_stream_in_use = false;
      } else {
         d_outgoing[p].d_stream_in_use = true;
         d_outgoing[p].d_stream =
            new MessageStream(bytes, MessageStream::Write);
         for (ITERATOR pack(d_send_set[p]); pack; pack++) {
            pack()->packStream(*d_outgoing[p].d_stream);
         }
#ifdef HAVE_MPI
         SAMRAI_MPI::updateOutgoingStatistics(1, bytes);
         int ierr = MPI_Isend(d_outgoing[p].d_stream->getBufferStart(),
                              d_outgoing[p].d_stream->getCurrentSize(),
                              MPI_BYTE,
                              p,
                              SCHEDULE_DATA_TAG,
                              SAMRAI_MPI::getCommunicator(),
                              &d_outgoing[p].d_request_id);
         TBOX_ASSERT(ierr == 0);
#endif
      }
   }
}

/*
*************************************************************************
*									*
* Perform all of the local memory-to-memory copies for this processor.	*
*									*
*************************************************************************
*/

void Schedule::performLocalCopies()
{
   SAMRAI_SETUP_TIMER_AND_SCOPE("tbox::Schedule::performLocalCopies()");
   for (ITERATOR local(d_local_set); local; local++) {
      local()->copyLocalData();
   }
}

/*
*************************************************************************
*									*
* Wait until a message arrives and then unpack the data.  Note that we	*
* want to unpack data is it comes in, not in processor order.		*
*									*
*************************************************************************
*/

void Schedule::processIncomingMessages()
{
   SAMRAI_SETUP_TIMER_AND_SCOPE("tbox::Schedule::processIncomingMessages()");
   for (int p = 0; p < d_nnodes; p++) {
      if (d_incoming[p].d_stream_in_use) {
         for (ITERATOR recv(d_recv_set[p]); recv; recv++) {
            recv()->unpackStream(*d_incoming[p].d_stream);
         }

         d_incoming[p].d_stream_in_use = false;
         d_incoming[p].d_stream.setNull();
      }
   }
}

/*
*************************************************************************
*									*
* Wait until all message sends have completed and deallocate buffers.	*
*									*
*************************************************************************
*/

void Schedule::deallocateSendBuffers()
{
   SAMRAI_SETUP_TIMER_AND_SCOPE("tbox::Schedule::deallocateSendBuffers()");
   for (int p = 0; p < d_nnodes; p++) {
      if (d_outgoing[p].d_stream_in_use) {
         d_outgoing[p].d_stream_in_use = false;
         d_outgoing[p].d_stream.setNull();
      }
   }
}

/*
*************************************************************************
*									*
* Print class data to the specified output stream.			*
*									*
*************************************************************************
*/

void Schedule::printClassData(std::ostream& stream) const
{
   stream << "Schedule::printClassData()" << std::endl;
   stream << "-------------------------------" << std::endl;

   stream << "Number of nodes: " << d_nnodes << std::endl;

   for (int s = 0; s < d_nnodes; s++) {
      stream << "Incoming Message Stream: " << s << std::endl;
      stream << "   bytes in stream: " << d_incoming[s].d_bytes_in_stream
              << std::endl;
      stream << "   must communicate bytes size: "
              << d_incoming[s].d_must_communicate_byte_size << std::endl;
      stream << "   stream in use: " << d_incoming[s].d_stream_in_use << std::endl;
      stream << "Outgoing Message Stream: " << s << std::endl;
      stream << "   bytes in stream: " << d_outgoing[s].d_bytes_in_stream
              << std::endl;
      stream << "   must communicate bytes size: "
              << d_outgoing[s].d_must_communicate_byte_size << std::endl;
      stream << "   stream in use: " << d_outgoing[s].d_stream_in_use << std::endl;
   }

   for (int ss = 0; ss < d_nnodes; ss++) {
      stream << "Send Set: " << ss << std::endl;
      for (ITERATOR send(d_send_set[ss]); send; send++) {
         send()->printClassData(stream);
      }
   }

   for (int rs = 0; rs < d_nnodes; rs++) {
      stream << "Receive Set: " << rs << std::endl;
      for (ITERATOR recv(d_recv_set[rs]); recv; recv++) {
         recv()->printClassData(stream);
      }
   }

   stream << "Local Set" << std::endl;
   for (ITERATOR local(d_local_set); local; local++) {
      local()->printClassData(stream);
   }
}

}
}

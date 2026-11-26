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
         d_peers[src_id].d_recv_set.addItem(transaction);
      } else if (my_id == src_id) {
         d_peers[dst_id].d_send_set.addItem(transaction);
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
         d_peers[src_id].d_recv_set.appendItem(transaction);
      } else if (my_id == src_id) {
         d_peers[dst_id].d_send_set.appendItem(transaction);
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
   for (PeerMap::iterator peer = d_peers.begin(); peer != d_peers.end();
        ++peer) {
      ScheduleMessageStream& incoming = peer->second.d_incoming;
      ScheduleMessageStream& outgoing = peer->second.d_outgoing;
      if (incoming.d_stream_in_use) {
         requests.push_back(incoming.d_request_id);
         incoming.d_request_id = MPI_REQUEST_NULL;
      }
      if (outgoing.d_stream_in_use) {
         requests.push_back(outgoing.d_request_id);
         outgoing.d_request_id = MPI_REQUEST_NULL;
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
   for (PeerMap::iterator peer = d_peers.begin(); peer != d_peers.end();
        ++peer) {
      ScheduleMessageStream& outgoing = peer->second.d_outgoing;
      int bytes = 0;
      outgoing.d_must_communicate_byte_size = false;
      for (ITERATOR s(peer->second.d_send_set); s; s++) {
         if (!s()->canEstimateIncomingMessageSize()) {
            outgoing.d_must_communicate_byte_size = true;
         }
         bytes += s()->computeOutgoingMessageSize();
      }
      outgoing.d_bytes_in_stream = bytes;
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
   for (PeerMap::iterator peer = d_peers.begin(); peer != d_peers.end();
        ++peer) {
      const int p = peer->first;
      ScheduleMessageStream& incoming = peer->second.d_incoming;
      int bytes = 0;
      bool needs_bytes = false;

      for (ITERATOR r(peer->second.d_recv_set); r && !needs_bytes; r++) {
         if (!r()->canEstimateIncomingMessageSize()) {
            needs_bytes = true;
         }
         bytes += r()->computeIncomingMessageSize();
      }

      if (!needs_bytes) {
         incoming.d_must_communicate_byte_size = false;
         incoming.d_bytes_in_stream            = bytes;
         incoming.d_stream_in_use              = false;
      } else {
         incoming.d_must_communicate_byte_size = true;
         incoming.d_stream_in_use              = true;
#ifdef HAVE_MPI
         SAMRAI_MPI::updateIncomingStatistics(1, sizeof(int));
         ierr = MPI_Irecv(&incoming.d_bytes_in_stream,
                          1,
                          MPI_INT,
                          p,
                          SCHEDULE_SIZE_TAG,
                          SAMRAI_MPI::getCommunicator(),
                          &incoming.d_request_id);
         TBOX_ASSERT(ierr == 0);
         requests.push_back(incoming.d_request_id);
#endif
      }
   }

   /*
    * Walk the send list and send out size messages if size data
    * is needed
    */

   for (PeerMap::iterator peer = d_peers.begin(); peer != d_peers.end();
        ++peer) {
      const int q = peer->first;
      ScheduleMessageStream& outgoing = peer->second.d_outgoing;
      if (outgoing.d_must_communicate_byte_size) {
#ifdef HAVE_MPI
         SAMRAI_MPI::updateOutgoingStatistics(1, sizeof(int));
         ierr = MPI_Isend(&outgoing.d_bytes_in_stream,
                          1,
                          MPI_INT,
                          q,
                          SCHEDULE_SIZE_TAG,
                          SAMRAI_MPI::getCommunicator(),
                          &outgoing.d_request_id);
         TBOX_ASSERT(ierr == 0);
         requests.push_back(outgoing.d_request_id);
#endif
      }
   }

#ifdef HAVE_MPI
   ierr = MPI_Waitall(static_cast<int>(requests.size()),
                      requests.data(),
                      MPI_STATUSES_IGNORE);
   TBOX_ASSERT(ierr == 0);
#endif

   for (PeerMap::iterator peer = d_peers.begin(); peer != d_peers.end();
        ++peer) {
      peer->second.d_incoming.d_stream_in_use = false;
      peer->second.d_outgoing.d_stream_in_use = false;
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
   for (PeerMap::iterator peer = d_peers.begin(); peer != d_peers.end();
        ++peer) {
      const int p = peer->first;
      ScheduleMessageStream& incoming = peer->second.d_incoming;
      const int bytes = incoming.d_bytes_in_stream;
      if (bytes == 0) {
         incoming.d_stream_in_use = false;
      } else {
         incoming.d_stream_in_use = true;
         incoming.d_stream = MessageStream(bytes, MessageStream::Read);
#ifdef HAVE_MPI
         SAMRAI_MPI::updateIncomingStatistics(1, bytes);
         int ierr = MPI_Irecv(incoming.d_stream.getBufferStart(),
                              bytes,
                              MPI_BYTE,
                              p,
                              SCHEDULE_DATA_TAG,
                              SAMRAI_MPI::getCommunicator(),
                              &incoming.d_request_id);
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
   for (PeerMap::iterator peer = d_peers.begin(); peer != d_peers.end();
        ++peer) {
      const int p = peer->first;
      ScheduleMessageStream& outgoing = peer->second.d_outgoing;
      const int bytes = outgoing.d_bytes_in_stream;
      if (bytes == 0) {
         outgoing.d_stream_in_use = false;
      } else {
         outgoing.d_stream_in_use = true;
         outgoing.d_stream = MessageStream(bytes, MessageStream::Write);
         for (ITERATOR pack(peer->second.d_send_set); pack; pack++) {
            pack()->packStream(outgoing.d_stream);
         }
#ifdef HAVE_MPI
         SAMRAI_MPI::updateOutgoingStatistics(1, bytes);
         int ierr = MPI_Isend(outgoing.d_stream.getBufferStart(),
                              outgoing.d_stream.getCurrentSize(),
                              MPI_BYTE,
                              p,
                              SCHEDULE_DATA_TAG,
                              SAMRAI_MPI::getCommunicator(),
                              &outgoing.d_request_id);
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
   for (PeerMap::iterator peer = d_peers.begin(); peer != d_peers.end();
        ++peer) {
      ScheduleMessageStream& incoming = peer->second.d_incoming;
      if (incoming.d_stream_in_use) {
         for (ITERATOR recv(peer->second.d_recv_set); recv; recv++) {
            recv()->unpackStream(incoming.d_stream);
         }

         incoming.d_stream_in_use = false;
         incoming.d_stream = {};
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
   for (PeerMap::iterator peer = d_peers.begin(); peer != d_peers.end();
        ++peer) {
      ScheduleMessageStream& outgoing = peer->second.d_outgoing;
      if (outgoing.d_stream_in_use) {
         outgoing.d_stream_in_use = false;
         outgoing.d_stream = {};
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

   stream << "Number of peers: " << d_peers.size() << std::endl;

   for (PeerMap::const_iterator peer = d_peers.begin(); peer != d_peers.end();
        ++peer) {
      const int s = peer->first;
      const ScheduleMessageStream& incoming = peer->second.d_incoming;
      const ScheduleMessageStream& outgoing = peer->second.d_outgoing;
      stream << "Incoming Message Stream: " << s << std::endl;
      stream << "   bytes in stream: " << incoming.d_bytes_in_stream
              << std::endl;
      stream << "   must communicate bytes size: "
              << incoming.d_must_communicate_byte_size << std::endl;
      stream << "   stream in use: " << incoming.d_stream_in_use << std::endl;
      stream << "Outgoing Message Stream: " << s << std::endl;
      stream << "   bytes in stream: " << outgoing.d_bytes_in_stream
              << std::endl;
      stream << "   must communicate bytes size: "
              << outgoing.d_must_communicate_byte_size << std::endl;
      stream << "   stream in use: " << outgoing.d_stream_in_use << std::endl;
   }

   for (PeerMap::const_iterator peer = d_peers.begin(); peer != d_peers.end();
        ++peer) {
      stream << "Send Set: " << peer->first << std::endl;
      for (ITERATOR send(peer->second.d_send_set); send; send++) {
         send()->printClassData(stream);
      }
   }

   for (PeerMap::const_iterator peer = d_peers.begin(); peer != d_peers.end();
        ++peer) {
      stream << "Receive Set: " << peer->first << std::endl;
      for (ITERATOR recv(peer->second.d_recv_set); recv; recv++) {
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

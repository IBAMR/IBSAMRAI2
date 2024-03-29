//
// File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-4-4/source/toolbox/timers/Statistic.C $
// Package:     SAMRAI toolbox
// Copyright:   (c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:    $LastChangedRevision: 2043 $
// Modified:    $LastChangedDate: 2008-03-12 09:14:32 -0700 (Wed, 12 Mar 2008) $
// Description: Class to record statistics during program execution.
//

#include "tbox/Statistic.h"

#include "tbox/Array.h"
#include "tbox/SAMRAI_MPI.h"
#include "tbox/Utilities.h"
#include "tbox/MathUtilities.h"

#define TBOX_STATISTIC_VERSION (1)

#ifdef DEBUG_NO_INLINE
#include "tbox/Statistic.I"
#endif

namespace SAMRAI {
   namespace tbox {

#ifndef ARRAY_INCREMENT
#define ARRAY_INCREMENT (100)
#endif


double Statistic::s_empty_seq_tag_entry = -99999999.;

/*
*************************************************************************
*                                                                       *
* Statistic constructor and destructor.                                 * 
*                                                                       *
*************************************************************************
*/

Statistic::Statistic(const std::string& name,
                               const std::string& stat_type,
                               int instance_id)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!name.empty());
   TBOX_ASSERT(!stat_type.empty());
   TBOX_ASSERT(instance_id > -1);
#endif 

   if ( !(stat_type == "PROC_STAT" || stat_type == "PATCH_STAT") ) {
      TBOX_ERROR("Invalid stat_type passed to Statistic constructor"
         << "   object name = " << name
         << ".\n   Valid stat types are `PROC_STAT' and `PATCH_STAT'" << std::endl);
   }

   d_object_name = name;
   d_instance_id = instance_id;

   if (stat_type == "PROC_STAT") {
      d_stat_type = PROC_STAT;
   }
   if (stat_type == "PATCH_STAT") {
      d_stat_type = PATCH_STAT;
   }

   d_seq_counter           = 0;
   d_total_patch_entries   = 0;
   d_proc_stat_array_size  = 0;
   d_patch_stat_array_size = 0;
   
}

Statistic::~Statistic()
{
   reset();
}

/*
*************************************************************************
*                                                                       *
* Utility functions to record statistic record data.                    *
*                                                                       *
*************************************************************************
*/

void Statistic::recordProcStat(double value,
                                    int seq_num) 
{
   if (d_stat_type != PROC_STAT) {
      TBOX_ERROR("Statistic::recordProcStat error...\n"
                 << "    Statistic type is `PATCH_STAT'" << std::endl);
   }

   /*
    * Resize array of processor stats, if necessary.
    */
   checkArraySizes(seq_num);

   /*  
    * The statistic maintains its own sequence counter but the user may
    * override the data maintained in the statistic by supplying a seq_num.
    * The seq_num argument for this method is -1 if the user does not
    * supply an alternative value.  The value is recorded using the following
    * logic:
    * 1) If seq_num < 0, increment counter and record a new value at seq_cnt. 
    * 2) If 0 < seq_num < seq_cnt, overwrite the previous record at 
    *    seq_num with the new value.
    * 3) If seq_cnt < seq_num, set seq_cnt = seq_num and add new record
    *    at seq_cnt.  
    */
   if (seq_num < 0) {
      d_proc_array[d_seq_counter].value    = value;
   } else if (seq_num < d_seq_counter) {
      d_proc_array[seq_num].value          = value;
   } else {      
      d_seq_counter                        = seq_num;
      d_proc_array[d_seq_counter].value    = value;
   }
   d_seq_counter++;
}

void Statistic::recordPatchStat(int patch_num,
                                     double value,
                                     int seq_num)
{
   if (d_stat_type != PATCH_STAT) {
      TBOX_ERROR("Statistic::recordPatchStat error...\n"
         << "    Statistic type is `PROC_STAT'" << std::endl);
   }

   /*
    * Resize array of processor stats, if necessary.
    */
   checkArraySizes(seq_num);

   /*  
    * The patch statistic differs from the processor statistic in that
    * each entry of the array of seq numbers contains a LIST of 
    * PatchStatRecord entries, one for each patch on the processor.  
    * The recording logic is thus slightly different than the Processor stat:
    *
    *   If seq_num < seq_counter {  
    *     - Check the entries in the list of records at array index seq_num.
    *       Iterate through the list and see if the patch_id of any record
    *       matches the specified patch_num.
    *
    *       If patch_num entry exists {
    *          - overwrite existing entry
    *       } else {
    *          - create a new entry and append to end of list
    *
    *   } 
    *
    *   If seq_num >= seq_counter 
    *      - create new entry and append to end of list at the seq_num
    *        array index.
    *      - set seq_counter = seq_num
    *   }
    */
   if (seq_num < d_seq_counter) {
      List<Statistic::PatchStatRecord>& records = 
         d_patch_array[seq_num].patch_records;
      bool found_patch_id = false;
      List<Statistic::PatchStatRecord>::Iterator ir(records);
      for ( ; ir; ir++ ) {
         if (ir().patch_id == patch_num) {
            ir().value      = value;
            found_patch_id  = true;
         }
      }
      if (!found_patch_id) {
         PatchStatRecord patchitem_record;
         patchitem_record.value    = value;
         patchitem_record.patch_id = patch_num;
         d_patch_array[seq_num].patch_records.appendItem(patchitem_record);
         d_total_patch_entries++;
      }
      
   }
   
   if (seq_num >= d_seq_counter) {
      PatchStatRecord patchitem_record;
      patchitem_record.value    = value;
      patchitem_record.patch_id = patch_num;
      d_patch_array[seq_num].patch_records.appendItem(patchitem_record);
      d_total_patch_entries++;
      d_seq_counter = seq_num+1;
   }

}

/*
*************************************************************************
*                                                                       *
* Utility function for communicating statistic data in parallel.        *
*                                                                       *
* Stream data size includes 4 ints (instance id, proc rank,             *
*                                   stat type, seq length).             *
*                                                                       *
* Additionally, data stream size includes space needed for statistic    *
* data values:                                                          *
*                                                                       *
*    o for processor stat, this is 1 double (value) for each seq entry. *
*                                                                       *
*    o for patch stat, this is 1 int (#patches) for each sequence       *
*      entry + 1 int (patch_id) + 1 double (value) for each patch       *
*      entry.                                                           *
*                                                                       *
*************************************************************************
*/

int Statistic::getDataStreamSize()
{
   int byte_size = AbstractStream::sizeofInt(4);
   if (d_stat_type == PROC_STAT) {
      byte_size += AbstractStream::sizeofDouble(d_seq_counter);
   } else { // d_stat_type == PATCH_STAT
      byte_size += AbstractStream::sizeofInt(d_seq_counter);
      byte_size += AbstractStream::sizeofInt(d_total_patch_entries);
      byte_size += AbstractStream::sizeofDouble(d_total_patch_entries);
   }
   return(byte_size);
}

void Statistic::packStream(AbstractStream& stream)
{
   if (SAMRAI_MPI::getRank() == 0) {
      TBOX_ERROR("Statistic::packStream error...\n"
         << "    Processor zero should not pack stat data" << std::endl);
   }

   int num_int     = 4;
   int num_double  = 0;
   if (d_stat_type == PROC_STAT) {
      num_double   = d_seq_counter;
   }
   if (d_stat_type == PATCH_STAT) {
      num_int      += d_seq_counter + d_total_patch_entries;
      num_double   = d_total_patch_entries;
   }
   Array<int> idata(num_int);
   Array<double> ddata(num_double);

   idata[0] = SAMRAI_MPI::getRank();
   idata[1] = d_instance_id;
   idata[2] = d_stat_type;
   idata[3] = d_seq_counter;

   int is = 0;
   if (d_stat_type == PROC_STAT) {

      for (is = 0; is < d_seq_counter; is++) {
         ddata[is] = d_proc_array[is].value;
      }

   } else {  // d_stat_type == PATCH_STAT

      int mark = 4 + d_seq_counter;
      int isr = 0;

      for (is = 0; is < d_seq_counter; is++) {
         List<Statistic::PatchStatRecord>& lrec = 
            d_patch_array[is].patch_records;
         idata[4+is] = lrec.getNumberOfItems();

         List<Statistic::PatchStatRecord>::Iterator ilr(lrec);
         for ( ; ilr; ilr++) { 
            idata[mark+isr] = ilr().patch_id;
            ddata[isr] = ilr().value;
            isr++;
         }
      }
   }

   if (num_int > 0) {
      stream.pack(idata.getPointer(), num_int);
   }
   if (num_double > 0) {
      stream.pack(ddata.getPointer(), num_double);
   }
   
}

void Statistic::unpackStream(AbstractStream& stream)
{
   if (SAMRAI_MPI::getRank() != 0) {
      TBOX_ERROR("Statistic::unpackStream error...\n"
         << "    Only processor zero should unpack stat data" << std::endl);
   }

   int src_rank, stat_id, stat_type, seq_len;
   
   stream >> src_rank;
   stream >> stat_id;
   stream >> stat_type;
   stream >> seq_len;

   if (src_rank == 0) {
      TBOX_ERROR("Statistic::unpackStream error...\n"
         << "     Processor zero should not send stat data" << std::endl);
   }
   if (stat_id != d_instance_id) {
      TBOX_ERROR("Statistic::unpackStream error...\n"
         << "    Incompatible statistic number ids" << std::endl);
   }
   if (stat_type != d_stat_type) {
      TBOX_ERROR("Statistic::unpackStream error...\n"
         << "    Incompatible statistic types" << std::endl);
   }

   int is;
   if (d_stat_type == PROC_STAT) {

      Array<double> ddata(seq_len);

      if (seq_len > 0) {
         stream.unpack(ddata.getPointer(), seq_len);
         for (is = 0; is < seq_len; is++) {
            recordProcStat(ddata[is], is);
         } 
      }
      

   } else { // d_stat_type == PATCH_STAT
 
      if (seq_len > 0) {
         Array<int> inum_patches_data(seq_len);
         stream.unpack(inum_patches_data.getPointer(), seq_len);

         int total_seq_items = 0;
         for (is = 0; is < seq_len; is++) {
            total_seq_items += inum_patches_data[is];
         }

         Array<int> ipatch_num_data(total_seq_items);
         Array<double> ddata(total_seq_items);

         stream.unpack(ipatch_num_data.getPointer(), total_seq_items);
         stream.unpack(ddata.getPointer(), total_seq_items);

         int isr = 0;
         for (is = 0; is < seq_len; is++) {
            for (int ipsr = 0; ipsr < inum_patches_data[is]; ipsr++) {
               recordPatchStat(ipatch_num_data[isr], ddata[isr], is);
               isr++;
            }
         }
      }
   }

}

void Statistic::printClassData(std::ostream& stream,
                                    int precision) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
      TBOX_ASSERT(precision > 0);
#endif

   stream.precision(precision);

   stream << "Local Data for " << getType() << " : " << getName() << std::endl;
   stream << "   Processor id = " << SAMRAI_MPI::getRank() << std::endl;
   stream << "   Instance id = " << getInstanceId() << std::endl;
   stream << "   Sequence length = " << getStatSequenceLength() << std::endl;

   int is = 0;
   if (d_stat_type == PROC_STAT) {
      for (is = 0; is < d_seq_counter; is++) {
         stream << "     sequence[" << is 
                << "] : value = " << d_proc_array[is].value << std::endl;
      }
   } else {
      for (is = 0; is < d_seq_counter; is++) {
         stream << "     sequence[" << is 
                << "]" << std::endl;

         List<Statistic::PatchStatRecord>::Iterator ilr(
            d_patch_array[is].patch_records);
         for ( ; ilr; ilr++) {
            stream << "        patch # = " << ilr().patch_id
                   << " : value = " << ilr().value << std::endl;
         }
      }
   }
    
}

void Statistic::checkArraySizes(int seq_num) 
{
   /*
    * Verify that seq_num is less than array size.  If so, drop through.
    * If not, resize and initialize the array.
    */
   int high_mark = tbox::MathUtilities<int>::Max(seq_num, d_seq_counter);

   if (d_stat_type == PROC_STAT) {

      if (high_mark >= d_proc_stat_array_size) {
         int old_array_size = d_proc_stat_array_size;
         d_proc_stat_array_size += ARRAY_INCREMENT;
         d_proc_array.resizeArray(d_proc_stat_array_size);
         for (int i = old_array_size; i < d_proc_stat_array_size; i++) {
            d_proc_array[i].value    = s_empty_seq_tag_entry;
         }
         
         
      }

   } else if (d_stat_type == PATCH_STAT) {

      if (high_mark >= d_patch_stat_array_size) {
         d_patch_stat_array_size += ARRAY_INCREMENT;
         d_patch_array.resizeArray(d_patch_stat_array_size);
      }

   }

}

void Statistic::putToDatabase(
   Pointer<Database> db)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!db.isNull());
#endif
   db->putInteger("TBOX_STATISTIC_VERSION",
                  TBOX_STATISTIC_VERSION);
   
   db->putString("object_name", d_object_name);
   db->putInteger("stat_type", d_stat_type);
   db->putInteger("instance_id", d_instance_id);
   db->putInteger("seq_counter", d_seq_counter);
   db->putInteger("total_patch_entries", d_total_patch_entries);
   db->putInteger("proc_stat_array_size", d_proc_stat_array_size);
   db->putInteger("patch_stat_array_size", d_patch_stat_array_size);
   
   int i;
   
   if (d_stat_type == PROC_STAT) {
      Array<double> ddata(d_seq_counter);
      for (i = 0; i < d_seq_counter; i++) {
         ddata[i] = d_proc_array[i].value;
      }
      
      if (d_seq_counter > 0) {        
         db->putDoubleArray("ddata", ddata);
      }
      
   }
   
   if (d_stat_type == PATCH_STAT) {
      Array<int> idata(d_seq_counter + d_total_patch_entries);
      Array<double> ddata(d_total_patch_entries);

      int il = 0;
      int mark = d_seq_counter;

      for (i = 0; i < d_seq_counter; i++) {
         List<Statistic::PatchStatRecord>& records = 
            d_patch_array[i].patch_records;
         idata[i] = records.getNumberOfItems();  // # patches at seq num
         List<Statistic::PatchStatRecord>::Iterator ir(records);
         for ( ; ir; ir++ ) {
            idata[mark + il] = ir().patch_id;
            ddata[il] = ir().value;
            il++;
         }
      }

      if (d_seq_counter > 0) {
         db->putIntegerArray("idata", idata);
         if (d_total_patch_entries > 0) {
            db->putDoubleArray("ddata", ddata);
         }
      }
   }
}

void Statistic::getFromRestart(
   Pointer<Database> db)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!db.isNull());
#endif
   int ver = db->getInteger("TBOX_STATISTIC_VERSION");
   if (ver != TBOX_STATISTIC_VERSION) {
      TBOX_ERROR("Restart file version different than class version.");
   }

   d_object_name = db->getString("object_name");
   d_stat_type   = db->getInteger("stat_type");
   d_instance_id = db->getInteger("instance_id");
   int seq_entries = db->getInteger("seq_counter");
   int total_patches = db->getInteger("total_patch_entries");
   d_proc_stat_array_size = db->getInteger("proc_stat_array_size");
   d_patch_stat_array_size = db->getInteger("patch_stat_array_size");
   
   d_proc_array.resizeArray(d_proc_stat_array_size);
   d_patch_array.resizeArray(d_patch_stat_array_size);
   
   int i;
   if (d_stat_type == PROC_STAT) {
      if (seq_entries > 0) {
         Array<double> ddata = db->getDoubleArray("ddata");
         for (i = 0; i < seq_entries; i++) {
            recordProcStat(ddata[i], i);
         } 
      }
   } 

   if (d_stat_type == PATCH_STAT) {
      if (seq_entries > 0) {
         Array<int> idata = db->getIntegerArray("idata");
         
         Array<int> inum_patches(seq_entries);
         for (i = 0; i < seq_entries; i++) {
            inum_patches[i] = idata[i];
         }

         if (total_patches > 0) {
            Array<double> ddata = db->getDoubleArray("ddata");

            int il = 0;
            int mark = seq_entries;
            for (i = 0; i < seq_entries; i++) {
               for (int ipsr = 0; ipsr < inum_patches[i]; ipsr++) {
                  int patch_num = idata[mark+il];
                  double val = ddata[il];
                  recordPatchStat(patch_num, val, i);
                  il++;
               }
            }
         }
      }
   }
}


}
}

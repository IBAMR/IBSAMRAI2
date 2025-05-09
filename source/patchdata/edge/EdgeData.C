//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-4-4/source/patchdata/edge/EdgeData.C $
// Package:	SAMRAI patch data
// Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 2043 $
// Modified:	$LastChangedDate: 2008-03-12 09:14:32 -0700 (Wed, 12 Mar 2008) $
// Description:	Templated edge centered patch data type
//

#ifndef included_pdat_EdgeData_C
#define included_pdat_EdgeData_C

#include "EdgeData.h"

#include "Box.h"
#include "BoxList.h"
#include "EdgeGeometry.h"
#include "EdgeOverlap.h"
#include "tbox/Arena.h"
#include "tbox/ArenaManager.h"
#include "tbox/Utilities.h"
#include "tbox/TimerManager.h"

#define PDAT_EDGEDATA_VERSION 1

#ifdef DEBUG_NO_INLINE
#include "EdgeData.I"
#endif
namespace SAMRAI {
    namespace pdat {

/*
*************************************************************************
*									*
* Constructor and destructor for edge data objects.  The constructor	*
* simply initializes data variables and sets up the array data.		*
*									*
*************************************************************************
*/

template<int DIM, class TYPE>
EdgeData<DIM,TYPE>::EdgeData(const hier::Box<DIM>& box,
                             int depth,
                             const hier::IntVector<DIM>& ghosts,
                             tbox::Pointer<tbox::Arena> pool)
:  hier::PatchData<DIM>(box, ghosts),
   d_depth(depth)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(depth > 0);
   TBOX_ASSERT(ghosts.min() >= 0);
#endif
   if (pool.isNull()) {
      pool = tbox::ArenaManager::getManager()->getStandardAllocator();
   }
   for (int d = 0; d < DIM; d++) {
      const hier::Box<DIM> edge_box = 
         EdgeGeometry<DIM>::toEdgeBox(this -> getGhostBox(), d);
      d_data[d].initializeArray(edge_box, depth, pool);
   }
}

template<int DIM, class TYPE>
EdgeData<DIM,TYPE>::~EdgeData()
{
}

/*
*************************************************************************
*									*
* The following are private and cannot be used, but they are defined	*
* here for compilers that require that every template declaration have	*
* a definition (a stupid requirement, if you ask me).			*
*									*
*************************************************************************
*/

template<int DIM, class TYPE>
EdgeData<DIM,TYPE>::EdgeData(const EdgeData<DIM,TYPE>& foo)
:  hier::PatchData<DIM>(foo.getBox(), foo.getGhostCellWidth())
{
   NULL_USE(foo);
}

template<int DIM, class TYPE>
void EdgeData<DIM,TYPE>::operator=(const EdgeData<DIM,TYPE>& foo)
{
   NULL_USE(foo);
}

/*
*************************************************************************
*									*
* Perform a fast copy between two edge centered arrays where their	*
* index spaces overlap.							*
*									*
*************************************************************************
*/

template<int DIM, class TYPE>
void EdgeData<DIM,TYPE>::copy(const hier::PatchData<DIM>& src)
{
   const EdgeData<DIM,TYPE> *t_src = 
      dynamic_cast<const EdgeData<DIM,TYPE> *>(&src);
   if (t_src == NULL) {
      src.copy2(*this);
   } else {
      for (int d = 0; d < DIM; d++) {
         const hier::Box<DIM> box = d_data[d].getBox() * t_src->d_data[d].getBox();
         if (!box.empty()) {
            d_data[d].copy(t_src->d_data[d], box);
         }
      }
   }
}

template<int DIM, class TYPE>
void EdgeData<DIM,TYPE>::copy2(hier::PatchData<DIM>& dst) const
{
   EdgeData<DIM,TYPE> *t_dst =
      dynamic_cast<EdgeData<DIM,TYPE> *>(&dst);
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(t_dst != NULL);
#endif
   for (int d = 0; d < DIM; d++) {
      const hier::Box<DIM> box = d_data[d].getBox() * t_dst->d_data[d].getBox();
      if (!box.empty()) {
         t_dst->d_data[d].copy(d_data[d], box);
      }
   }
}

/*
*************************************************************************
*									*
* Copy data from the source into the destination according to the	*
* overlap descriptor.							*
*									*
*************************************************************************
*/

template<int DIM, class TYPE>
void EdgeData<DIM,TYPE>::copy(const hier::PatchData<DIM>& src,
                                const hier::BoxOverlap<DIM>& overlap)
{
   const EdgeData<DIM,TYPE> *t_src =
      dynamic_cast<const EdgeData<DIM,TYPE> *>(&src);

   const EdgeOverlap<DIM> *t_overlap =
      dynamic_cast<const EdgeOverlap<DIM> *>(&overlap);

   if ((t_src == NULL) || (t_overlap == NULL)) {
      src.copy2(*this, overlap);
   } else {
      const hier::IntVector<DIM>& src_offset = t_overlap->getSourceOffset();
      for (int d = 0; d < DIM; d++) {
         const hier::BoxList<DIM>& box_list = t_overlap->getDestinationBoxList(d);
         d_data[d].copy(t_src->d_data[d], box_list, src_offset);
      }
   }
}

template<int DIM, class TYPE>
void EdgeData<DIM,TYPE>::copy2(hier::PatchData<DIM>& dst,
                                 const hier::BoxOverlap<DIM>& overlap) const
{
   EdgeData<DIM,TYPE> *t_dst =
      dynamic_cast<EdgeData<DIM,TYPE> *>(&dst);
   const EdgeOverlap<DIM> *t_overlap =
      dynamic_cast<const EdgeOverlap<DIM> *>(&overlap);

#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(t_dst != NULL);
   TBOX_ASSERT(t_overlap != NULL);
#endif
   const hier::IntVector<DIM>& src_offset = t_overlap->getSourceOffset();
   for (int d = 0; d < DIM; d++) {
      const hier::BoxList<DIM>& box_list = t_overlap->getDestinationBoxList(d);
      t_dst->d_data[d].copy(d_data[d], box_list, src_offset);
   }
}

/*
*************************************************************************
*									*
* Perform a fast copy from an edge data object to this edge data        *
* object at the specified depths, where their index spaces overlap.     *
*									*
*************************************************************************
*/

template<int DIM, class TYPE>
void EdgeData<DIM,TYPE>::copyDepth(int dst_depth,
				     const EdgeData<DIM,TYPE>& src,
				     int src_depth)
{
   for (int d = 0; d < DIM; d++) {
      const hier::Box<DIM> box = d_data[d].getBox() * src.d_data[d].getBox();
      if (!box.empty()) {
	 d_data[d].copyDepth(dst_depth, src.d_data[d], src_depth, box);
      }
   }
}

/*
*************************************************************************
*									*
* Calculate the buffer space needed to pack/unpack messages on the box	*
* region using the overlap descriptor.					*
*									*
*************************************************************************
*/

template<int DIM, class TYPE>
bool EdgeData<DIM,TYPE>::canEstimateStreamSizeFromBox() const
{
   return(ArrayData<DIM,TYPE>::canEstimateStreamSizeFromBox());
}

template<int DIM, class TYPE>
int EdgeData<DIM,TYPE>::getDataStreamSize(
   const hier::BoxOverlap<DIM>& overlap) const
{
   const EdgeOverlap<DIM> *t_overlap =
      dynamic_cast<const EdgeOverlap<DIM> *>(&overlap);

#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(t_overlap != NULL);
#endif

   const hier::IntVector<DIM>& offset = t_overlap->getSourceOffset();

   int size = 0;
   for (int d = 0; d < DIM; d++) {
      size += d_data[d].getDataStreamSize(t_overlap->getDestinationBoxList(d),
                                          offset);
   }
   return(size);
}

/*
*************************************************************************
*									*
* Pack/unpack data into/out of the message streams using the index	*
* space in the overlap descriptor.					*
*									*
*************************************************************************
*/

template<int DIM, class TYPE>
void EdgeData<DIM,TYPE>::packStream(tbox::AbstractStream& stream,
                                      const hier::BoxOverlap<DIM>& overlap) const
{
   SAMRAI_SETUP_TIMER_AND_SCOPE("pdat::EdgeData::packStream()");
   const EdgeOverlap<DIM> *t_overlap =
      dynamic_cast<const EdgeOverlap<DIM> *>(&overlap);

#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(t_overlap != NULL);
#endif

   const hier::IntVector<DIM>& offset = t_overlap->getSourceOffset();
   for (int d = 0; d < DIM; d++) {
      const hier::BoxList<DIM>& boxes = t_overlap->getDestinationBoxList(d);
      if (boxes.getNumberOfItems() > 0) {
         d_data[d].packStream(stream, boxes, offset);
      }
   }
}

template<int DIM, class TYPE>
void EdgeData<DIM,TYPE>::unpackStream(tbox::AbstractStream& stream,
                                      const hier::BoxOverlap<DIM>& overlap)
{
   SAMRAI_SETUP_TIMER_AND_SCOPE("pdat::EdgeData::unpackStream()");
   const EdgeOverlap<DIM> *t_overlap =
      dynamic_cast<const EdgeOverlap<DIM> *>(&overlap);
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(t_overlap != NULL);
#endif

   const hier::IntVector<DIM>& offset = t_overlap->getSourceOffset();
   for (int d = 0; d < DIM; d++) {
      const hier::BoxList<DIM>& boxes = t_overlap->getDestinationBoxList(d);
      if (boxes.getNumberOfItems() > 0) {
         d_data[d].unpackStream(stream, boxes, offset);
      }
   }
}

/*
*************************************************************************
*									*
* Calculate the amount of memory space needed to represent the data	*
* for a  edge centered grid.						*
*									*
*************************************************************************
*/

template<int DIM, class TYPE>
size_t EdgeData<DIM,TYPE>::getSizeOfData(const hier::Box<DIM>& box,
                                         int depth,
                                         const hier::IntVector<DIM>& ghosts)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(depth > 0);
#endif
   size_t size = 0;
   const hier::Box<DIM> ghost_box = hier::Box<DIM>::grow(box, ghosts);
   for (int d = 0; d < DIM; d++) {
      const hier::Box<DIM> edge_box = EdgeGeometry<DIM>::toEdgeBox(ghost_box, d);
      size += ArrayData<DIM,TYPE>::getSizeOfData(edge_box, depth);
   }
   return(size);
}

/*
*************************************************************************
*									*
* Fill the edge centered box with the given value.			*
*									*
*************************************************************************
*/

template<int DIM, class TYPE>
void EdgeData<DIM,TYPE>::fill(const TYPE& t, int d)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT((d >= 0) && (d < d_depth));
#endif
   for (int i = 0; i < DIM; i++) {
      d_data[i].fill(t, d);
   }
}

template<int DIM, class TYPE>
void EdgeData<DIM,TYPE>::fill(const TYPE& t,
                                const hier::Box<DIM>& box,
                                int d)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT((d >= 0) && (d < d_depth));
#endif
   for (int i = 0; i < DIM; i++) {
      d_data[i].fill(t, EdgeGeometry<DIM>::toEdgeBox(box, i), d);
   }
}

template<int DIM, class TYPE>
void EdgeData<DIM,TYPE>::fillAll(const TYPE& t)
{
   for (int i = 0; i < DIM; i++) {
      d_data[i].fillAll(t);
   }
}

template<int DIM, class TYPE>
void EdgeData<DIM,TYPE>::fillAll(const TYPE& t, const hier::Box<DIM>& box)
{
   for (int i = 0; i < DIM; i++) {
      d_data[i].fillAll(t, EdgeGeometry<DIM>::toEdgeBox(box, i));
   }
}

/*
*************************************************************************
*                                                                       *
* Print edge centered data.  Note:  makes call to specialized printAxis *
* routine in EdgeDataSpecialized.C                                      *
*                                                                       *
*************************************************************************
*/

template<int DIM, class TYPE>
void EdgeData<DIM,TYPE>::print(const hier::Box<DIM>& box, 
                               std::ostream& os, 
                               int prec) const
{
   for (int axis = 0; axis < DIM; axis++) {
      os << "Array axis = " << axis << std::endl;
      printAxis(axis, box, os, prec);
   }
}

template<int DIM, class TYPE>
void EdgeData<DIM,TYPE>::print(const hier::Box<DIM>& box,
                               int depth,
                               std::ostream& os,
                               int prec) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT((depth >= 0) && (depth < d_depth));
#endif
   for (int axis = 0; axis < DIM; axis++) {
      os << "Array axis = " << axis << std::endl;
      printAxis(axis, box, depth, os, prec);
   }
}

template<int DIM, class TYPE>
void EdgeData<DIM,TYPE>::printAxis(int axis,
                                   const hier::Box<DIM>& box,
                                   std::ostream& os,
                                   int prec) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT((axis >= 0) && (axis < DIM));
#endif
   for (int d = 0; d < d_depth; d++) {
      os << "Array depth = " << d << std::endl;
      printAxis(axis, box, d, os, prec);
   }
}

template <int DIM, class TYPE>
void EdgeData<DIM,TYPE>::printAxis(int axis, 
                                   const hier::Box<DIM>& box,
                                   int depth, 
                                   std::ostream& os,
                                   int prec) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT((depth >= 0) && (depth < d_depth));
   TBOX_ASSERT((axis >= 0) && (axis < DIM));
#endif
   os.precision(prec);
   for (EdgeIterator<DIM> i(box, axis); i; i++) {
      os << "array" << i() << " = " << d_data[axis](i(),depth) 
         << std::endl << std::flush;
   } 
}

/*
*************************************************************************
*                                                                       *
* Checks that class version and restart file version are equal.  If so, *
* reads in the d_depth data member to the database.  Then tells		*
* d_data to read itself in from the database.				*
*                                                                       *
*************************************************************************
*/

template<int DIM, class TYPE>
void EdgeData<DIM,TYPE>::getSpecializedFromDatabase(
   tbox::Pointer<tbox::Database> database)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!database.isNull());
#endif
  
   int ver = database->getInteger("PDAT_EDGEDATA_VERSION");
   if (ver != PDAT_EDGEDATA_VERSION) {
      TBOX_ERROR("EdgeData<DIM>::getSpecializedFromDatabase error...\n"
          << " : Restart file version different than class version" << std::endl);
   }

   d_depth = database->getInteger("d_depth");

   tbox::Pointer<tbox::Database> array_database;
   for (int i = 0; i < DIM; i++) {
      std::string array_name = "d_data" + tbox::Utilities::intToString(i);
      array_database = database->getDatabase(array_name);
      (d_data[i]).getFromDatabase(array_database);
   }
}

/*
*************************************************************************
*                                                                       *
* Write out the class version number, d_depth data member to the	*
* database.  Then tells d_data to write itself to the database.		*
*                                                                       *
*************************************************************************
*/

template<int DIM, class TYPE>
void EdgeData<DIM,TYPE>::putSpecializedToDatabase(
   tbox::Pointer<tbox::Database> database)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!database.isNull());
#endif

   database->putInteger("PDAT_EDGEDATA_VERSION", PDAT_EDGEDATA_VERSION);

   database->putInteger("d_depth", d_depth);

   tbox::Pointer<tbox::Database> array_database;
   for (int i = 0; i < DIM; i++) {
      std::string array_name = "d_data" + tbox::Utilities::intToString(i);
      array_database = database->putDatabase(array_name);
      (d_data[i]).putToDatabase(array_database);
   }
}


}
}

#endif

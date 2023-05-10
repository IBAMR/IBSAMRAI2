//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-4-4/source/patchdata/edge/EdgeData.h $
// Package:	SAMRAI patch data
// Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 2132 $
// Modified:	$LastChangedDate: 2008-04-14 14:51:47 -0700 (Mon, 14 Apr 2008) $
// Description:	Templated edge centered patch data type
//

#ifndef included_pdat_EdgeData
#define included_pdat_EdgeData

#include "SAMRAI_config.h"

#ifndef included_iostream
#define included_iostream
#include <iostream>
#endif

#include "PatchData.h"
#include "ArrayData.h"
#include "EdgeIndex.h"
#include "EdgeIterator.h"
#include "tbox/Arena.h"
#include "tbox/Complex.h"
#include "tbox/PIO.h"
#include "tbox/Pointer.h"

namespace SAMRAI {
    namespace pdat {

/*!
 * @brief Class EdgeData<DIM> provides an implementation for data defined
 * at cell edges on AMR patches.  It is derived from the hier::PatchData
 * interface common to all SAMRAI patch data types.  Given a CELL-centered
 * AMR index space box, an edge data object represents data of some template 
 * TYPE and depth on the edges of the cells in the box.  Here, depth
 * indicates the number of data values at each edge index location.  The 
 * EdgeGeometry class provides the translation between the standard SAMRAI 
 * cell-centered AMR index space and edge-centered data.
 *
 * Edge data is stored in DIM arrays, each of which contains the 
 * data for the edges tangent to a corresponding coordinate direction.  
 * Within each array, data is stored in (i,...,k,d) order, where i,...,k 
 * indicates a spatial index and the d indicates the component depth at 
 * that locaion.  Memory allocation is in column-major ordering (e.g., Fortran  
 * style) so that the leftmost index runs fastest in memory.  
 * For example, a three-dimensional edge data object created over a
 * CELL-centered AMR index space [l0:u0,l1:u1,l2:u2] allocates three data 
 * arrays dimensioned as follows: 
 * \verbatim

   axis 0 
     [ l0 : u0 ,
       l1 : u1+1 ,
       l2 : u2+1 , d ]   ,

   axis 1 
     [ l0 : u0+1 ,
       l1 : u1 ,
       l2 : u2+1 , d ]   ,

   axis 2 
     [ l0 : u0+1 ,
       l1 : u1+1 ,
       l2 : u2 , d ]   ,
 
 * \endverbatim
 * Here the axis directions 0, 1, 2 can be thought of as the x, y, and z 
 * edge directions, respectively, and d is the depth index (i.e., number 
 * of values at each face index location).  Other spatial dimensions are 
 * represented similarly.
 *  
 * The data type TYPE must define a default constructor (that takes no
 * arguments) and also the assignment operator.
 *
 * @see pdat::ArrayData
 * @see hier::PatchData
 * @see pdat::EdgeDataFactory
 * @see pdat::EdgeIndex
 * @see pdat::EdgeIterator
 * @see pdat::EdgeGeometry
 */

template<int DIM, class TYPE>
class EdgeData : public hier::PatchData<DIM>
{
public:
   /*!
    * @brief Calculate the amount of memory needed to represent edge- 
    * centered data over a CELL-centered AMR index space box.
    *
    * This function assumes that the amount of memory
    * needed for TYPE is sizeof(TYPE).  If this is not the case, then a
    * specialized function must be defined.
    *
    * @param box const Box reference describing the interior of the
    *            standard CELL-centered index box over which the
    *            edge data object will be created.
    * @param depth gives the number of components for each
    *              spatial location in the array.
    * @param ghosts const IntVector reference indicating the width
    *              of the ghost cell region around the box over which
    *              the edge data will be allocated.
    */
   static size_t getSizeOfData(const hier::Box<DIM>& box,
                               int depth,
                               const hier::IntVector<DIM>& ghosts);

   /*!
    * @brief The constructor for an edge data object.
    *
    * @param box const Box reference describing the interior of the
    *            standard CELL-centered index box over which the
    *            edge data object will be created.
    * @param depth gives the number of components for each
    *              spatial location in the array.
    * @param ghosts const IntVector reference indicating the width
    *              of the ghost cell region around the box over which
    *              the edge data will be allocated.
    * @param pool memory arena.  If not given, then the
    *             standard arena is used.
    */
   EdgeData(const hier::Box<DIM>& box,
            int depth,
            const hier::IntVector<DIM>& ghosts,
            tbox::Pointer<tbox::Arena> pool = tbox::Pointer<tbox::Arena>(NULL));

   /*!
    * @brief The virtual destructor for an edge data object.
    */
   virtual ~EdgeData();

   /*!
    * @brief Return the depth (e.g., the number of components in each spatial
    * location) of the array.
    */
   int getDepth() const;

   /*!
    * @brief Get a pointer to the beginning of a particular axis and
    * depth component of the edge centered array.
    */
   TYPE* getPointer(int axis, int depth = 0);

   /*!
    * @brief Get a const pointer to the beginning of a particular axis and
    * depth component of the edge centered array.
    */
   const TYPE* getPointer(int axis, int depth = 0) const;

   /*!
    * @brief Return a reference to the data entry corresponding
    * to a given edge index and depth.
    */
   TYPE& operator()(const EdgeIndex<DIM>& i, int depth = 0);

   /*!
    * @brief Return a const reference to the data entry corresponding
    * to a given edge index and depth.
    */
   const TYPE& operator()(const EdgeIndex<DIM>& i, int depth = 0) const;

   /*!
    * @brief Return a reference to the array data object for the
    * given axis of the edge centered data object.
    */
   ArrayData<DIM,TYPE>& getArrayData(int axis);

   /*!
    * @brief Return a const reference to the array data object for the
    * given axis of the edge centered data object.
    */
   const ArrayData<DIM,TYPE>& getArrayData(int axis) const;

   /*!
    * @brief A fast copy from source to destination (i.e., this)
    * patch data object.
    *
    * Data is copied where there is overlap in the underlying index space.
    * The copy is performed on the interior plus the ghost cell width (for
    * both the source and destination).  Currently, source data must be
    * an EdgeData of the same DIM and TYPE.  If not, then an unrecoverable
    * error results.
    */
   virtual void copy(const hier::PatchData<DIM>& src);

   /*!
    * @brief A fast copy from source (i.e., this) to destination
    * patch data object.
    *
    * Data is copied where there is overlap in the underlying index space.
    * The copy is performed on the interior plus the ghost cell width (for
    * both the source and destination).  Currently, destination data must be
    * an EdgeData of the same DIM and TYPE.  If not, then an unrecoverable
    * error results.
    */
   virtual void copy2(hier::PatchData<DIM>& dst) const;

   /*!
    * @brief Copy data from source to destination (i.e., this)
    * patch data object on the given overlap.
    *
    * Currently, source data must be EdgeData of the same DIM and TYPE 
    * and the overlap must be a EdgeOverlap of the same DIM. If not, 
    * then an unrecoverable error results.
    */
   virtual void copy(const hier::PatchData<DIM>& src,
                     const hier::BoxOverlap<DIM>& overlap);

   /*!
    * @brief Copy data from source (i.e., this) to destination
    * patch data object on the given overlap.
    *
    * Currently, destination data must be EdgeData of the same DIM and TYPE
    * and the overlap must be a EdgeOverlap of the same DIM.  If not, 
    * then an unrecoverable error results.
    */
   virtual void copy2(hier::PatchData<DIM>& dst,
                      const hier::BoxOverlap<DIM>& overlap) const;

   /*!
    * @brief Copy data from source to destination (i.e., this)
    * patch data object on the given CELL-centered AMR index box.
    */
   void copyOnBox(const EdgeData<DIM,TYPE>& src, 
                  const hier::Box<DIM>& box);
   
   /*!
    * @brief Fast copy (i.e., source and this edge data objects are
    * defined over the same box) to this destination edge data object 
    * from the given source edge data object at the specified depths.
    */
   void copyDepth(int dst_depth,
		  const EdgeData<DIM,TYPE>& src,
		  int src_depth);

   /*!
    * @brief Return true if the patch data object can estimate the
    * stream size required to fit its data using only index
    * space information (i.e., a box).
    *
    * This routine is defined for the standard types (bool, char,
    * double, float, int, and dcomplex).
    */
   virtual bool canEstimateStreamSizeFromBox() const;

   /*!
    * @brief Return the number of bytes needed to stream the data
    * in this patch data object lying in the specified box overlap
    * region.
    *
    * This routine is defined for the standard types (bool, char,
    * double, float, int, and dcomplex). 
    */
   virtual int getDataStreamSize(const hier::BoxOverlap<DIM>& overlap) const;

   /*!
    * @brief Pack data in this patch data object lying in the specified
    * box overlap region into the stream.  The overlap must be an
    * EdgeOverlap of the same DIM.
    */
   virtual void packStream(tbox::AbstractStream& stream,
                           const hier::BoxOverlap<DIM>& overlap) const;
 
   /*!
    * @brief Unpack data from stream into this patch data object over
    * the specified box overlap region. The overlap must be an
    * EdgeOverlap of the same DIM.
    */
   virtual void unpackStream(tbox::AbstractStream& stream,
                             const hier::BoxOverlap<DIM>& overlap);

   /*!
    * @brief Fill all values at depth d with the value t.
    */
   void fill(const TYPE& t,
             int d = 0);
 
   /*!
    * @brief Fill all values at depth d within the box with the value t.
    */
   void fill(const TYPE& t,
             const hier::Box<DIM>& box,
             int d = 0);
 
   /*!
    * @brief Fill all depth components with value t.
    */
   void fillAll(const TYPE& t);
 
   /*!
    * @brief Fill all depth components within the box with value t.
    */
   void fillAll(const TYPE& t,
                const hier::Box<DIM>& box);

   /*!
    * @brief Print all edge data values residing in the specified box.
    * If the depth of the array is greater than one, all depths are printed.
    *
    * @param box  const reference to box over whioch to print data. Note box
    *        is assumed to reside in standard cell-centered index space
    *        and will be converted to edge index space.
    * @param os   reference to output stream.
    * @param prec integer precision for printing floating point numbers
    *        (i.e., TYPE = float, double, or dcomplex). The default
    *        is 12 decimal places for double and complex floating point numbers,
    *        and the default is 6 decimal places floats.  For other types, this
    *        value is ignored.
    */
   void print(const hier::Box<DIM>& box, 
              std::ostream& os = tbox::plog, 
              int prec = 12) const;

   /*!
    * @brief Print all edge data values at the given array depth in
    * the specified box.
    *
    * @param box  const reference to box over whioch to print data. Note box
    *        is assumed to reside in standard cell-centered index space
    *        and will be converted to edge index space.
    * @param depth integer depth component, must satisfy
    *              0 <= depth < actual depth of data array
    * @param os   reference to output stream.
    * @param prec integer precision for printing floating point numbers
    *        (i.e., TYPE = float, double, or dcomplex). The default
    *        is 12 decimal places for double and complex floating point numbers,
    *        and the default is 6 decimal places floats.  For other types, this
    *        value is ignored.
    */
   void print(const hier::Box<DIM>& box, 
              int depth, 
              std::ostream& os = tbox::plog,
              int prec = 12) const;

   /*!
    * @brief Print all edge centered data values for specified axis
    * residing in the specified box.  If the depth of the data is 
    * greater than one, all depths are printed.   
    *
    * @param axis  integer edge data coordinate axis,
    *              must satisfy 0 <= axis < DIM
    * @param box  const reference to box over whioch to print data. Note box
    *        is assumed to reside in standard cell-centered index space
    *        and will be converted to edge index space.
    * @param os    reference to output stream.
    * @param prec integer precision for printing floating point numbers
    *        (i.e., TYPE = float, double, or dcomplex). The default
    *        is 12 decimal places for double and complex floating point numbers,
    *        and the default is 6 decimal places floats.  For other types, this
    *        value is ignored.
    */
   void printAxis(int axis,
                  const hier::Box<DIM>& box,
                  std::ostream& os = tbox::plog,
                  int prec = 12) const;

   /*!
    * @brief Print all edge centered data values for specified axis
    * residing in the specified box.  If the depth of the data is
    * greater than one, all depths are printed.
    *
    * @param axis  integer edge data coordinate axis,
    *              must satisfy 0 <= axis < DIM
    * @param box  const reference to box over whioch to print data. Note box
    *        is assumed to reside in standard cell-centered index space
    *        and will be converted to edge index space.
    * @param depth integer depth component, must satisfy
    *              0 <= depth < actual depth of data array
    * @param os    reference to output stream.
    * @param prec integer precision for printing floating point numbers
    *        (i.e., TYPE = float, double, or dcomplex). The default
    *        is 12 decimal places for double and complex floating point numbers,
    *        and the default is 6 decimal places floats.  For other types, this
    *        value is ignored.
    */
   void printAxis(int axis,
                  const hier::Box<DIM>& box,
                  int depth,
                  std::ostream& os = tbox::plog,
                  int prec = 12) const;  

   /*!
    * Check that class version and restart file version are equal.  If so,
    * read data members from the database.
    *
    * Assertions: database must be non-null pointer.
    */
   virtual void getSpecializedFromDatabase(
           tbox::Pointer<tbox::Database> database);

   /*!
    * Write out the class version number and other data members to
    * the database.
    *
    * Assertions: database must be non-null pointer.
    */
   virtual void putSpecializedToDatabase(
           tbox::Pointer<tbox::Database> database);

   /*!
    * The edge iterator iterates over the elements on one axis of an edge
    * centered box geometry.  This typedef is a convenience for
    * using the EdgeIterator<DIM> class.
    */
   typedef EdgeIterator<DIM> Iterator;

private:
   EdgeData(const EdgeData<DIM,TYPE>&);	// not implemented
   void operator=(const EdgeData<DIM,TYPE>&);		// not implemented

   int d_depth;
   ArrayData<DIM,TYPE> d_data[DIM];

};

}
}

#ifndef DEBUG_NO_INLINE
#include "EdgeData.I"
#endif
#endif

#ifdef INCLUDE_TEMPLATE_IMPLEMENTATION
#include "EdgeData.C"
#endif

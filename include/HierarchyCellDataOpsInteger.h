//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-4-4/source/mathops/cell/HierarchyCellDataOpsInteger.h $
// Package:     SAMRAI mathops
// Copyright:   (c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:    $LastChangedRevision: 2132 $
// Modified:    $LastChangedDate: 2008-04-14 14:51:47 -0700 (Mon, 14 Apr 2008) $
// Description: Operations for integer cell data on multiple levels.
//

#ifndef included_math_HierarchyCellDataOpsInteger
#define included_math_HierarchyCellDataOpsInteger

#include "SAMRAI_config.h"
#ifndef included_iostream
#define included_iostream
#include <iostream>
#endif
#include "HierarchyDataOpsInteger.h"
#include "PatchCellDataOpsInteger.h"
#include "Box.h"
#include "PatchHierarchy.h"
#include "tbox/Pointer.h"

namespace SAMRAI {
    namespace math {

/**
 * Class HierarchyCellDataOpsInteger<DIM> provides a collection of 
 * operations that manipulate integer cell-centered patch data components over 
 * multiple levels in an AMR hierarchy.  It is derived from the abstract
 * base class HierarchyDataOpsInteger<DIM> which defines the interface to 
 * similar operations for cell-centered, face-centered, node-centered patch
 * data objects where the data is of type integer.  The operations include 
 * basic arithmetic and some ordering operations.  On each patch, the
 * operations are performed by the PatchCellDataOpsInteger<DIM> data member.
 *
 * The patch hierarchy and set of levels within that hierarcy over which the
 * operations will be performed are set in the constructor.  However, note
 * that the constructor accepts default arguments for the coarsest and finest
 * level numbers.  If the level numbers are not specified when calling the
 * constructor the levels which exist in the hierarchy will be assumed in
 * all operations.  The hierarchy and levels may be changed at any time using
 * the proper member functions.
 *
 * Note that, when it makes sense, an operation accept a boolean argument 
 * which indicates whether the operation should be performed on all of the 
 * data or just those data elements corresponding to the patch interiors.  
 * If no boolean argument is provided, the default behavior is to treat only 
 * the patch interiors.  Also, a similar set of operations for real (double
 * and float) and complex cell-centered data is provided in the classes 
 * HierarchyCellDataOpsReal<DIM> and HierarchyCellDataOpsComplex<DIM>, 
 * respectively.
 * 
 * @see math::PatchCellDataOpsInteger
 */

template<int DIM>
class HierarchyCellDataOpsInteger : public HierarchyDataOpsInteger<DIM>
{
public:
   /**
    * The constructor for the HierarchyCellDataOpsInteger<DIM> class sets
    * the default patch hierarchy and coarsest and finest patch levels
    * in that hierarchy over which operations will be performed.  The
    * hierarchy and operations may be reset using the member functions
    * setPatchHierarchy() and resetLevels() below.  If no level number
    * arguments are given here, the levels over which the operations will
    * be performed are those already existing in the hierarchy.  If the 
    * hierarchy level configuration changes, the operations must be explicitly
    * reset by calling the resetLevels() function.
    */
   HierarchyCellDataOpsInteger(
      tbox::Pointer< hier::PatchHierarchy<DIM> > hierarchy,
      const int coarsest_level = -1,
      const int finest_level = -1);

   /**
    * Virtual destructor for the HierarchyCellDataOpsInteger<DIM> class.
    */
   virtual ~HierarchyCellDataOpsInteger<DIM>();

   /**
    * Reset patch hierarchy over which operations occur.
    */
   void setPatchHierarchy(tbox::Pointer< hier::PatchHierarchy<DIM> > hierarchy);

   /**
    * Reset range of patch levels over which operations occur.
    * The levels must exist in the hierarchy or an assertion will result. 
    */
   void resetLevels(const int coarsest_level, 
                    const int finest_level);

   /**
    * Return const pointer to patch hierarchy associated with operations.
    */
   const tbox::Pointer< hier::PatchHierarchy<DIM> > getPatchHierarchy() const;

   /**
    * Return the total number of data values for the component on the set 
    * of hierarchy levels.  If the boolean argument is true, the number of
    * elements will be summed over patch interiors.  If the boolean argument
    * is false, all elements will be counted (including ghost values) 
    * over all patches.
    */
   int numberOfEntries(const int data_id,
                       const bool interior_only = true) const;

   /**
    * Copy source data to destination data.
    */
   void copyData(const int dst_id,
                 const int src_id,
                 const bool interior_only = true) const;

   /**
    * Swap data pointers (i.e., storage) between two data components.
    */
   void swapData(const int data1_id,
                 const int data2_id) const;

   /**
    * Print data over multiple levels to specified output stream.
    */
   void printData(const int data_id,
                  std::ostream& s,
                  const bool interior_only = true) const;

   /**
    * Set data component to given scalar.
    */
   void setToScalar(const int data_id,
                    const int& alpha,
                    const bool interior_only = true) const;

   /**
    * Set destination to source multiplied by given scalar, pointwise. 
    */
   void scale(const int dst_id,
              const int& alpha,
              const int src_id,
              const bool interior_only = true) const;

   /**
    * Add scalar to each entry in source data and set destination to result.
    */
   void addScalar(const int dst_id,
                  const int src_id,
                  const int& alpha,
                  const bool interior_only = true) const;

   /**
    * Set destination to sum of two source components, pointwise.
    */
   void add(const int dst_id,
            const int src1_id,
            const int src2_id,
            const bool interior_only = true) const;

   /**
    * Subtract second source component from first source component pointwise 
    * and set destination data component to result.
    */
   void subtract(const int dst_id,
                 const int src1_id,
                 const int src2_id,
                 const bool interior_only = true) const;

   /**
    * Set destination component to product of two source components, pointwise.
    */
   void multiply(const int dst_id,
                 const int src1_id,
                 const int src2_id,
                 const bool interior_only = true) const;

   /**
    * Divide first data component by second source component pointwise
    * and set destination data component to result.
    */
   void divide(const int dst_id,
               const int src1_id,
               const int src2_id,
               const bool interior_only = true) const;

   /**
    * Set each entry of destination component to reciprocal of corresponding 
    * source data component entry.
    */
   void reciprocal(const int dst_id,
                   const int src_id,
                   const bool interior_only = true) const;

   /**
    * Set \f$d = \alpha s_1 + \beta s_2\f$, where \f$d\f$ is the destination patch
    * data component and \f$s_1, s_2\f$ are the first and second source components,
    * respectively.  Here \f$\alpha, \beta\f$ are scalar values.
    */
   void linearSum(const int dst_id,
                  const int& alpha,
                  const int src1_id,
                  const int& beta,
                  const int src2_id,
                  const bool interior_only = true) const;

   /**
    * Set \f$d = \alpha s_1 + s_2\f$, where \f$d\f$ is the destination patch data
    * component and \f$s_1, s_2\f$ are the first and second source components,
    * respectively.  Here \f$\alpha\f$ is a scalar. 
    */
   void axpy(const int dst_id,
             const int& alpha,
             const int src1_id,
             const int src2_id,
             const bool interior_only = true) const;

   /**
    * Set \f$d = \alpha s_1 - s_2\f$, where \f$d\f$ is the destination patch data
    * component and \f$s_1, s_2\f$ are the first and second source components,
    * respectively.  Here \f$\alpha\f$ is a scalar.
    */
   void axmy(const int dst_id,
             const int& alpha,
             const int src1_id,
             const int src2_id,
             const bool interior_only = true) const;

   /**
    * Set destination data to absolute value of source data, pointwise.
    */
   void abs(const int dst_id,
            const int src_id,
            const bool interior_only = true) const;

   /**
    * Return minimum data value over all patches in the collection of levels.
    */
   int min(const int data_id,
           const bool interior_only = true) const;

   /**
    * Return maximum data value over all patches in the collection of levels.
    */
   int max(const int data_id,
           const bool interior_only = true) const;

   /**
    * Set data entries to random values.  See the operations in the 
    * array data operation classes for details on the generation of 
    * the random values.
    */
   void setRandomValues(const int data_id,
                        const int& width,
                        const int& low,
                        const bool interior_only = true) const;

private:
   // The following are not implemented
   HierarchyCellDataOpsInteger(const HierarchyCellDataOpsInteger<DIM>&);
   void operator=(const HierarchyCellDataOpsInteger<DIM>&);

   tbox::Pointer< hier::PatchHierarchy<DIM> > d_hierarchy;
   int  d_coarsest_level;
   int  d_finest_level; 
   PatchCellDataOpsInteger<DIM> d_patch_ops;

};

}
}
#endif

#ifdef INCLUDE_TEMPLATE_IMPLEMENTATION
#include "HierarchyCellDataOpsInteger.C"
#endif

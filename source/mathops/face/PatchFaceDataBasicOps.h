//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-4-4/source/mathops/face/PatchFaceDataBasicOps.h $
// Package:	SAMRAI mathops
// Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 2132 $
// Modified:	$LastChangedDate: 2008-04-14 14:51:47 -0700 (Mon, 14 Apr 2008) $
// Description:	Basic templated face-centered patch data operations.
//

#ifndef included_math_PatchFaceDataBasicOps
#define included_math_PatchFaceDataBasicOps

#include "SAMRAI_config.h"
#ifndef included_iostream
#define included_iostream
#include <iostream>
#endif
#include "tbox/Complex.h"
#include "ArrayDataBasicOps.h"
#include "Box.h"
#include "FaceData.h"
#include "tbox/Pointer.h"

namespace SAMRAI {
    namespace math {

/**
 * Class PatchFaceDataBasicOps<DIM> provides access to a collection 
 * of basic numerical operations that may be applied to numerical face-
 * centered patch data.  These operations include simple arithmetic 
 * operations as well as min and max, etc.  The primary intent of this
 * class is to provide the interface to these standard operations for
 * an PatchFaceDataOps<DIM> object which provides access to a complete set
 * of operations that may be used to manipulate face-centered patch data
 * objects.   Each member function accepts a box argument indicating the 
 * region of index space on which the operation should be performed.  The 
 * operation will be performed on the intersection of this box and those 
 * boxes corresponding to the patch data objects involved.
 *
 * These operations typically apply only to the numerical standard built-in
 * types, such as double, float, and int, and the complex type (which may or
 * may not be a built-in type depending on the C++ compiler).  Thus, this
 * templated class should only be used to instantiate objects with those
 * types as the template parameter.  None of the operations are implemented
 * for any other type.
 *
 * @see math::ArrayDataBasicOps
 */

template<int DIM, class TYPE>
class PatchFaceDataBasicOps
{
public:
   /** 
    * Empty constructor and destructor.
    */
   PatchFaceDataBasicOps();

   virtual ~PatchFaceDataBasicOps<DIM,TYPE>();

   /**
    * Set dst = alpha * src, elementwise.
    */
   void scale(tbox::Pointer< pdat::FaceData<DIM,TYPE> >& dst,
              const TYPE& alpha,
              const tbox::Pointer< pdat::FaceData<DIM,TYPE> >& src,
              const hier::Box<DIM>& box) const;

   /**
    * Set dst = src + alpha, elementwise.
    */
   void addScalar(tbox::Pointer< pdat::FaceData<DIM,TYPE> >& dst,
                  const tbox::Pointer< pdat::FaceData<DIM,TYPE> >& src,
                  const TYPE& alpha,
                  const hier::Box<DIM>& box) const;
  
   /**
    * Set dst = src1 + src2, elementwise.
    */
   void add(tbox::Pointer< pdat::FaceData<DIM,TYPE> >& dst,
            const tbox::Pointer< pdat::FaceData<DIM,TYPE> >& src1,
            const tbox::Pointer< pdat::FaceData<DIM,TYPE> >& src2,
            const hier::Box<DIM>& box) const;

   /**
    * Set dst = src1 - src2, elementwise.
    */
   void subtract(tbox::Pointer< pdat::FaceData<DIM,TYPE> >& dst,
                 const tbox::Pointer< pdat::FaceData<DIM,TYPE> >& src1,
                 const tbox::Pointer< pdat::FaceData<DIM,TYPE> >& src2,
                 const hier::Box<DIM>& box) const;

   /**
    * Set dst = src1 * src2, elementwise.
    */
   void multiply(tbox::Pointer< pdat::FaceData<DIM,TYPE> >& dst,
                 const tbox::Pointer< pdat::FaceData<DIM,TYPE> >& src1,
                 const tbox::Pointer< pdat::FaceData<DIM,TYPE> >& src2,
                 const hier::Box<DIM>& box) const;

   /**
    * Set dst = src1 / src2, elementwise.  No check for division by zero.
    */
   void divide(tbox::Pointer< pdat::FaceData<DIM,TYPE> >& dst,
               const tbox::Pointer< pdat::FaceData<DIM,TYPE> >& src1,
               const tbox::Pointer< pdat::FaceData<DIM,TYPE> >& src2,
               const hier::Box<DIM>& box) const;

   /**
    * Set dst = 1 / src, elementwise.  No check for division by zero.
    */
   void reciprocal(tbox::Pointer< pdat::FaceData<DIM,TYPE> >& dst,
                   const tbox::Pointer< pdat::FaceData<DIM,TYPE> >& src,
                   const hier::Box<DIM>& box) const;

   /**
    * Set dst = alpha * src1 + beta * src2, elementwise.
    */
   void linearSum(tbox::Pointer< pdat::FaceData<DIM,TYPE> >& dst,
                  const TYPE& alpha,
                  const tbox::Pointer< pdat::FaceData<DIM,TYPE> >& src1,
                  const TYPE& beta,
                  const tbox::Pointer< pdat::FaceData<DIM,TYPE> >& src2,
                  const hier::Box<DIM>& box) const;

   /**
    * Set dst = alpha * src1 + src2, elementwise.
    */
   void axpy(tbox::Pointer< pdat::FaceData<DIM,TYPE> >& dst,
             const TYPE& alpha,
             const tbox::Pointer< pdat::FaceData<DIM,TYPE> >& src1,
             const tbox::Pointer< pdat::FaceData<DIM,TYPE> >& src2,
             const hier::Box<DIM>& box) const;

   /**
    * Set dst = alpha * src1 - src2, elementwise.
    */ 
   void axmy(tbox::Pointer< pdat::FaceData<DIM,TYPE> >& dst,
             const TYPE& alpha,
             const tbox::Pointer< pdat::FaceData<DIM,TYPE> >& src1,
             const tbox::Pointer< pdat::FaceData<DIM,TYPE> >& src2,
             const hier::Box<DIM>& box) const;

   /**
    * Return the minimum patch data component entry  When the data is
    * complex, the result is the data element with the smallest norm.
    */
   TYPE min(const tbox::Pointer< pdat::FaceData<DIM,TYPE> >& data,
            const hier::Box<DIM>& box) const;

   /**
    * Return the maximum patch data component entry  When the data is
    * complex, the result is the data element with the largest norm.
    */
   TYPE max(const tbox::Pointer< pdat::FaceData<DIM,TYPE> >& data,
            const hier::Box<DIM>& box) const;

   /**
    * Set patch data to random values.  See the operations in the
    * ArrayDataBasicOps<DIM> class for details on the generation
    * of the random values for each data type.
    */
   void setRandomValues(tbox::Pointer< pdat::FaceData<DIM,TYPE> >& dst,
                        const TYPE& width,
                        const TYPE& low,
                        const hier::Box<DIM>& box) const;

private:
   // The following are not implemented:
   PatchFaceDataBasicOps(
      const PatchFaceDataBasicOps<DIM,TYPE>&);
   void operator=(const PatchFaceDataBasicOps<DIM,TYPE>&);

   ArrayDataBasicOps<DIM,TYPE> d_array_ops;
};

}
}
#endif

#ifdef INCLUDE_TEMPLATE_IMPLEMENTATION
#include "PatchFaceDataBasicOps.C"
#endif

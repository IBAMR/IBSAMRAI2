//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-4-4/source/mathops/face/PatchFaceDataNormOpsReal.h $
// Package:	SAMRAI mathops
// Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 2132 $
// Modified:	$LastChangedDate: 2008-04-14 14:51:47 -0700 (Mon, 14 Apr 2008) $
// Description:	Templated norm operations for real face-centered data.
//

#ifndef included_math_PatchFaceDataNormOpsReal
#define included_math_PatchFaceDataNormOpsReal

#include "SAMRAI_config.h"
#ifndef included_iostream
#define included_iostream
#include <iostream>
#endif
#include "ArrayDataNormOpsReal.h"
#include "Box.h"
#include "FaceData.h"
#include "tbox/Pointer.h"

namespace SAMRAI {
    namespace math {

#ifndef NULL
#define NULL (0)
#endif

/**
 * Class PatchFaceDataNormOpsReal<DIM> provides a collection of common 
 * norm operations that may be applied to real (double or float)
 * numerical face-centered patch data.  The primary intent of this class is 
 * to define part of the interface for an PatchFaceDataOpsReal<DIM> object 
 * which provides access operations that may be used to manipulate 
 * face-centered patch data.  Each member function accepts a box argument 
 * indicating the region of index space on which the operation should be 
 * performed.  The operation will be performed on the intersection of this 
 * box and those boxes corresponding to the patch data objects.  Also, each 
 * operation allows an additional face-centered patch data object to be used
 * to represent a control volume that weights the contribution of each data 
 * entry in the given norm calculation.  Note that the control volume patch 
 * data must be of type double and have face-centered geometry (i.e., the 
 * same as the data itself).  The use of control volumes is important when 
 * these operations are used in vector kernels where the data resides over 
 * multiple levels of spatial resolution in an AMR hierarchy.  If the control 
 * volume is not given in the function call, it will be ignored in the 
 * calculation.  Also, note that the depth of the control volume patch data
 * object must be either 1 or be equal to the depth of the other data objects.
 * 
 * These operations typically apply only to the numerical standard built-in
 * types, such as double, float, and the complex type (which may or may not 
 * be a built-in type depending on the C++ compiler).  This templated 
 * class should only be used to instantiate objects with double or float as 
 * the template parameter.  Note that a similar set of norm operations is 
 * implemented for complex patch data in the class 
 * PatchFaceDataNormOpsComplex<DIM>.
 *
 * @see math::ArrayDataNormOpsReal
 */

template<int DIM, class TYPE>
class PatchFaceDataNormOpsReal
{
public:
   /** 
    * Empty constructor and destructor.
    */
   PatchFaceDataNormOpsReal();

   virtual ~PatchFaceDataNormOpsReal();

   /**
    * Return the number of data values for the face-centered data object
    * in the given box.  Note that it is assumed that the box refers to
    * the cell-centered index space corresponding to the patch hierarchy.
    */
   int numberOfEntries(const tbox::Pointer< pdat::FaceData<DIM,TYPE> >& data,
                       const hier::Box<DIM>& box) const;

   /**
    * Return sum of control volume entries for the face-centered data object.
    */
   double sumControlVolumes(
      const tbox::Pointer< pdat::FaceData<DIM,TYPE> >& data, 
      const tbox::Pointer< pdat::FaceData<DIM,double> >& cvol, 
      const hier::Box<DIM>& box) const;

   /**
    * Set destination component to absolute value of source component.
    * That is, each destination entry is set to \f$d_i = \| s_i \|\f$.
    */
   void abs(tbox::Pointer< pdat::FaceData<DIM,TYPE> >& dst,
            const tbox::Pointer< pdat::FaceData<DIM,TYPE> >& src,
            const hier::Box<DIM>& box) const;

   /**
    * Return discrete \f$L_1\f$-norm of the data using the control volume to
    * weight the contribution of each data entry to the sum.  That is, the
    * return value is the sum \f$\sum_i ( \| data_i \| cvol_i )\f$.  If the 
    * control volume is NULL, the return value is \f$\sum_i ( \| data_i \| )\f$.
    */
   double L1Norm(
      const tbox::Pointer< pdat::FaceData<DIM,TYPE> >& data, 
      const hier::Box<DIM>& box, 
      const tbox::Pointer< pdat::FaceData<DIM,double> > cvol = NULL) const;

   /**
    * Return discrete \f$L_2\f$-norm of the data using the control volume to
    * weight the contribution of each data entry to the sum.  That is, the
    * return value is the sum \f$\sqrt{ \sum_i ( (data_i)^2 cvol_i ) }\f$.
    * If the control volume is NULL, the return value is 
    * \f$\sqrt{ \sum_i ( (data_i)^2 cvol_i ) }\f$. 
    */
   double L2Norm(
      const tbox::Pointer< pdat::FaceData<DIM,TYPE> >& data, 
      const hier::Box<DIM>& box, 
      const tbox::Pointer< pdat::FaceData<DIM,double> > cvol = NULL) const;

   /**
    * Return discrete weighted \f$L_2\f$-norm of the data using the control
    * volume to weight the contribution of the data and weight entries to
    * the sum.  That is, the return value is the sum \f$\sqrt{ \sum_i (
    * (data_i * weight_i)^2 cvol_i ) }\f$.  If the control volume is NULL,
    * the return value is \f$\sqrt{ \sum_i ( (data_i * weight_i)^2 ) }\f$.
    */
   double weightedL2Norm(
      const tbox::Pointer< pdat::FaceData<DIM,TYPE> >& data,
      const tbox::Pointer< pdat::FaceData<DIM,TYPE> >& weight,
      const hier::Box<DIM>& box,
      const tbox::Pointer< pdat::FaceData<DIM,double> > cvol = NULL) const;

   /**
    * Return discrete root mean squared norm of the data.  If the control
    * volume is not NULL, the return value is the \f$L_2\f$-norm divided by
    * the square root of the sum of the control volumes.  Otherwise, the
    * return value is the \f$L_2\f$-norm divided by the square root of the 
    * number of data entries.
    */
   double RMSNorm(
      const tbox::Pointer< pdat::FaceData<DIM,TYPE> >& data,
      const hier::Box<DIM>& box,
      const tbox::Pointer< pdat::FaceData<DIM,double> > cvol = NULL) const;

   /**
    * Return discrete weighted root mean squared norm of the data.  If the 
    * control volume is not NULL, the return value is the weighted \f$L_2\f$-norm 
    * divided by the square root of the sum of the control volumes.  Otherwise, 
    * the return value is the weighted \f$L_2\f$-norm divided by the square root 
    * of the number of data entries.
    */
   double weightedRMSNorm(
      const tbox::Pointer< pdat::FaceData<DIM,TYPE> >& data,
      const tbox::Pointer< pdat::FaceData<DIM,TYPE> >& weight,
      const hier::Box<DIM>& box,
      const tbox::Pointer< pdat::FaceData<DIM,double> > cvol = NULL) const;

   /**
    * Return the \f$\max\f$-norm of the data using the control volume to weight
    * the contribution of each data entry to the maximum.  That is, the return
    * value is \f$\max_i ( \| data_i \| )\f$, where the max is over the data
    * elements where \f$cvol_i > 0\f$.  If the control volume is NULL, it is
    * ignored during the computation of the maximum.
    */
   double maxNorm(
      const tbox::Pointer< pdat::FaceData<DIM,TYPE> >& data,
      const hier::Box<DIM>& box,
      const tbox::Pointer< pdat::FaceData<DIM,double> > cvol = NULL) const;

   /**
    * Return the dot product of the two data arrays using the control volume
    * to weight the contribution of each product to the sum.  That is, the
    * return value is the sum \f$\sum_i ( data1_i * data2_i * cvol_i )\f$.
    * If the control volume is NULL, it is ignored during the summation.
    */
   TYPE dot(
      const tbox::Pointer< pdat::FaceData<DIM,TYPE> >& data1,
      const tbox::Pointer< pdat::FaceData<DIM,TYPE> >& data2,
      const hier::Box<DIM>& box,
      const tbox::Pointer< pdat::FaceData<DIM,double> > cvol = NULL) const;

   /**
    * Return the integral of the function represented by the data array.
    * The return value is the sum \f$\sum_i ( data_i * vol_i )\f$.
    */
   TYPE integral(
      const tbox::Pointer< pdat::FaceData<DIM,TYPE> >& data,
      const hier::Box<DIM>& box,
      const tbox::Pointer< pdat::FaceData<DIM,double> > vol) const;

private:
   // The following are not implemented:
   PatchFaceDataNormOpsReal(
      const PatchFaceDataNormOpsReal<DIM,TYPE>&);
   void operator=(const PatchFaceDataNormOpsReal<DIM,TYPE>&);

   ArrayDataNormOpsReal<DIM,TYPE> d_array_ops;
};

}
}
#endif

#ifdef INCLUDE_TEMPLATE_IMPLEMENTATION
#include "PatchFaceDataNormOpsReal.C"
#endif

//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-4-4/source/patchdata/operators/time_interpolate/face/FaceFloatLinearTimeInterpolateOp.h $
// Package:	SAMRAI patchdata
// Copyright:   (c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 2132 $
// Modified:	$LastChangedDate: 2008-04-14 14:51:47 -0700 (Mon, 14 Apr 2008) $
// Description: Linear time interp operator for face-centered float patch data.
//

#ifndef included_pdat_FaceFloatLinearTimeInterpolateOp
#define included_pdat_FaceFloatLinearTimeInterpolateOp

#include "SAMRAI_config.h"
#include "tbox/Pointer.h"
#ifndef included_String
#include <string>
#define included_String
#endif
#include "TimeInterpolateOperator.h"

namespace SAMRAI {
    namespace pdat {

/**
 * Class FaceFloatLinearTimeInterpolateOp<DIM> implements standard
 * linear time interpolation for face-centered float patch data. 
 * It is derived from the xfer::TimeInterpolateOperator<DIM> base class.
 * The interpolation uses FORTRAN numerical routines.
 *
 * The findCoarsenOperator() operator function returns true if the input
 * variable is a face-centered float, and the string is 
 * "STD_LINEAR_TIME_INTERPOLATE".
 * 
 * @see xfer::TimeInterpolateOperator
 */

template<int DIM> class FaceFloatLinearTimeInterpolateOp 
: public xfer::TimeInterpolateOperator<DIM>
{
public:
   /**
    * Uninteresting default constructor.
    */
   FaceFloatLinearTimeInterpolateOp();

   /**
    * Uninteresting virtual destructor.
    */
   virtual ~FaceFloatLinearTimeInterpolateOp<DIM>();

   /**
    * Return true if the variable and name string match the standard
    * face-centered float interpolation; otherwise, return false.
    */
   bool findTimeInterpolateOperator(const tbox::Pointer< hier::Variable<DIM> >& var,
                                    const std::string &op_name) const;

   /**
    * Perform linear time interpolation between two face-centered float
    * patch data sources and place result in the destination patch data.
    * Time interpolation is performed on the intersection of the destination 
    * patch data and the input box.  The time to which data is interpolated 
    * is provided by the destination data.
    */
   void timeInterpolate(hier::PatchData<DIM>& dst_data,
                        const hier::Box<DIM>& where,
                        const hier::PatchData<DIM>& src_data_old,
                        const hier::PatchData<DIM>& src_data_new) const;

private:

};

}
}
#endif

#ifdef INCLUDE_TEMPLATE_IMPLEMENTATION
#include "FaceFloatLinearTimeInterpolateOp.C"
#endif

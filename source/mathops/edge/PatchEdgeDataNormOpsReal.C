//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-4-4/source/mathops/edge/PatchEdgeDataNormOpsReal.C $
// Package:	SAMRAI mathops
// Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 2249 $
// Modified:	$LastChangedDate: 2008-07-03 08:17:20 -0700 (Thu, 03 Jul 2008) $
// Description:	Templated norm operations for real edge-centered patch data.
//

#ifndef included_math_PatchEdgeDataNormOpsReal_C
#define included_math_PatchEdgeDataNormOpsReal_C

#include "PatchEdgeDataNormOpsReal.h"

#include "tbox/MathUtilities.h"
#include "tbox/Utilities.h"
#include "EdgeGeometry.h"

namespace SAMRAI {
    namespace math {

template<int DIM, class TYPE>
PatchEdgeDataNormOpsReal<DIM,TYPE>::PatchEdgeDataNormOpsReal()
{
}

template<int DIM, class TYPE>
PatchEdgeDataNormOpsReal<DIM,TYPE>::~PatchEdgeDataNormOpsReal()
{
}

/*
*************************************************************************
*                                                                       *
* The const constructor and assignment operator are not actually used   *
* but are defined here for compilers that require an implementation for *
* every declaration.                                                    *
*                                                                       *
*************************************************************************
*/

template<int DIM, class TYPE>
PatchEdgeDataNormOpsReal<DIM,TYPE>::PatchEdgeDataNormOpsReal(
   const PatchEdgeDataNormOpsReal<DIM,TYPE>& foo)
{
   NULL_USE(foo); 
}

template<int DIM, class TYPE>
void PatchEdgeDataNormOpsReal<DIM,TYPE>::operator=(
   const PatchEdgeDataNormOpsReal<DIM,TYPE>& foo)
{
   NULL_USE(foo); 
}

/*
*************************************************************************
*                                                                       *
* Compute the number of data entries on a patch in the given box.       *
*                                                                       *
*************************************************************************
*/

template<int DIM, class TYPE>
int PatchEdgeDataNormOpsReal<DIM,TYPE>::numberOfEntries(
   const tbox::Pointer< pdat::EdgeData<DIM,TYPE> >& data,
   const hier::Box<DIM>& box) const
{
   int retval = 0;
   const hier::Box<DIM> ibox = box * data->getGhostBox();
   for (int d = 0; d < DIM; d++) {
      const hier::Box<DIM> dbox = pdat::EdgeGeometry<DIM>::toEdgeBox(ibox, d);
      retval += (dbox.size() * data->getDepth());
   }
   return( retval );
}

/*
*************************************************************************
*                                                                       *
* Templated norm operations for real edge-centered data.                *
*                                                                       *
*************************************************************************
*/

template<int DIM, class TYPE>
double PatchEdgeDataNormOpsReal<DIM,TYPE>::sumControlVolumes(
   const tbox::Pointer< pdat::EdgeData<DIM,TYPE> >& data,
   const tbox::Pointer< pdat::EdgeData<DIM,double> >& cvol,
   const hier::Box<DIM>& box) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!data.isNull() && !cvol.isNull());
#endif
   double retval = 0.0;
   for (int d = 0; d < DIM; d++) {
      const hier::Box<DIM> edge_box = pdat::EdgeGeometry<DIM>::toEdgeBox(box, d);
      retval += d_array_ops.sumControlVolumes(data->getArrayData(d),
                                              cvol->getArrayData(d),
                                              edge_box);
   }
   return( retval );
}

template<int DIM, class TYPE>
void PatchEdgeDataNormOpsReal<DIM,TYPE>::abs(
   tbox::Pointer< pdat::EdgeData<DIM,TYPE> >& dst,
   const tbox::Pointer< pdat::EdgeData<DIM,TYPE> >& src,
   const hier::Box<DIM>& box) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!dst.isNull() && !src.isNull());
#endif
   for (int d = 0; d < DIM; d++) {
      const hier::Box<DIM> edge_box = pdat::EdgeGeometry<DIM>::toEdgeBox(box, d);
      d_array_ops.abs(dst->getArrayData(d),
                      src->getArrayData(d),
                      edge_box);
   }
}

template<int DIM, class TYPE>
double PatchEdgeDataNormOpsReal<DIM,TYPE>::L1Norm(
   const tbox::Pointer< pdat::EdgeData<DIM,TYPE> >& data,
   const hier::Box<DIM>& box,
   const tbox::Pointer< pdat::EdgeData<DIM,double> > cvol) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!data.isNull());
#endif
   double retval = 0.0;
   if (cvol.isNull()) {
      for (int d = 0; d < DIM; d++) {
         const hier::Box<DIM> edge_box = pdat::EdgeGeometry<DIM>::toEdgeBox(box, d);
         retval += d_array_ops.L1Norm(data->getArrayData(d), edge_box);
      }
   } else {
      for (int d = 0; d < DIM; d++) {
         const hier::Box<DIM> edge_box = pdat::EdgeGeometry<DIM>::toEdgeBox(box, d);
         retval += d_array_ops.L1NormWithControlVolume(data->getArrayData(d),
                                                       cvol->getArrayData(d),
                                                       edge_box);
      }
   }
   return( retval );
}

template<int DIM, class TYPE>
double PatchEdgeDataNormOpsReal<DIM,TYPE>::L2Norm(
   const tbox::Pointer< pdat::EdgeData<DIM,TYPE> >& data,
   const hier::Box<DIM>& box,
   const tbox::Pointer< pdat::EdgeData<DIM,double> > cvol) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!data.isNull());
#endif
   double retval = 0.0;
   if (cvol.isNull()) {
      for (int d = 0; d < DIM; d++) {
         const hier::Box<DIM> edge_box = pdat::EdgeGeometry<DIM>::toEdgeBox(box, d);
         double aval = d_array_ops.L2Norm(data->getArrayData(d), edge_box);
         retval += aval * aval;
      }
   } else {
      for (int d = 0; d < DIM; d++) {
         const hier::Box<DIM> edge_box = pdat::EdgeGeometry<DIM>::toEdgeBox(box, d);
         double aval = d_array_ops.L2NormWithControlVolume(
                                   data->getArrayData(d), 
                                   cvol->getArrayData(d),
                                   edge_box);
         retval += aval * aval;
      }
   }
   return( sqrt(retval) );
}

template<int DIM, class TYPE>
double PatchEdgeDataNormOpsReal<DIM,TYPE>::weightedL2Norm(
   const tbox::Pointer< pdat::EdgeData<DIM,TYPE> >& data,
   const tbox::Pointer< pdat::EdgeData<DIM,TYPE> >& weight,
   const hier::Box<DIM>& box,
   const tbox::Pointer< pdat::EdgeData<DIM,double> > cvol) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!data.isNull() && !weight.isNull());
#endif
   double retval = 0.0;
   if (cvol.isNull()) {
      for (int d = 0; d < DIM; d++) {
         const hier::Box<DIM> edge_box = pdat::EdgeGeometry<DIM>::toEdgeBox(box, d);
         double aval = d_array_ops.weightedL2Norm(data->getArrayData(d),
                                                  weight->getArrayData(d),
                                                  edge_box);
         retval += aval * aval;
      }
   } else {
      for (int d = 0; d < DIM; d++) {
         const hier::Box<DIM> edge_box = pdat::EdgeGeometry<DIM>::toEdgeBox(box, d);
         double aval = d_array_ops.weightedL2NormWithControlVolume(
                                   data->getArrayData(d),
                                   weight->getArrayData(d),
                                   cvol->getArrayData(d),
                                   edge_box);
         retval += aval * aval;
      }
   }
   return( sqrt(retval) );
}

template<int DIM, class TYPE>
double PatchEdgeDataNormOpsReal<DIM,TYPE>::RMSNorm(
   const tbox::Pointer< pdat::EdgeData<DIM,TYPE> >& data,
   const hier::Box<DIM>& box,
   const tbox::Pointer< pdat::EdgeData<DIM,double> > cvol) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!data.isNull());
#endif
   double retval = L2Norm(data, box, cvol);
   if (cvol.isNull()) {
      retval /= sqrt( (double)numberOfEntries(data, box) );
   } else {
      retval /= sqrt( sumControlVolumes(data, cvol, box) );
   }
   return( retval );
}

template<int DIM, class TYPE>
double PatchEdgeDataNormOpsReal<DIM,TYPE>::weightedRMSNorm(
   const tbox::Pointer< pdat::EdgeData<DIM,TYPE> >& data,
   const tbox::Pointer< pdat::EdgeData<DIM,TYPE> >& weight,
   const hier::Box<DIM>& box,
   const tbox::Pointer< pdat::EdgeData<DIM,double> > cvol) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!data.isNull() && !weight.isNull());
#endif
   double retval = weightedL2Norm(data, weight, box, cvol);
   if (cvol.isNull()) {
      retval /= sqrt( (double)numberOfEntries(data, box) );
   } else {
      retval /= sqrt( sumControlVolumes(data, cvol, box) );
   }
   return( retval );
}

template<int DIM, class TYPE>
double PatchEdgeDataNormOpsReal<DIM,TYPE>::maxNorm(
   const tbox::Pointer< pdat::EdgeData<DIM,TYPE> >& data,
   const hier::Box<DIM>& box,
   const tbox::Pointer< pdat::EdgeData<DIM,double> > cvol) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!data.isNull());
#endif
   double retval = 0.0;
   if (cvol.isNull()) {
      for (int d = 0; d < DIM; d++) {
         const hier::Box<DIM> edge_box = 
            pdat::EdgeGeometry<DIM>::toEdgeBox(box, d);
         retval = tbox::MathUtilities<double>::Max(retval,
                     d_array_ops.maxNorm(data->getArrayData(d), edge_box) );
      }
   } else {
      for (int d = 0; d < DIM; d++) {
         const hier::Box<DIM> edge_box = 
            pdat::EdgeGeometry<DIM>::toEdgeBox(box, d);
         retval = tbox::MathUtilities<double>::Max(retval,
                     d_array_ops.maxNormWithControlVolume(
                        data->getArrayData(d), 
                        cvol->getArrayData(d),
                        edge_box));
      }
   }
   return( retval );
}

template<int DIM, class TYPE>
TYPE PatchEdgeDataNormOpsReal<DIM,TYPE>::dot(
   const tbox::Pointer< pdat::EdgeData<DIM,TYPE> >& data1,
   const tbox::Pointer< pdat::EdgeData<DIM,TYPE> >& data2,
   const hier::Box<DIM>& box,
   const tbox::Pointer< pdat::EdgeData<DIM,double> > cvol) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!data1.isNull() && !data2.isNull());
#endif
   TYPE retval = 0.0;
   if (cvol.isNull()) {
      for (int d = 0; d < DIM; d++) {
         const hier::Box<DIM> edge_box = pdat::EdgeGeometry<DIM>::toEdgeBox(box, d);
         retval += d_array_ops.dot(data1->getArrayData(d),
                                   data2->getArrayData(d),
                                   edge_box);
      }
   } else {
      for (int d = 0; d < DIM; d++) {
         const hier::Box<DIM> edge_box = pdat::EdgeGeometry<DIM>::toEdgeBox(box, d);
         retval += d_array_ops.dotWithControlVolume(
                               data1->getArrayData(d),
                               data2->getArrayData(d),
                               cvol->getArrayData(d),
                               edge_box);
      }
   }
   return( retval );
}

template<int DIM, class TYPE>
TYPE PatchEdgeDataNormOpsReal<DIM,TYPE>::integral(
   const tbox::Pointer< pdat::EdgeData<DIM,TYPE> >& data,
   const hier::Box<DIM>& box,
   const tbox::Pointer< pdat::EdgeData<DIM,double> > vol) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!data.isNull());
#endif
   TYPE retval = 0.0;

   for (int d = 0; d < DIM; d++) {
      const hier::Box<DIM> side_box = pdat::EdgeGeometry<DIM>::toEdgeBox(box, d);
      retval += d_array_ops.integral(
                            data->getArrayData(d),
                            vol->getArrayData(d),
                            side_box);
   }

   return( retval );
}

}
}
#endif

//
// File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-4-4/source/mathops/array/ArrayDataNormOpsComplex.C $
// Package:     SAMRAI mathops
// Copyright:   (c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:    $LastChangedRevision: 1917 $
// Modified:    $LastChangedDate: 2008-01-25 13:28:01 -0800 (Fri, 25 Jan 2008) $
// Description: Norm operations for complex data arrays.
//

#ifndef included_math_ArrayDataNormOpsComplex_C
#define included_math_ArrayDataNormOpsComplex_C

#include "ArrayDataNormOpsComplex.h"

#include "tbox/MathUtilities.h"
#ifdef DEBUG_CHECK_ASSERTIONS
#include "tbox/Utilities.h"
#endif

namespace SAMRAI {
    namespace math {

template<int DIM>  ArrayDataNormOpsComplex<DIM>::ArrayDataNormOpsComplex()
{
}

template<int DIM>  ArrayDataNormOpsComplex<DIM>::~ArrayDataNormOpsComplex()
{
}

/*
*************************************************************************
*                                                                       *
* Norm operations for complex array data.				*
*                                                                       *
*************************************************************************
*/

template<int DIM> void ArrayDataNormOpsComplex<DIM>::abs(
   pdat::ArrayData<DIM,double>& dst,
   const pdat::ArrayData<DIM,dcomplex>& src,
   const hier::Box<DIM>& box) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(dst.getDepth() == src.getDepth());
#endif
   const hier::Box<DIM> dst_box = dst.getBox();
   const hier::Box<DIM> src_box = src.getBox();
   const hier::Box<DIM> ibox = box * dst_box * src_box;

   if (!ibox.empty()) {

      int box_w[DIM];
      int dst_w[DIM];
      int src_w[DIM];
      int dim_counter[DIM];
      for (int i = 0; i < DIM; i++) {
         box_w[i] = ibox.numberCells(i);
         dst_w[i] = dst_box.numberCells(i);
         src_w[i] = src_box.numberCells(i);
         dim_counter[i] = 0;
      }

      const int dst_offset = dst.getOffset();
      const int src_offset = src.getOffset();

      const int num_d0_blocks = ibox.size() / box_w[0];

      int dst_begin = dst_box.offset(ibox.lower());
      int src_begin = src_box.offset(ibox.lower());

      double* dd         = dst.getPointer();
      const dcomplex* sd = src.getPointer();

      const int ddepth = dst.getDepth();
      for (int d = 0; d < ddepth; d++) {

         int dst_counter = dst_begin;
         int src_counter = src_begin;

         int dst_b[DIM];
         int src_b[DIM];
         for (int nd = 0; nd < DIM; nd++) {
            dst_b[nd] = dst_counter;
            src_b[nd] = src_counter;
         }

         for (int nb = 0; nb < num_d0_blocks; nb++) {

            for (int i0 = 0; i0 < box_w[0]; i0++) {
               dd[dst_counter+i0] = sqrt(norm(sd[src_counter+i0]));
            }

            int dim_jump = 0;

            for (int j = 1; j < DIM; j++) {
               if (dim_counter[j] < box_w[j]-1) {
                  ++dim_counter[j];
                  dim_jump = j;
                  break;
               } else {
                  dim_counter[j] = 0;
               }
            }

            if (dim_jump > 0) {
               int dst_step = 1;
               int src_step = 1;
               for (int k = 0; k < dim_jump; k++) {
                  dst_step *= dst_w[k];
                  src_step *= src_w[k];
               }
               dst_counter = dst_b[dim_jump-1] + dst_step;
               src_counter = src_b[dim_jump-1] + src_step;

               for (int m = 0; m < dim_jump; m++) {
                  dst_b[m] = dst_counter;
                  src_b[m] = src_counter;
               }
            }
         }

         dst_begin += dst_offset;
         src_begin += src_offset;

      }
   }
}

template<int DIM> double ArrayDataNormOpsComplex<DIM>::sumControlVolumes(
   const pdat::ArrayData<DIM,dcomplex>& data,
   const pdat::ArrayData<DIM,double>& cvol,
   const hier::Box<DIM>& box) const
{
   double sum = 0.0;

   const hier::Box<DIM> d_box  = data.getBox();
   const hier::Box<DIM> cv_box = cvol.getBox();
   const hier::Box<DIM> ibox = box * d_box * cv_box;

   if (!ibox.empty()) {

      int box_w[DIM];
      int cv_w[DIM];
      int dim_counter[DIM];
      for (int i = 0; i < DIM; i++) {
         box_w[i] = ibox.numberCells(i);
         cv_w[i] = cv_box.numberCells(i);
         dim_counter[i] = 0;
      }

      const int cv_offset = cvol.getOffset();

      int cv_begin = cv_box.offset(ibox.lower());

      const int num_d0_blocks = ibox.size() / box_w[0];

      const int ddepth = cvol.getDepth();
#ifdef DEBUG_CHECK_ASSERTIONS
      TBOX_ASSERT((ddepth == data.getDepth()) || (ddepth == 1));
#endif

      const double* cvd = cvol.getPointer();

      for (int d = 0; d < ddepth; d++) {

         int cv_counter = cv_begin;

         int cv_b[DIM];
         for (int nd = 0; nd < DIM; nd++) {
            cv_b[nd] = cv_counter;
         }

         for (int nb = 0; nb < num_d0_blocks; nb++) {

            for (int i0 = 0; i0 < box_w[0]; i0++) {
               sum += cvd[cv_counter+i0];
            }
            int dim_jump = 0;

            for (int j = 1; j < DIM; j++) {
               if (dim_counter[j] < box_w[j]-1) {
                  ++dim_counter[j];
                  dim_jump = j;
                  break;
               } else {
                  dim_counter[j] = 0;
               }
            }
            if (dim_jump > 0) {
               int cv_step = 1;
               for (int k = 0; k < dim_jump; k++) {
                  cv_step *= cv_w[k];
               }
               cv_counter = cv_b[dim_jump-1] + cv_step;

               for (int m = 0; m < dim_jump; m++) {
                  cv_b[m] = cv_counter;
               }
            }
         }
         cv_begin += cv_offset;
      }

      if (ddepth != data.getDepth()) sum *= data.getDepth();

   }

   return(sum);
}

template<int DIM> double ArrayDataNormOpsComplex<DIM>::L1NormWithControlVolume(
   const pdat::ArrayData<DIM,dcomplex>& data,
   const pdat::ArrayData<DIM,double>& cvol,
   const hier::Box<DIM>& box) const
{
   double l1norm = 0.0;

   const hier::Box<DIM> d_box = data.getBox();
   const hier::Box<DIM> cv_box = cvol.getBox();
   const hier::Box<DIM> ibox = box * d_box * cv_box;

   if (!ibox.empty()) {
      const int ddepth  = data.getDepth();
      const int cvdepth = cvol.getDepth();
#ifdef DEBUG_CHECK_ASSERTIONS
      TBOX_ASSERT((ddepth == cvdepth) || (cvdepth == 1));
#endif

      int box_w[DIM];
      int d_w[DIM];
      int cv_w[DIM];
      int dim_counter[DIM];
      for (int i = 0; i < DIM; i++) {
         box_w[i] = ibox.numberCells(i);
         d_w[i] = d_box.numberCells(i);
         cv_w[i] = cv_box.numberCells(i);
         dim_counter[i] = 0;
      }

      const int d_offset = data.getOffset();
      const int cv_offset = ( (cvdepth == 1) ? 0 : cvol.getOffset() );

      int d_begin = d_box.offset(ibox.lower());
      int cv_begin = cv_box.offset(ibox.lower());

      const int num_d0_blocks = ibox.size() / box_w[0];

      const dcomplex* dd = data.getPointer();
      const double* cvd  = cvol.getPointer();

      for (int d = 0; d < ddepth; d++) {

         int d_counter = d_begin;
         int cv_counter = cv_begin;

         int d_b[DIM];
         int cv_b[DIM];
         for (int nd = 0; nd < DIM; nd++) {
            d_b[nd] = d_counter;
            cv_b[nd] = cv_counter;
         }

         for (int nb = 0; nb < num_d0_blocks; nb++) {

            for (int i0 = 0; i0 < box_w[0]; i0++) {
               l1norm += sqrt(norm(dd[d_counter+i0])) * cvd[cv_counter+i0];
            }
            int dim_jump = 0;

            for (int j = 1; j < DIM; j++) {
               if (dim_counter[j] < box_w[j]-1) {
                  ++dim_counter[j];
                  dim_jump = j;
                  break;
               } else {
                  dim_counter[j] = 0;
               }
            }

            if (dim_jump > 0) {
               int d_step = 1;
               int cv_step = 1;
               for (int k = 0; k < dim_jump; k++) {
                  d_step *= d_w[k];
                  cv_step *= cv_w[k];
               }
               d_counter = d_b[dim_jump-1] + d_step;
               cv_counter = cv_b[dim_jump-1] + cv_step;

               for (int m = 0; m < dim_jump; m++) {
                  d_b[m] = d_counter;
                  cv_b[m] = cv_counter;
               }
            }
         }

         d_begin += d_offset;
         cv_begin += cv_offset;
      }
   }

   return(l1norm);
}

template<int DIM> double ArrayDataNormOpsComplex<DIM>::L1Norm(
   const pdat::ArrayData<DIM,dcomplex>& data,
   const hier::Box<DIM>& box) const
{
   double l1norm = 0.0;

   const hier::Box<DIM> d_box = data.getBox();
   const hier::Box<DIM> ibox = box * d_box;

   if (!ibox.empty()) {

      int box_w[DIM];
      int d_w[DIM];
      int dim_counter[DIM];
      for (int i = 0; i < DIM; i++) {
         box_w[i] = ibox.numberCells(i);
         d_w[i] = d_box.numberCells(i);
         dim_counter[i] = 0;
      }

      const int d_offset = data.getOffset();

      int d_begin = d_box.offset(ibox.lower());

      const int num_d0_blocks = ibox.size() / box_w[0];

      const int ddepth = data.getDepth();

      const dcomplex* dd    = data.getPointer();

      for (int d = 0; d < ddepth; d++) {

         int d_counter = d_begin;

         int d_b[DIM];
         for (int nd = 0; nd < DIM; nd++) {
            d_b[nd] = d_counter;
         }

         for (int nb = 0; nb < num_d0_blocks; nb++) {

            for (int i0 = 0; i0 < box_w[0]; i0++) {
               l1norm += sqrt(norm(dd[d_counter+i0]));
            }
            int dim_jump = 0;

            for (int j = 1; j < DIM; j++) {
               if (dim_counter[j] < box_w[j]-1) {
                  ++dim_counter[j];
                  dim_jump = j;
                  break;
               } else {
                  dim_counter[j] = 0;
               }
            }

            if (dim_jump > 0) {
               int d_step = 1;
               for (int k = 0; k < dim_jump; k++) {
                  d_step *= d_w[k];
               }
               d_counter = d_b[dim_jump-1] + d_step;

               for (int m = 0; m < dim_jump; m++) {
                  d_b[m] = d_counter;
               }
            }
         }

         d_begin += d_offset;
      }
   }

   return(l1norm);
}

template<int DIM> dcomplex ArrayDataNormOpsComplex<DIM>::dotWithControlVolume(
   const pdat::ArrayData<DIM,dcomplex>& data1,
   const pdat::ArrayData<DIM,dcomplex>& data2,
   const pdat::ArrayData<DIM,double>& cvol,
   const hier::Box<DIM>& box) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(data1.getDepth() == data2.getDepth());
#endif

   dcomplex dprod = dcomplex(0.0, 0.0);

   const hier::Box<DIM> d1_box = data1.getBox();
   const hier::Box<DIM> d2_box = data2.getBox();
   const hier::Box<DIM> cv_box = cvol.getBox();
   const hier::Box<DIM> ibox = box * d1_box * d2_box * cv_box;

   if (!ibox.empty()) {
      const int d1depth = data1.getDepth();
      const int cvdepth = cvol.getDepth();
#ifdef DEBUG_CHECK_ASSERTIONS
      TBOX_ASSERT(d1depth == data2.getDepth());
      TBOX_ASSERT((d1depth == cvdepth) || (cvdepth == 1));
#endif

      int box_w[DIM];
      int d1_w[DIM];
      int d2_w[DIM];
      int cv_w[DIM];
      int dim_counter[DIM];
      for (int i = 0; i < DIM; i++) {
         box_w[i] = ibox.numberCells(i);
         d1_w[i] = d1_box.numberCells(i);
         d2_w[i] = d2_box.numberCells(i);
         cv_w[i] = cv_box.numberCells(i);
         dim_counter[i] = 0;
      }

      const int d1_offset = data1.getOffset();
      const int d2_offset = data2.getOffset();
      const int cv_offset = ( (cvdepth == 1) ? 0 : cvol.getOffset() );

      const int num_d0_blocks = ibox.size() / box_w[0];

      int d1_begin = d1_box.offset(ibox.lower());
      int d2_begin = d2_box.offset(ibox.lower());
      int cv_begin = cv_box.offset(ibox.lower());

      const dcomplex* dd1   = data1.getPointer();
      const dcomplex* dd2   = data2.getPointer();
      const double* cvd = cvol.getPointer();

      for (int d = 0; d < d1depth; d++) {

         int d1_counter = d1_begin;
         int d2_counter = d2_begin;
         int cv_counter = cv_begin;

         int d1_b[DIM];
         int d2_b[DIM];
         int cv_b[DIM];
         for (int nd = 0; nd < DIM; nd++) {
            d1_b[nd] = d1_counter;
            d2_b[nd] = d2_counter;
            cv_b[nd] = cv_counter;
         }

         for (int nb = 0; nb < num_d0_blocks; nb++) {

            for (int i0 = 0; i0 < box_w[0]; i0++) {
               dprod += dd1[d1_counter+i0] *
                        conj(dd2[d2_counter+i0]) * cvd[cv_counter+i0];
            }
            int dim_jump = 0;

            for (int j = 1; j < DIM; j++) {
               if (dim_counter[j] < box_w[j]-1) {
                  ++dim_counter[j];
                  dim_jump = j;
                  break;
               } else {
                  dim_counter[j] = 0;
               }
            }

            if (dim_jump > 0) {
               int d1_step = 1;
               int d2_step = 1;
               int cv_step = 1;
               for (int k = 0; k < dim_jump; k++) {
                  d1_step *= d1_w[k];
                  d2_step *= d2_w[k];
                  cv_step *= cv_w[k];
               }
               d1_counter = d1_b[dim_jump-1] + d1_step;
               d2_counter = d2_b[dim_jump-1] + d2_step;
               cv_counter = cv_b[dim_jump-1] + cv_step;

               for (int m = 0; m < dim_jump; m++) {
                  d1_b[m] = d1_counter;
                  d2_b[m] = d2_counter;
                  cv_b[m] = cv_counter;
               }
            }
         }

         d1_begin += d1_offset;
         d2_begin += d2_offset;
         cv_begin += cv_offset;

      }

   }

   return(dprod);
}

template<int DIM> dcomplex ArrayDataNormOpsComplex<DIM>::dot(
   const pdat::ArrayData<DIM,dcomplex>& data1,
   const pdat::ArrayData<DIM,dcomplex>& data2,
   const hier::Box<DIM>& box) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(data1.getDepth() == data2.getDepth());
#endif

   dcomplex dprod = dcomplex(0.0, 0.0);

   const hier::Box<DIM> d1_box = data1.getBox();
   const hier::Box<DIM> d2_box = data2.getBox();
   const hier::Box<DIM> ibox = box * d1_box * d2_box;

   if (!ibox.empty()) {
      const int d1depth = data1.getDepth();
#ifdef DEBUG_CHECK_ASSERTIONS
      TBOX_ASSERT(d1depth == data2.getDepth());
#endif

      int box_w[DIM];
      int d1_w[DIM];
      int d2_w[DIM];
      int dim_counter[DIM];
      for (int i = 0; i < DIM; i++) {
         box_w[i] = ibox.numberCells(i);
         d1_w[i] = d1_box.numberCells(i);
         d2_w[i] = d2_box.numberCells(i);
         dim_counter[i] = 0;
      }

      const int d1_offset = data1.getOffset();
      const int d2_offset = data2.getOffset();

      const int num_d0_blocks = ibox.size() / box_w[0];

      int d1_begin = d1_box.offset(ibox.lower());
      int d2_begin = d2_box.offset(ibox.lower());

      const dcomplex* dd1   = data1.getPointer();
      const dcomplex* dd2   = data2.getPointer();

      for (int d = 0; d < d1depth; d++) {

         int d1_counter = d1_begin;
         int d2_counter = d2_begin;

         int d1_b[DIM];
         int d2_b[DIM];
         for (int nd = 0; nd < DIM; nd++) {
            d1_b[nd] = d1_counter;
            d2_b[nd] = d2_counter;
         }

         for (int nb = 0; nb < num_d0_blocks; nb++) {

            for (int i0 = 0; i0 < box_w[0]; i0++) {
               dprod += dd1[d1_counter+i0] * conj(dd2[d2_counter+i0]);
            }
            int dim_jump = 0;

            for (int j = 1; j < DIM; j++) {
               if (dim_counter[j] < box_w[j]-1) {
                  ++dim_counter[j];
                  dim_jump = j;
                  break;
               } else {
                  dim_counter[j] = 0;
               }
            }

            if (dim_jump > 0) {
               int d1_step = 1;
               int d2_step = 1;
               for (int k = 0; k < dim_jump; k++) {
                  d1_step *= d1_w[k];
                  d2_step *= d2_w[k];
               }
               d1_counter = d1_b[dim_jump-1] + d1_step;
               d2_counter = d2_b[dim_jump-1] + d2_step;

               for (int m = 0; m < dim_jump; m++) {
                  d1_b[m] = d1_counter;
                  d2_b[m] = d2_counter;
               }
            }
         }

         d1_begin += d1_offset;
         d2_begin += d2_offset;

      }

   }

   return(dprod);
}

template<int DIM> dcomplex ArrayDataNormOpsComplex<DIM>::integral(
   const pdat::ArrayData<DIM,dcomplex>& data,
   const pdat::ArrayData<DIM,double>& vol,
   const hier::Box<DIM>& box) const
{

   dcomplex integral = dcomplex(0.0,0.0);

   const hier::Box<DIM> d_box = data.getBox();
   const hier::Box<DIM> v_box = vol.getBox();
   const hier::Box<DIM> ibox = box * d_box * v_box;

   if (!ibox.empty()) {
      const int ddepth = data.getDepth();
      const int vdepth = vol.getDepth();
#ifdef DEBUG_CHECK_ASSERTIONS
      TBOX_ASSERT((ddepth == vdepth) || (vdepth == 1));
#endif

      int box_w[DIM];
      int d_w[DIM];
      int v_w[DIM];
      int dim_counter[DIM];
      for (int i = 0; i < DIM; i++) {
         box_w[i] = ibox.numberCells(i);
         d_w[i] = d_box.numberCells(i);
         v_w[i] = v_box.numberCells(i);
         dim_counter[i] = 0;
      }

      const int d_offset = data.getOffset();
      const int v_offset = ( (vdepth == 1) ? 0 : vol.getOffset() );

      const int num_d0_blocks = ibox.size() / box_w[0];

      int d_begin = d_box.offset(ibox.lower());
      int v_begin = v_box.offset(ibox.lower());

      const dcomplex* dd    = data.getPointer();
      const double* vd = vol.getPointer();

      for (int d = 0; d < ddepth; d++) {

         int d_counter = d_begin;
         int v_counter = v_begin;

         int d_b[DIM];
         int v_b[DIM];
         for (int nd = 0; nd < DIM; nd++) {
            d_b[nd] = d_counter;
            v_b[nd] = v_counter;
         }

         for (int nb = 0; nb < num_d0_blocks; nb++) {

            for (int i0 = 0; i0 < box_w[0]; i0++) {
               integral += dd[d_counter+i0] * vd[v_counter+i0];
            }
            int dim_jump = 0;

            for (int j = 1; j < DIM; j++) {
               if (dim_counter[j] < box_w[j]-1) {
                  ++dim_counter[j];
                  dim_jump = j;
                  break;
               } else {
                  dim_counter[j] = 0;
               }
            }

            if (dim_jump > 0) {
               int d_step = 1;
               int v_step = 1;
               for (int k = 0; k < dim_jump; k++) {
                  d_step *= d_w[k];
                  v_step *= v_w[k];
               }
               d_counter = d_b[dim_jump-1] + d_step;
               v_counter = v_b[dim_jump-1] + v_step;

               for (int m = 0; m < dim_jump; m++) {
                  d_b[m] = d_counter;
                  v_b[m] = v_counter;
               }
            }
         }

         d_begin += d_offset;
         v_begin += v_offset;
      }
   }

   return(integral);
}

template<int DIM> double ArrayDataNormOpsComplex<DIM>::L2NormWithControlVolume(
   const pdat::ArrayData<DIM,dcomplex>& data,
   const pdat::ArrayData<DIM,double>& cvol,
   const hier::Box<DIM>& box) const
{
   return( sqrt(real( dotWithControlVolume(data, data, cvol, box) )) );
}

template<int DIM> double ArrayDataNormOpsComplex<DIM>::L2Norm(
   const pdat::ArrayData<DIM,dcomplex>& data,
   const hier::Box<DIM>& box) const
{
   return( sqrt(real( dot(data, data, box) )) );
}

template<int DIM> double ArrayDataNormOpsComplex<DIM>::weightedL2NormWithControlVolume(
   const pdat::ArrayData<DIM,dcomplex>& data,
   const pdat::ArrayData<DIM,dcomplex>& wgt,
   const pdat::ArrayData<DIM,double>& cvol,
   const hier::Box<DIM>& box) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(data.getDepth() == wgt.getDepth());
#endif

   double wl2norm = 0.0;

   const hier::Box<DIM> d_box  = data.getBox();
   const hier::Box<DIM> w_box  = wgt.getBox();
   const hier::Box<DIM> cv_box = cvol.getBox();
   const hier::Box<DIM> ibox = box * d_box * w_box * cv_box;

   if (!ibox.empty()) {
      const int ddepth = data.getDepth();
      const int cvdepth = cvol.getDepth();
#ifdef DEBUG_CHECK_ASSERTIONS
      TBOX_ASSERT((ddepth == cvdepth) || (cvdepth == 1));
#endif

      int box_w[DIM];
      int w_w[DIM];
      int d_w[DIM];
      int cv_w[DIM];
      int dim_counter[DIM];
      for (int i = 0; i < DIM; i++) {
         box_w[i] = ibox.numberCells(i);
         w_w[i] = w_box.numberCells(i);
         d_w[i] = d_box.numberCells(i);
         cv_w[i] = cv_box.numberCells(i);
         dim_counter[i] = 0;
      }

      const int d_offset = data.getOffset();
      const int w_offset = wgt.getOffset();
      const int cv_offset = ( (cvdepth == 1) ? 0 : cvol.getOffset() );

      const int num_d0_blocks = ibox.size() / box_w[0];

      int d_begin = d_box.offset(ibox.lower());
      int w_begin = w_box.offset(ibox.lower());
      int cv_begin = cv_box.offset(ibox.lower());

      const dcomplex* dd    = data.getPointer();
      const dcomplex* wd    = wgt.getPointer();
      const double* cvd = cvol.getPointer();

      for (int d = 0; d < ddepth; d++) {

         int d_counter = d_begin;
         int w_counter = w_begin;
         int cv_counter = cv_begin;

         int d_b[DIM];
         int w_b[DIM];
         int cv_b[DIM];
         for (int nd = 0; nd < DIM; nd++) {
            d_b[nd] = d_counter;
            w_b[nd] = w_counter;
            cv_b[nd] = cv_counter;
         }

         for (int nb = 0; nb < num_d0_blocks; nb++) {

            for (int i0 = 0; i0 < box_w[0]; i0++) {
               dcomplex val = dd[d_counter+i0] * wd[w_counter+i0];
               wl2norm += real( val * conj(val) ) * cvd[cv_counter+i0];
            }
            int dim_jump = 0;

            for (int j = 1; j < DIM; j++) {
               if (dim_counter[j] < box_w[j]-1) {
                  ++dim_counter[j];
                  dim_jump = j;
                  break;
               } else {
                  dim_counter[j] = 0;
               }
            }

            if (dim_jump > 0) {
               int d_step = 1;
               int w_step = 1;
               int cv_step = 1;
               for (int k = 0; k < dim_jump; k++) {
                  d_step *= d_w[k];
                  w_step *= w_w[k];
                  cv_step *= cv_w[k];
               }
               d_counter = d_b[dim_jump-1] + d_step;
               w_counter = w_b[dim_jump-1] + w_step;
               cv_counter = cv_b[dim_jump-1] + cv_step;

               for (int m = 0; m < dim_jump; m++) {
                  d_b[m] = d_counter;
                  w_b[m] = w_counter;
                  cv_b[m] = cv_counter;
               }
            }
         }

         d_begin += d_offset;
         w_begin += w_offset;
         cv_begin += cv_offset;
      }
   }

   return(sqrt(wl2norm));
}

template<int DIM> double ArrayDataNormOpsComplex<DIM>::weightedL2Norm(
   const pdat::ArrayData<DIM,dcomplex>& data,
   const pdat::ArrayData<DIM,dcomplex>& wgt,
   const hier::Box<DIM>& box) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(data.getDepth() == wgt.getDepth());
#endif

   double wl2norm = 0.0;

   const hier::Box<DIM> d_box  = data.getBox();
   const hier::Box<DIM> w_box  = wgt.getBox();
   const hier::Box<DIM> ibox = box * d_box * w_box;

   if (!ibox.empty()) {

      int box_w[DIM];
      int w_w[DIM];
      int d_w[DIM];
      int dim_counter[DIM];
      for (int i = 0; i < DIM; i++) {
         box_w[i] = ibox.numberCells(i);
         w_w[i] = w_box.numberCells(i);
         d_w[i] = d_box.numberCells(i);
         dim_counter[i] = 0;
      }

      const int d_offset = data.getOffset();
      const int w_offset = wgt.getOffset();

      const int num_d0_blocks = ibox.size() / box_w[0];

      int d_begin = d_box.offset(ibox.lower());
      int w_begin = w_box.offset(ibox.lower());

      const dcomplex* dd    = data.getPointer();
      const dcomplex* wd    = wgt.getPointer();

      const int ddepth = data.getDepth();

      for (int d = 0; d < ddepth; d++) {

         int d_counter = d_begin;
         int w_counter = w_begin;

         int d_b[DIM];
         int w_b[DIM];
         for (int nd = 0; nd < DIM; nd++) {
            d_b[nd] = d_counter;
            w_b[nd] = w_counter;
         }

         for (int nb = 0; nb < num_d0_blocks; nb++) {
            for (int i0 = 0; i0 < box_w[0]; i0++) {
               dcomplex val = dd[d_counter+i0] * wd[w_counter+i0];
               wl2norm += real(val * conj(val));
            }
            int dim_jump = 0;

            for (int j = 1; j < DIM; j++) {
               if (dim_counter[j] < box_w[j]-1) {
                  ++dim_counter[j];
                  dim_jump = j;
                  break;
               } else {
                  dim_counter[j] = 0;
               }
            }

            if (dim_jump > 0) {
               int d_step = 1;
               int w_step = 1;
               for (int k = 0; k < dim_jump; k++) {
                  d_step *= d_w[k];
                  w_step *= w_w[k];
               }
               d_counter = d_b[dim_jump-1] + d_step;
               w_counter = w_b[dim_jump-1] + w_step;

               for (int m = 0; m < dim_jump; m++) {
                  d_b[m] = d_counter;
                  w_b[m] = w_counter;
               }
            }
         }

         d_begin += d_offset;
         w_begin += w_offset;
      }
   }

   return(sqrt(wl2norm));
}

template<int DIM> double ArrayDataNormOpsComplex<DIM>::maxNormWithControlVolume(
   const pdat::ArrayData<DIM,dcomplex>& data,
   const pdat::ArrayData<DIM,double>& cvol,
   const hier::Box<DIM>& box) const
{
   double maxnorm = 0.0;

   const hier::Box<DIM> d_box = data.getBox();
   const hier::Box<DIM> cv_box = cvol.getBox();
   const hier::Box<DIM> ibox = box * d_box * cv_box;

   if (!ibox.empty()) {
      const int ddepth  = data.getDepth();
      const int cvdepth = cvol.getDepth();
#ifdef DEBUG_CHECK_ASSERTIONS
      TBOX_ASSERT((ddepth == cvdepth) || (cvdepth == 1));
#endif

      int box_w[DIM];
      int d_w[DIM];
      int cv_w[DIM];
      int dim_counter[DIM];
      for (int i = 0; i < DIM; i++) {
         box_w[i] = ibox.numberCells(i);
         d_w[i] = d_box.numberCells(i);
         cv_w[i] = cv_box.numberCells(i);
         dim_counter[i] = 0;
      }

      const int d_offset = data.getOffset();
      const int cv_offset = ( (cvdepth == 1) ? 0 : cvol.getOffset() );

      const int num_d0_blocks = ibox.size() / box_w[0];

      int d_begin = d_box.offset(ibox.lower());
      int cv_begin = cv_box.offset(ibox.lower());

      const dcomplex* dd    = data.getPointer();
      const double* cvd = cvol.getPointer();

      for (int d = 0; d < ddepth; d++) {

         int d_counter = d_begin;
         int cv_counter = cv_begin;

         int d_b[DIM];
         int cv_b[DIM];
         for (int nd = 0; nd < DIM; nd++) {
            d_b[nd] = d_counter;
            cv_b[nd] = cv_counter;
         }

         for (int nb = 0; nb < num_d0_blocks; nb++) {

            for (int i0 = 0; i0 < box_w[0]; i0++) {
               if (cvd[cv_counter+i0] > 0.0) {
                  maxnorm = tbox::MathUtilities<double>::Max( maxnorm,
                                                  norm(dd[d_counter+i0]) );
               }
            }
            int dim_jump = 0;

            for (int j = 1; j < DIM; j++) {
               if (dim_counter[j] < box_w[j]-1) {
                  ++dim_counter[j];
                  dim_jump = j;
                  break;
               } else {
                  dim_counter[j] = 0;
               }
            }

            if (dim_jump > 0) {
               int d_step = 1;
               int cv_step = 1;
               for (int k = 0; k < dim_jump; k++) {
                  d_step *= d_w[k];
                  cv_step *= cv_w[k];
               }
               d_counter = d_b[dim_jump-1] + d_step;
               cv_counter = cv_b[dim_jump-1] + cv_step;

               for (int m = 0; m < dim_jump; m++) {
                  d_b[m] = d_counter;
                  cv_b[m] = cv_counter;
               }
            }
         }

         d_begin += d_offset;
         cv_begin += cv_offset;
      }
   }

   return( sqrt(maxnorm) );
}

template<int DIM> double ArrayDataNormOpsComplex<DIM>::maxNorm(
   const pdat::ArrayData<DIM,dcomplex>& data,
   const hier::Box<DIM>& box) const
{
   double maxnorm = 0.0;

   const hier::Box<DIM> d_box = data.getBox();
   const hier::Box<DIM> ibox = box * d_box;

   if (!ibox.empty()) {

      int box_w[DIM];
      int d_w[DIM];
      int dim_counter[DIM];
      for (int i = 0; i < DIM; i++) {
         box_w[i] = ibox.numberCells(i);
         d_w[i] = d_box.numberCells(i);
         dim_counter[i] = 0;
      }

      const int d_offset = data.getOffset();

      const int num_d0_blocks = ibox.size() / box_w[0];

      int d_begin = d_box.offset(ibox.lower());

      const dcomplex* dd    = data.getPointer();

      const int ddepth = data.getDepth();
      for (int d = 0; d < ddepth; d++) {

         int d_counter = d_begin;

         int d_b[DIM];
         for (int nd = 0; nd < DIM; nd++) {
            d_b[nd] = d_counter;
         }

         for (int nb = 0; nb < num_d0_blocks; nb++) {

            for (int i0 = 0; i0 < box_w[0]; i0++) {
               maxnorm = tbox::MathUtilities<double>::Max( maxnorm,
                                               norm(dd[d_counter+i0]) );
            }
            int dim_jump = 0;

            for (int j = 1; j < DIM; j++) {
               if (dim_counter[j] < box_w[j]-1) {
                  ++dim_counter[j];
                  dim_jump = j;
                  break;
               } else {
                  dim_counter[j] = 0;
               }
            }

            if (dim_jump > 0) {
               int d_step = 1;
               for (int k = 0; k < dim_jump; k++) {
                  d_step *= d_w[k];
               }
               d_counter = d_b[dim_jump-1] + d_step;

               for (int m = 0; m < dim_jump; m++) {
                  d_b[m] = d_counter;
               }
            }
         }

         d_begin += d_offset;
      }
   }

   return( sqrt(maxnorm) );
}

}
}
#endif

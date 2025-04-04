c
c
c  File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-4-4/source/patchdata/fortran/pdat_m4arrdim2d.i $
c  Package:     SAMRAI patchdata
c  Copyright:   (c) 1997-2008 Lawrence Livermore National Security, LLC
c  Revision:    $LastChangedRevision: 1917 $
c  Modified:    $LastChangedDate: 2008-01-25 13:28:01 -0800 (Fri, 25 Jan 2008) $
c  Description: m4 include file for dimensioning 2d arrays in FORTRAN routines.
c
c
      subroutine mathops_array_dotwithcv2d(
     &     data1,
     &     data1_ilower0,data1_iupper0,
     &     data1_ilower1,data1_iupper1,
     &     data1_depth,
     &     data2,
     &     data2_ilower0,data2_iupper0,
     &     data2_ilower1,data2_iupper1,
     &     data2_depth,
     &     cv,
     &     cv_ilower0,cv_iupper0,
     &     cv_ilower1,cv_iupper1,
     &     ilower0,iupper0,
     &     ilower1,iupper1,
     &     dprod)
c
      implicit none
c
c     Input.
c
      integer data1_ilower0,data1_iupper0
      integer data1_ilower1,data1_iupper1
      integer data1_depth

      integer data2_ilower0,data2_iupper0
      integer data2_ilower1,data2_iupper1
      integer data2_depth

      integer cv_ilower0,cv_iupper0
      integer cv_ilower1,cv_iupper1

      integer ilower0,iupper0
      integer ilower1,iupper1

      double precision data1(data1_ilower0:data1_iupper0,
     &          data1_ilower1:data1_iupper1,0:data1_depth-1)
      double precision data2(data2_ilower0:data2_iupper0,
     &          data2_ilower1:data2_iupper1,0:data2_depth-1)
      double precision cv(cv_ilower0:cv_iupper0,
     &          cv_ilower1:cv_iupper1)
c
c     Input/Output.
c
      double precision dprod
c
c     Local variables.
c
      integer i0,i1,d
c
c     Compute the dot product with control volume.
c
      dprod = 0.d0

      do d = 0,MIN(data1_depth-1,data2_depth-1)
         do i1 = ilower1,iupper1
            do i0 = ilower0,iupper0
               dprod = dprod + data1(i0,i1,d) *
     &                         data2(i0,i1,d) *
     &                         cv(i0,i1)
            enddo
         enddo
      enddo
c
      return
      end
c

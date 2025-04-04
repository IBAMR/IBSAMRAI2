c
define(REAL,`double precision')dnl
define(INTEGER,`integer')dnl
include(SAMRAI_FORTDIR/pdat_m4arrdim3d.i)dnl
c
      subroutine mathops_array_dotwithcv3d(
     &     data1,
     &     data1_ilower0,data1_iupper0,
     &     data1_ilower1,data1_iupper1,
     &     data1_ilower2,data1_iupper2,
     &     data1_depth,
     &     data2,
     &     data2_ilower0,data2_iupper0,
     &     data2_ilower1,data2_iupper1,
     &     data2_ilower2,data2_iupper2,
     &     data2_depth,
     &     cv,
     &     cv_ilower0,cv_iupper0,
     &     cv_ilower1,cv_iupper1,
     &     cv_ilower2,cv_iupper2,
     &     ilower0,iupper0,
     &     ilower1,iupper1,
     &     ilower2,iupper2,
     &     dprod)
c
      implicit none
c
c     Input.
c
      INTEGER data1_ilower0,data1_iupper0
      INTEGER data1_ilower1,data1_iupper1
      INTEGER data1_ilower2,data1_iupper2
      INTEGER data1_depth

      INTEGER data2_ilower0,data2_iupper0
      INTEGER data2_ilower1,data2_iupper1
      INTEGER data2_ilower2,data2_iupper2
      INTEGER data2_depth

      INTEGER cv_ilower0,cv_iupper0
      INTEGER cv_ilower1,cv_iupper1
      INTEGER cv_ilower2,cv_iupper2

      INTEGER ilower0,iupper0
      INTEGER ilower1,iupper1
      INTEGER ilower2,iupper2

      REAL data1(CELL3d(data1_ilower,data1_iupper,0),0:data1_depth-1)
      REAL data2(CELL3d(data2_ilower,data2_iupper,0),0:data2_depth-1)
      double precision cv(CELL3d(cv_ilower,cv_iupper,0))
c
c     Input/Output.
c
      REAL dprod
c
c     Local variables.
c
      INTEGER i0,i1,i2,d
c
c     Compute the dot product with control volume.
c
      dprod = 0.d0

      do d = 0,MIN(data1_depth-1,data2_depth-1)
         do i2 = ilower2,iupper2
            do i1 = ilower1,iupper1
               do i0 = ilower0,iupper0
                  dprod = dprod + data1(i0,i1,i2,d) *
     &                            data2(i0,i1,i2,d) *
     &                            cv(i0,i1,i2)
               enddo
            enddo
         enddo
      enddo
c
      return
      end
c

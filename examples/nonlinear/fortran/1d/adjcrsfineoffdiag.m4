define(NDIM,1)dnl
include(SAMRAI_FORTDIR/pdat_m4arrdim1d.i)dnl

      subroutine adjcrsfineoffdiag1d(
     & flo0, fhi0,
     & clo0, chi0,
     & direction, side, 
     & offdiag0 )

      implicit none

      integer clo0

      integer chi0

      integer flo0

      integer fhi0

      integer direction
      integer side

      double precision offdiag0(FACE1d(flo,fhi,0))

      integer i

      double precision factor
      parameter      ( factor=2.0d0/3.0d0)

      if ( direction .eq. 0 ) then
         i = clo0-side+1
         offdiag0(i) = offdiag0(i)*factor
      endif

      return
      end

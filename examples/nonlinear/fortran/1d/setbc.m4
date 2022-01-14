define(NDIM,1)dnl
include(SAMRAI_FORTDIR/pdat_m4arrdim1d.i)dnl

      subroutine setbc(
     & lo0, hi0, ghostcells,
     & u,
     & bdrySegLo, bdrySegHi,
     & bdrySide )

c  Set boundary conditions in modified Bratu problem.

      implicit none

      integer lo0

      integer hi0

      integer ghostcells

      integer bdrySide

      integer bdrySegHi(0:NDIM-1)
      integer bdrySegLo(0:NDIM-1)

      double precision u(CELL1d(lo,hi,ghostcells))

      integer ihi
      integer ilo

      double precision zero
      parameter      ( zero=0.0d0 )

      if (bdrySide .eq. 0) then

         ilo = bdrySegLo(0)+1
         u(ilo-1) = zero

      else if (bdrySide .eq. 1) then

         ihi = bdrySegHi(0)-1
         u(ihi+1) = zero

      end if

      return 
      end

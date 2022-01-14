define(NDIM,1)dnl
include(SAMRAI_FORTDIR/pdat_m4arrdim1d.i)dnl

      subroutine evalFaceFluxes(
     & lo0, hi0, ghostcells,
     & dew,
     & u,
     & dx,
     & gew )

c  Evaluate face-centered fluxes in div-grad operator.

      implicit none

      integer lo0

      integer hi0

      integer ghostcells

      double precision dx(0:NDIM-1)

      double precision dew(SIDE1d(lo,hi,0))
      double precision u(CELL1d(lo,hi,ghostcells))
      double precision gew(SIDE1d(lo,hi,0))

      integer i

      do i = lo0, hi0+1
         gew(i) = dew(i)*(u(i) - u(i-1))/dx(0)
      end do      

      return
      end

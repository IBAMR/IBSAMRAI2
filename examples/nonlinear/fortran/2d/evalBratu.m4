define(NDIM,2)dnl
include(SAMRAI_FORTDIR/pdat_m4arrdim2d.i)dnl

      subroutine evalbratu(
     & lo0, hi0, lo1, hi1, ghostcells,
     & gew, gns,
     & f, lexpu,
     & v,
     & u,
     & dx, dt,
     & r )

c  Evaluate modified bratu problem at u by assembling fluxes
c  (in gew, gns), sources (in f and lexpu), and previous
c  time step (in v).

      implicit none

      integer lo0
      integer lo1

      integer hi0
      integer hi1

      integer ghostcells

      double precision dt

      double precision dx(0:NDIM-1)

      double precision gew(SIDE2d0(lo,hi,0))
      double precision gns(SIDE2d1(lo,hi,0))
      double precision f(CELL2d(lo,hi,0))
      double precision lexpu(CELL2d(lo,hi,0))
      double precision v(CELL2d(lo,hi,0))
      double precision u(CELL2d(lo,hi,ghostcells))
      double precision r(CELL2d(lo,hi,0))
     
      integer i
      integer j

      do j = lo1, hi1
         do i = lo0, hi0
            r(i,j) = u(i,j) - v(i,j)
     &               - dt*((gew(i+1,j) - gew(i,j))/dx(0) +
     &                     (gns(i,j+1) - gns(i,j))/dx(1) +
     &                     lexpu(i,j) + f(i,j))
         end do
      end do

      return
      end

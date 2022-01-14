define(NDIM,3)dnl
include(SAMRAI_FORTDIR/pdat_m4arrdim3d.i)dnl

      subroutine evalSource(
     & lo0, hi0, lo1, hi1, lo2, hi2,
     & lambda,
     & xlo, xhi, dx,
     & t,
     & f )

c  Evaluate source term in modified bratu problem.

      implicit none

      integer lo0
      integer lo1
      integer lo2

      integer hi0
      integer hi1
      integer hi2

      double precision lambda
      double precision t

      double precision xlo(0:NDIM-1)
      double precision xhi(0:NDIM-1)
      double precision dx(0:NDIM-1)

      double precision f(CELL3d(lo,hi,0))

      integer i
      integer j
      integer k

      double precision vol
      double precision xi, xterm
      double precision yj, yterm
      double precision zk, zterm

      double precision half,       one,       two
      parameter      ( half=0.5d0, one=1.0d0, two=2.0d0 )

      vol = dx(0)*dx(1)*dx(2)
      zk = xlo(2) + dx(2)*half
      do k = lo2, hi2
         zterm = zk*(one - zk)
         yj = xlo(1) + dx(1)*half
         do j = lo1, hi1
            yterm = yj*(one - yj)
            xi = xlo(0) + dx(0)*half
            do i = lo0, hi0
               xterm = xi*(one - xi)
               f(i,j,k) = xterm*yterm*zterm
     &                  + t*(two*(xterm*zterm+xterm*yterm+yterm*zterm))
     &                  - lambda*exp(t*xterm*yterm*zterm)
               xi = xi + dx(0)
            end do
            yj = yj + dx(1)
         end do
         zk = zk + dx(2)
      end do

      return
      end

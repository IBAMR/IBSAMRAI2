c
c  File:        algs_upfluxsum2d.m4
c  Package:     SAMRAI algorithms
c  Copyright:   (c) 1997-2008 Lawrence Livermore National Security, LLC
c  Revision:    $LastChangedRevision: 1917 $
c  Modified:    $LastChangedDate: 2008-01-25 13:28:01 -0800 (Fri, 25 Jan 2008) $
c  Description: F77 routines for updating 2d flux sums from fluxes.
c
define(SAMRAI_FORTDIR,../../patchdata/fortran)dnl
define(REAL,`double precision')dnl
include(SAMRAI_FORTDIR/pdat_m4arrdim2d.i)dnl
c
c***********************************************************************
c Add flux integrals to fluxsums
c***********************************************************************
c
define(upfluxsumface_case_2d,`dnl
      subroutine upfluxsumface2d$1(
     &  ilo0,ilo1,ihi0,ihi1,
     &  flxgc0,flxgc1,
     &  iface,
     &  flux,fluxsum)
c***********************************************************************
      implicit none
c
      integer
     &  ilo0,ilo1,ihi0,ihi1,
     &  flxgc0,flxgc1,
     &  iface
      REAL
     &  flux(SIDE2d$1VECG(ilo,ihi,flxgc)),
     &  fluxsum(OUTERSIDE2d$1(ilo,ihi,0))
      integer ie$1,ic$2
c
c***********************************************************************
c
      if (iface.eq.0) then
        ie$1 = ilo$1
      else 
        ie$1 = ihi$1+1
      endif 
 
      do ic$2=ilo$2,ihi$2
         fluxsum(ic$2)=fluxsum(ic$2)+flux(ie$1,ic$2)
      enddo
c
      return
      end
')dnl
define(upfluxsumside_case_2d,`dnl
      subroutine upfluxsumside2d$1(
     &  ilo0,ilo1,ihi0,ihi1,
     &  flxgc0,flxgc1,
     &  iface,
     &  flux,fluxsum)
c***********************************************************************
      implicit none
c
      integer
     &  ilo0,ilo1,ihi0,ihi1,
     &  flxgc0,flxgc1,
     &  iface
      REAL
     &  flux(SIDE2d$1VECG(ilo,ihi,flxgc)),
     &  fluxsum(OUTERSIDE2d$1(ilo,ihi,0))
      integer ic$1,ic$2
c
c***********************************************************************
c
      if (iface.eq.0) then
        ic$1 = ilo$1
      else 
        ic$1 = ihi$1+1
      endif 
 
      do ic$2=ilo$2,ihi$2
         fluxsum(ic$2)=fluxsum(ic$2)+flux(ic0,ic1)
      enddo
c
      return
      end
')dnl
upfluxsumface_case_2d(0,1)dnl
c
upfluxsumface_case_2d(1,0)dnl
c
upfluxsumside_case_2d(0,1)dnl
c
upfluxsumside_case_2d(1,0)dnl
c

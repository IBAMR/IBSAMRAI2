c
c  File:        algs_upfluxsum3d.m4
c  Package:     SAMRAI algorithms
c  Copyright:   (c) 1997-2008 Lawrence Livermore National Security, LLC
c  Revision:    $LastChangedRevision: 1917 $
c  Modified:    $LastChangedDate: 2008-01-25 13:28:01 -0800 (Fri, 25 Jan 2008) $
c  Description: F77 routines for updating 3d flux sums from fluxes.
c
c
c  File:        pdat_m4arrdim3d.i
c  Package:     SAMRAI patchdata
c  Copyright:   (c) 1997-2008 Lawrence Livermore National Security, LLC
c  Revision:    $LastChangedRevision: 1917 $
c  Modified:    $LastChangedDate: 2008-01-25 13:28:01 -0800 (Fri, 25 Jan 2008) $
c  Description: m4 include file for dimensioning 3d arrays in FORTRAN routines.
c
c
c***********************************************************************
c Add flux integrals to fluxsums
c***********************************************************************
c
c
c
      subroutine upfluxsumface3d0(
     &  ilo0,ilo1,ilo2,ihi0,ihi1,ihi2,
     &  flxgc0,flxgc1,flxgc2,
     &  iface,
     &  flux,fluxsum)
c***********************************************************************
      implicit none
c
      integer
     &  ilo0,ilo1,ilo2,ihi0,ihi1,ihi2,
     &  flxgc0,flxgc1,flxgc2,
     &  iface
      double precision
     &  flux(ilo0-flxgc0:ihi0+1+flxgc0,
     &          ilo1-flxgc1:ihi1+flxgc1,
     &          ilo2-flxgc2:ihi2+flxgc2),
     &  fluxsum(ilo1:ihi1,
     &          ilo2:ihi2)
      integer ie0,ic1,ic2
c
c***********************************************************************
c
      if (iface.eq.0) then
        ie0 = ilo0
      else
        ie0 = ihi0+1
      endif 
 
      do ic2=ilo2,ihi2
         do ic1=ilo1,ihi1
            fluxsum(ic1,ic2)=fluxsum(ic1,ic2)+flux(ie0,ic1,ic2)
         enddo
      enddo
c
      return
      end
c
      subroutine upfluxsumface3d1(
     &  ilo0,ilo1,ilo2,ihi0,ihi1,ihi2,
     &  flxgc0,flxgc1,flxgc2,
     &  iface,
     &  flux,fluxsum)
c***********************************************************************
      implicit none
c
      integer
     &  ilo0,ilo1,ilo2,ihi0,ihi1,ihi2,
     &  flxgc0,flxgc1,flxgc2,
     &  iface
      double precision
     &  flux(ilo0-flxgc0:ihi0+flxgc0,
     &          ilo1-flxgc1:ihi1+1+flxgc1,
     &          ilo2-flxgc2:ihi2+flxgc2),
     &  fluxsum(ilo0:ihi0,
     &          ilo2:ihi2)
      integer ie1,ic2,ic0
c
c***********************************************************************
c
      if (iface.eq.0) then
        ie1 = ilo1
      else
        ie1 = ihi1+1
      endif 
 
      do ic0=ilo0,ihi0
         do ic2=ilo2,ihi2
            fluxsum(ic2,ic0)=fluxsum(ic2,ic0)+flux(ie1,ic2,ic0)
         enddo
      enddo
c
      return
      end
c
      subroutine upfluxsumface3d2(
     &  ilo0,ilo1,ilo2,ihi0,ihi1,ihi2,
     &  flxgc0,flxgc1,flxgc2,
     &  iface,
     &  flux,fluxsum)
c***********************************************************************
      implicit none
c
      integer
     &  ilo0,ilo1,ilo2,ihi0,ihi1,ihi2,
     &  flxgc0,flxgc1,flxgc2,
     &  iface
      double precision
     &  flux(ilo0-flxgc0:ihi0+flxgc0,
     &          ilo1-flxgc1:ihi1+flxgc1,
     &          ilo2-flxgc2:ihi2+1+flxgc2),
     &  fluxsum(ilo0:ihi0,
     &          ilo1:ihi1)
      integer ie2,ic0,ic1
c
c***********************************************************************
c
      if (iface.eq.0) then
        ie2 = ilo2
      else
        ie2 = ihi2+1
      endif 
 
      do ic1=ilo1,ihi1
         do ic0=ilo0,ihi0
            fluxsum(ic0,ic1)=fluxsum(ic0,ic1)+flux(ie2,ic0,ic1)
         enddo
      enddo
c
      return
      end
c
      subroutine upfluxsumside3d0(
     &  ilo0,ilo1,ilo2,ihi0,ihi1,ihi2,
     &  flxgc0,flxgc1,flxgc2,
     &  iside,
     &  flux,fluxsum)
c***********************************************************************
      implicit none
c
      integer
     &  ilo0,ilo1,ilo2,ihi0,ihi1,ihi2,
     &  flxgc0,flxgc1,flxgc2,
     &  iside
      double precision
     &  flux(ilo0-flxgc0:ihi0+1+flxgc0,
     &          ilo1-flxgc1:ihi1+flxgc1,
     &          ilo2-flxgc2:ihi2+flxgc2),
     &  fluxsum(ilo1:ihi1,
     &          ilo2:ihi2)
      integer ic0,ic1,ic2
c
c***********************************************************************
c
      if (iside.eq.0) then
        ic0 = ilo0
      else
        ic0 = ihi0+1
      endif 
 
      do ic2=ilo2,ihi2
         do ic1=ilo1,ihi1
            fluxsum(ic1,ic2)=fluxsum(ic1,ic2)+flux(ic0,ic1,ic2)
         enddo
      enddo
c
      return
      end
c
      subroutine upfluxsumside3d1(
     &  ilo0,ilo1,ilo2,ihi0,ihi1,ihi2,
     &  flxgc0,flxgc1,flxgc2,
     &  iside,
     &  flux,fluxsum)
c***********************************************************************
      implicit none
c
      integer
     &  ilo0,ilo1,ilo2,ihi0,ihi1,ihi2,
     &  flxgc0,flxgc1,flxgc2,
     &  iside
      double precision
     &  flux(ilo0-flxgc0:ihi0+flxgc0,
     &          ilo1-flxgc1:ihi1+1+flxgc1,
     &          ilo2-flxgc2:ihi2+flxgc2),
     &  fluxsum(ilo0:ihi0,
     &          ilo2:ihi2)
      integer ic1,ic0,ic2
c
c***********************************************************************
c
      if (iside.eq.0) then
        ic1 = ilo1
      else
        ic1 = ihi1+1
      endif 
 
      do ic2=ilo2,ihi2
         do ic0=ilo0,ihi0
            fluxsum(ic0,ic2)=fluxsum(ic0,ic2)+flux(ic0,ic1,ic2)
         enddo
      enddo
c
      return
      end
c
      subroutine upfluxsumside3d2(
     &  ilo0,ilo1,ilo2,ihi0,ihi1,ihi2,
     &  flxgc0,flxgc1,flxgc2,
     &  iside,
     &  flux,fluxsum)
c***********************************************************************
      implicit none
c
      integer
     &  ilo0,ilo1,ilo2,ihi0,ihi1,ihi2,
     &  flxgc0,flxgc1,flxgc2,
     &  iside
      double precision
     &  flux(ilo0-flxgc0:ihi0+flxgc0,
     &          ilo1-flxgc1:ihi1+flxgc1,
     &          ilo2-flxgc2:ihi2+1+flxgc2),
     &  fluxsum(ilo0:ihi0,
     &          ilo1:ihi1)
      integer ic2,ic0,ic1
c
c***********************************************************************
c
      if (iside.eq.0) then
        ic2 = ilo2
      else
        ic2 = ihi2+1
      endif 
 
      do ic1=ilo1,ihi1
         do ic0=ilo0,ihi0
            fluxsum(ic0,ic1)=fluxsum(ic0,ic1)+flux(ic0,ic1,ic2)
         enddo
      enddo
c
      return
      end

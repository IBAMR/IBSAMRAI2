c
c  File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-4-4/source/apputils/boundary/appu_cartbdry2d.f $
c  Package:     SAMRAI application utilities
c  Copyright:   (c) 1997-2008 Lawrence Livermore National Security, LLC
c  Revision:    $LastChangedRevision: 1917 $
c  Modified:    $LastChangedDate: 2008-01-25 13:28:01 -0800 (Fri, 25 Jan 2008) $
c  Description: F77 routines for Cartesian 2d boundary conditions
c
c
c  File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-4-4/source/apputils/boundary/appu_cartbdry2d.f $
c  Package:     SAMRAI patchdata
c  Copyright:   (c) 1997-2008 Lawrence Livermore National Security, LLC
c  Revision:    $LastChangedRevision: 1917 $
c  Modified:    $LastChangedDate: 2008-01-25 13:28:01 -0800 (Fri, 25 Jan 2008) $
c  Description: m4 include file for dimensioning 2d arrays in FORTRAN routines.
c
c
c general boundary condition cases
c
c
c node boundary conditions
c
c
c edge boundary conditions
c
c***********************************************************************
c***********************************************************************
      subroutine getcartedgebdry2d(
     &  ifirst0,ilast0,ifirst1,ilast1,
     &  ibeg0,iend0,ibeg1,iend1,
     &  ngc0,ngc1,
     &  dx,
     &  bdry_loc,
     &  bdry_cond,
     &  edge_values,
     &  arrdata,
     &  arrdepth)
c***********************************************************************
      implicit none
c
c  File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-4-4/source/apputils/boundary/appu_cartbdry2d.f $
c  Package:     SAMRAI application utilities
c  Copyright:   (c) 1997-2008 Lawrence Livermore National Security, LLC
c  Revision:    $LastChangedRevision: 1917 $
c  Modified:    $LastChangedDate: 2008-01-25 13:28:01 -0800 (Fri, 25 Jan 2008) $
c  Description: m4 include file for 2d boundary constant common blocks
c
      common/cartbdrylocparams2d/
     &  XLEFT,XRIGHT,YLEFT,YRIGHT,
     &  X0Y0,X1Y0,X0Y1,X1Y1
      integer
     &  XLEFT,XRIGHT,YLEFT,YRIGHT,
     &  X0Y0,X1Y0,X0Y1,X1Y1
c
c
      common/cartbdrycondparams2d/
     &  FLOW,XFLOW,YFLOW,
     &  REFLECT,XREFLECT,YREFLECT,
     &  DIRICHLET,XDIRICHLET,YDIRICHLET,
     &  NEUMANN,XNEUMANN,YNEUMANN
      integer
     &  FLOW,XFLOW,YFLOW,
     &  REFLECT,XREFLECT,YREFLECT,
     &  DIRICHLET,XDIRICHLET,YDIRICHLET,
     &  NEUMANN,XNEUMANN,YNEUMANN
c***********************************************************************
      integer ifirst0,ilast0,ifirst1,ilast1
      integer ibeg0,iend0,ibeg1,iend1
      integer ngc0,ngc1
c
      double precision dx(0:2-1)
c
      integer bdry_loc,
     &        arrdepth, 
     &        bdry_cond
      double precision
     &        edge_values(0:arrdepth-1,0:2*2-1)
c
      double precision
     &  arrdata(ifirst0-ngc0:ilast0+ngc0,
     &          ifirst1-ngc1:ilast1+ngc1,0:arrdepth-1)
c
      integer ic1,ic0,ict0,ict1
      integer k
      integer ipivot0,ipivot1
      integer edge_loc
      double precision    dirsign0,dirsign1
c
c***********************************************************************
c   bdry_loc   index     bdry_loc  index   
c       0     (-1, 0)       2     (0,-1)  
c       1     ( 1, 0)       3     (0, 1)  
c***********************************************************************
c***********************************************************************
      edge_loc = bdry_loc
      if ((bdry_loc.eq.XLEFT).or.(bdry_loc.eq.XRIGHT)) then
         if (bdry_loc.eq.XLEFT) then
c                             x0 boundary
            ipivot0  = ifirst0
            dirsign0 = -1.d0
         else
c                             x1 boundary
            ipivot0  = ilast0
            dirsign0 = 1.d0
         endif
           if (bdry_cond.eq.DIRICHLET) then
             do ic1=ibeg1,iend1
               do ic0=ibeg0,iend0
         do k=0,arrdepth-1
           arrdata(ic0,ic1,k)=edge_values(k,edge_loc)
         enddo
               enddo
             enddo
           else if (bdry_cond.eq.NEUMANN) then
             do ic1=ibeg1,iend1
               ict1 = ic1
               do ic0=ibeg0,iend0
                 ict0 = ipivot0
         do k=0,arrdepth-1
           arrdata(ic0,ic1,k)=arrdata(ict0,ict1,k) + dirsign0
     &       * edge_values(k,edge_loc) * dble(ic0-ict0) * dx(0)
         enddo
               enddo
             enddo
           else if (bdry_cond.eq.FLOW) then
             do ic1=ibeg1,iend1
               ict1 = ic1
               do ic0=ibeg0,iend0
                 ict0 = ipivot0
         do k=0,arrdepth-1
           arrdata(ic0,ic1,k)=arrdata(ict0,ict1,k)
         enddo
               enddo
             enddo
           else if (bdry_cond.eq.REFLECT) then
             do ic1=ibeg1,iend1
               ict1 = ic1
               do ic0=ibeg0,iend0
                 ict0 = ipivot0
         do k=0,arrdepth-1
           arrdata(ic0,ic1,k)=arrdata(ict0,ict1,k)
         enddo
         arrdata(ic0,ic1,0)= -arrdata(ic0,ic1,0)
               enddo
             enddo
          else
             write(6,*) "INVALID EDGE bdry_cond in getcartedgebdry2d"
             write(6,*) "bdry_loc = ",bdry_loc,
     &                  "bdry_cond = ",bdry_cond
          endif
      else if ((bdry_loc.eq.YLEFT).or.(bdry_loc.eq.YRIGHT)) then
         if (bdry_loc.eq.YLEFT) then
c                             y0 boundary
            ipivot1  = ifirst1
            dirsign1 = -1.d0
         else
c                             y1 boundary
            ipivot1  = ilast1
            dirsign1 = 1.d0
         endif
           if (bdry_cond.eq.DIRICHLET) then
             do ic1=ibeg1,iend1
               do ic0=ibeg0,iend0
         do k=0,arrdepth-1
           arrdata(ic0,ic1,k)=edge_values(k,edge_loc)
         enddo
               enddo
             enddo
           else if (bdry_cond.eq.NEUMANN) then
             do ic1=ibeg1,iend1
               ict1 = ipivot1
               do ic0=ibeg0,iend0
                 ict0 = ic0
         do k=0,arrdepth-1
           arrdata(ic0,ic1,k)=arrdata(ict0,ict1,k) + dirsign1
     &       * edge_values(k,edge_loc) * dble(ic1-ict1) * dx(1)
         enddo
               enddo
             enddo
           else if (bdry_cond.eq.FLOW) then
             do ic1=ibeg1,iend1
               ict1 = ipivot1
               do ic0=ibeg0,iend0
                 ict0 = ic0
         do k=0,arrdepth-1
           arrdata(ic0,ic1,k)=arrdata(ict0,ict1,k)
         enddo
               enddo
             enddo
           else if (bdry_cond.eq.REFLECT) then
             do ic1=ibeg1,iend1
               ict1 = ipivot1
               do ic0=ibeg0,iend0
                 ict0 = ic0
         do k=0,arrdepth-1
           arrdata(ic0,ic1,k)=arrdata(ict0,ict1,k)
         enddo
         arrdata(ic0,ic1,1)= -arrdata(ic0,ic1,1)
               enddo
             enddo
          else
             write(6,*) "INVALID EDGE bdry_cond in getcartedgebdry2d"
             write(6,*) "bdry_loc = ",bdry_loc,
     &                  "bdry_cond = ",bdry_cond
          endif
      endif
c
      return
      end
c***********************************************************************
c***********************************************************************
      subroutine getcartnodebdry2d(
     &  ifirst0,ilast0,ifirst1,ilast1,
     &  ibeg0,iend0,ibeg1,iend1,
     &  ngc0,ngc1,
     &  dx,
     &  bdry_loc,
     &  bdry_cond,
     &  edge_values,
     &  arrdata,
     &  arrdepth)
c***********************************************************************
      implicit none
c
c  File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-4-4/source/apputils/boundary/appu_cartbdry2d.f $
c  Package:     SAMRAI application utilities
c  Copyright:   (c) 1997-2008 Lawrence Livermore National Security, LLC
c  Revision:    $LastChangedRevision: 1917 $
c  Modified:    $LastChangedDate: 2008-01-25 13:28:01 -0800 (Fri, 25 Jan 2008) $
c  Description: m4 include file for 2d boundary constant common blocks
c
      common/cartbdrylocparams2d/
     &  XLEFT,XRIGHT,YLEFT,YRIGHT,
     &  X0Y0,X1Y0,X0Y1,X1Y1
      integer
     &  XLEFT,XRIGHT,YLEFT,YRIGHT,
     &  X0Y0,X1Y0,X0Y1,X1Y1
c
c
      common/cartbdrycondparams2d/
     &  FLOW,XFLOW,YFLOW,
     &  REFLECT,XREFLECT,YREFLECT,
     &  DIRICHLET,XDIRICHLET,YDIRICHLET,
     &  NEUMANN,XNEUMANN,YNEUMANN
      integer
     &  FLOW,XFLOW,YFLOW,
     &  REFLECT,XREFLECT,YREFLECT,
     &  DIRICHLET,XDIRICHLET,YDIRICHLET,
     &  NEUMANN,XNEUMANN,YNEUMANN
c***********************************************************************
      integer ifirst0,ilast0,ifirst1,ilast1
      integer ibeg0,iend0,ibeg1,iend1
      integer ngc0,ngc1
c
      double precision dx(0:2-1)
c
      integer bdry_loc,
     &        arrdepth,
     &        bdry_cond
      double precision
     &        edge_values(0:arrdepth-1,0:2*2-1)
c
      double precision
     &  arrdata(ifirst0-ngc0:ilast0+ngc0,
     &          ifirst1-ngc1:ilast1+ngc1,0:arrdepth-1)
c
      integer ic1,ic0,ict0,ict1
      integer k
      integer ipivot0,ipivot1
      integer edge_loc0,edge_loc1
      double precision    dirsign0,dirsign1
c
c***********************************************************************
c***********************************************************************
c    bdry_loc   index 
c       0     (-1,-1)
c       1     ( 1,-1)
c       2     (-1, 1)
c       3     ( 1, 1)
c***********************************************************************
c***********************************************************************
      if (bdry_loc.eq.X0Y0) then
         edge_loc0 = XLEFT
         ipivot0   = ifirst0
         dirsign0  = -1.d0
         edge_loc1 = YLEFT
         ipivot1   = ifirst1
         dirsign1  = -1.d0
      else if (bdry_loc.eq.X1Y0) then
         edge_loc0 = XRIGHT
         ipivot0   = ilast0
         dirsign0  = 1.d0
         edge_loc1 = YLEFT
         ipivot1   = ifirst1
         dirsign1  = -1.d0
      else if (bdry_loc.eq.X0Y1) then
         edge_loc0 = XLEFT
         ipivot0   = ifirst0
         dirsign0  = -1.d0
         edge_loc1 = YRIGHT
         ipivot1   = ilast1
         dirsign1  = 1.d0
      else if (bdry_loc.eq.X1Y1) then
         edge_loc0 = XRIGHT
         ipivot0   = ilast0
         dirsign0  = 1.d0
         edge_loc1 = YRIGHT
         ipivot1   = ilast1
         dirsign1  = 1.d0
      endif
           if (bdry_cond.eq.XDIRICHLET) then
             do ic1=ibeg1,iend1
               do ic0=ibeg0,iend0
         do k=0,arrdepth-1
           arrdata(ic0,ic1,k)=edge_values(k,edge_loc0)
         enddo
               enddo
             enddo
           else if (bdry_cond.eq.YDIRICHLET) then
             do ic1=ibeg1,iend1
               do ic0=ibeg0,iend0
         do k=0,arrdepth-1
           arrdata(ic0,ic1,k)=edge_values(k,edge_loc1)
         enddo
               enddo
             enddo
           else if (bdry_cond.eq.XNEUMANN) then
             do ic1=ibeg1,iend1
               ict1 = ic1
                do ic0=ibeg0,iend0
                 ict0 = ipivot0
         do k=0,arrdepth-1
           arrdata(ic0,ic1,k)=arrdata(ict0,ict1,k) + dirsign0
     &       * edge_values(k,edge_loc0) * dble(ic0-ict0) * dx(0)
         enddo
               enddo
             enddo
           else if (bdry_cond.eq.YNEUMANN) then
             do ic1=ibeg1,iend1
               ict1 = ipivot1
                do ic0=ibeg0,iend0
                 ict0 = ic0
         do k=0,arrdepth-1
           arrdata(ic0,ic1,k)=arrdata(ict0,ict1,k) + dirsign1
     &       * edge_values(k,edge_loc1) * dble(ic1-ict1) * dx(1)
         enddo
               enddo
             enddo
           else if (bdry_cond.eq.XFLOW) then
             do ic1=ibeg1,iend1
               ict1 = ic1
               do ic0=ibeg0,iend0
                 ict0 = ipivot0
         do k=0,arrdepth-1
           arrdata(ic0,ic1,k)=arrdata(ict0,ict1,k)
         enddo
                enddo
             enddo
           else if (bdry_cond.eq.YFLOW) then
             do ic1=ibeg1,iend1
               ict1 = ipivot1
               do ic0=ibeg0,iend0
                 ict0 = ic0
         do k=0,arrdepth-1
           arrdata(ic0,ic1,k)=arrdata(ict0,ict1,k)
         enddo
                enddo
             enddo
           else if (bdry_cond.eq.XREFLECT) then
             do ic1=ibeg1,iend1
               ict1 = ic1
               do ic0=ibeg0,iend0
                 ict0 = ipivot0
         do k=0,arrdepth-1
           arrdata(ic0,ic1,k)=arrdata(ict0,ict1,k)
         enddo
         arrdata(ic0,ic1,0)= -arrdata(ic0,ic1,0)
               enddo
             enddo
           else if (bdry_cond.eq.YREFLECT) then
             do ic1=ibeg1,iend1
               ict1 = ipivot1
               do ic0=ibeg0,iend0
                 ict0 = ic0
         do k=0,arrdepth-1
           arrdata(ic0,ic1,k)=arrdata(ict0,ict1,k)
         enddo
         arrdata(ic0,ic1,1)= -arrdata(ic0,ic1,1)
               enddo
             enddo
           else
             write(6,*) "INVALID NODE bdry_cond in getcartnodebdry2d"
             write(6,*) "bdry_loc = ",bdry_loc,
     &                  "bdry_cond = ",bdry_cond
           endif
c
      return
      end

c***********************************************************************
c***********************************************************************
      subroutine stufcartbdryloc2d(
     &  XLEFTin, XRIGHTin, YLEFTin, YRIGHTin,
     &  X0Y0in, X1Y0in, X0Y1in, X1Y1in)
      implicit none
      integer
     &  XLEFTin, XRIGHTin, YLEFTin, YRIGHTin,
     &  X0Y0in, X1Y0in, X0Y1in, X1Y1in
c
c  File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-4-4/source/apputils/boundary/appu_cartbdry2d.f $
c  Package:     SAMRAI application utilities
c  Copyright:   (c) 1997-2008 Lawrence Livermore National Security, LLC
c  Revision:    $LastChangedRevision: 1917 $
c  Modified:    $LastChangedDate: 2008-01-25 13:28:01 -0800 (Fri, 25 Jan 2008) $
c  Description: m4 include file for 2d boundary constant common blocks
c
      common/cartbdrylocparams2d/
     &  XLEFT,XRIGHT,YLEFT,YRIGHT,
     &  X0Y0,X1Y0,X0Y1,X1Y1
      integer
     &  XLEFT,XRIGHT,YLEFT,YRIGHT,
     &  X0Y0,X1Y0,X0Y1,X1Y1
c
c
      common/cartbdrycondparams2d/
     &  FLOW,XFLOW,YFLOW,
     &  REFLECT,XREFLECT,YREFLECT,
     &  DIRICHLET,XDIRICHLET,YDIRICHLET,
     &  NEUMANN,XNEUMANN,YNEUMANN
      integer
     &  FLOW,XFLOW,YFLOW,
     &  REFLECT,XREFLECT,YREFLECT,
     &  DIRICHLET,XDIRICHLET,YDIRICHLET,
     &  NEUMANN,XNEUMANN,YNEUMANN

c 2d edges
        XLEFT=XLEFTin
        XRIGHT=XRIGHTin
        YLEFT=YLEFTin
        YRIGHT=YRIGHTin

c 2d nodes
        X0Y0=X0Y0in
        X1Y0=X1Y0in
        X0Y1=X0Y1in
        X1Y1=X1Y1in

      return
      end
c***********************************************************************
c***********************************************************************
      subroutine stufcartbdrycond2d(
     &  FLOWin, XFLOWin, YFLOWin, 
     &  REFLECTin, XREFLECTin, YREFLECTin,
     &  DIRICHLETin, XDIRICHLETin, YDIRICHLETin,
     &  NEUMANNin, XNEUMANNin, YNEUMANNin)
      implicit none
      integer
     &  FLOWin, XFLOWin, YFLOWin, 
     &  REFLECTin, XREFLECTin, YREFLECTin,
     &  DIRICHLETin, XDIRICHLETin, YDIRICHLETin,
     &  NEUMANNin, XNEUMANNin, YNEUMANNin
c
c  File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-4-4/source/apputils/boundary/appu_cartbdry2d.f $
c  Package:     SAMRAI application utilities
c  Copyright:   (c) 1997-2008 Lawrence Livermore National Security, LLC
c  Revision:    $LastChangedRevision: 1917 $
c  Modified:    $LastChangedDate: 2008-01-25 13:28:01 -0800 (Fri, 25 Jan 2008) $
c  Description: m4 include file for 2d boundary constant common blocks
c
      common/cartbdrylocparams2d/
     &  XLEFT,XRIGHT,YLEFT,YRIGHT,
     &  X0Y0,X1Y0,X0Y1,X1Y1
      integer
     &  XLEFT,XRIGHT,YLEFT,YRIGHT,
     &  X0Y0,X1Y0,X0Y1,X1Y1
c
c
      common/cartbdrycondparams2d/
     &  FLOW,XFLOW,YFLOW,
     &  REFLECT,XREFLECT,YREFLECT,
     &  DIRICHLET,XDIRICHLET,YDIRICHLET,
     &  NEUMANN,XNEUMANN,YNEUMANN
      integer
     &  FLOW,XFLOW,YFLOW,
     &  REFLECT,XREFLECT,YREFLECT,
     &  DIRICHLET,XDIRICHLET,YDIRICHLET,
     &  NEUMANN,XNEUMANN,YNEUMANN

        FLOW=FLOWin
        XFLOW=XFLOWin
        YFLOW=YFLOWin

        REFLECT=REFLECTin
        XREFLECT=XREFLECTin
        YREFLECT=YREFLECTin

        DIRICHLET=DIRICHLETin
        XDIRICHLET=XDIRICHLETin
        YDIRICHLET=YDIRICHLETin

        NEUMANN=NEUMANNin
        XNEUMANN=XNEUMANNin
        YNEUMANN=YNEUMANNin

      return
      end

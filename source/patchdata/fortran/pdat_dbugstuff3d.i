c
c  File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-4-4/source/patchdata/fortran/pdat_dbugstuff3d.i $
c  Package:     SAMRAI patchdata
c  Copyright:   (c) 1997-2008 Lawrence Livermore National Security, LLC
c  Revision:    $LastChangedRevision: 1917 $
c  Modified:    $LastChangedDate: 2008-01-25 13:28:01 -0800 (Fri, 25 Jan 2008) $
c  Description: m4 include file for 3d patchdata debugging routines
c
define(SAMRAI_FORTDIR,.)dnl
include(SAMRAI_FORTDIR/pdat_m4arrdim3d.i)dnl
c
define(pdat_debug_subroutine_head_3d,`dnl
     &  fi0,la0,fi1,la1,fi2,la2,ng,
     &  ibeg0,iend0,ibeg1,iend1,ibeg2,iend2,array)
c     =============================================================
      implicit none
      integer fi0,la0,fi1,la1,fi2,la2,ng,
     &   ibeg0,iend0,ibeg1,iend1,ibeg2,iend2
')dnl
c
define(pdat_debug_body_3d,`dnl

      do $3=$6
         do $2=$5
            do $1=$4 
               write(6,*) "array[",$1,",",$2,",",$3,"] = ",
     &         array($1,$2,$3)
            enddo
         enddo
      enddo

      call flush(6)
      return
      end
')dnl
define(pdat_debug_cell_3d,`dnl
pdat_debug_subroutine_head_3d()dnl
      $1
     &  array(CELL3d(fi,la,ng))
      integer ic0,ic1,ic2
c     =============================================================
pdat_debug_body_3d(`ic0',`ic1',`ic2',`ibeg0,iend0',`ibeg1,iend1',`ibeg2,iend2')dnl
')dnl
define(pdat_debug_face_3d,`dnl
pdat_debug_subroutine_head_3d()dnl
      $1
     &  array(FACE3d$2(fi,la,ng))
      integer ie$2,ic$3,ic$4
c     =============================================================
pdat_debug_body_3d(`ie$2',`ic$3',`ic$4',`ibeg$2,iend$2+1',`ibeg$3,iend$3',`ibeg$4,iend$4')dnl
')dnl
define(pdat_debug_node_3d,`dnl
pdat_debug_subroutine_head_3d()dnl
      $1
     &  array(NODE3d(fi,la,ng))
      integer ie0,ie1,ie2
c     =============================================================
pdat_debug_body_3d(`ie0',`ie1',`ie2',`ibeg0,iend0+1',`ibeg1,iend1+1',`ibeg2,iend2+1')dnl
')dnl
define(pdat_debug_outerface_3d,`dnl
pdat_debug_subroutine_head_3d()dnl
      $1
     &  array(OUTERFACE3d$2(fi,la,ng))
      integer ic$3,ic$4  
c     =============================================================

      do ic$4=ibeg$4,iend$4
         do ic$3=ibeg$3,iend$3
            write(6,*) "array[",ic$3,",",ic$4,"] = ",
     &      array(ic$3,ic$4)
         enddo
      enddo
 
      call flush(6)
      return
      end
')dnl
define(pdat_debug_outerside_3d,`dnl
pdat_debug_subroutine_head_3d()dnl
      $1
     &  array(OUTERSIDE3d$2(fi,la,ng))
      integer ic$3,ic$4  
c     =============================================================

      do ic$4=ibeg$4,iend$4
         do ic$3=ibeg$3,iend$3
            write(6,*) "array[",ic$3,",",ic$4,"] = ",
     &      array(ic$3,ic$4)
         enddo
      enddo
 
      call flush(6)
      return
      end
')dnl
define(pdat_debug_side_3d,`dnl
pdat_debug_subroutine_head_3d()dnl
      $1
     &  array(SIDE3d$2(fi,la,ng))
      integer ie0,ic1,ic2
c     =============================================================
ifelse($2,`0',`pdat_debug_body_3d(`ie0',`ic1',`ic2',`ibeg0,iend0+1',`ibeg1,iend1',`ibeg2,iend2')dnl
',`')dnl
ifelse($2,`1',`pdat_debug_body_3d(`ie0',`ic1',`ic2',`ibeg0,iend0',`ibeg1,iend1+1',`ibeg2,iend2')dnl
',`')dnl
ifelse($2,`2',`pdat_debug_body_3d(`ie0',`ic1',`ic2',`ibeg0,iend0',`ibeg1,iend1',`ibeg2,iend2+1')dnl
',`')dnl
')
define(pdat_debug_edge_3d,`dnl
pdat_debug_subroutine_head_3d()dnl
      $1
     &  array(EDGE3d$2(fi,la,ng))
      integer ie0,ic1,ic2
c     =============================================================
ifelse($2,`0',`pdat_debug_body_3d(`ie0',`ic1',`ic2',`ibeg0,iend0',`ibeg1,iend1+1',`ibeg2,iend2+1')dnl
',`')dnl
ifelse($2,`1',`pdat_debug_body_3d(`ie0',`ic1',`ic2',`ibeg0,iend0+1',`ibeg1,iend1',`ibeg2,iend2+1')dnl
',`')dnl
ifelse($2,`2',`pdat_debug_body_3d(`ie0',`ic1',`ic2',`ibeg0,iend0+1',`ibeg1,iend1+1',`ibeg2,iend2')dnl
',`')dnl
')

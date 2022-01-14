c
c  File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-4-4/source/patchdata/operators/constant/fortran/pdat_m4concoarsenops2d.i $
c  Package:     SAMRAI patchdata
c  Copyright:   (c) 1997-2008 Lawrence Livermore National Security, LLC
c  Revision:    $LastChangedRevision: 1917 $
c  Modified:    $LastChangedDate: 2008-01-25 13:28:01 -0800 (Fri, 25 Jan 2008) $
c  Description: m4 include file for 2d Cartesian coarsen operators
c
define(NDIM,2)dnl
define(SAMRAI_FORTDIR,../../../../patchdata/fortran)dnl
include(SAMRAI_FORTDIR/pdat_m4arrdim2d.i)dnl
define(con_coarsen_op_subroutine_head_2d,`dnl
     &  ifirstc0,ifirstc1,ilastc0,ilastc1,
     &  filo0,filo1,fihi0,fihi1,
     &  cilo0,cilo1,cihi0,cihi1,
     &  ratio,
     &  arrayf,arrayc)
c***********************************************************************
      implicit none
      double precision zero
      parameter (zero=0.d0)
c
      integer
     &  ifirstc0,ifirstc1,ilastc0,ilastc1,
     &  cilo0,cilo1,cihi0,cihi1,
     &  filo0,filo1,fihi0,fihi1
      integer ratio(0:NDIM-1)
')dnl
c
define(conavg_node_body_2d,`dnl
c
c***********************************************************************
c
      do ie1=ifirstc1,ilastc1+1
         if1=ie1*ratio(1)
         do ie0=ifirstc0,ilastc0+1
            arrayc(ie0,ie1)=arrayf(ie0*ratio(0),if1)
         enddo
      enddo
c
      return
      end
')dnl
c
define(conavg_op_node_2d,`dnl
con_coarsen_op_subroutine_head_2d()dnl
      $1
     &  arrayf(NODE2d(filo,fihi,0)),
     &  arrayc(NODE2d(cilo,cihi,0))
      integer ie0,ie1,if1
conavg_node_body_2d()dnl
')dnl
c
define(conavg_op_outernode_2d0,`dnl
con_coarsen_op_subroutine_head_2d()dnl
      $1
     &  arrayf(OUTERNODE2d0(filo,fihi,0)),
     &  arrayc(OUTERNODE2d0(cilo,cihi,0))
      integer ic1,if1
c
c***********************************************************************
c
      do ic1=ifirstc1+1,ilastc1
         if1=ic1*ratio(1)
         arrayc(ic1)=arrayf(if1)
      enddo

      return
      end
')dnl
c

define(conavg_op_outernode_2d1,`dnl 
con_coarsen_op_subroutine_head_2d()dnl 
      $1
     &  arrayf(OUTERNODE2d1(filo,fihi,0)),
     &  arrayc(OUTERNODE2d1(cilo,cihi,0))
      integer ic0,if0
c
c***********************************************************************
c
      do ic0=ifirstc0,ilastc0+1
         if0=ic0*ratio(0)
         arrayc(ic0)=arrayf(if0)
      enddo

      return
      end
')dnl


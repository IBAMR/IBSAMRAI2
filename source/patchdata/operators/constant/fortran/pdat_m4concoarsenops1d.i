c
c  File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-4-4/source/patchdata/operators/constant/fortran/pdat_m4concoarsenops1d.i $
c  Package:     SAMRAI patchdata
c  Copyright:   (c) 1997-2008 Lawrence Livermore National Security, LLC
c  Revision:    $LastChangedRevision: 1917 $
c  Modified:    $LastChangedDate: 2008-01-25 13:28:01 -0800 (Fri, 25 Jan 2008) $
c  Description: m4 include file for 1d constant coarsen operators
c
define(NDIM,1)dnl
define(SAMRAI_FORTDIR,../../../../patchdata/fortran)dnl
include(SAMRAI_FORTDIR/pdat_m4arrdim1d.i)dnl
define(con_coarsen_op_subroutine_head_1d,`dnl
     &  ifirstc0,ilastc0,
     &  filo0,fihi0,
     &  cilo0,cihi0,
     &  ratio,
     &  arrayf,arrayc)
c***********************************************************************
      implicit none
      double precision zero
      parameter (zero=0.d0)
c
      integer
     &  ifirstc0,ilastc0,
     &  filo0,fihi0,
     &  cilo0,cihi0
      integer ratio(0:NDIM-1)
')dnl
c
define(con_injection_body_1d,`dnl
c
c***********************************************************************
c
      do $1=$2
         arrayc($1)=arrayf($1*ratio(0))
      enddo
c
      return
      end
')dnl
c
c
define(conavg_op_node_1d,`dnl
con_coarsen_op_subroutine_head_1d()dnl
      $1
     &  arrayf(NODE1d(filo,fihi,0)),
     &  arrayc(NODE1d(cilo,cihi,0))
      integer ie0
con_injection_body_1d(`ie0',`ifirstc0,ilastc0+1')dnl
')dnl
c

c
define(conavg_op_outernode_1d,`dnl
con_coarsen_op_subroutine_head_1d()dnl
      $1
     &  arrayf(OUTERNODE1d(filo,fihi,0)),
     &  arrayc(OUTERNODE1d(cilo,cihi,0))
c
c***********************************************************************
c
      arrayc(1)=arrayf(1)

      return
      end
')dnl


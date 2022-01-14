c
c  File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-4-4/source/patchdata/operators/time_interpolate/fortran/pdat_m4lintimeintops2d.i $
c  Package:     SAMRAI patchdata
c  Copyright:   (c) 1997-2008 Lawrence Livermore National Security, LLC
c  Revision:    $LastChangedRevision: 1917 $
c  Modified:    $LastChangedDate: 2008-01-25 13:28:01 -0800 (Fri, 25 Jan 2008) $
c  Description: m4 include file for 2d std linear time interpolation operators.
c
define(SAMRAI_FORTDIR,../../../../patchdata/fortran)dnl
include(SAMRAI_FORTDIR/pdat_m4arrdim2d.i)dnl
define(lin_time_int_subroutine_head_2d,`dnl
     &  ifirst0,ifirst1,ilast0,ilast1,
     &  oilo0,oilo1,oihi0,oihi1,
     &  nilo0,nilo1,nihi0,nihi1,
     &  dilo0,dilo1,dihi0,dihi1,
     &  tfrac,
     &  arrayold,arraynew,
     &  arraydst)
c***********************************************************************
      implicit none
      double precision one
      parameter (one=1.d0)
c
      integer
     &  ifirst0,ifirst1,ilast0,ilast1,
     &  oilo0,oilo1,oihi0,oihi1,
     &  nilo0,nilo1,nihi0,nihi1,
     &  dilo0,dilo1,dihi0,dihi1
      double precision
     &  tfrac, oldfrac
')dnl
c
define(lin_time_int_body_2d,`dnl
c
c***********************************************************************
c
      oldfrac=one-tfrac

      do $2=$4
         do $1=$3
            arraydst($1,$2)=arrayold($1,$2)*oldfrac
     &                       +arraynew($1,$2)*tfrac
         enddo
      enddo
c
      return
      end
')dnl
c
define(lin_time_int_op_cell_2d,`dnl
lin_time_int_subroutine_head_2d()dnl
      $1
     &  arrayold(CELL2d(oilo,oihi,0)),
     &  arraynew(CELL2d(nilo,nihi,0)),
     &  arraydst(CELL2d(dilo,dihi,0))
      integer ic0,ic1
lin_time_int_body_2d(`ic0',`ic1',`ifirst0,ilast0',`ifirst1,ilast1')dnl
')dnl
c
define(lin_time_int_op_edge_2d,`dnl
lin_time_int_subroutine_head_2d()dnl
      $1
     &  arrayold(EDGE2d$2(oilo,oihi,0)),
     &  arraynew(EDGE2d$2(nilo,nihi,0)),
     &  arraydst(EDGE2d$2(dilo,dihi,0))
      integer ie0,ie1
lin_time_int_body_2d(`ie0',`ie1',`ifirst0,ilast0+$2',`ifirst1,ilast1+$3')dnl
')dnl
c
define(lin_time_int_op_face_2d,`dnl
lin_time_int_subroutine_head_2d()dnl
      $1
     &  arrayold(FACE2d$2(oilo,oihi,0)),
     &  arraynew(FACE2d$2(nilo,nihi,0)),
     &  arraydst(FACE2d$2(dilo,dihi,0))
      integer ie$2,ic$3
lin_time_int_body_2d(`ie$2',`ic$3',`ifirst$2,ilast$2+1',`ifirst$3,ilast$3')dnl
')dnl
c
define(lin_time_int_op_node_2d,`dnl
lin_time_int_subroutine_head_2d()dnl
      $1
     &  arrayold(NODE2d(oilo,oihi,0)),
     &  arraynew(NODE2d(nilo,nihi,0)),
     &  arraydst(NODE2d(dilo,dihi,0))
      integer ie0,ie1
lin_time_int_body_2d(`ie0',`ie1',`ifirst0,ilast0+1',`ifirst1,ilast1+1')dnl
')dnl
c
define(lin_time_int_op_outerface_2d,`dnl
lin_time_int_subroutine_head_2d()dnl
      $1
     &  arrayold(OUTERFACE2d$2(oilo,oihi,0)),
     &  arraynew(OUTERFACE2d$2(nilo,nihi,0)),
     &  arraydst(OUTERFACE2d$2(dilo,dihi,0))
      integer ic$3 
c
c***********************************************************************
c
      oldfrac=one-tfrac

      do ic$3=ifirst$3,ilast$3
         arraydst(ic$3)=arrayold(ic$3)*oldfrac
     &                +arraynew(ic$3)*tfrac
      enddo
c
      return
      end
')dnl
c
define(lin_time_int_op_outerside_2d,`dnl
lin_time_int_subroutine_head_2d()dnl
      $1
     &  arrayold(OUTERSIDE2d$2(oilo,oihi,0)),
     &  arraynew(OUTERSIDE2d$2(nilo,nihi,0)),
     &  arraydst(OUTERSIDE2d$2(dilo,dihi,0))
      integer ic$3 
c
c***********************************************************************
c
      oldfrac=one-tfrac

      do ic$3=ifirst$3,ilast$3
         arraydst(ic$3)=arrayold(ic$3)*oldfrac
     &                +arraynew(ic$3)*tfrac
      enddo
c
      return
      end
')dnl
c
define(lin_time_int_op_side_2d,`dnl
lin_time_int_subroutine_head_2d()dnl
      $1
     &  arrayold(SIDE2d$2(oilo,oihi,0)),
     &  arraynew(SIDE2d$2(nilo,nihi,0)),
     &  arraydst(SIDE2d$2(dilo,dihi,0))
      integer ie0,ie1
lin_time_int_body_2d(`ie0',`ie1',`ifirst0,ilast0+$3',`ifirst1,ilast1+$2')dnl
')dnl

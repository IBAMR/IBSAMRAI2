c
c  File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-4-4/source/geometry/cartesian/operators/fortran/geom_cartcoarsen3d.m4 $
c  Package:     SAMRAI geometry
c  Copyright:   (c) 1997-2008 Lawrence Livermore National Security, LLC
c  Revision:    $LastChangedRevision: 1917 $
c  Modified:    $LastChangedDate: 2008-01-25 13:28:01 -0800 (Fri, 25 Jan 2008) $
c  Description: FORTRAN routines for spatial coarsening of 3d patch data
c               on a regular Cartesian mesh.
c
include(geom_m4cartcoarsenops3d.i)dnl
c
c***********************************************************************
c Weighted averaging for 3d cell-centered double data
c***********************************************************************
c
      subroutine cartwgtavgcelldoub3d(
cart_wgtavg_op_cell_3d(`double precision')dnl
c
c***********************************************************************
c Weighted averaging for 3d cell-centered float data
c***********************************************************************
c
      subroutine cartwgtavgcellflot3d(
cart_wgtavg_op_cell_3d(`real')dnl
c
c***********************************************************************
c Weighted averaging for 3d cell-centered complex data
c***********************************************************************
c
      subroutine cartwgtavgcellcplx3d(
cart_wgtavg_op_cell_3d(`double complex')dnl
c
c***********************************************************************
c Weighted averaging for 3d edge-centered double data
c***********************************************************************
c
      subroutine cartwgtavgedgedoub3d0(
cart_wgtavg_op_edge_3d(`double precision',0,1,2)dnl
c
      subroutine cartwgtavgedgedoub3d1(
cart_wgtavg_op_edge_3d(`double precision',1,2,0)dnl
c
      subroutine cartwgtavgedgedoub3d2(
cart_wgtavg_op_edge_3d(`double precision',2,0,1)dnl
c
c***********************************************************************
c Weighted averaging for 3d edge-centered float data
c***********************************************************************
c
      subroutine cartwgtavgedgeflot3d0(
cart_wgtavg_op_edge_3d(`real',0,1,2)dnl
c
      subroutine cartwgtavgedgeflot3d1(
cart_wgtavg_op_edge_3d(`real',1,2,0)dnl
c
      subroutine cartwgtavgedgeflot3d2(
cart_wgtavg_op_edge_3d(`real',2,0,1)dnl
c
c***********************************************************************
c Weighted averaging for 3d edge-centered complex data
c***********************************************************************
c
      subroutine cartwgtavgedgecplx3d0(
cart_wgtavg_op_edge_3d(`double complex',0,1,2)dnl
c
      subroutine cartwgtavgedgecplx3d1(
cart_wgtavg_op_edge_3d(`double complex',1,2,0)dnl
c
      subroutine cartwgtavgedgecplx3d2(
cart_wgtavg_op_edge_3d(`double complex',2,0,1)dnl
c
c***********************************************************************
c Weighted averaging for 3d face-centered double data
c***********************************************************************
c
      subroutine cartwgtavgfacedoub3d0(
cart_wgtavg_op_face_3d(`double precision',0,1,2)dnl
c
      subroutine cartwgtavgfacedoub3d1(
cart_wgtavg_op_face_3d(`double precision',1,2,0)dnl
c
      subroutine cartwgtavgfacedoub3d2(
cart_wgtavg_op_face_3d(`double precision',2,0,1)dnl
c
c***********************************************************************
c Weighted averaging for 3d face-centered float data
c***********************************************************************
c
      subroutine cartwgtavgfaceflot3d0(
cart_wgtavg_op_face_3d(`real',0,1,2)dnl
c
      subroutine cartwgtavgfaceflot3d1(
cart_wgtavg_op_face_3d(`real',1,2,0)dnl
c
      subroutine cartwgtavgfaceflot3d2(
cart_wgtavg_op_face_3d(`real',2,0,1)dnl
c
c***********************************************************************
c Weighted averaging for 3d face-centered complex data
c***********************************************************************
c
      subroutine cartwgtavgfacecplx3d0(
cart_wgtavg_op_face_3d(`double complex',0,1,2)dnl
c
      subroutine cartwgtavgfacecplx3d1(
cart_wgtavg_op_face_3d(`double complex',1,2,0)dnl
c
      subroutine cartwgtavgfacecplx3d2(
cart_wgtavg_op_face_3d(`double complex',2,0,1)dnl
c
c***********************************************************************
c Weighted averaging for 3d outerface double data
c***********************************************************************
c
      subroutine cartwgtavgoutfacedoub3d0(
cart_wgtavg_op_outerface_3d(`double precision',0,1,2)dnl
c
      subroutine cartwgtavgoutfacedoub3d1(
cart_wgtavg_op_outerface_3d(`double precision',1,2,0)dnl
c
      subroutine cartwgtavgoutfacedoub3d2(
cart_wgtavg_op_outerface_3d(`double precision',2,0,1)dnl
c
c***********************************************************************
c Weighted averaging for 3d outerface float data
c***********************************************************************
c
      subroutine cartwgtavgoutfaceflot3d0(
cart_wgtavg_op_outerface_3d(`real',0,1,2)dnl
c
      subroutine cartwgtavgoutfaceflot3d1(
cart_wgtavg_op_outerface_3d(`real',1,2,0)dnl
c
      subroutine cartwgtavgoutfaceflot3d2(
cart_wgtavg_op_outerface_3d(`real',2,0,1)dnl
c
c***********************************************************************
c Weighted averaging for 3d outerface complex data
c***********************************************************************
c
      subroutine cartwgtavgoutfacecplx3d0(
cart_wgtavg_op_outerface_3d(`double complex',0,1,2)dnl
c
      subroutine cartwgtavgoutfacecplx3d1(
cart_wgtavg_op_outerface_3d(`double complex',1,2,0)dnl
c
      subroutine cartwgtavgoutfacecplx3d2(
cart_wgtavg_op_outerface_3d(`double complex',2,0,1)dnl
c
c***********************************************************************
c Weighted averaging for 3d outerside double data
c***********************************************************************
c
      subroutine cartwgtavgoutsidedoub3d0(
cart_wgtavg_op_outerside_3d(`double precision',0,1,2)dnl
c
      subroutine cartwgtavgoutsidedoub3d1(
cart_wgtavg_op_outerside_3d(`double precision',1,0,2)dnl
c
      subroutine cartwgtavgoutsidedoub3d2(
cart_wgtavg_op_outerside_3d(`double precision',2,0,1)dnl
c
c***********************************************************************
c Weighted averaging for 3d side-centered double data
c***********************************************************************
c
      subroutine cartwgtavgsidedoub3d0(
cart_wgtavg_op_side_3d(`double precision',0,1,2)dnl
c
      subroutine cartwgtavgsidedoub3d1(
cart_wgtavg_op_side_3d(`double precision',1,2,0)dnl
c
      subroutine cartwgtavgsidedoub3d2(
cart_wgtavg_op_side_3d(`double precision',2,0,1)dnl
c
c***********************************************************************
c Weighted averaging for 3d side-centered float data
c***********************************************************************
c
      subroutine cartwgtavgsideflot3d0(
cart_wgtavg_op_side_3d(`real',0,1,2)dnl
c
      subroutine cartwgtavgsideflot3d1(
cart_wgtavg_op_side_3d(`real',1,2,0)dnl
c
      subroutine cartwgtavgsideflot3d2(
cart_wgtavg_op_side_3d(`real',2,0,1)dnl
c
c***********************************************************************
c Weighted averaging for 3d side-centered complex data
c***********************************************************************
c
      subroutine cartwgtavgsidecplx3d0(
cart_wgtavg_op_side_3d(`double complex',0,1,2)dnl
c
      subroutine cartwgtavgsidecplx3d1(
cart_wgtavg_op_side_3d(`double complex',1,2,0)dnl
c
      subroutine cartwgtavgsidecplx3d2(
cart_wgtavg_op_side_3d(`double complex',2,0,1)dnl
c

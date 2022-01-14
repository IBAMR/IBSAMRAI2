c
c  File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-4-4/source/geometry/cartesian/operators/fortran/geom_cartcoarsen2d.m4 $
c  Package:     SAMRAI geometry
c  Copyright:   (c) 1997-2008 Lawrence Livermore National Security, LLC
c  Revision:    $LastChangedRevision: 1917 $
c  Modified:    $LastChangedDate: 2008-01-25 13:28:01 -0800 (Fri, 25 Jan 2008) $
c  Description: FORTRAN routines for spatial coarsening of 2d patch data
c               on a regular Cartesian mesh.
c
include(geom_m4cartcoarsenops2d.i)dnl
c
c***********************************************************************
c Weighted averaging for 2d cell-centered double data
c***********************************************************************
c
      subroutine cartwgtavgcelldoub2d(
cart_wgtavg_op_cell_2d(`double precision')dnl
c
c***********************************************************************
c Weighted averaging for 2d cell-centered float data
c***********************************************************************
c
      subroutine cartwgtavgcellflot2d(
cart_wgtavg_op_cell_2d(`real')dnl
c
c***********************************************************************
c Weighted averaging for 2d cell-centered complex data
c***********************************************************************
c
      subroutine cartwgtavgcellcplx2d(
cart_wgtavg_op_cell_2d(`double complex')dnl
c
c***********************************************************************
c Weighted averaging for 2d edge-centered double data
c***********************************************************************
c
      subroutine cartwgtavgedgedoub2d0(
cart_wgtavg_op_edge_2d(`double precision',0,1)dnl
c
      subroutine cartwgtavgedgedoub2d1(
cart_wgtavg_op_edge_2d(`double precision',1,0)dnl
c
c***********************************************************************
c Weighted averaging for 2d edge-centered float data
c***********************************************************************
c
      subroutine cartwgtavgedgeflot2d0(
cart_wgtavg_op_edge_2d(`real',0,1)dnl
c
      subroutine cartwgtavgedgeflot2d1(
cart_wgtavg_op_edge_2d(`real',1,0)dnl
c
c***********************************************************************
c Weighted averaging for 2d edge-centered complex data
c***********************************************************************
c
      subroutine cartwgtavgedgecplx2d0(
cart_wgtavg_op_edge_2d(`double complex',0,1)dnl
c
      subroutine cartwgtavgedgecplx2d1(
cart_wgtavg_op_edge_2d(`double complex',1,0)dnl
c
c***********************************************************************
c Weighted averaging for 2d face-centered double data
c***********************************************************************
c
      subroutine cartwgtavgfacedoub2d0(
cart_wgtavg_op_face_2d(`double precision',0,1)dnl
c
      subroutine cartwgtavgfacedoub2d1(
cart_wgtavg_op_face_2d(`double precision',1,0)dnl
c
c***********************************************************************
c Weighted averaging for 2d face-centered float data
c***********************************************************************
c
      subroutine cartwgtavgfaceflot2d0(
cart_wgtavg_op_face_2d(`real',0,1)dnl
c
      subroutine cartwgtavgfaceflot2d1(
cart_wgtavg_op_face_2d(`real',1,0)dnl
c
c***********************************************************************
c Weighted averaging for 2d face-centered complex data
c***********************************************************************
c
      subroutine cartwgtavgfacecplx2d0(
cart_wgtavg_op_face_2d(`double complex',0,1)dnl
c
      subroutine cartwgtavgfacecplx2d1(
cart_wgtavg_op_face_2d(`double complex',1,0)dnl
c
c***********************************************************************
c Weighted averaging for 2d outerface double data
c***********************************************************************
c
      subroutine cartwgtavgoutfacedoub2d0(
cart_wgtavg_op_outerface_2d(`double precision',0,1)dnl
c
      subroutine cartwgtavgoutfacedoub2d1(
cart_wgtavg_op_outerface_2d(`double precision',1,0)dnl
c
c***********************************************************************
c Weighted averaging for 2d outerface float data
c***********************************************************************
c
      subroutine cartwgtavgoutfaceflot2d0(
cart_wgtavg_op_outerface_2d(`real',0,1)dnl
c
      subroutine cartwgtavgoutfaceflot2d1(
cart_wgtavg_op_outerface_2d(`real',1,0)dnl
c
c***********************************************************************
c Weighted averaging for 2d outerface complex data
c***********************************************************************
c
      subroutine cartwgtavgoutfacecplx2d0(
cart_wgtavg_op_outerface_2d(`double complex',0,1)dnl
c
      subroutine cartwgtavgoutfacecplx2d1(
cart_wgtavg_op_outerface_2d(`double complex',1,0)dnl
c
c***********************************************************************
c Weighted averaging for 2d outerside double data
c***********************************************************************
c
      subroutine cartwgtavgoutsidedoub2d0(
cart_wgtavg_op_outerside_2d(`double precision',0,1)dnl
c
      subroutine cartwgtavgoutsidedoub2d1(
cart_wgtavg_op_outerside_2d(`double precision',1,0)dnl
c
c***********************************************************************
c Weighted averaging for 2d side-centered double data
c***********************************************************************
c
      subroutine cartwgtavgsidedoub2d0(
cart_wgtavg_op_side_2d(`double precision',0,1)dnl
c
      subroutine cartwgtavgsidedoub2d1(
cart_wgtavg_op_side_2d(`double precision',1,0)dnl
c
c***********************************************************************
c Weighted averaging for 2d side-centered float data
c***********************************************************************
c
      subroutine cartwgtavgsideflot2d0(
cart_wgtavg_op_side_2d(`real',0,1)dnl
c
      subroutine cartwgtavgsideflot2d1(
cart_wgtavg_op_side_2d(`real',1,0)dnl
c
c***********************************************************************
c Weighted averaging for 2d side-centered complex data
c***********************************************************************
c
      subroutine cartwgtavgsidecplx2d0(
cart_wgtavg_op_side_2d(`double complex',0,1)dnl
c
      subroutine cartwgtavgsidecplx2d1(
cart_wgtavg_op_side_2d(`double complex',1,0)dnl
c

c
c  File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-4-4/source/patchdata/operators/constant/fortran/pdat_conrefine3d.m4 $
c  Package:     SAMRAI patchdata
c  Copyright:   (c) 1997-2008 Lawrence Livermore National Security, LLC
c  Revision:    $LastChangedRevision: 1917 $
c  Modified:    $LastChangedDate: 2008-01-25 13:28:01 -0800 (Fri, 25 Jan 2008) $
c  Description: FORTRAN routines for spatial refining of 3d patch data
c               on a regular Cartesian mesh.
c
include(pdat_m4conrefineops3d.i)dnl
c
c***********************************************************************
c Constant interpolation for 3d cell-centered double data
c***********************************************************************
c
      subroutine conrefcelldoub3d(
conref_op_cell_3d(`double precision')dnl
c
c***********************************************************************
c Constant interpolation for 3d cell-centered float data
c***********************************************************************
c
      subroutine conrefcellflot3d(
conref_op_cell_3d(`real')dnl
c
c***********************************************************************
c Constant interpolation for 3d cell-centered complex data
c***********************************************************************
c
      subroutine conrefcellcplx3d(
conref_op_cell_3d(`double complex')dnl
c
c***********************************************************************
c Constant interpolation for 3d cell-centered integer data
c***********************************************************************
c
      subroutine conrefcellintg3d(
conref_op_cell_3d(`integer')dnl
c
c***********************************************************************
c Constant interpolation for 3d edge-centered double data
c***********************************************************************
c
      subroutine conrefedgedoub3d0(
conref_op_edge_3d(`double precision',0,1,2)dnl
c
      subroutine conrefedgedoub3d1(
conref_op_edge_3d(`double precision',1,2,0)dnl
c
      subroutine conrefedgedoub3d2(
conref_op_edge_3d(`double precision',2,0,1)dnl
c
c***********************************************************************
c Constant interpolation for 3d edge-centered float data
c***********************************************************************
c
      subroutine conrefedgeflot3d0(
conref_op_edge_3d(`real',0,1,2)dnl
c
      subroutine conrefedgeflot3d1(
conref_op_edge_3d(`real',1,2,0)dnl
c
      subroutine conrefedgeflot3d2(
conref_op_edge_3d(`real',2,0,1)dnl
c
c***********************************************************************
c Constant interpolation for 3d edge-centered complex data
c***********************************************************************
c
      subroutine conrefedgecplx3d0(
conref_op_edge_3d(`double complex',0,1,2)dnl
c
      subroutine conrefedgecplx3d1(
conref_op_edge_3d(`double complex',1,2,0)dnl
c
      subroutine conrefedgecplx3d2(
conref_op_edge_3d(`double complex',2,0,1)dnl
c
c***********************************************************************
c Constant interpolation for 3d edge-centered integer data
c***********************************************************************
c
      subroutine conrefedgeintg3d0(
conref_op_edge_3d(`integer',0,1,2)dnl
c
      subroutine conrefedgeintg3d1(
conref_op_edge_3d(`integer',1,2,0)dnl
c
      subroutine conrefedgeintg3d2(
conref_op_edge_3d(`integer',2,0,1)dnl
c
c***********************************************************************
c Constant interpolation for 3d face-centered double data
c***********************************************************************
c
      subroutine conreffacedoub3d0(
conref_op_face_3d(`double precision',0,1,2)dnl
c
      subroutine conreffacedoub3d1(
conref_op_face_3d(`double precision',1,2,0)dnl
c
      subroutine conreffacedoub3d2(
conref_op_face_3d(`double precision',2,0,1)dnl
c
c***********************************************************************
c Constant interpolation for 3d face-centered float data
c***********************************************************************
c
      subroutine conreffaceflot3d0(
conref_op_face_3d(`real',0,1,2)dnl
c
      subroutine conreffaceflot3d1(
conref_op_face_3d(`real',1,2,0)dnl
c
      subroutine conreffaceflot3d2(
conref_op_face_3d(`real',2,0,1)dnl
c
c***********************************************************************
c Constant interpolation for 3d face-centered complex data
c***********************************************************************
c
      subroutine conreffacecplx3d0(
conref_op_face_3d(`double complex',0,1,2)dnl
c
      subroutine conreffacecplx3d1(
conref_op_face_3d(`double complex',1,2,0)dnl
c
      subroutine conreffacecplx3d2(
conref_op_face_3d(`double complex',2,0,1)dnl
c
c***********************************************************************
c Constant interpolation for 3d face-centered integer data
c***********************************************************************
c
      subroutine conreffaceintg3d0(
conref_op_face_3d(`integer',0,1,2)dnl
c
      subroutine conreffaceintg3d1(
conref_op_face_3d(`integer',1,2,0)dnl
c
      subroutine conreffaceintg3d2(
conref_op_face_3d(`integer',2,0,1)dnl
c
c***********************************************************************
c Constant interpolation for 3d outerface double data
c***********************************************************************
c
      subroutine conrefoutfacedoub3d0(
conref_op_outerface_3d(`double precision',0,1,2)dnl
c
      subroutine conrefoutfacedoub3d1(
conref_op_outerface_3d(`double precision',1,2,0)dnl
c
      subroutine conrefoutfacedoub3d2(
conref_op_outerface_3d(`double precision',2,0,1)dnl
c
c***********************************************************************
c Constant interpolation for 3d outerface float data
c***********************************************************************
c
      subroutine conrefoutfaceflot3d0(
conref_op_outerface_3d(`real',0,1,2)dnl
c
      subroutine conrefoutfaceflot3d1(
conref_op_outerface_3d(`real',1,2,0)dnl
c
      subroutine conrefoutfaceflot3d2(
conref_op_outerface_3d(`real',2,0,1)dnl
c
c***********************************************************************
c Constant interpolation for 3d outerface complex data
c***********************************************************************
c
      subroutine conrefoutfacecplx3d0(
conref_op_outerface_3d(`double complex',0,1,2)dnl
c
      subroutine conrefoutfacecplx3d1(
conref_op_outerface_3d(`double complex',1,2,0)dnl
c
      subroutine conrefoutfacecplx3d2(
conref_op_outerface_3d(`double complex',2,0,1)dnl
c
c***********************************************************************
c Constant interpolation for 3d outerface integer data
c***********************************************************************
c
      subroutine conrefoutfaceintg3d0(
conref_op_outerface_3d(`integer',0,1,2)dnl
c
      subroutine conrefoutfaceintg3d1(
conref_op_outerface_3d(`integer',1,2,0)dnl
c
      subroutine conrefoutfaceintg3d2(
conref_op_outerface_3d(`integer',2,0,1)dnl
c
c***********************************************************************
c Constant interpolation for 3d side-centered double data
c***********************************************************************
c
      subroutine conrefsidedoub3d0(
conref_op_side_3d(`double precision',0,1,2)dnl
c
      subroutine conrefsidedoub3d1(
conref_op_side_3d(`double precision',1,2,0)dnl
c
      subroutine conrefsidedoub3d2(
conref_op_side_3d(`double precision',2,0,1)dnl
c
c***********************************************************************
c Constant interpolation for 3d side-centered float data
c***********************************************************************
c
      subroutine conrefsideflot3d0(
conref_op_side_3d(`real',0,1,2)dnl
c
      subroutine conrefsideflot3d1(
conref_op_side_3d(`real',1,2,0)dnl
c
      subroutine conrefsideflot3d2(
conref_op_side_3d(`real',2,0,1)dnl
c
c***********************************************************************
c Constant interpolation for 3d side-centered complex data
c***********************************************************************
c
      subroutine conrefsidecplx3d0(
conref_op_side_3d(`double complex',0,1,2)dnl
c
      subroutine conrefsidecplx3d1(
conref_op_side_3d(`double complex',1,2,0)dnl
c
      subroutine conrefsidecplx3d2(
conref_op_side_3d(`double complex',2,0,1)dnl
c
c***********************************************************************
c Constant interpolation for 3d side-centered integer data
c***********************************************************************
c
      subroutine conrefsideintg3d0(
conref_op_side_3d(`integer',0,1,2)dnl
c
      subroutine conrefsideintg3d1(
conref_op_side_3d(`integer',1,2,0)dnl
c
      subroutine conrefsideintg3d2(
conref_op_side_3d(`integer',2,0,1)dnl
c

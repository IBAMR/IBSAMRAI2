c
c  File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-4-4/source/patchdata/fortran/pdat_dbugfort1d.m4 $
c  Package:     SAMRAI patchdata
c  Copyright:   (c) 1997-2008 Lawrence Livermore National Security, LLC
c  Revision:    $LastChangedRevision: 1917 $
c  Modified:    $LastChangedDate: 2008-01-25 13:28:01 -0800 (Fri, 25 Jan 2008) $
c  Description: F77 routines for debugging 1d patch data types.
c
include(pdat_dbugstuff1d.i)dnl
c
c***********************************************************************
c Debugging routines for 1d cell-centered data
c***********************************************************************
c
      subroutine dbugcelldoub1d(
pdat_debug_cell_1d(`double precision')dnl
c
      subroutine dbugcellflot1d(
pdat_debug_cell_1d(`real')dnl
c
      subroutine dbugcellcplx1d(
pdat_debug_cell_1d(`double complex')dnl
c
c
c***********************************************************************
c Debugging routines for 1d face-centered data
c***********************************************************************
c
      subroutine dbugfacedoub1d(
pdat_debug_face_1d(`double precision')dnl
c
      subroutine dbugfaceflot1d(
pdat_debug_face_1d(`real')dnl
c
      subroutine dbugfacecplx1d(
pdat_debug_face_1d(`double complex')dnl
c
c
c***********************************************************************
c Debugging routines for 1d node-centered data
c***********************************************************************
c
      subroutine dbugnodedoub1d(
pdat_debug_node_1d(`double precision')dnl
c
      subroutine dbugnodeflot1d(
pdat_debug_node_1d(`real')dnl
c
      subroutine dbugnodecplx1d(
pdat_debug_node_1d(`double complex')dnl
c
c
c***********************************************************************
c Debugging routines for 1d outerface data
c***********************************************************************
c
      subroutine dbugoutfacedoub1d(
pdat_debug_outerface_1d(`double precision')dnl
c
      subroutine dbugoutfaceflot1d(
pdat_debug_outerface_1d(`real')dnl
c
      subroutine dbugoutfacecplx1d(
pdat_debug_outerface_1d(`double complex')dnl
c
c
c***********************************************************************
c Debugging routines for 1d outerside data
c***********************************************************************
c
      subroutine dbugoutsidedoub1d(
pdat_debug_outerside_1d(`double precision')dnl
c
      subroutine dbugoutsideflot1d(
pdat_debug_outerside_1d(`real')dnl
c
      subroutine dbugoutsidecplx1d(
pdat_debug_outerside_1d(`double complex')dnl
c
c
c***********************************************************************
c Debugging routines for 1d side-centered data
c***********************************************************************
c
      subroutine dbugsidedoub1d(
pdat_debug_side_1d(`double precision')dnl
c
      subroutine dbugsideflot1d(
pdat_debug_side_1d(`real')dnl
c
      subroutine dbugsidecplx1d(
pdat_debug_side_1d(`double complex')dnl
c
c
c***********************************************************************
c Debugging routines for 1d edge-centered data
c***********************************************************************
c
      subroutine dbugedgedoub1d(
pdat_debug_edge_1d(`double precision')dnl
c
      subroutine dbugedgeflot1d(
pdat_debug_edge_1d(`real')dnl
c
      subroutine dbugedgecplx1d(
pdat_debug_edge_1d(`double complex')dnl
c

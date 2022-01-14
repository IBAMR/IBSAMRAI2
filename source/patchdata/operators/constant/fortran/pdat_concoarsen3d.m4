c
c  File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-4-4/source/patchdata/operators/constant/fortran/pdat_concoarsen3d.m4 $
c  Package:     SAMRAI patchdata
c  Copyright:   (c) 1997-2008 Lawrence Livermore National Security, LLC
c  Revision:    $LastChangedRevision: 1917 $
c  Modified:    $LastChangedDate: 2008-01-25 13:28:01 -0800 (Fri, 25 Jan 2008) $
c  Description: FORTRAN routines for spatial coarsening of 3d patch data
c               on a regular Cartesian mesh.
c
include(pdat_m4concoarsenops3d.i)dnl
c
c***********************************************************************
c Constant coarsening for 3d node-centered double data
c***********************************************************************
c
      subroutine conavgnodedoub3d(
conavg_op_node_3d(`double precision')dnl
c
c***********************************************************************
c Constant coarsening for 3d node-centered float data 
c***********************************************************************
c
      subroutine conavgnodeflot3d(
conavg_op_node_3d(`real')dnl
c
c***********************************************************************
c Constant coarsening for 3d node-centered complex data
c***********************************************************************
c
      subroutine conavgnodecplx3d(
conavg_op_node_3d(`double complex')dnl
c
c***********************************************************************
c Constant coarsening for 3d node-centered integer data
c***********************************************************************
c
      subroutine conavgnodeintg3d(
conavg_op_node_3d(`integer')dnl
c

c***********************************************************************
c Constant coarsening for 3d outernode-centered double data
c***********************************************************************
c
      subroutine conavgouternodedoub3d0(
conavg_op_outernode_3d0(`double precision')dnl

      subroutine conavgouternodedoub3d1(
conavg_op_outernode_3d1(`double precision')dnl

      subroutine conavgouternodedoub3d2(
conavg_op_outernode_3d2(`double precision')dnl

c***********************************************************************
c Constant coarsening for 3d outernode-centered float data
c***********************************************************************
c
      subroutine conavgouternodeflot3d0(
conavg_op_outernode_3d0(`real')dnl

      subroutine conavgouternodeflot3d1(
conavg_op_outernode_3d1(`real')dnl

      subroutine conavgouternodeflot3d2(
conavg_op_outernode_3d2(`real')dnl

c
c***********************************************************************
c Constant coarsening for 3d outernode-centered complex data
c***********************************************************************
c
      subroutine conavgouternodecplx3d0(
conavg_op_outernode_3d0(`complex')dnl

      subroutine conavgouternodecplx3d1(
conavg_op_outernode_3d1(`complex')dnl

      subroutine conavgouternodecplx3d2(
conavg_op_outernode_3d2(`complex')dnl


c
c***********************************************************************
c Constant coarsening for 3d outernode-centered integer data
c***********************************************************************
c
      subroutine conavgouternodeint3d0(
conavg_op_outernode_3d0(`integer')dnl

      subroutine conavgouternodeint3d1(
conavg_op_outernode_3d1(`integer')dnl

      subroutine conavgouternodeint3d2(
conavg_op_outernode_3d2(`integer')dnl




c
c  File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-4-4/source/patchdata/operators/constant/fortran/pdat_concoarsen2d.m4 $
c  Package:     SAMRAI patchdata
c  Copyright:   (c) 1997-2008 Lawrence Livermore National Security, LLC
c  Revision:    $LastChangedRevision: 1917 $
c  Modified:    $LastChangedDate: 2008-01-25 13:28:01 -0800 (Fri, 25 Jan 2008) $
c  Description: FORTRAN routines for spatial coarsening of 2d patch data
c               on a regular Cartesian mesh.
c
include(pdat_m4concoarsenops2d.i)dnl
c
c***********************************************************************
c Constant coarsening for 2d node-centered double data
c***********************************************************************
c
      subroutine conavgnodedoub2d(
conavg_op_node_2d(`double precision')dnl
c
c***********************************************************************
c Constant coarsening for 2d node-centered float data 
c***********************************************************************
c
      subroutine conavgnodeflot2d(
conavg_op_node_2d(`real')dnl
c
c***********************************************************************
c Constant coarsening for 2d node-centered complex data
c***********************************************************************
c
      subroutine conavgnodecplx2d(
conavg_op_node_2d(`double complex')dnl
c
c***********************************************************************
c Constant coarsening for 2d node-centered integer data
c***********************************************************************
c
      subroutine conavgnodeintg2d(
conavg_op_node_2d(`integer')dnl
c
c

c***********************************************************************
c Constant coarsening for 2d outernode-centered double data
c***********************************************************************
c
      subroutine conavgouternodedoub2d0(
conavg_op_outernode_2d0(`double precision')dnl

      subroutine conavgouternodedoub2d1(
conavg_op_outernode_2d1(`double precision')dnl

c***********************************************************************
c Constant coarsening for 2d outernode-centered float data
c***********************************************************************
c
      subroutine conavgouternodeflot2d0(
conavg_op_outernode_2d0(`real')dnl

      subroutine conavgouternodeflot2d1(
conavg_op_outernode_2d1(`real')dnl

c
c***********************************************************************
c Constant coarsening for 2d outernode-centered complex data
c***********************************************************************
c
      subroutine conavgouternodecplx2d0(
conavg_op_outernode_2d0(`complex')dnl

      subroutine conavgouternodecplx2d1(
conavg_op_outernode_2d1(`complex')dnl

c
c***********************************************************************
c Constant coarsening for 2d outernode-centered integer data
c***********************************************************************
c
      subroutine conavgouternodeint2d0(
conavg_op_outernode_2d0(`integer')dnl

      subroutine conavgouternodeint2d1(
conavg_op_outernode_2d1(`integer')dnl

c


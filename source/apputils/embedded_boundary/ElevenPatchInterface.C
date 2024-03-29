//
// File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-4-4/source/apputils/embedded_boundary/ElevenPatchInterface.C $
// Package:     SAMRAI 
//              Structured Adaptive Mesh Refinement Applications Infrastructure
// Copyright:   (c) 1997-2008 Lawrence Livermore National Security, LLC
// Release:     $Name:  $
// Revision:    $LastChangedRevision: 1917 $
// Modified:    $LastChangedDate: 2008-01-25 13:28:01 -0800 (Fri, 25 Jan 2008) $
// Description: SAMRAI interface to Eleven library
//              
//
// THIS CLASS IS CURRENTLY EMPTY BECAUSE INTERFACES TO ELEVEN ARE
// STILL BEING DEVELOPED.  PLEASE CONTACT SAMRAI DEVELOPERS IF YOU
// ARE INTERESTED IN USING ELEVEN INTERFACES.
// 
#ifndef included_ElevenPatchInterface_C
#define included_ElevenPatchInterface_C

#include "ElevenPatchInterface.h"

#include "CartesianPatchGeometry.h"
#include "CellData.h"
#include "NodeData.h"
#include "CutCell.h"
#include "EmbeddedBoundaryDefines.h"
#include "IndexData.h"
#include "tbox/Utilities.h"



namespace SAMRAI {
  namespace appu {

/*
*******************************************************************
*
*  Constructor
*
*******************************************************************
*/
template<int DIM> 
ElevenPatchInterface<DIM>::ElevenPatchInterface(
   const std::string& object_name,
   tbox::Pointer<tbox::Database> input_db)
{
  d_object_name = object_name;

  getFromInput(input_db);
}

/*
*******************************************************************
*
*  Destructor
*
*******************************************************************
*/
template<int DIM> 
ElevenPatchInterface<DIM>::~ElevenPatchInterface()
{
  
}

/*
*******************************************************************
*
*  Return eleven patch information
*
*******************************************************************
*/
template<int DIM>
void ElevenPatchInterface<DIM>::calculateCutCellInfo(
   tbox::Pointer<hier::Patch<DIM> >& patch,
   const int cell_flag_data_id,
   const int cell_vol_data_id,
   const int node_flag_data_id,
   const int cutcell_data_id)
{
   NULL_USE(patch);
   NULL_USE(cell_flag_data_id);
   NULL_USE(cell_vol_data_id);
   NULL_USE(node_flag_data_id);
   NULL_USE(cutcell_data_id);
   
   TBOX_ERROR(d_object_name << ":Unable to use ELEVEN."
              << "\nPlease contact SAMRAI developers to find out more"
              << "\ninformation." << std::endl);
}

/*
*******************************************************************
*
*  Check whether a particular node is inside the shape
*
*******************************************************************
*/
template<int DIM>
bool ElevenPatchInterface<DIM>::isInside(
   const double* xyz) const
{
   NULL_USE(xyz);

   TBOX_ERROR(d_object_name << ":Unable to use ELEVEN."
              << "\nPlease contact SAMRAI developers to find out more"
              << "\ninformation." << std::endl);

  return true;
}

/*
*******************************************************************
*
*  Mark a SET of nodes on a patch as being inside or outside the 
*  shape. 
*
*     nx[DIM] is the dimensions of the patch 
*     dx[DIM] is the delta x in each direction of the patch 
*     origin[DIM] is the location of the lower left node of the 
*                 patch
*     inout[nx[0]*nx[1]*nx[2]] is the array holding whether the
*                 nodes are inside or outside (1=inside, 0=outside) 
*
*******************************************************************
*/
template<int DIM>
void ElevenPatchInterface<DIM>::isInside(
   const int* nx,
   const double* dx,
   const double* origin,
   int* inout) const
{
   NULL_USE(nx);
   NULL_USE(dx);
   NULL_USE(origin);
   NULL_USE(inout);
   
   TBOX_ERROR(d_object_name << ":Unable to use ELEVEN."
              << "\nPlease contact SAMRAI developers to find out more"
              << "\ninformation." << std::endl);
}


/*
*******************************************************************
*
*  Read info from input
*
*******************************************************************
*/
template<int DIM>
void ElevenPatchInterface<DIM>::getFromInput(
   tbox::Pointer<tbox::Database> db)
{
   NULL_USE(db);

#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!db.isNull());
#endif

   TBOX_ERROR(d_object_name << ":Unable to use ELEVEN."
              << "\nPlease contact SAMRAI developers to find out more"
              << "\ninformation." << std::endl);

}

/*
*******************************************************************
*
*  Dump (to supplied os) class information
*
*******************************************************************
*/
template<int DIM>
void ElevenPatchInterface<DIM>::printClassData(
   std::ostream& os) const
{
   os << "d_object_name = " << d_object_name << std::endl;
}


}
}
#endif


//
// File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-4-4/examples/communication/OuternodeDataTest.h $
// Package:     SAMRAI tests
// Copyright:   (c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:    $LastChangedRevision: 2132 $
// Modified:    $LastChangedDate: 2008-04-14 14:51:47 -0700 (Mon, 14 Apr 2008) $
// Description: AMR communication tests for node-centered patch data
//

#ifndef included_pdat_OuternodeDataTest
#define included_pdat_OuternodeDataTest

#include "SAMRAI_config.h"

#include "tbox/Array.h"
#include "Box.h"
#include "CartesianGridGeometry.h"
#include "CartesianPatchGeometry.h"
#include "tbox/Database.h"
#include "IntVector.h"
#include "NodeData.h"
#include "OuternodeData.h"
#include "Patch.h"
#include "PatchDataTestStrategy.h"
#include "tbox/Pointer.h"
#ifndef included_String
#include <string>
using namespace std;
#define included_String
#endif
#include "Variable.h"

namespace SAMRAI {

class CommTester;

/**
 * Class OuternodeDataTest provides routines to test communication operations
 * for node-centered patch data on an AMR patch hierarchy.
 *
 * Required input keys and data types:
 * 

 
 *
 *   Double values that define linear function initial data to test refine
 *   operations (Ax + By + Cz + D = f(x,y,z), where f(x,y,z) is the value 
 *   assigned to each array value at initialization and against which 
 *   linear interpolation is tested: 
 *  
 *    Acoef, Dcoef always required.
 *    If (NDIM > 1), Bcoef is needed.
 *    If (NDIM > 2), Ccoef is needed.
 *
 * 


 *
 * See PatchDataTestStrategy header file comments for variable and
 * refinement input data description.
 */

class OuternodeDataTest : public PatchDataTestStrategy
{
public:
  /**
   * The constructor initializes variable data arrays to zero length.
   */
   OuternodeDataTest(const string& object_name,
                tbox::Pointer<tbox::Database> main_input_db,
                bool do_refine,
                bool do_coarsen,
                const string& refine_option);

   /**
    * Virtual destructor for OuternodeDataTest.
    */
   ~OuternodeDataTest();

   /**
    * User-supplied boundary conditions.  Note that we do not implement
    * user-defined coarsen and refine operations.
    */
   virtual void setPhysicalBoundaryConditions(hier::Patch<NDIM>& patch,
                                              const double time,
                                              const hier::IntVector<NDIM>&) const;

   /**
    * This function is called from the CommTester constructor.  Its
    * purpose is to register variables used in the patch data test
    * and appropriate communication parameters (ghost node widths,
    * coarsen/refine operations) with the CommTester object, which
    * manages the variable storage.
    */
   void registerVariables(CommTester* commtest);

   /**
    * Function for setting data on new patch in hierarchy.
    *
    * @param src_or_dst Flag set to 's' for source or 'd' for destination
    *        to indicate variables to set data for.
    */
   virtual void initializeDataOnPatch(hier::Patch<NDIM>& patch,
                                      const tbox::Pointer<hier::PatchHierarchy<NDIM> > hierarchy,
                                      int level_number,
				      char src_or_dst);

   /**
    * Function for checking results of communication operations.
    */
   bool verifyResults(hier::Patch<NDIM>& patch,
                      const tbox::Pointer<hier::PatchHierarchy<NDIM> > hierarchy,
                      int level_number);

private:
   /**
    * Function for reading test data from input file.
    */
   void readTestInput(tbox::Pointer<tbox::Database> db);

   /**
    * Set linear function data for testing interpolation
    */
   void setLinearData(tbox::Pointer< pdat::OuternodeData<NDIM,double> > data,
                      const hier::Box<NDIM>& box,
                      hier::Patch<NDIM>& patch) const;

   /**
    * Set linear function data for testing interpolation
    */
   void setLinearData(tbox::Pointer< pdat::NodeData<NDIM,double> > data,
                      const hier::Box<NDIM>& box,
                      hier::Patch<NDIM>& patch) const;

   void checkPatchInteriorData(
      const tbox::Pointer< pdat::OuternodeData<NDIM,double> >& data,
      const hier::Box<NDIM>& interior,
      const tbox::Pointer<geom::CartesianPatchGeometry<NDIM> >& pgeom) const;

   /*
    * Object string identifier for error reporting
    */
   string d_object_name;

   /*
    * Data members specific to this node data test.
    */
   tbox::Pointer<geom::CartesianGridGeometry<NDIM> > d_cart_grid_geometry;

   /*
    * Data members specific to this node data test.
    */
   double d_Acoef;
   double d_Bcoef;
   double d_Ccoef;
   double d_Dcoef;

   bool d_do_refine;
   bool d_do_coarsen;
   string d_refine_option;
   int d_finest_level_number;

   tbox::Array< tbox::Pointer<hier::Variable<NDIM> > > d_variables_src;
   tbox::Array< tbox::Pointer<hier::Variable<NDIM> > > d_variables_dst;

};

}
#endif
